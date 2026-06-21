#!/usr/bin/env python3
"""
Create SAIGE group files for gene set (rare variant) analysis from a VEP-annotated VCF.

This script processes a VCF file (annotated with VEP) and a list of genes to create
group files suitable for SAIGE rare variant analysis. Variants are clustered into
annotation tiers using one of three schemes:

  IMPACT   - VEP IMPACT (HIGH/MODERATE/LOW/MODIFIER)
  LOFTEE   - LOFTEE-aware consequences (pLoF/missense_LC/synonymous)
  ENSEMBLE - ensemble missense/splice tiers
             (pLoF/damMis/weakMis/neutralMis/synonymous)

Output (one file): SAIGE group format, two lines per gene:
  {gene} var  v1 v2 v3 ...
  {gene} anno a1 a2 a3 ...
"""

import os
import argparse
import sys
import gzip
import re

# --- IMPACT Mappings ---
IMPACT_SCORES = {
    "HIGH": 4,
    "MODERATE": 3,
    "LOW": 2,
    "MODIFIER": 1,
    "": 0,
}
SCORE_TO_IMPACT = {v: k for k, v in IMPACT_SCORES.items() if k != ""}

# --- LOFTEE Mappings ---
LOFTEE_SCORES = {
    "pLoF": 3,
    "missense_LC": 2,
    "synonymous": 1,
    "": 0,
}
SCORE_TO_LOFTEE = {v: k for k, v in LOFTEE_SCORES.items() if k != ""}

# --- ENSEMBLE Mappings ---
ENSEMBLE_SCORES = {
    "pLoF": 5,
    "damMis": 4,
    "weakMis": 3,
    "neutralMis": 2,
    "synonymous": 1,
    "": 0,
}
SCORE_TO_ENSEMBLE = {v: k for k, v in ENSEMBLE_SCORES.items() if k != ""}

REVERSE_MAPS = {
    'IMPACT': SCORE_TO_IMPACT,
    'LOFTEE': SCORE_TO_LOFTEE,
    'ENSEMBLE': SCORE_TO_ENSEMBLE,
}

PLOF_CONSEQUENCES = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant"}
MISSENSE_CONSEQUENCES = {"missense_variant", "start_lost", "stop_lost", "inframe_insertion", "inframe_deletion", "protein_altering_variant"}
SYNONYMOUS_CONSEQUENCES = {"synonymous_variant", "stop_retained_variant"}
CLINVAR_PATHOGENIC = {"pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic"}

ENSEMBLE_FIELD_CANDIDATES = {
    "spliceai_max": ["SpliceAI_DS_max"],
    "spliceai_ds": [
        "SpliceAI_pred_DS_AG",
        "SpliceAI_pred_DS_AL",
        "SpliceAI_pred_DS_DG",
        "SpliceAI_pred_DS_DL",
        "SpliceAI_DS_AG",
        "SpliceAI_DS_AL",
        "SpliceAI_DS_DG",
        "SpliceAI_DS_DL",
    ],
    "revel": ["REVEL", "REVEL_score"],
    "primateai": ["PrimateAI3D_score", "PrimateAI_score"],
    "alphamissense": ["AlphaMissense_score"],
    "alphamissense_class": ["am_class", "AlphaMissense_pred", "AlphaMissense_class"],
    "clinvar_sig": ["ClinVar_CLNSIG", "CLNSIG", "CLIN_SIG"],
    "clinvar_revstat": ["ClinVar_CLNREVSTAT", "CLNREVSTAT", "CLINREVSTAT"],
}

AUTOSOMES = set(str(c) for c in range(1, 23))


def _find_idx(csq_format, names):
    for name in names:
        if name in csq_format:
            return csq_format.index(name)
    return None


def _find_indices(csq_format, names):
    return [csq_format.index(name) for name in names if name in csq_format]


def _csq_value(annos, idx):
    if idx is None or len(annos) <= idx:
        return ""
    return annos[idx]


def _iter_plugin_values(value):
    for token in re.split(r"[&,;]", value):
        token = token.strip()
        if token and token != ".":
            yield token


def _max_float(values):
    best = None
    for value in values:
        for token in _iter_plugin_values(value):
            try:
                score = float(token)
            except ValueError:
                continue
            if best is None or score > best:
                best = score
    return best


def _max_csq_float(annos, indices):
    return _max_float(_csq_value(annos, idx) for idx in indices)


def _clinvar_is_pathogenic(value):
    for token in _iter_plugin_values(value.replace(" ", "_")):
        token = token.lower()
        if token in CLINVAR_PATHOGENIC:
            return True
    return False


def _alphamissense_is_likely_pathogenic(value):
    return any(token.lower() == "likely_pathogenic" for token in _iter_plugin_values(value))


def _warn_missing_ensemble_fields(ensemble_indices):
    if not ensemble_indices["spliceai_indices"]:
        print("Warning: ENSEMBLE mode could not find SpliceAI DS fields.", file=sys.stderr)
    if ensemble_indices["revel_idx"] is None:
        print("Warning: ENSEMBLE mode could not find REVEL score field.", file=sys.stderr)
    if ensemble_indices["primateai_idx"] is None:
        print("Warning: ENSEMBLE mode could not find PrimateAI score field.", file=sys.stderr)
    if ensemble_indices["alphamissense_idx"] is None and ensemble_indices["alphamissense_class_idx"] is None:
        print("Warning: ENSEMBLE mode could not find AlphaMissense score/class field.", file=sys.stderr)
    if ensemble_indices["clinvar_sig_idx"] is None:
        print("Warning: ENSEMBLE mode could not find ClinVar significance field.", file=sys.stderr)


def parse_vcf_header(vcf_handle, anno_type):
    """
    Parse VCF header to find the CSQ field indices needed for the given annotation type.
    """
    csq_format = None

    for line in vcf_handle:
        if line.startswith("#CHROM"):
            break

        if line.startswith("##INFO=<ID=CSQ"):
            try:
                description = line.split('Description="')[1].split('"')[0]
                if "Format: " in description:
                    format_str = description.split("Format: ")[1].strip()
                    csq_format = format_str.split('|')
            except IndexError:
                continue

    if not csq_format:
        print("Error: Could not find CSQ Format in VCF header. "
              "Ensure VCF is VEP-annotated.", file=sys.stderr)
        sys.exit(1)

    try:
        symbol_idx = csq_format.index("SYMBOL")
    except ValueError:
        sys.exit("Error: Could not find SYMBOL in CSQ format.")

    impact_idx = None
    consequence_idx = None
    lof_idx = None
    ensemble_indices = {}

    if anno_type == 'IMPACT':
        try:
            impact_idx = csq_format.index("IMPACT")
        except ValueError:
            sys.exit("Error: Could not find IMPACT in CSQ format.")

    if anno_type in {'LOFTEE', 'ENSEMBLE'}:
        try:
            consequence_idx = csq_format.index("Consequence")
        except ValueError:
            sys.exit(f"Error: Could not find Consequence in CSQ format (required for {anno_type} mode).")

        if "LoF" in csq_format:
            lof_idx = csq_format.index("LoF")

    if anno_type == 'ENSEMBLE':
        impact_idx = _find_idx(csq_format, ["IMPACT"])
        spliceai_max_idx = _find_idx(csq_format, ENSEMBLE_FIELD_CANDIDATES["spliceai_max"])
        spliceai_indices = _find_indices(csq_format, ENSEMBLE_FIELD_CANDIDATES["spliceai_ds"])
        if spliceai_max_idx is not None:
            spliceai_indices.append(spliceai_max_idx)

        ensemble_indices = {
            "spliceai_indices": spliceai_indices,
            "revel_idx": _find_idx(csq_format, ENSEMBLE_FIELD_CANDIDATES["revel"]),
            "primateai_idx": _find_idx(csq_format, ENSEMBLE_FIELD_CANDIDATES["primateai"]),
            "alphamissense_idx": _find_idx(csq_format, ENSEMBLE_FIELD_CANDIDATES["alphamissense"]),
            "alphamissense_class_idx": _find_idx(csq_format, ENSEMBLE_FIELD_CANDIDATES["alphamissense_class"]),
            "clinvar_sig_idx": _find_idx(csq_format, ENSEMBLE_FIELD_CANDIDATES["clinvar_sig"]),
            "clinvar_revstat_idx": _find_idx(csq_format, ENSEMBLE_FIELD_CANDIDATES["clinvar_revstat"]),
        }
        _warn_missing_ensemble_fields(ensemble_indices)

    try:
        gnomad_af_idx = csq_format.index("gnomADe_AF")
    except ValueError:
        gnomad_af_idx = None

    return (
        symbol_idx, impact_idx, consequence_idx, lof_idx,
        gnomad_af_idx, ensemble_indices
    )


def _spliceai_tier(score):
    """Map a SpliceAI delta score to an ENSEMBLE tier using the author-recommended
    0.2 (high recall) / 0.5 (recommended) / 0.8 (high precision) cutoffs."""
    if score is None:
        return 0
    if score >= 0.8:
        return ENSEMBLE_SCORES["pLoF"]
    if score >= 0.5:
        return ENSEMBLE_SCORES["damMis"]
    if score >= 0.2:
        return ENSEMBLE_SCORES["weakMis"]
    return 0


def classify_ensemble_annotation(annos, impact_idx, consequence_idx, lof_idx, ensemble_indices):
    consequences = set(_csq_value(annos, consequence_idx).split('&'))
    impact = _csq_value(annos, impact_idx)
    lof_val = _csq_value(annos, lof_idx)
    spliceai_score = _max_csq_float(annos, ensemble_indices["spliceai_indices"])
    clinvar_sig = _csq_value(annos, ensemble_indices["clinvar_sig_idx"])

    splice_score = _spliceai_tier(spliceai_score)
    is_loftee_hc_plof = lof_val == "HC" and bool(consequences.intersection(PLOF_CONSEQUENCES))
    is_clinvar_pathogenic = _clinvar_is_pathogenic(clinvar_sig)

    # pLoF: confident LoF annotation, ClinVar P/LP, or high-precision splice disruption
    if is_loftee_hc_plof or is_clinvar_pathogenic or splice_score == ENSEMBLE_SCORES["pLoF"]:
        return "pLoF", ENSEMBLE_SCORES["pLoF"]

    is_loftee_lc = lof_val == "LC"
    is_missense_candidate = (
        impact == "MODERATE"
        or bool(consequences.intersection(MISSENSE_CONSEQUENCES))
        or is_loftee_lc
    )
    missense_score = 0
    if is_missense_candidate:
        vote_count = 0
        revel_score = _max_csq_float(annos, [ensemble_indices["revel_idx"]])
        primateai_score = _max_csq_float(annos, [ensemble_indices["primateai_idx"]])
        alphamissense_score = _max_csq_float(annos, [ensemble_indices["alphamissense_idx"]])
        alphamissense_class = _csq_value(annos, ensemble_indices["alphamissense_class_idx"])

        if revel_score is not None and revel_score >= 0.75:
            vote_count += 1
        if primateai_score is not None and primateai_score >= 0.75:
            vote_count += 1
        if (alphamissense_score is not None and alphamissense_score >= 0.564) \
                or _alphamissense_is_likely_pathogenic(alphamissense_class):
            vote_count += 1

        if vote_count >= 2:
            missense_score = ENSEMBLE_SCORES["damMis"]
        elif vote_count == 1 or is_loftee_lc:
            missense_score = ENSEMBLE_SCORES["weakMis"]
        else:
            missense_score = ENSEMBLE_SCORES["neutralMis"]

    # Splice and missense evidence feed the same tiers; take the stronger.
    # This applies splice scoring to ANY consequence (synonymous/intronic included)
    # and lets splice disruption upgrade a benign-scoring missense.
    best = max(splice_score, missense_score)
    if best > 0:
        return SCORE_TO_ENSEMBLE[best], best

    # True synonymous control: synonymous consequence with no predicted splice effect.
    if consequences.intersection(SYNONYMOUS_CONSEQUENCES):
        if spliceai_score is not None and spliceai_score < 0.2:
            return "synonymous", ENSEMBLE_SCORES["synonymous"]

    return "", 0


def classify_transcript(annos, anno_type, impact_idx, consequence_idx, lof_idx,
                        ensemble_indices):
    """
    Classify a single CSQ transcript annotation and return its tier score
    (0 means "no tier / excluded").
    """
    if anno_type == 'IMPACT':
        impact = _csq_value(annos, impact_idx)
        return IMPACT_SCORES.get(impact, 0)

    if anno_type == 'LOFTEE':
        consequences = set(_csq_value(annos, consequence_idx).split('&'))
        lof_val = _csq_value(annos, lof_idx)
        category = ""

        # pLoF (downgrade Low-Confidence pLoFs to missense_LC)
        if consequences.intersection(PLOF_CONSEQUENCES):
            if lof_idx is not None:
                if lof_val == 'HC':
                    category = "pLoF"
                elif lof_val == 'LC':
                    category = "missense_LC"
                else:
                    category = "pLoF"  # LOFTEE field present but empty for a pLoF
            else:
                category = "pLoF"

        if not category and consequences.intersection(MISSENSE_CONSEQUENCES):
            category = "missense_LC"

        if not category and consequences.intersection(SYNONYMOUS_CONSEQUENCES):
            category = "synonymous"

        return LOFTEE_SCORES.get(category, 0)

    # ENSEMBLE
    _category, score = classify_ensemble_annotation(
        annos, impact_idx, consequence_idx, lof_idx, ensemble_indices
    )
    return score


def _chrom_sort_rank(chrom):
    try:
        return int(chrom)
    except ValueError:
        return {"X": 23, "Y": 24, "MT": 25, "M": 25}.get(chrom, 99)


def _variant_sort_key(var_id):
    parts = var_id.split(":")
    return (_chrom_sort_rank(parts[0]), int(parts[1]), parts[2], parts[3])


def get_gene_variants(vcf_file, gene_set, anno_type, autosomes_only=False,
                      gnomad_only=False, variant_allowlist=None):
    """
    Reads the VCF line by line, classifies each variant per gene, and groups
    variants by gene.

    Returns:
        gene_data: Dict { gene_name: { var_id: max_tier_score } }
    """
    gene_data = {gene: {} for gene in gene_set}
    n_excluded_no_gnomad = 0
    n_excluded_not_in_list = 0

    opener = gzip.open if vcf_file.endswith(".gz") else open

    print(f"Parsing VCF file: {vcf_file} using {anno_type} annotations...",
          file=sys.stderr)

    with opener(vcf_file, 'rt') as f:
        (
            symbol_idx, impact_idx, consequence_idx, lof_idx,
            gnomad_af_idx, ensemble_indices
        ) = parse_vcf_header(f, anno_type)

        if gnomad_only and gnomad_af_idx is None:
            sys.exit("Error: --gnomad-only specified but gnomADe_AF not found in CSQ.")

        # Reset file pointer to beginning for data reading
        f.seek(0)

        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split('\t')

            # Standard VCF columns: 0=CHROM, 1=POS, 2=ID, 3=REF, 4=ALT, 7=INFO
            chrom = parts[0].replace("chr", "")

            if autosomes_only and chrom not in AUTOSOMES:
                continue

            pos = parts[1]
            ref = parts[3]
            alt = parts[4]
            info = parts[7]

            # Construct Variant ID (chr:pos:ref:alt)
            var_id = f"{chrom}:{pos}:{ref}:{alt}"

            # Parse AC and CSQ from INFO
            ac = None
            csq_str = None
            for field in info.split(';'):
                if field.startswith("CSQ="):
                    csq_str = field[4:]
                elif gnomad_only and field.startswith("AC="):
                    ac = int(field[3:].split(',')[0])

            if not csq_str:
                continue

            # CSQ can have multiple transcripts separated by comma
            in_gnomad = False
            variant_in_genes = False
            pending_gene_impacts = []

            for trans in csq_str.split(','):
                annos = trans.split('|')

                if len(annos) <= symbol_idx:
                    continue

                # Check gnomADe_AF for presence in gnomAD
                if gnomad_af_idx is not None and not in_gnomad:
                    if len(annos) > gnomad_af_idx and annos[gnomad_af_idx]:
                        in_gnomad = True

                gene = annos[symbol_idx]
                if gene not in gene_set:
                    continue

                score = classify_transcript(
                    annos, anno_type, impact_idx, consequence_idx,
                    lof_idx, ensemble_indices
                )
                if score > 0:
                    variant_in_genes = True
                    pending_gene_impacts.append((gene, score))

            # If gnomad-only: keep singletons, exclude non-singleton non-gnomAD variants
            if gnomad_only and not in_gnomad and variant_in_genes:
                if ac != 1:
                    n_excluded_no_gnomad += 1
                    continue

            # Filter to allowlist if provided
            if variant_allowlist is not None and var_id not in variant_allowlist:
                if variant_in_genes:
                    n_excluded_not_in_list += 1
                continue

            # Commit gene-variant annotations (keep the most severe tier per variant)
            for gene, score in pending_gene_impacts:
                if var_id not in gene_data[gene] or score > gene_data[gene][var_id]:
                    gene_data[gene][var_id] = score

    if gnomad_only:
        print(f"Excluded {n_excluded_no_gnomad} non-singleton variants not found in gnomAD.",
              file=sys.stderr)
    if variant_allowlist is not None:
        print(f"Excluded {n_excluded_not_in_list} variants not in the supplied variant list.",
              file=sys.stderr)

    return gene_data


def write_output(gene_data, anno_type, output_file):
    """
    Writes the gene data to the output file in SAIGE group format.
    """
    reverse_map = REVERSE_MAPS[anno_type]
    genes_processed = 0
    total_genes = len(gene_data)
    total_variants = 0

    print(f"Writing results to: {output_file}", file=sys.stderr)

    with open(output_file, "w") as f:
        for gene, variants in gene_data.items():
            if not variants:
                continue

            genes_processed += 1
            sorted_vars = sorted(variants.keys(), key=_variant_sort_key)
            annos_list = [reverse_map[variants[v]] for v in sorted_vars]
            total_variants += len(sorted_vars)

            f.write(f"{gene} var {' '.join(sorted_vars)}\n")
            f.write(f"{gene} anno {' '.join(annos_list)}\n")

    print(f"\nResults:", file=sys.stderr)
    print(f"  Genes with variants: {genes_processed}/{total_genes}", file=sys.stderr)
    print(f"  Total variants: {total_variants}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Create SAIGE group files from a VEP-annotated VCF (using INFO/CSQ)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '-g', '--gene-list',
        required=True,
        help='Path to gene list file (one gene per line)'
    )

    parser.add_argument(
        '-v', '--vcf-file',
        required=True,
        help='Path to VCF file (can be .vcf or .vcf.gz)'
    )

    parser.add_argument(
        '-o', '--output-folder',
        required=True,
        help='Output directory for SAIGE group files'
    )

    parser.add_argument(
        '-n', '--vep-name',
        default='cohort',
        help='Name identifier for the cohort (default: cohort)'
    )

    parser.add_argument(
        '--anno-type', choices=['IMPACT', 'LOFTEE', 'ENSEMBLE'], default='IMPACT',
        help='Choose whether to cluster variants by VEP IMPACT (HIGH/MODERATE/LOW/MODIFIER), '
             'LOFTEE consequences (pLoF/missense_LC/synonymous), or ENSEMBLE tiers '
             '(pLoF/damMis/weakMis/neutralMis/synonymous). Default: IMPACT'
    )

    parser.add_argument(
        '--autosomes-only', action='store_true',
        help='Keep only autosomes (chr1-22), excluding X, Y, MT, etc.'
    )

    parser.add_argument(
        '--gnomad-only', action='store_true',
        help='Exclude non-singleton variants not found in gnomAD'
    )

    parser.add_argument(
        '-l', '--variant-list',
        default=None,
        help='Path to variant ID list file (one ID per line, format chrom:pos:ref:alt); '
             'only these variants will be included'
    )

    args = parser.parse_args()

    # Validation
    if not os.path.exists(args.gene_list):
        sys.exit(f"Error: Gene list file not found: {args.gene_list}")
    if not os.path.exists(args.vcf_file):
        sys.exit(f"Error: VCF file not found: {args.vcf_file}")
    if args.variant_list and not os.path.exists(args.variant_list):
        sys.exit(f"Error: Variant list file not found: {args.variant_list}")
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    # Get Gene List Name
    gene_list_name = os.path.splitext(os.path.basename(args.gene_list))[0]
    if gene_list_name.endswith('.txt') or gene_list_name.endswith('.tsv'):
        gene_list_name = os.path.splitext(gene_list_name)[0]

    output_filename = f"{gene_list_name}_{args.vep_name}_saige_group.txt"
    output_path = os.path.join(args.output_folder, output_filename)

    # Read Gene List
    with open(args.gene_list, 'r') as f:
        genes = set(line.strip() for line in f if line.strip())

    print(f"Processing {len(genes)} genes from list.", file=sys.stderr)

    # Read variant allowlist if provided
    variant_allowlist = None
    if args.variant_list:
        with open(args.variant_list, 'r') as f:
            variant_allowlist = set(line.strip() for line in f if line.strip())
        print(f"Loaded {len(variant_allowlist)} variants from allowlist.", file=sys.stderr)

    # Process
    gene_data = get_gene_variants(
        args.vcf_file, genes, args.anno_type,
        autosomes_only=args.autosomes_only,
        gnomad_only=args.gnomad_only,
        variant_allowlist=variant_allowlist)
    write_output(gene_data, args.anno_type, output_path)


if __name__ == "__main__":
    main()
