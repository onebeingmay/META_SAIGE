#!/usr/bin/env python3
"""Merge SAIGE group files from multiple cohorts.

A SAIGE group file has two whitespace-separated lines per gene:
  {gene} var  v1 v2 v3 ...
  {gene} anno a1 a2 a3 ...
where the i-th annotation describes the i-th variant.

Merging unions the variant set per gene across cohorts. When two cohorts assign
different annotations to the same variant in the same gene, the conflict is
resolved by annotation priority within a recognized scheme (IMPACT / LOFTEE /
ENSEMBLE), falling back to first-seen for unrecognized or mixed schemes.
"""
import argparse, os, sys

IMPACT_RANKS = {
    "HIGH": 4,
    "MODERATE": 3,
    "LOW": 2,
    "MODIFIER": 1,
}

LOFTEE_RANKS = {
    "pLoF": 3,
    "missense_LC": 2,
    "synonymous": 1,
}

ENSEMBLE_RANKS = {
    "pLoF": 5,
    "damMis": 4,
    "weakMis": 3,
    "neutralMis": 2,
    "synonymous": 1,
}

ANNOTATION_RANKS = {
    "IMPACT": IMPACT_RANKS,
    "LOFTEE": LOFTEE_RANKS,
    "ENSEMBLE": ENSEMBLE_RANKS,
}

def _chrom_rank(c):
    try: return int(c)
    except ValueError: return {"X":23,"Y":24,"MT":25,"M":25}.get(c,99)

def _var_key(v):
    p = v.split(":")
    return (_chrom_rank(p[0]), int(p[1]), p[2], p[3])

def read_saige_group(path):
    """Parse a SAIGE group file -> { gene: { var_id: annotation } }.

    Pairs the 'var' and 'anno' lines for each gene by position. Genes missing
    one of the two lines, or with mismatched lengths, are reported and skipped.
    """
    gene_vars, gene_annos = {}, {}
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) < 2:
                continue
            gene, tag, values = p[0], p[1], p[2:]
            if tag == "var":
                gene_vars[gene] = values
            elif tag == "anno":
                gene_annos[gene] = values

    d = {}
    for gene, variants in gene_vars.items():
        annos = gene_annos.get(gene)
        if annos is None:
            print(f"Warning: {gene} in {os.path.basename(path)} has a 'var' line "
                  f"but no 'anno' line; skipping.", file=sys.stderr)
            continue
        if len(annos) != len(variants):
            print(f"Warning: {gene} in {os.path.basename(path)} has "
                  f"{len(variants)} variants but {len(annos)} annotations; skipping.",
                  file=sys.stderr)
            continue
        var_map = d.setdefault(gene, {})
        for v, a in zip(variants, annos):
            # Keep the higher-priority annotation if a variant is listed twice.
            if v not in var_map:
                var_map[v] = a
            else:
                var_map[v], _scheme, _mode = _resolve_annotation(var_map[v], a)

    for gene in gene_annos:
        if gene not in gene_vars:
            print(f"Warning: {gene} in {os.path.basename(path)} has an 'anno' line "
                  f"but no 'var' line; skipping.", file=sys.stderr)
    return d

def _annotation_scheme(annotation):
    for scheme, ranks in ANNOTATION_RANKS.items():
        if annotation in ranks:
            return scheme
    return None

def _ranking_for_pair(existing, incoming):
    matches = []
    for scheme, ranks in ANNOTATION_RANKS.items():
        if existing in ranks and incoming in ranks:
            matches.append((scheme, ranks))
    if not matches:
        return None, None
    if len(matches) == 1:
        return matches[0]
    return "/".join(scheme for scheme, _ranks in matches), matches[0][1]

def _resolve_annotation(existing, incoming):
    existing_scheme = _annotation_scheme(existing)
    incoming_scheme = _annotation_scheme(incoming)
    pair_scheme, ranks = _ranking_for_pair(existing, incoming)

    if pair_scheme:
        if ranks[incoming] > ranks[existing]:
            return incoming, pair_scheme, "ranked"
        return existing, pair_scheme, "ranked"

    return existing, existing_scheme or incoming_scheme or "UNKNOWN", "first_seen"

def merge_groups(group_list, names):
    """Merge per-cohort { gene: { var: anno } } dicts into one, resolving
    annotation conflicts and unioning variants per gene."""
    merged, disc, resolved = {}, 0, 0
    for groups, name in zip(group_list, names):
        for gene, var_map in groups.items():
            target = merged.setdefault(gene, {})
            for var_id, annotation in var_map.items():
                if var_id not in target:
                    target[var_id] = annotation
                elif target[var_id] != annotation:
                    disc += 1
                    earlier = target[var_id]
                    kept, scheme, mode = _resolve_annotation(earlier, annotation)
                    if mode == "ranked":
                        resolved += 1
                        reason = f"Keeping {kept} by {scheme} priority."
                    else:
                        reason = f"Keeping {earlier} (first seen; unrecognized or mixed annotation scheme)."
                    target[var_id] = kept
                    if disc <= 20:
                        print(
                            f"Warning: {var_id} in {gene}: {annotation} ({name}) vs "
                            f"{earlier} (earlier). {reason}",
                            file=sys.stderr,
                        )
    if disc > 20: print(f"... {disc-20} more discrepancies suppressed.", file=sys.stderr)
    if disc:
        print(
            f"Total annotation discrepancies: {disc} "
            f"({resolved} resolved by annotation priority, {disc - resolved} kept as first seen).",
            file=sys.stderr,
        )
    return merged

def _gene_sort_key(gene, var_map):
    return min(_var_key(v) for v in var_map)

def sorted_genes(merged):
    return sorted(merged, key=lambda g: (_gene_sort_key(g, merged[g]), g))

def write_saige_group(merged, out):
    genes_written, total_variants = 0, 0
    with open(out, 'w') as f:
        for gene in sorted_genes(merged):
            var_map = merged[gene]
            if not var_map:
                continue
            sorted_vars = sorted(var_map.keys(), key=_var_key)
            annos = [var_map[v] for v in sorted_vars]
            f.write(f"{gene} var {' '.join(sorted_vars)}\n")
            f.write(f"{gene} anno {' '.join(annos)}\n")
            genes_written += 1
            total_variants += len(sorted_vars)
    return genes_written, total_variants

def main():
    ap = argparse.ArgumentParser(description="Merge SAIGE group files from multiple cohorts")
    ap.add_argument('-g', '--group-files', nargs='+', required=True,
                    help='Two or more SAIGE group files to merge')
    ap.add_argument('-o', '--output-file', required=True,
                    help='Path to write the merged SAIGE group file')
    args = ap.parse_args()

    for f in args.group_files:
        if not os.path.exists(f): sys.exit(f"Error: {f} not found")

    out_dir = os.path.dirname(args.output_file)
    if out_dir: os.makedirs(out_dir, exist_ok=True)

    names = [os.path.basename(f).replace('_saige_group.txt', '') for f in args.group_files]
    print(f"Reading {len(args.group_files)} cohorts...", file=sys.stderr)
    group_list = [read_saige_group(f) for f in args.group_files]
    merged = merge_groups(group_list, names)
    ng, nv = write_saige_group(merged, args.output_file)
    print(f"Unique genes: {ng}, total variants: {nv}", file=sys.stderr)
    print(f"Wrote {ng} genes to {args.output_file}", file=sys.stderr)

if __name__ == "__main__":
    main()
