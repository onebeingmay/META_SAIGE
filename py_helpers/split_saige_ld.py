#!/usr/bin/env python
"""Split a combined SAIGE step3_LDmat.R LD matrix into per-gene files for META_SAIGE.

Background
----------
SAIGE `step3_LDmat.R` (as shipped in wzhou88/saige:1.5.0) writes ONE combined
sparse LD matrix per chromosome, NOT one file per gene:

    <prefix>.LDmat.txt       concatenated per-gene sparse triplet blocks "i j x".
                             Indices are LOCAL 0-based *within each gene block*
                             (they reset to 0 at the start of every gene).
    <prefix>.index.txt       one line per gene: "<start> <end> <gene>" giving the
                             0-based, inclusive LINE range of that gene's block
                             inside <prefix>.LDmat.txt.
    <prefix>.marker_info.txt per-gene marker info (Index column is already local).

META_SAIGE's loader (load_cohort in R/MetaSAIGE.R) instead expects ONE LD file
per gene, found at  paste0(gene_file_prefix, gene, ".txt")  -- e.g. the bundled
example uses cohort1_chr_7_GCK.txt, cohort1_chr_7_AASS.txt, ...
(README: "Step3 generates a sparse LD matrix for each gene").

This script bridges that format gap. For each gene it copies the gene's block of
lines out of the combined LDmat verbatim -- because the block indices are already
local 0-based, no re-indexing is needed -- into "<out_prefix><gene>.txt".

Usage
-----
    python split_saige_ld.py \
        --ldmat  /path/hmf_AS_step3.22.LDmat.txt \
        --index  /path/hmf_AS_step3.22.index.txt \
        --out-prefix /path/ld_split/hmf_AS_step3.22_

Then point META_SAIGE's --gene_file_prefix at the same out-prefix:
    --gene_file_prefix /path/ld_split/hmf_AS_step3.22_
"""
import argparse
import os
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ldmat", required=True,
                    help="combined <prefix>.LDmat.txt from SAIGE step3_LDmat.R")
    ap.add_argument("--index", required=True,
                    help="<prefix>.index.txt: 'start end gene' line ranges into --ldmat")
    ap.add_argument("--out-prefix", required=True, dest="out_prefix",
                    help="output prefix; writes <out_prefix><gene>.txt per gene")
    args = ap.parse_args()

    out_dir = os.path.dirname(args.out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # The LDmat for one chromosome is small (a few MB at most); read it once.
    with open(args.ldmat) as fh:
        lines = fh.readlines()
    n_lines = len(lines)

    n_written = 0
    n_empty = 0
    with open(args.index) as idx:
        for lineno, raw in enumerate(idx, start=1):
            raw = raw.strip()
            if not raw:
                continue
            parts = raw.split()
            if len(parts) < 3:
                sys.exit("ERROR: %s line %d is malformed (expected 'start end gene'): %r"
                         % (args.index, lineno, raw))
            try:
                start, end = int(parts[0]), int(parts[1])
            except ValueError:
                sys.exit("ERROR: %s line %d has non-integer range: %r"
                         % (args.index, lineno, raw))
            # Gene symbols never contain whitespace; join defensively just in case.
            gene = "_".join(parts[2:])

            if end < start:
                # Empty block: SAIGE step3 writes "<start> <start-1> <gene>" (zero
                # lines) for a gene with no qualifying variants in this cohort. That
                # gene is also absent from marker_info, so META_SAIGE skips it; we
                # still emit an empty file so the per-gene path exists explicitly.
                block = []
                n_empty += 1
            elif start < 0 or end >= n_lines:
                sys.exit("ERROR: gene %s range [%d, %d] is out of bounds for %s "
                         "(%d lines). LDmat and index.txt are mismatched."
                         % (gene, start, end, args.ldmat, n_lines))
            else:
                block = lines[start:end + 1]

            with open("%s%s.txt" % (args.out_prefix, gene), "w") as out:
                out.writelines(block)
            n_written += 1

    msg = "[split_saige_ld] wrote %d per-gene LD files with prefix '%s'" % (
        n_written, args.out_prefix)
    if n_empty:
        msg += " (%d empty)" % n_empty
    print(msg)


if __name__ == "__main__":
    main()
