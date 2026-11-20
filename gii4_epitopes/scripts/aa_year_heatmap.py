#!/usr/bin/env python3
"""Tier A1 from Parra 2017 / Tohma 2017 methods: aa-difference x year-gap heatmap.

For every pair of strains in a genotype, count the amino-acid differences and
the absolute difference in isolation year, then tabulate counts in a matrix.
This is the plot Parra et al. used to show static vs evolving genotypes and the
one Tohma et al. used for GII.2: synonymous nt change accumulates while aa
change does not.

Input: a TSV with columns accession, genotype, year (and an AA FASTA whose ids
match accessions). The parsed S3 Table has genotype/year; the bundled reference
set supplies the VP1 proteins. Output: one TSV per genotype (heatmap counts)
plus a summary of max aa difference within each genotype.

    python dev/gii4_epitopes/scripts/aa_year_heatmap.py \
        inputs/S3.tsv <vp1_aa_fasta> results/heatmaps/
"""
import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS, read_fasta, write_tsv  # noqa: E402


def main(s3, fasta, outdir):
    meta = {}
    with open(s3, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            acc = r["accession"].split(".")[0]
            try:
                year = int(r["year"])
            except (ValueError, TypeError):
                continue
            meta[acc] = (r["genotype"], year)
    seqs = {k.split(".")[0]: v for k, v in read_fasta(fasta).items()}
    by_geno = {}
    for acc, (geno, year) in meta.items():
        if acc in seqs:
            by_geno.setdefault(geno, []).append((acc, year, seqs[acc]))

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    summary = []
    for geno, items in sorted(by_geno.items()):
        counts = Counter()
        max_aa = 0
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                _, y1, s1 = items[i]
                _, y2, s2 = items[j]
                n = min(len(s1), len(s2))
                d = sum(a != b for a, b in zip(s1[:n], s2[:n]))
                max_aa = max(max_aa, d)
                counts[(d, abs(y1 - y2))] += 1
        write_tsv(out / f"heatmap_{geno.replace('.', '_')}.tsv",
                  [{"aa_difference": d, "year_gap": y, "n_pairs": n}
                   for (d, y), n in sorted(counts.items())],
                  ["aa_difference", "year_gap", "n_pairs"])
        summary.append({"genotype": geno, "n_sequences": len(items),
                        "max_aa_difference": max_aa,
                        "static_under_5pct": int(max_aa <= 0.05 * len(items[0][2]))})
    write_tsv(out / "summary.tsv", summary,
              ["genotype", "n_sequences", "max_aa_difference", "static_under_5pct"])
    print(f"wrote {len(summary)} genotype heatmaps + summary to {out}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
