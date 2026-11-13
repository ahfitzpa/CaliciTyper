#!/usr/bin/env python3
"""Rebuild calicityper/data/putative_novel_bounds.tsv from the bundled references.

One row per reference cluster (genotype and genogroup, VP1 and RdRp): the mean
and SD of distances between its members and the 2x SD bound (mean + 2 SD) the
putative-novel-type screen compares a query's average distance against. Clusters with
fewer than three references borrow the pooled SD (`sd_source`), shown on the
`*pooled*` rows. VP1 is complete amino acid sequence, RdRp the nucleotide
typing region. Rerun whenever the reference sets change.

    python scripts/build_putative_novel_bounds.py
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from calicityper import putative_novel  # noqa: E402

rows = putative_novel.build_bounds()
cols = ["region", "level", "cluster", "distance", "n_members", "n_pairs", "alphabet",
        "mean", "sd", "sd_source", "bound"]
with open(putative_novel.BOUNDS_FILE, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for r in rows:
        w.writerow({k: (f"{r[k]:.4f}" if isinstance(r[k], float) else r[k]) for k in cols})
print(f"wrote {len(rows)} rows to {putative_novel.BOUNDS_FILE}")
for r in rows:
    if r["cluster"] == "*pooled*":
        print(f"  {r['distance']:<9} {r['region']} {r['level']:<9} pooled: mean {r['mean']:.3f} "
              f"sd {r['sd']:.3f} -> bound {r['bound']:.3f}")
