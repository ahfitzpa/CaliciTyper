#!/usr/bin/env python3
"""Rebuild calicityper/data/reference_accessions.tsv from the files on disk.

The accession list ships with the reference deposit so a reader can see which
sequences the calls are made against without unpacking the archive. It is
documentation: no code reads it, which is exactly why it had drifted.

Two defects it had, both fixed here by generating it rather than maintaining
it. It listed 44 duplicate rows, because a strain in a per-genogroup .orf is
also in all.orf and both were walked. And it labelled the region `VP1` or
`RdRp` for norovirus and sapovirus alike, so `VP1 / GI` named two unrelated
reference sets; a `virus` column now separates them.

    python scripts/build_reference_accessions.py [--apply]
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from calicityper import config as cg          # noqa: E402

OUT = cg.DATA / "reference_accessions.tsv"
COLUMNS = ["virus", "region", "genogroup", "type", "accession", "nt"]


def from_refsets():
    """The typing reference sets: one row per (virus, region, accession).

    all.orf is skipped: it is the same records pooled, and walking it as well
    is what produced the duplicate rows.
    """
    rows = []
    sets = [("norovirus", "VP1", cg.REF_DIR["orf2"]),
            ("norovirus", "RdRp", cg.REF_DIR["orf1"]),
            ("sapovirus", "VP1", cg.SAPOVIRUS_REF_DIR["orf2"]),
            ("sapovirus", "RdRp", cg.SAPOVIRUS_REF_DIR["orf1"])]
    for virus, region, d in sets:
        d = Path(d)
        if not d.exists():
            continue
        for orf in sorted(d.glob("*.orf")):
            if orf.stem == "all":
                continue
            for r in SeqIO.parse(orf, "fasta"):
                parts = r.id.split("|")
                group, label, acc = (parts if len(parts) == 3
                                     else (orf.stem, parts[0], parts[-1]))
                rows.append(dict(virus=virus, region=region, genogroup=group,
                                 type=label, accession=acc,
                                 nt=len(str(r.seq).replace("-", ""))))
    return rows


def from_blastdb():
    """The genogroup assignment database, which is neither VP1 nor RdRp."""
    fa = cg.DATA / "database" / "refseq_2.fasta"
    if not fa.exists():
        return []
    rows = []
    for r in SeqIO.parse(fa, "fasta"):
        parts = r.description.split("|")
        group = parts[0].split(maxsplit=1)[-1] if len(parts) == 1 else parts[0]
        acc = r.id.split("|")[-1]
        rows.append(dict(virus="caliciviridae", region="genogroup_blastdb",
                         genogroup=group, type="", accession=acc,
                         nt=len(str(r.seq).replace("-", ""))))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--apply", action="store_true", help="write the file")
    a = ap.parse_args()

    rows = from_refsets() + from_blastdb()
    seen, unique = set(), []
    for r in rows:
        k = tuple(r[c] for c in COLUMNS)
        if k in seen:
            continue
        seen.add(k)
        unique.append(r)
    unique.sort(key=lambda r: (r["virus"], r["region"], r["genogroup"],
                               r["type"], r["accession"]))

    dropped = len(rows) - len(unique)
    old = sum(1 for _ in open(OUT)) - 1 if OUT.exists() else 0
    print(f"{len(unique)} rows ({old} before); {dropped} duplicate row(s) dropped")
    for virus in sorted({r["virus"] for r in unique}):
        n = sum(1 for r in unique if r["virus"] == virus)
        regions = sorted({r["region"] for r in unique if r["virus"] == virus})
        print(f"  {virus:14s} {n:4d} rows over {', '.join(regions)}")

    if not a.apply:
        print("\nreport only; pass --apply to write")
        return 0
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(unique)
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
