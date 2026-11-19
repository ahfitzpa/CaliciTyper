#!/usr/bin/env python3
"""Cross-check the Parra 2017 S3/S4 accessions against the bundled reference set.

Reports, per genotype, how many of the paper's accessions are present in the
installed reference set and writes a reviewable TSV. It does NOT assign cluster
labels (A/B/C, C/D): those are only in the paper's figures, not in the tables.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import INPUTS, RESULTS, check_data_dir, read_tsv, write_tsv  # noqa: E402


def bundled():
    out = {}
    data = check_data_dir()
    path = data / "reference_accessions.tsv"
    if not path.exists():
        raise SystemExit(f"reference accessions not found at {path}")
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        # The bundled file is region/genogroup/type/accession; the HPC copy's
        # S3-derived accessions.tsv may be accession/genotype/... . Accept both.
        if r.get("region", "VP1") != "VP1":
            continue
        acc = r.get("accession", "")
        label = r.get("type") or r.get("genotype") or ""
        if acc:
            out[acc.split(".")[0]] = label
    return out


def main():
    refs = bundled()
    s4 = read_tsv(INPUTS / "parra2017_s4_gii6_gii12_gii17.tsv",
                  required=["genotype", "accession"])
    rows = []
    for r in s4:
        acc = r["accession"].split(".")[0]
        in_ref = acc in refs
        rows.append({"paper_genotype": r["genotype"], "accession": r["accession"],
                     "strain": r["strain"], "in_bundled_VP1_set": int(in_ref),
                     "bundled_label": refs.get(acc, ""),
                     "cluster_label": "NOT_IN_TABLE"})
    write_tsv(RESULTS / "parra_s4_bundled_crosscheck.tsv", rows,
              ["paper_genotype", "accession", "strain", "in_bundled_VP1_set",
               "bundled_label", "cluster_label"])
    for g in ("GII.6", "GII.12", "GII.17"):
        sub = [r for r in rows if r["paper_genotype"] == g]
        n = sum(r["in_bundled_VP1_set"] for r in sub)
        print(f"{g}: {n}/{len(sub)} paper accessions in bundled VP1 set")
    print("Cluster labels A/B/C and C/D are in Parra Figs 2A/S2, not in S4 table; left unresolved.")


if __name__ == "__main__":
    main()
