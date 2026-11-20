#!/usr/bin/env python3
"""Extract the two Parra 2017 DOCX supplementary tables into TSV inputs.

    python dev/gii4_epitopes/scripts/extract_parra_tables.py <S3.docx> <S4.docx>

S3 (ppat.1006136.s009.docx): ORF2 sequences used to build the immunotype tree:
    idx, accession, host, genotype, strain, year, country
S4 (ppat.1006136.s010.docx): GII.6/GII.12/GII.17 sequences:
    genotype, idx, accession, strain

The DOCX files are not committed (publisher supplementary material); the
extracted TSVs are. Licence: PLoS Pathog 2017 is CC0 1.0.
"""
import csv
import sys
from pathlib import Path

from docx import Document

HERE = Path(__file__).resolve().parents[1]
INPUTS = HERE / "inputs"


def as_rows(table):
    out = []
    for r in table.rows[1:]:
        cells = [c.text.strip() for c in r.cells]
        if len(cells) >= 2 and cells[1]:
            out.append(cells)
    return out


def main(s3, s4):
    d = Document(s3)
    rows3 = [r for r in as_rows(d.tables[0]) if len(r) >= 7]
    p3 = INPUTS / "parra2017_s3_immunotype_sequences.tsv"
    with open(p3, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["idx", "accession", "host", "genotype", "strain", "year", "country"])
        w.writerows([r[:7] for r in rows3])
    print(f"wrote {p3} ({len(rows3)} rows, {len(set(r[3] for r in rows3))} genotypes)")

    d = Document(s4)
    out4 = []
    for label, table in zip(("GII.6", "GII.12", "GII.17"), d.tables):
        for r in as_rows(table):
            out4.append([label, r[0], r[1], r[2] if len(r) > 2 else ""])
    p4 = INPUTS / "parra2017_s4_gii6_gii12_gii17.tsv"
    with open(p4, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["genotype", "idx", "accession", "strain"])
        w.writerows(out4)
    print(f"wrote {p4} ({len(out4)} rows)")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2])
