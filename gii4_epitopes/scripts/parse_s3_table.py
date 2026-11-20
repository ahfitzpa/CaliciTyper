#!/usr/bin/env python3
"""Extract the accession list from Parra et al. 2017 S3 Table.

S3 Table is "Norovirus ORF2 sequences used in the construction of trees
representing the clustering of immunotypes" -- that is, the exact input to
Figure 5, whose per-genotype immunotype assignments are otherwise only
readable off the image. Rebuilding the tree from these accessions is the way
to recover the assignments with their own method rather than approximating
it with ours.

    python scripts/parse_s3_table.py            # writes inputs/parra_s3/accessions.tsv
"""
from __future__ import annotations

import csv
import html
import re
import sys
import zipfile

from common import INPUTS

SRC = INPUTS / "parra_s3" / "S3_Table.docx"
OUT = INPUTS / "parra_s3" / "accessions.tsv"
COLS = ["accession", "host", "genotype", "strain", "year", "country"]


def rows_from_docx(path):
    xml = zipfile.ZipFile(path).read("word/document.xml").decode("utf8")
    for tr in re.findall(r"<w:tr[ >].*?</w:tr>", xml, re.S):
        cells = [html.unescape("".join(re.findall(r"<w:t[^>]*>([^<]*)</w:t>", tc))).strip()
                 for tc in re.findall(r"<w:tc[ >].*?</w:tc>", tr, re.S)]
        if cells:
            yield cells


def main():
    if not SRC.exists():
        sys.exit(f"missing {SRC}\n  download it from "
                 "https://journals.plos.org/plospathogens/article/file"
                 "?type=supplementary&id=10.1371/journal.ppat.1006136.s009")
    acc_re = re.compile(r"^(?:[A-Z]{1,2}\d{5,6}|NC_\d{6})(?:\.\d+)?$")
    out, seen = [], set()
    for cells in rows_from_docx(SRC):
        # the body rows carry a leading row number, the header does not
        body = cells[1:] if cells and cells[0].isdigit() else cells
        if len(body) < 3 or not acc_re.match(body[0]):
            continue
        acc = body[0].split(".")[0]
        if acc in seen:
            continue
        seen.add(acc)
        out.append(dict(zip(COLS, (body + [""] * 6)[:6])))

    if not out:
        sys.exit("no accession rows parsed; the table layout may have changed")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        w.writerows(out)

    from collections import Counter
    by_gt = Counter(r["genotype"] for r in out)
    print(f"{len(out)} accessions over {len(by_gt)} genotypes -> {OUT}")
    print("  " + ", ".join(f"{g}:{n}" for g, n in sorted(by_gt.items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
