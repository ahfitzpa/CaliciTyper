#!/usr/bin/env python3
"""Rebuild the sapovirus taxonomy CSVs from the published supplementary table.

Source: validation/12864_2024_11048_MOESM1_ESM.xlsx -- Supplemental Table 1 of
Zhao et al. 2024, BMC Genomics 25:1131, "Representative strains and the range of
pairwise distances for all genogroups/P-groups and major genotypes/P types".

The CSVs were first transcribed from an image of that table, and the
transcription put six accessions in the wrong row or changed a character. This
script replaces that with a deterministic read of the spreadsheet, so the
groupings are the published ones by construction.

Four accessions are corrected against NCBI because the PAPER's value does not
exist or is not a sapovirus (see PAPER_TYPOS). One, KC309412 for VP1 genogroup
GXXVI, has no resolution: it is a 366 nt Pseudomonas aeruginosa 16S record. It
is kept in the CSV, flagged, and excluded from any reference set.

    python scripts/build_sapovirus_taxonomy.py
"""
from __future__ import annotations

import csv
import html
import re
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
XLSX = ROOT / "validation" / "12864_2024_11048_MOESM1_ESM.xlsx"
VP1_CSV = ROOT / "validation" / "sapovirus_vp1_taxonomy.csv"
RDRP_CSV = ROOT / "validation" / "sapovirus_rdrp_taxonomy.csv"

# Accession in the table -> (accession to use, why). Each was checked with NCBI
# esummary: the published string is either absent or not a calicivirus, and the
# replacement is the sapovirus record the strain name points to.
PAPER_TYPOS = {
    "L504335": ("LC504335", "L504335 does not resolve; LC504335 is Sapovirus GI Hu/GI.1/JPN/2001/HI-90"),
    "KC215898": ("LC215898", "KC215898 is a scorpion 18S record; LC215898 is Sapovirus GXI/HgTa2-2-2"),
    "KC469052": ("LC469052", "KC469052 is an uncultured bacterium clone; LC469052 is Porcine sapovirus Pig/HgYa2018"),
    "LC215696": ("LC215896", "LC215696 is a whitefly COI record; LC215896 is Sapovirus GX/HgTa2"),
}
# Accessions that cannot be used, with the reason recorded in the CSV.
UNRESOLVED = {
    "KC309412": "does not resolve to a calicivirus (Pseudomonas aeruginosa 16S, 366 nt); "
                "published table appears to be in error",
}


def read_sheet(path):
    """[{column letter: value}] for every row, without an xlsx library."""
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        x = z.read("xl/sharedStrings.xml").decode("utf8")
        shared = [html.unescape(re.sub(r"<[^>]+>", "", m))
                  for m in re.findall(r"<si>(.*?)</si>", x, re.S)]
    sheet = z.read("xl/worksheets/sheet1.xml").decode("utf8")
    rows = []
    for rm in re.finditer(r"<row[^>]*>(.*?)</row>", sheet, re.S):
        cells = {}
        for cm in re.finditer(r'<c r="([A-Z]+)\d+"([^>]*)>(?:<v>([^<]*)</v>|<is>(.*?)</is>)?</c>',
                              rm.group(1), re.S):
            colname, attrs, v, inline = cm.groups()
            if inline is not None:
                val = html.unescape(re.sub(r"<[^>]+>", "", inline))
            elif v is None:
                val = ""
            elif 't="s"' in attrs:
                val = shared[int(v)]
            else:
                val = v
            cells[colname] = val.strip()
        rows.append(cells)
    return rows


def tidy_accession(acc):
    """NC 000940 / NC027026 -> NC_000940 / NC_027026, and drop any version."""
    acc = acc.replace(" ", "_").strip()
    acc = re.sub(r"^(NC)[_ ]?(\d)", r"\1_\2", acc)
    return acc.split(".")[0]


def tidy_range(text):
    """'0.32400000000000001' -> '0.324'; '1.2E-2' -> '0.012'; '0.-0.086' -> '0-0.086'."""
    t = (text or "").strip()
    if not t or t.upper() == "NA":
        return "NA"
    t = t.replace("0.-", "0-")
    parts = re.split(r"\s*-\s*(?=[\d.])", t) if t.count("-") and "E-" not in t.upper() else [t]
    out = []
    for p in parts:
        try:
            out.append(f"{float(p):g}")
        except ValueError:
            out.append(p)
    return "-".join(out)


def parse_block(rows, label_col, acc_col, range_col, start_row=4):
    """The table repeats a label once and then lists its strains on the rows
    below, so a blank label means 'same group as above'."""
    groups, current = {}, None
    for r in rows[start_row:]:
        label = r.get(label_col, "").strip()
        acc = r.get(acc_col, "").strip()
        rng = r.get(range_col, "").strip()
        if label:
            current = label
            groups.setdefault(current, {"accs": [], "range": tidy_range(rng)})
        if current is None:
            continue
        if rng and groups[current]["range"] == "NA":
            groups[current]["range"] = tidy_range(rng)
        if acc:
            groups[current]["accs"].append(acc)
    # genogroup/P-group rows carry a distance range but no strains of their own
    return {k: v for k, v in groups.items() if v["accs"]}


HEADER = """\
# Sapovirus {what} reference strains, from Supplemental Table 1 of
# Zhao et al. 2024, BMC Genomics 25:1131 (dual VP1/RdRp typing proposal).
#
# Generated by scripts/build_sapovirus_taxonomy.py directly from
# validation/{xlsx}, so the groupings are the
# published ones. Do not edit by hand; rerun the script.
#
# `note` is non-empty only where the published accession had to be corrected
# against NCBI, or could not be resolved at all.
"""


def write_csv(path, what, key, groups):
    corrected = unresolved = 0
    with open(path, "w", newline="") as fh:
        fh.write(HEADER.format(what=what, xlsx=XLSX.name))
        w = csv.writer(fh)
        w.writerow([key, "accession", "pairwise_distance_range", "note"])
        for label in groups:
            for acc in groups[label]["accs"]:
                bare = tidy_accession(acc)
                note = ""
                if bare in PAPER_TYPOS:
                    bare, why = PAPER_TYPOS[bare]
                    note = f"corrected from published table: {why}"
                    corrected += 1
                elif bare in UNRESOLVED:
                    note = f"UNRESOLVED: {UNRESOLVED[bare]}"
                    unresolved += 1
                w.writerow([label, bare, groups[label]["range"], note])
    n = sum(len(v["accs"]) for v in groups.values())
    print(f"{path.name}: {len(groups)} labels, {n} accessions "
          f"({corrected} corrected, {unresolved} unresolved)")


def main():
    if not XLSX.exists():
        sys.exit(f"missing source spreadsheet: {XLSX}")
    rows = read_sheet(XLSX)
    write_csv(VP1_CSV, "VP1 genotype", "genotype", parse_block(rows, "A", "B", "C"))
    write_csv(RDRP_CSV, "RdRp P-type", "ptype", parse_block(rows, "E", "F", "G"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
