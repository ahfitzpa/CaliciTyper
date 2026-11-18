#!/usr/bin/env python3
"""Convert Zhao et al. (2024) Supplemental Table 1 (.xlsx) into a single
machine-readable CSV of unique accessions with VP1 and/or RdRp (P-type)
designations and pairwise distance ranges.

Usage:
    python3 convert_sapovirus_table.py \
        --input /home/antigravity/12864_2024_11048_MOESM1_ESM.xlsx \
        --output /home/antigravity/sapovirus_vp1_ptype_distances.csv
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import OrderedDict
from pathlib import Path
import openpyxl


KNOWN_TABLE_ANOMALIES = {
    "L504335": "Published typo in RdRp GI.P1 (does not resolve in NCBI; corresponds to LC504335)",
    "KC215898": "Published typo in VP1 GXI (scorpion 18S in NCBI; corresponds to LC215898 in RdRp GXI.P)",
    "KC469052": "Published typo in VP1 GXI (bacterium clone in NCBI; corresponds to LC469052 in RdRp GXI.P)",
    "LC215696": "Published typo in VP1 GX (whitefly COI in NCBI; corresponds to LC215896 in RdRp GX.P)",
    "KC309412": "Published accession in VP1 GXXVI resolves to Pseudomonas aeruginosa 16S in NCBI",
}


def clean_cell(val) -> str:
    if val is None:
        return ""
    return re.sub(r"\s+", " ", str(val)).strip()


def normalize_accession(raw: str) -> str:
    """Normalize whitespace, RefSeq NC_ underscore, and strip .version suffix."""
    s = raw.strip().replace(" ", "_")
    s = re.sub(r"^(NC)[_]?(\d)", r"\1_\2", s)
    return s.split(".")[0]


def tidy_range(raw: str) -> tuple[str, str]:
    """Return (cleaned_range_str, extra_note).
    Handles Excel float precision ('0.32400000000000001'), '0.-0.086', and text notes.
    """
    t = raw.strip()
    if not t:
        return ("", "")
    if t.upper() == "NA":
        return ("NA", "")
    # Detect prose notes placed in the distance column (e.g. row 125)
    if re.search(r"[a-zA-Z]{3,}", t):
        return ("", t)
    note = ""
    if "0.-" in t:
        note = f"Fixed range typo '{t}' -> '{t.replace('0.-', '0-')}'"
        t = t.replace("0.-", "0-")
    parts = (
        re.split(r"\s*-\s*(?=[\d.])", t)
        if ("-" in t and "E-" not in t.upper())
        else [t]
    )
    formatted = []
    for p in parts:
        val = float(p)
        # Format cleanly up to 6 decimal places without trailing zeros, keeping 0 as '0'
        s = f"{val:.6f}".rstrip("0").rstrip(".")
        formatted.append(s if s else "0")
    return ("-".join(formatted), note)


def parse_min_max(rng: str) -> tuple[str, str]:
    """Extract numeric (min, max) strings for float columns; ('', '') if NA or empty."""
    if not rng or rng == "NA" or ";" in rng:
        return ("", "")
    if "-" in rng:
        lo, hi = rng.split("-", 1)
        return (lo, hi)
    return (rng, rng)


def parse_side(
    ws,
    label_col: int,
    acc_col: int,
    range_col: int,
    is_ptype: bool,
    start_row: int = 4,
    end_row: int = 164,
) -> list[dict]:
    """Parse either the VP1 (cols 1-3) or RdRp/P-type (cols 5-7) block."""
    group_ranges: dict[str, str] = {}
    cur_group = ""
    cur_type = ""
    cur_type_range = ""
    records: list[dict] = []

    for r in range(start_row, end_row + 1):
        lbl = clean_cell(ws.cell(r, label_col).value)
        acc = clean_cell(ws.cell(r, acc_col).value)
        rng_raw = clean_cell(ws.cell(r, range_col).value)
        rng, rng_note = tidy_range(rng_raw)

        if lbl and not acc:
            # Summary header row for a genogroup / P-group (e.g. GI, GII, GI.P, GV.P)
            cur_group = lbl
            group_ranges[cur_group] = rng
            cur_type = ""
            cur_type_range = ""
            continue

        if lbl and acc:
            cur_type = lbl
            cur_type_range = rng if rng else "NA"
            if is_ptype:
                m = re.match(r"^([A-Z0-9]+\.P)", lbl)
                cur_group = m.group(1) if m else lbl
            else:
                cur_group = lbl.split(".")[0]

        if acc:
            notes = []
            if rng_note:
                notes.append(rng_note)
            if not lbl and rng and rng != cur_type_range:
                notes.append(f"Extra distance cell value: {rng}")
            norm = normalize_accession(acc)
            if norm != acc:
                notes.append(f"Normalized accession from '{acc}'")
            records.append(
                {
                    "row": r,
                    "raw_acc": acc,
                    "accession": norm,
                    "group": cur_group,
                    "group_range": group_ranges.get(cur_group, ""),
                    "type": cur_type,
                    "type_range": cur_type_range,
                    "notes": notes,
                }
            )

    return records


def collapse_entries(entries: list[dict], region_label: str) -> dict:
    """Collapse 1+ rows for the same accession within a single region."""
    raw_accs = list(dict.fromkeys(e["raw_acc"] for e in entries))
    groups = list(dict.fromkeys(e["group"] for e in entries))
    group_ranges = [e["group_range"] for e in entries]
    # Keep distinct non-empty group ranges aligned if unique
    uniq_group_ranges = list(dict.fromkeys(gr for gr in group_ranges if gr))

    types = [e["type"] for e in entries]
    type_ranges = [e["type_range"] for e in entries]
    notes: list[str] = []
    for e in entries:
        notes.extend(e["notes"])

    if len(entries) > 1:
        dup_desc = ", ".join(
            f"{e['type']} ({e['type_range']}, row {e['row']})" for e in entries
        )
        notes.append(f"Duplicate {region_label} listing in source table: {dup_desc}")
        rng_joined = ";".join(type_ranges)
        d_min, d_max = ("", "")
    else:
        rng_joined = type_ranges[0]
        d_min, d_max = parse_min_max(rng_joined)

    return {
        "raw_acc": ";".join(raw_accs),
        "group": ";".join(groups),
        "group_range": ";".join(uniq_group_ranges),
        "type": ";".join(types),
        "range": rng_joined,
        "min": d_min,
        "max": d_max,
        "notes": notes,
    }


def build_unified_table(xlsx_path: Path, out_csv: Path) -> None:
    wb = openpyxl.load_workbook(xlsx_path, data_only=True)
    ws = wb.active

    vp1_records = parse_side(ws, 1, 2, 3, is_ptype=False)
    ptype_records = parse_side(ws, 5, 6, 7, is_ptype=True)

    # Preserve first-seen accession order (VP1 order, then RdRp-only in sheet order)
    by_acc: OrderedDict[str, dict[str, list[dict]]] = OrderedDict()
    for rec in vp1_records:
        by_acc.setdefault(rec["accession"], {"vp1": [], "ptype": []})["vp1"].append(rec)
    for rec in ptype_records:
        by_acc.setdefault(rec["accession"], {"vp1": [], "ptype": []})["ptype"].append(rec)

    fieldnames = [
        "accession",
        "vp1_raw_accession",
        "vp1_genogroup",
        "vp1_genogroup_distance_range",
        "vp1_genotype",
        "vp1_distance_range",
        "vp1_distance_min",
        "vp1_distance_max",
        "ptype_raw_accession",
        "p_group",
        "p_group_distance_range",
        "p_type",
        "ptype_distance_range",
        "ptype_distance_min",
        "ptype_distance_max",
        "regions_present",
        "notes",
    ]

    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for acc, buckets in by_acc.items():
            vp1 = collapse_entries(buckets["vp1"], "VP1") if buckets["vp1"] else None
            pt = collapse_entries(buckets["ptype"], "P-type") if buckets["ptype"] else None

            if vp1 and pt:
                presence = "both"
            elif vp1:
                presence = "vp1_only"
            else:
                presence = "ptype_only"

            all_notes: list[str] = []
            if acc in KNOWN_TABLE_ANOMALIES:
                all_notes.append(KNOWN_TABLE_ANOMALIES[acc])
            if vp1:
                all_notes.extend(vp1["notes"])
            if pt:
                all_notes.extend(pt["notes"])

            writer.writerow(
                {
                    "accession": acc,
                    "vp1_raw_accession": vp1["raw_acc"] if vp1 else "",
                    "vp1_genogroup": vp1["group"] if vp1 else "",
                    "vp1_genogroup_distance_range": vp1["group_range"] if vp1 else "",
                    "vp1_genotype": vp1["type"] if vp1 else "",
                    "vp1_distance_range": vp1["range"] if vp1 else "",
                    "vp1_distance_min": vp1["min"] if vp1 else "",
                    "vp1_distance_max": vp1["max"] if vp1 else "",
                    "ptype_raw_accession": pt["raw_acc"] if pt else "",
                    "p_group": pt["group"] if pt else "",
                    "p_group_distance_range": pt["group_range"] if pt else "",
                    "p_type": pt["type"] if pt else "",
                    "ptype_distance_range": pt["range"] if pt else "",
                    "ptype_distance_min": pt["min"] if pt else "",
                    "ptype_distance_max": pt["max"] if pt else "",
                    "regions_present": presence,
                    "notes": " | ".join(dict.fromkeys(all_notes)),
                }
            )

    print(
        f"Wrote {len(by_acc)} unique accessions "
        f"(from {len(vp1_records)} VP1 entries and {len(ptype_records)} P-type entries) "
        f"to {out_csv}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("/home/antigravity/12864_2024_11048_MOESM1_ESM.xlsx"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("/home/antigravity/sapovirus_vp1_ptype_distances.csv"),
    )
    args = parser.parse_args()
    build_unified_table(args.input, args.output)


if __name__ == "__main__":
    main()
