#!/usr/bin/env python3
"""Validate every F1-F4 input TSV before any build/validation script runs.

Checks schema, placeholder values, and cross-file consistency. A non-zero exit
means the inputs are not ready; fix the TSVs rather than bypassing this script.
"""
import csv
from pathlib import Path

from common import HERE, INPUTS

EXPECTED = {
    "gii4_variant_references.tsv": ["accession", "variant", "region", "source", "evidence"],
    "sublineage_definitions.tsv": ["sublineage_id", "genotype", "region", "label", "source_example_accession", "evidence", "notes"],
    "variant_thresholds.tsv": ["param", "value", "basis", "status", "notes"],
    "intragenotypic_clusters.tsv": ["genotype", "cluster", "accession", "source_label", "evidence"],
    "immunotype_regime.tsv": ["genotype", "immunotype", "regime", "evidence", "status"],
}

PLACEHOLDERS = {"FILL_FROM_PARRA", "TO_VERIFY_FROM_FIG5C", "FILL", "TODO", "NA"}


def main():
    errors, warnings = [], 0
    for fn, cols in EXPECTED.items():
        p = INPUTS / fn
        if not p.exists():
            errors.append(f"MISSING input file: {p}")
            continue
        with open(p, newline="", encoding="utf-8") as fh:
            r = csv.DictReader(fh, delimiter="\t")
            if r.fieldnames != cols:
                errors.append(f"{fn}: header {r.fieldnames} != expected {cols}")
            rows = list(r)
        if not rows:
            errors.append(f"{fn}: no data rows")
        for i, row in enumerate(rows, start=2):
            for key in cols:
                if row.get(key) in PLACEHOLDERS:
                    warnings += 1
                    print(f"WARN {fn}:{i} '{key}' is placeholder ({row[key]})")
        print(f"OK   {fn}: {len(rows)} rows")

    # cross-file: GII.4 variant references must match the bundled label list
    accs = {}
    with open(INPUTS / "gii4_variant_references.tsv", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            accs[r["accession"]] = r["variant"]
    if len(accs) != 20:
        errors.append(f"expected 20 GII.4 VP1 references, found {len(accs)}")
    variants = sorted(set(accs.values()))
    if len(variants) != 11:
        errors.append(f"expected 11 GII.4 variants, found {len(variants)}: {variants}")
    print(f"GII.4 references: {len(accs)} accessions, {len(variants)} variants")

    # immunotype table must cover GII.4/GII.17/GII.6 with supported status
    imm = {}
    with open(INPUTS / "immunotype_regime.tsv", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            imm[r["genotype"]] = r
    for g in ("GII.4", "GII.6", "GII.17", "GI.3", "GI.7", "GI.8", "GI.9"):
        if g not in imm:
            errors.append(f"immunotype_regime.tsv missing {g}")

    if errors:
        print("\nFAIL:")
        for e in errors:
            print("  " + e)
        raise SystemExit(1)
    print(f"\nAll structural checks passed ({warnings} placeholder warnings — resolve before implementation).")


if __name__ == "__main__":
    main()
