#!/usr/bin/env python3
"""F2: assemble the reviewable sub-lineage input table from supplied definitions.

Combines inputs/sublineage_definitions.tsv with the GII.4 variant references
and emits a single candidates table with a 'ready' flag. Rows whose source
accession is a FILL placeholder are marked not ready. This is an INPUT
preparation step; the definition of a divergent sub-lineage and its threshold
must still be validated before any flag is implemented.
"""
from common import INPUTS, RESULTS, read_tsv, write_tsv

BAD = {"FILL_FROM_PARRA", "FILL", "TODO", "NA"}


def main():
    defs = read_tsv(INPUTS / "sublineage_definitions.tsv",
                    required=["sublineage_id", "genotype", "region", "label", "evidence"])
    refs = read_tsv(INPUTS / "gii4_variant_references.tsv",
                    required=["accession", "variant"])
    rows = []
    for r in defs:
        acc = r["source_example_accession"]
        rows.append({
            "sublineage_id": r["sublineage_id"], "genotype": r["genotype"],
            "region": r["region"], "label": r["label"],
            "source_example_accession": acc,
            "ready": int(acc not in BAD),
            "evidence": r["evidence"], "notes": r.get("notes", ""),
        })
    for r in refs:
        rows.append({
            "sublineage_id": r["variant"], "genotype": "GII.4", "region": "VP1",
            "label": "named GII.4 variant", "source_example_accession": r["accession"],
            "ready": 1, "evidence": r["evidence"], "notes": "bundled variant label; not an immunotype",
        })
    write_tsv(RESULTS / "f2_sublineage_candidates.tsv", rows,
              ["sublineage_id", "genotype", "region", "label",
               "source_example_accession", "ready", "evidence", "notes"])
    print(f"{len(rows)} candidate sub-lineages, {sum(r['ready'] for r in rows)} with a real example accession.")
    print("Threshold rule UNVALIDATED — do not add a divergent-sublineage flag to the tool yet.")


if __name__ == "__main__":
    main()
