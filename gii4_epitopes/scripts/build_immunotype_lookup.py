#!/usr/bin/env python3
"""F4: build a versioned immunotype/regime lookup JSON from inputs/immunotype_regime.tsv.

The lookup is a STATIC ANNOTATION of genotypes, never a query-inferred property.
Only rows whose status is SUPPORTED or SUPPORTED_WITH_CAVEAT are emitted by
default; TO_VERIFY rows are listed as unverified and must be curated from
Parra 2017 Fig 5C before the lookup is complete. The output carries provenance
(source DOI, manifest hash) and a schema version, and is intended to be added
to the reference set and manifest — not to code.
"""
import json
from datetime import datetime, timezone

from common import INPUTS, RESULTS, read_tsv, provenance

ACCEPTED = {"SUPPORTED", "SUPPORTED_WITH_CAVEAT"}


def main():
    rows = read_tsv(INPUTS / "immunotype_regime.tsv",
                    required=["genotype", "immunotype", "regime", "evidence", "status"])
    verified, unverified = [], []
    for r in rows:
        (verified if r["status"] in ACCEPTED else unverified).append(r)

    doc = {
        "schema": {"name": "calicityper.immunotype_regime", "version": "0.1.0-draft"},
        "source": {
            "primary": "Parra et al. 2017, PLoS Pathog 13(1):e1006136",
            "doi": "10.1371/journal.ppat.1006136",
            "note": "Immunotypes A-L grouped at ~>=20% VP1 aa difference between groups. "
                    "GII.17 is 'Static' in Table 1 with a candidate emerging lineage; "
                    "it is annotated static_emerging_candidate, not evolving.",
        },
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "reference_set": provenance(),
        "status": "DRAFT — not approved for inclusion in the tool",
        "verified": {r["genotype"]: {"immunotype": r["immunotype"],
                                     "regime": r["regime"], "evidence": r["evidence"]}
                     for r in verified},
        "unverified": [r["genotype"] for r in unverified],
    }
    out = RESULTS / "immunotype_regime.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=2), encoding="utf-8")
    print(f"wrote {out}: {len(verified)} verified, {len(unverified)} unverified genotypes")
    if unverified:
        print("Unverified (curate from Parra Fig 5C before use): " + ", ".join(r["genotype"] for r in unverified))


if __name__ == "__main__":
    main()
