#!/usr/bin/env python3
"""Cross-check our sapovirus tables against the published tables they derive from.

Our VP1 and RdRp taxonomy tables are generated from Supplemental Table 1 of
Zhao et al. 2024. The same paper publishes two further tables that assign the
same strains independently:

  Table 2  representatives per VP1 genotype and P-type, with the dual-typing
           designation            -> validation/zhao_table2_dual_types.csv
  pig/bat  current and new VP1 genogroup, and the RdRp P-GROUP, for 21 pig and
           bat strains            -> validation/zhao_genogroup_reassignment.csv

Where a strain appears in more than one of these the assignments must agree.
Where they do not, either the published tables disagree with each other or our
parse of the supplement is wrong, and either way it needs resolving rather
than being carried forward silently.

These tables are also the reference standard for sapovirus P-type, which no
other tool provides: the Norovirus Typing Tool does not assign sapovirus
P-types at all.

    python scripts/check_against_zhao_tables.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
V = ROOT / "validation"

# Deviations we have adjudicated: the two published tables disagree and the
# supplement wins, so our table is deliberately not what the other table says.
# These are errors in the paper that we cannot fix, only record, so they are
# reported as accepted rather than failing the check. Anything NOT listed here
# is a real problem.
ACCEPTED = {
    ("MH003837", "VP1 genogroup"): (
        "the pig/bat table says GXV, the supplement says GNA4. GNA4 is a genotype name "
        "coined in this paper and MH003837 is its only member, so were it GXV the name "
        "would never have been coined and the paper's own count of 53 VP1 genotypes "
        "would be 52. We follow the supplement."),
    ("MH003837", "RdRp P-group"): (
        "follows the VP1 decision above: GNA4.P is this strain's own P-group."),
}


def read(path, skip_retired=True):
    rows = []
    for r in csv.DictReader(l for l in open(path) if not l.startswith("#")):
        note = (r.get("note") or "")
        if skip_retired and note.startswith(("UNRESOLVED", "RETIRED")):
            continue
        rows.append(r)
    return rows


def ours():
    vp1, rdrp, retired = {}, {}, {}
    for r in read(V / "sapovirus_vp1_taxonomy.csv"):
        vp1[r["accession"].split(".")[0]] = r["genotype"]
    for r in read(V / "sapovirus_vp1_taxonomy.csv", skip_retired=False):
        if (r.get("note") or "").startswith("RETIRED"):
            retired.setdefault(r["accession"].split(".")[0], []).append(r["genotype"])
    for r in read(V / "sapovirus_rdrp_taxonomy.csv"):
        rdrp[r["accession"].split(".")[0]] = r["ptype"]
    return vp1, rdrp, retired


def main():
    vp1, rdrp, retired = ours()
    print(f"our tables: {len(vp1)} VP1 and {len(rdrp)} RdRp assignments; "
          f"{len(retired)} strains also carry a retired VP1 label\n")
    problems, confirmed, missing, accepted = [], [], [], []

    t2 = V / "zhao_table2_dual_types.csv"
    if t2.exists():
        print("=" * 78)
        print("Table 2 -- representatives: VP1 genotype and P-type")
        print("=" * 78)
        for r in read(t2):
            a = r["accession"].split(".")[0]
            for field, mine, label in (("genotype", vp1.get(a), "VP1 genotype"),
                                       ("ptype", rdrp.get(a), "P-type")):
                want = (r.get(field) or "").strip()
                if not want:
                    continue
                if mine is None:
                    missing.append((a, label, want))
                elif mine != want:
                    if (a, label) in ACCEPTED:
                        accepted.append((a, label, want, mine))
                        continue
                    print(f"  {a:11s} {label:13s} published {want:12s} ours {mine:12s}  MISMATCH")
                    problems.append((a, label, want, mine))
                elif a in retired:
                    confirmed.append((a, label, want, retired[a]))

    gg = V / "zhao_genogroup_reassignment.csv"
    if gg.exists():
        print("\n" + "=" * 78)
        print("Pig/bat table -- VP1 genogroup and RdRp P-group")
        print("=" * 78)
        for r in read(gg):
            a = r["accession"].split(".")[0]
            want_g, want_pg = r["vp1_genogroup"].strip(), r["rdrp_group"].strip()
            mine_gt, mine_pt = vp1.get(a), rdrp.get(a)
            # our genotype carries its genogroup as a prefix; this table gives the
            # genogroup alone, so compare at genogroup level
            mine_g = mine_gt.split(".")[0] if mine_gt else None
            mine_pg = (mine_pt.split(".P")[0] + ".P") if mine_pt else None
            if mine_g is None:
                missing.append((a, "VP1 genogroup", want_g))
            elif mine_g != want_g:
                if (a, "VP1 genogroup") in ACCEPTED:
                    accepted.append((a, "VP1 genogroup", want_g, mine_g))
                else:
                    print(f"  {a:11s} VP1 genogroup published {want_g:9s} ours {mine_g:9s}  MISMATCH")
                    problems.append((a, "VP1 genogroup", want_g, mine_g))
            elif a in retired:
                confirmed.append((a, "VP1 genogroup", want_g, retired[a]))
            if mine_pg is None:
                missing.append((a, "RdRp P-group", want_pg))
            elif mine_pg != want_pg:
                if (a, "RdRp P-group") in ACCEPTED:
                    accepted.append((a, "RdRp P-group", want_pg, mine_pg))
                else:
                    print(f"  {a:11s} RdRp P-group  published {want_pg:9s} ours {mine_pg:9s}  MISMATCH")
                    problems.append((a, "RdRp P-group", want_pg, mine_pg))

    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    if confirmed:
        seen = set()
        print(f"\nINDEPENDENTLY CONFIRMED -- strains the supplement lists twice, which we")
        print("resolved on distance evidence alone, and a second published table agrees:")
        for a, label, want, old in confirmed:
            if (a, label) in seen:
                continue
            seen.add((a, label))
            print(f"  {a:11s} {label:13s} -> {want:10s}  (retired: {', '.join(sorted(set(old)))})")
    if accepted:
        print(f"\nACCEPTED DEVIATIONS ({len(accepted)}) -- the paper's two tables disagree and we")
        print("follow the supplement; these are errors in the publication, recorded not fixed:")
        for a, label, want, mine in accepted:
            print(f"  {a:11s} {label:13s} that table {want:10s} ours {mine:10s}")
            print(f"              {ACCEPTED[(a, label)]}")
    if missing:
        print(f"\nNOT HELD in our reference sets ({len(missing)}):")
        for a, label, want in missing[:12]:
            print(f"  {a:11s} {label:13s} published {want}")
        if len(missing) > 12:
            print(f"  ... and {len(missing) - 12} more")
    if problems:
        print(f"\nMISMATCHES ({len(problems)}) -- our tables disagree with a published table:")
        for a, label, want, mine in problems:
            print(f"  {a:11s} {label:13s} published {want:12s} ours {mine}")
        print("\nEach must be resolved against Supplemental Table 1 before release.")
    else:
        print("\nno mismatches")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
