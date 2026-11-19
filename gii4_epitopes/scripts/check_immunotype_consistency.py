#!/usr/bin/env python3
"""Sanity-check the immunotype table against the rule Parra et al. state.

Parra et al. 2017 group genotypes into immunotypes A-L at "an approximate
cut-off value of >=20% aa differences between groups". The per-genotype
assignments appear only in their Figure 5, which is an image; the text names
just a few (immunotype G is GII.4 and GII.20; GII.17 is J; GII.6 is used as
an example of H). The rest of inputs/immunotype_regime.tsv was drafted
without a source and is marked TO_VERIFY_FROM_FIG5C.

This does NOT reconstruct the assignments. It measures mean VP1 amino-acid
p-distance between genotypes in the bundled reference set and reports where
the drafted table places two genotypes together that are further apart than
the stated cut-off.

Read the output carefully, because it contains its own disclaimer: this
measure does not reproduce the paper's stated anchors either. GII.4 and
GII.20 are one immunotype by the text, and come out 26.8% apart here. So the
quantity measured here is not the quantity Parra measured -- different VP1
span, different reference strains, or a different distance -- and a
disagreement is therefore not evidence that the table is wrong.

One conclusion does survive that caveat: a pair from different genogroups,
tens of percent apart, cannot share an immunotype under any reading.

    python scripts/check_immunotype_consistency.py
"""
from __future__ import annotations

import itertools
import sys

import numpy as np

from common import HERE, ROOT, INPUTS, read_tsv     # noqa: F401
sys.path.insert(0, str(ROOT))
from calicityper import putative_novel as nv         # noqa: E402

CUTOFF = 0.20
ANCHORS = [("GII.4", "GII.20", "G"), ("GI.3", "GI.7", "C"),
           ("GI.3", "GI.8", "C"), ("GI.3", "GI.9", "C")]


def genotype_sequences():
    recs, _ = nv.vp1_aa_alignment()
    seqs = dict(recs)
    out = {}
    for acc, _gg, ty, _s in nv.vp1_references():
        if acc in seqs:
            out.setdefault(str(ty).split("_")[0], []).append(seqs[acc])
    return out


def pdist(a, b, min_sites=100):
    ok = [(x, y) for x, y in zip(a, b) if x not in "-X" and y not in "-X"]
    return (sum(x != y for x, y in ok) / len(ok)) if len(ok) >= min_sites else float("nan")


def mean_between(gt, a, b):
    v = [pdist(x, y) for x in gt.get(a, []) for y in gt.get(b, [])]
    v = [x for x in v if x == x]
    return float(np.mean(v)) if v else None


def main():
    gt = genotype_sequences()
    rows = read_tsv(INPUTS / "immunotype_regime.tsv")
    print(f"{len(gt)} genotypes with aligned amino acids; "
          f"{len(rows)} rows in the table\n")

    print("published anchors, as a control on the measure:")
    reproduced = 0
    for a, b, grp in ANCHORS:
        d = mean_between(gt, a, b)
        if d is None:
            print(f"  {a:7s} {b:7s} immunotype {grp}: not comparable")
            continue
        ok = d < CUTOFF
        reproduced += ok
        print(f"  {a:7s} {b:7s} immunotype {grp}: {d * 100:5.1f}%"
              f"  {'consistent' if ok else 'NOT reproduced by this measure'}")
    if reproduced < len(ANCHORS):
        print(f"\n  {len(ANCHORS) - reproduced} of {len(ANCHORS)} published anchors are not")
        print("  reproduced, so this measure is not the one the paper used. Treat the")
        print("  table below as a sanity check for impossibilities only, not as")
        print("  evidence for or against any particular assignment.\n")

    by_group = {}
    for r in rows:
        if r.get("immunotype"):
            by_group.setdefault(r["immunotype"], []).append(r["genotype"])

    print("pairs the table groups together, by distance:")
    impossible = []
    for grp, members in sorted(by_group.items()):
        for a, b in itertools.combinations(sorted(members), 2):
            d = mean_between(gt, a, b)
            if d is None:
                continue
            cross = a.split(".")[0] != b.split(".")[0]
            flag = ""
            if cross:
                flag = "  <- different genogroups; impossible under any reading"
                impossible.append((grp, a, b, d))
            print(f"  {grp}  {a:7s} {b:7s} {d * 100:5.1f}%{flag}")

    if impossible:
        print(f"\n{len(impossible)} grouping(s) place different genogroups in one")
        print("immunotype. Those rows should be cleared, not merely flagged.")
        return 1
    print("\nno cross-genogroup groupings")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
