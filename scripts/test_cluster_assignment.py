#!/usr/bin/env python3
"""Nearest-neighbour assignment against cluster-mean assignment for sapovirus.

Zhao et al. state the criterion as the MEAN distance between a sequence and
its nearest cluster(s), and derive the cut-offs as the mean +/- 2SD (RdRp) or
2.5SD (VP1) of pooled pairwise-distance peaks. Two consequences our
implementation currently ignores:

 1. A cut-off set at mean - 2SD of the between-type peak has, by construction,
    about 2.5% of between-type PAIRS below it. Asking whether the MINIMUM
    distance to another type clears the cut-off therefore tests the extreme
    tail and will fail for any reasonably sized set -- which is what we saw
    (1/19 for RdRp).
 2. Assigning a query to the type of its single nearest reference is not the
    published criterion. Where clusters are close, the nearest sequence can
    belong to a neighbouring type even when the query sits squarely inside
    its own cluster.

This compares, by leave-one-out over the reference sets:
    NN   assign the type of the nearest single reference (current behaviour)
    MEAN assign the type with the lowest mean distance to its members
and reports, for each, accuracy on types with a partner and the behaviour on
singletons, whose only correct outcome once removed is "untypeable".

    python scripts/test_cluster_assignment.py
"""
from __future__ import annotations

import csv
import itertools
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from Bio import AlignIO, SeqIO

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from calicityper import config as cg, tree   # noqa: E402

REGIONS = {"RdRp": ("orf1", 0.266), "VP1": ("orf2", 0.161)}


def p_matrix(codes):
    n = len(codes)
    D = np.full((n, n), np.nan)
    for i, j in itertools.combinations(range(n), 2):
        a, b = codes[i], codes[j]
        both = (a < 4) & (b < 4)
        m = int(both.sum())
        if m >= 50:
            D[i, j] = D[j, i] = float((a[both] != b[both]).sum()) / m
    return D


def main():
    for region, (key, cut) in REGIONS.items():
        d = cg.SAPOVIRUS_REF_DIR[key]
        res = {m: Counter() for m in ("NN", "MEAN")}
        singles = {m: Counter() for m in ("NN", "MEAN")}
        wrong = defaultdict(list)
        for orf in sorted(d.glob("*.orf")):
            if orf.stem == "all":
                continue
            types = {r.id.split("|")[2]: r.id.split("|")[1] for r in SeqIO.parse(orf, "fasta")}
            aln = AlignIO.read(orf.with_suffix(".aln"), "clustal")
            names = [r.id for r in aln]
            if len(names) < 2:
                continue
            D = p_matrix(tree.encode([str(r.seq).upper() for r in aln]))
            lab = [types[n] for n in names]
            counts = Counter(lab)
            for q in range(len(names)):
                truth = lab[q]
                singleton = counts[truth] == 1
                others = [i for i in range(len(names)) if i != q]
                if not others:
                    continue
                # --- nearest neighbour
                vals = [(D[q, i], lab[i]) for i in others if np.isfinite(D[q, i])]
                if not vals:
                    continue
                nn_d, nn_t = min(vals)
                nn_call = nn_t if nn_d < cut else "untypeable"
                # --- cluster mean
                means = []
                for t in set(lab[i] for i in others):
                    v = [D[q, i] for i in others if lab[i] == t and np.isfinite(D[q, i])]
                    if v:
                        means.append((float(np.mean(v)), t))
                mean_d, mean_t = min(means)
                mean_call = mean_t if mean_d < cut else "untypeable"
                for m, call in (("NN", nn_call), ("MEAN", mean_call)):
                    tgt = singles[m] if singleton else res[m]
                    if singleton:
                        tgt["correct" if call == "untypeable" else "given a known label"] += 1
                    else:
                        tgt["correct" if call == truth else
                            ("untypeable" if call == "untypeable" else "wrong")] += 1
                        if m == "MEAN" and call not in (truth, "untypeable"):
                            wrong["MEAN"].append((names[q], truth, call, round(mean_d, 3)))
                        if m == "NN" and call not in (truth, "untypeable"):
                            wrong["NN"].append((names[q], truth, call, round(nn_d, 3)))
        print(f"\n{'=' * 78}\n{region} (cut-off {cut})\n{'=' * 78}")
        for m in ("NN", "MEAN"):
            tot = sum(res[m].values())
            lab = "nearest reference (current)" if m == "NN" else "lowest mean distance to a cluster"
            print(f"  {lab}")
            print(f"    types with a partner (n={tot}): "
                  + ", ".join(f"{k} {v} ({100*v/tot:.0f}%)" for k, v in sorted(res[m].items())))
            st = sum(singles[m].values())
            if st:
                print(f"    singletons (n={st}): "
                      + ", ".join(f"{k} {v}" for k, v in sorted(singles[m].items())))
        for m in ("NN", "MEAN"):
            if wrong[m]:
                print(f"  wrong calls under {m}:")
                for a, t, c, v in wrong[m][:8]:
                    print(f"    {a:12s} {t:10s} -> {c:10s} (d={v})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
