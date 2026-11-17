#!/usr/bin/env python3
"""Does the choice of distance metric explain the sapovirus RdRp cut-off failure?

Zhao et al. computed pairwise distances in MEGA 7.0.26 "with default
parameters". MEGA's default nucleotide distance is Maximum Composite
Likelihood, a Tamura-Nei-based estimate that corrects for multiple
substitutions. We apply their cut-offs to an UNCORRECTED p-distance.

The two agree at low divergence and diverge as distance grows, which predicts
exactly the pattern observed: our within-type distances (low divergence)
reproduce the published ranges 18/20, while our between-type distances (high
divergence) fall below the published between-type threshold.

This compares both metrics against the published expectation:
    within a type   <= 0.201 (RdRp strain peak)   / <= 0.133 (VP1)
    between types   >= 0.266 (RdRp genotype peak) / >= 0.161 (VP1)

    python scripts/test_distance_metric.py
"""
from __future__ import annotations

import csv
import itertools
import sys
from pathlib import Path

import numpy as np
from Bio import AlignIO, SeqIO

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from calicityper import config as cg, tree   # noqa: E402

# (region key, taxonomy csv, label column, within-peak top, between-peak floor)
REGIONS = {
    "RdRp": ("orf1", ROOT / "validation" / "sapovirus_rdrp_taxonomy.csv", "ptype", 0.201, 0.266),
    "VP1": ("orf2", ROOT / "validation" / "sapovirus_vp1_taxonomy.csv", "genotype", 0.133, 0.161),
}


def published_ranges(path, key):
    out = {}
    for r in csv.DictReader(l for l in open(path) if not l.startswith("#")):
        if r["note"].startswith(("UNRESOLVED", "RETIRED")):
            continue
        t = r["pairwise_distance_range"].strip()
        if t and t != "NA":
            lo, _, hi = t.partition("-")
            out[r[key]] = (float(lo), float(hi or lo))
    return out


def p_matrix(codes):
    """Uncorrected p-distance, pairwise deletion -- what the tool uses today."""
    n = len(codes)
    D = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        a, b = codes[i], codes[j]
        both = (a < 4) & (b < 4)
        m = int(both.sum())
        D[i, j] = D[j, i] = (float((a[both] != b[both]).sum()) / m) if m >= 50 else np.nan
    return D


def main():
    for region, (key_, tax, key, within_top, between_floor) in REGIONS.items():
        pub = published_ranges(tax, key)
        d = cg.SAPOVIRUS_REF_DIR[key_]
        rows = []
        print(f"\n{'=' * 86}\n{region}: p-distance (current) against Tamura-Nei (MEGA-like correction)\n{'=' * 86}")
        for orf in sorted(d.glob("*.orf")):
            if orf.stem == "all":
                continue
            types = {r.id.split("|")[2]: r.id.split("|")[1] for r in SeqIO.parse(orf, "fasta")}
            aln = AlignIO.read(orf.with_suffix(".aln"), "clustal")
            names = [r.id for r in aln]
            if len(names) < 2:
                continue
            codes = tree.encode([str(r.seq).upper() for r in aln])
            Dp = p_matrix(codes)
            Dt = np.asarray(tree.tn93(codes), dtype=float)
            lab = [types[n] for n in names]
            for t in sorted(set(lab)):
                idx = [i for i, x in enumerate(lab) if x == t]
                oth = [i for i, x in enumerate(lab) if x != t]
                def agg(D, fn, pairs):
                    vals = [D[i, j] for i, j in pairs if np.isfinite(D[i, j])]
                    return fn(vals) if vals else np.nan
                win = list(itertools.combinations(idx, 2))
                bet = [(i, j) for i in idx for j in oth]
                rows.append(dict(
                    group=orf.stem, type=t, n=len(idx),
                    p_within=agg(Dp, max, win), t_within=agg(Dt, max, win),
                    p_between=agg(Dp, min, bet), t_between=agg(Dt, min, bet),
                    pub_hi=pub.get(t, (np.nan, np.nan))[1]))
        # ---- scoring
        def score(rows, wk, bk):
            w = [r for r in rows if np.isfinite(r[wk]) and np.isfinite(r["pub_hi"])]
            w_ok = sum(r[wk] <= r["pub_hi"] + 0.005 for r in w)
            wp = [r for r in rows if np.isfinite(r[wk])]
            wp_ok = sum(r[wk] <= within_top for r in wp)
            b = [r for r in rows if np.isfinite(r[bk])]
            b_ok = sum(r[bk] >= between_floor for r in b)
            return (w_ok, len(w)), (wp_ok, len(wp)), (b_ok, len(b))
        print(f"{'metric':>12} {'within published range':>24} {'within peak <=' + str(within_top):>22} "
              f"{'separated >=' + str(between_floor):>22}")
        for metric, wk, bk in (("p-distance", "p_within", "p_between"),
                               ("Tamura-Nei", "t_within", "t_between")):
            (a, A), (c, C), (e, E) = score(rows, wk, bk)
            print(f"{metric:>12} {a:12d}/{A:<11d} {c:12d}/{C:<9d} {e:12d}/{E:<9d}")
        print(f"\n  types where Tamura-Nei lifts the nearest other type above {between_floor}:")
        fixed = [r for r in rows if np.isfinite(r["p_between"]) and np.isfinite(r["t_between"])
                 and r["p_between"] < between_floor <= r["t_between"]]
        for r in sorted(fixed, key=lambda r: -r["t_between"])[:12]:
            print(f"    {r['type']:10s} p {r['p_between']:.3f} -> TN {r['t_between']:.3f}")
        print(f"  ({len(fixed)} types)")
        still = [r for r in rows if np.isfinite(r["t_between"]) and r["t_between"] < between_floor]
        if still:
            print(f"  still below after correction:")
            for r in sorted(still, key=lambda r: r["t_between"])[:8]:
                print(f"    {r['type']:10s} p {r['p_between']:.3f} -> TN {r['t_between']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
