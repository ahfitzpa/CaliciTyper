#!/usr/bin/env python3
"""Calibrate the acceptance thresholds against every dataset that has a truth.

The tool carries two tunable decision thresholds, and neither was set from
data:

  * bootstrap support >= 70% to accept a norovirus type call. 70% is a
    convention inherited from phylogenetics, not a measured operating point.
  * the 2 x SD bounds of the putative-novel screen, taken from the criterion
    Chhabra et al. describe.

This sweeps the support threshold over three datasets that differ in what
"truth" means, because the answer should not rest on one of them:

  1. norovirus amplicons   true genotype, known by construction
                           (Fitzpatrick et al. 2022, Zenodo 7271067)
  2. reference leave-one-out  true genotype, the published label of a
                           reference strain re-typed against the others
  3. the Caliciviridae panel  the Norovirus Typing Tool's own calls, which are
                           a reference standard rather than a ground truth

For each it reports precision, recall and F1 across thresholds, the error
taxonomy, and a per-class breakdown, so that a threshold change (or a decision
to leave it alone) can be argued from the numbers.

    python scripts/calibrate_thresholds.py --out-dir results/calibration
"""
from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
HOME = Path("~/scratch").expanduser()
RUNS = HOME / "calicityper_runs"

THRESHOLDS = [0, 30, 40, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95]


def load_amplicons():
    f = ROOT / "results" / "benchmark_zenodo" / "zenodo_comparison.csv"
    if not f.exists():
        return None
    d = pd.read_csv(f).fillna("")
    d = d[d.expected_genotype.astype(str) != ""].copy()
    return pd.DataFrame(dict(truth=d.expected_genotype.astype(str),
                             call=d.called_genotype.astype(str),
                             support=pd.to_numeric(d.support, errors="coerce"),
                             region="orf2", dataset="norovirus amplicons"))


def load_loo():
    f = RUNS / "benchmarking_wip_0.7.0" / "results" / "loo_0.7.0.csv"
    if not f.exists():
        return None
    d = pd.read_csv(f).fillna("")
    d = d[d.length.astype(str).isin(["0", ""]) | True]            # all lengths
    return pd.DataFrame(dict(truth=d.truth.astype(str), call=d.call.astype(str),
                             support=pd.to_numeric(d.support, errors="coerce"),
                             region=d.region.astype(str), dataset="reference leave-one-out"))


def load_panel():
    ours = RUNS / "rivm_genogroup_dev" / "all_typed.tsv"
    rivm = ROOT / "validation" / "results_rivm.csv"
    if not (ours.exists() and rivm.exists()):
        return None
    import validate_vs_rivm as V
    o = pd.read_csv(ours, sep="\t", dtype=str).fillna("")
    r = pd.read_csv(rivm, dtype=str).fillna("").set_index("name")
    rows = []
    for _, q in o.iterrows():
        if q["query"] not in r.index:
            continue
        for region, ourcol, rivmcol in (("orf2", "orf2_type", "capsid type"),
                                        ("orf1", "orf1_type", "polymerase type")):
            truth = V.genotype(r.at[q["query"], rivmcol])
            call = V.genotype(q.get(ourcol, ""))
            if truth in ("", "unassigned"):
                continue
            rows.append(dict(truth=truth, call=call,
                             support=pd.to_numeric(q.get(f"{region}_support", ""), errors="coerce"),
                             region=region, dataset="Caliciviridae panel (comparator labels)"))
    return pd.DataFrame(rows)


def sweep(d, name, out, fh):
    d = d.dropna(subset=["support"]).copy()
    d["typed"] = ~d.call.isin(["", "unassigned", "untypeable", "nan"])
    d["ok"] = d.typed & (d.call == d.truth)
    typed = d[d.typed]
    if not len(typed):
        return
    total_correct = int(d.ok.sum())
    print(f"\n{'=' * 78}\n{name}: {len(d)} comparisons, {len(typed)} typed, "
          f"{total_correct} correct\n{'=' * 78}")
    print(f"{'thresh':>6} {'accepted':>8} {'correct':>7} {'wrong':>5} "
          f"{'precision':>9} {'recall':>7} {'F1':>7}")
    best = None
    for t in THRESHOLDS:
        a = typed[typed.support >= t]
        if not len(a):
            continue
        p = float(a.ok.mean())
        rec = float(a.ok.sum()) / max(1, total_correct)
        f1 = 2 * p * rec / (p + rec) if (p + rec) else 0.0
        mark = "  <- current" if t == 70 else ""
        if best is None or f1 > best[1]:
            best = (t, f1)
        print(f"{t:>6} {len(a):>8} {int(a.ok.sum()):>7} {int((~a.ok).sum()):>5} "
              f"{p:>9.4f} {rec:>7.4f} {f1:>7.4f}{mark}")
        fh.writerow(dict(dataset=name, threshold=t, accepted=len(a),
                         correct=int(a.ok.sum()), wrong=int((~a.ok).sum()),
                         precision=round(p, 5), recall=round(rec, 5), f1=round(f1, 5)))
    print(f"best F1 at {best[0]}%")
    # smallest threshold that reaches perfect precision
    perfect = [t for t in THRESHOLDS
               if len(typed[typed.support >= t]) and typed[typed.support >= t].ok.all()]
    if perfect:
        t0 = min(perfect)
        kept = len(typed[typed.support >= t0])
        print(f"first threshold with no wrong call: {t0}% (keeps {kept}/{len(typed)} calls)")
    else:
        print("no threshold eliminates every wrong call")
    wrong = typed[~typed.ok]
    if len(wrong):
        sup = sorted(int(s) for s in wrong.support.dropna())
        print(f"support of wrong calls: {sup[:25]}{' ...' if len(sup) > 25 else ''}")
        print(f"support of correct calls: median {typed[typed.ok].support.median():.0f}, "
              f"5th pct {typed[typed.ok].support.quantile(0.05):.0f}")
        pairs = collections.Counter(zip(wrong.truth, wrong.call)).most_common(8)
        print("commonest confusions (truth -> call: n):")
        for (a, b), n in pairs:
            print(f"    {a:12s} -> {b:12s} {n}")
    for region, sub in d.groupby("region"):
        st = sub[sub.typed]
        if len(st):
            print(f"  [{region}] typed {len(st)}, correct {int(st.ok.sum())} "
                  f"({100*st.ok.mean():.1f}%)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out-dir", default="results/calibration")
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    datasets = [("norovirus amplicons (true genotype)", load_amplicons()),
                ("reference leave-one-out (true genotype)", load_loo()),
                ("Caliciviridae panel (comparator labels)", load_panel())]
    with open(out / "calibration_sweep.csv", "w", newline="") as f:
        fh = csv.DictWriter(f, fieldnames=["dataset", "threshold", "accepted", "correct",
                                           "wrong", "precision", "recall", "f1"])
        fh.writeheader()
        for name, d in datasets:
            if d is None or not len(d):
                print(f"\nSKIP {name}: inputs missing")
                continue
            sweep(d, name, out, fh)
    print(f"\nwrote {out}/calibration_sweep.csv")
    print("\nNote: the panel's labels are another tool's calls, so 'wrong' there means")
    print("disagreement, not error. Only the first two datasets carry a ground truth.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
