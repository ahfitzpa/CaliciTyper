#!/usr/bin/env python3
"""Sweep the recombination status rules over the Zhao Table 3 panel.

The question this answers: if the evidence rules are weakened, is there a
setting that raises F1 **without** giving up precision? Precision is the
property worth protecting -- a tool that announces recombination where there
is none is worse than one that stays quiet -- so a setting that trades it for
recall is not an improvement however much F1 rises.

Nothing is refitted here. The scoring already ran; this re-applies
`assign_status` to the stored per-event evidence (number of supporting
methods, phylogenetic support, parental distinctness) under different rules,
which is cheap and exactly reproduces what the tool would have reported.

The panel is 12 strains Zhao et al. report as recombinant plus JN420370, the
sea lion strain they report as NOT recombinant, which is the only negative.
One negative cannot measure a false-positive rate; it can only show that a
rule has stopped being conservative. Treat precision here as a tripwire, not
an estimate, and do not tune on 13 strains.

    python scripts/sweep_recomb_rules.py
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from calicityper.recomb import STATUSES, _status_rank    # noqa: E402

PANEL = ROOT / "analysis" / "benchmark_outputs" / "recomb_table3"
ORDER = list(STATUSES[:5])

# accession -> is a recombinant according to Zhao et al. Table 3
TRUTH = {
    "JN420370": False,          # sea lion: reported NOT recombinant
    "GQ261222": True, "MN794218": True, "MK291482": True, "MZ488271": True,
    "MZ488272": True, "MG012446": True, "KP067444": True, "MK296750": True,
    "MH541043": True, "MH541048": True, "AB455793": True, "LC504410": True,
}


def events():
    """{accession: [(n_methods, phylo_score, distinctness, event_score)]}.

    The event score is carried through because the strain-level result is the
    best event under `(status_rank, score)` -- status first, score second --
    which is how `recomb.py` now orders them. It used to order by score alone,
    which let a better-scoring weak event hide a strong one; the ordering is
    imported from `recomb.py` rather than restated here so the two cannot
    drift. The score itself does not depend on the status rules, so the stored
    value stays valid under every variant swept below.
    """
    out = {}
    for acc in sorted(TRUTH):
        p = PANEL / acc / "recomb" / "recombination.json"
        if not p.exists():
            print(f"  missing {p}", file=sys.stderr)
            continue
        ev = []
        for res in json.loads(p.read_text())["results"]:
            for e in res["recombination"]["events"]:
                ev.append((len(e["supporting_methods"]),
                           float(e.get("phylogenetic_detail", {}).get("score", 0.0)),
                           float(e.get("parental_distinctness", 0.0)),
                           float(e.get("score", 0.0)),
                           e.get("status", "none")))
        out[acc] = ev
    return out


def status(n, phylo, distinct, *, strong_methods, likely_methods,
           possible_methods, phylo_cut, strong_phylo, strong_distinct,
           single_method_phylo):
    """`recomb.assign_status`, with every rule exposed as a parameter."""
    phylo_ok, distinct_ok = phylo >= phylo_cut, distinct >= 0.5
    if (n >= strong_methods and (phylo_ok or not strong_phylo)
            and (distinct_ok or not strong_distinct)):
        return "strong"
    if n >= likely_methods and phylo_ok:
        return "likely"
    if n >= possible_methods or (n >= 1 and phylo >= single_method_phylo and distinct_ok):
        return "possible"
    if n >= 1:
        return "weak"
    return "none"


def score(ev, call_at, **rules):
    """Precision/recall/F1 when a strain counts as recombinant at `call_at`."""
    want = ORDER.index(call_at)
    tp = fp = fn = 0
    for acc, evs in ev.items():
        if not evs:
            top = "none"
        else:
            scored = [(status(n, ph, di, **rules), sc) for n, ph, di, sc, _ in evs]
            top = max(scored, key=lambda t: (_status_rank(t[0]), t[1]))[0]
        called = ORDER.index(top) >= want
        if called and TRUTH[acc]:
            tp += 1
        elif called:
            fp += 1
        elif TRUTH[acc]:
            fn += 1
    prec = tp / (tp + fp) if tp + fp else float("nan")
    rec = tp / (tp + fn) if tp + fn else float("nan")
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    return prec, rec, f1, tp, fp, fn


def main():
    ev = events()
    n_ev = sum(len(v) for v in ev.values())
    base = dict(strong_methods=3, likely_methods=3, possible_methods=2,
                phylo_cut=0.5, strong_phylo=True, strong_distinct=True,
                single_method_phylo=0.7)
    # Refuse to report a sweep from a reimplementation that does not reproduce
    # the tool. The check is per EVENT, against the status the tool stored on
    # that event, which tests `assign_status` without assuming anything about
    # how events are then ranked -- the aggregation itself is imported, not
    # restated. An earlier version checked only the aggregate and was wrong on
    # 3 of 13 strains, every one of them in the flattering direction.
    bad = [(acc, stored, status(n, ph, di, **base))
           for acc, evs in ev.items() for n, ph, di, _sc, stored in evs
           if status(n, ph, di, **base) != stored]
    if bad:
        raise SystemExit(f"recompute disagrees with the tool on {len(bad)} event(s), "
                         f"e.g. {bad[:3]}")
    print(f"recompute reproduces the tool's status for all {n_ev} stored events "
          f"across {len(ev)} strains\n")

    print(f"{len(ev)} strains, {n_ev} scored events "
          f"({sum(TRUTH.values())} recombinant, {sum(not v for v in TRUTH.values())} not)\n")

    print("=== as shipped ===")
    for call_at in ("strong", "likely", "possible"):
        p, r, f, tp, fp, fn = score(ev, call_at, **base)
        print(f"  call at {call_at:9s} P={p:.3f} R={r:.3f} F1={f:.3f}  "
              f"(TP {tp}, FP {fp}, FN {fn})")

    print("\n=== sweep: every variant that keeps precision at 1.000 ===")
    rows = []
    grid = itertools.product(
        ("strong", "likely"),          # where the positive call is made
        (2, 3),                        # strong_min_methods
        (2, 3),                        # likely_min_methods
        (0.3, 0.4, 0.5, 0.6),          # phylo_cut
        (True, False),                 # strong_require_phylo
        (True, False),                 # strong_require_distinct
    )
    for call_at, sm, lm, pc, sp, sd in grid:
        rules = dict(base, strong_methods=sm, likely_methods=lm, phylo_cut=pc,
                     strong_phylo=sp, strong_distinct=sd)
        p, r, f, tp, fp, fn = score(ev, call_at, **rules)
        rows.append((f, r, p, call_at, sm, lm, pc, sp, sd, tp, fp, fn))
    rows.sort(reverse=True)

    print(f"{'F1':>6} {'recall':>7} {'prec':>6}  call_at    strong_n likely_n "
          f"phylo_cut  req_phylo req_distinct   TP FP FN")
    shown = 0
    for f, r, p, call_at, sm, lm, pc, sp, sd, tp, fp, fn in rows:
        if p < 1.0:
            continue
        print(f"{f:6.3f} {r:7.3f} {p:6.3f}  {call_at:9s} {sm:8d} {lm:8d} "
              f"{pc:9.1f}  {str(sp):9s} {str(sd):12s} {tp:2d} {fp:2d} {fn:2d}")
        shown += 1
        if shown >= 14:
            print("  ...")
            break

    print("\n=== the best F1 at any precision, for contrast ===")
    f, r, p, call_at, sm, lm, pc, sp, sd, tp, fp, fn = rows[0]
    print(f"  F1={f:.3f} P={p:.3f} R={r:.3f}  call_at={call_at} strong_n={sm} "
          f"likely_n={lm} phylo_cut={pc} req_phylo={sp} req_distinct={sd}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
