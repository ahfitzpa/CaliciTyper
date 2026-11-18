#!/usr/bin/env python3
"""Re-rank the deposited Table 3 events after the ordering fix, without rescanning.

`recomb.py` used to report the status of the highest-scoring event; it now
reports the best event by (status rank, score). That changes only how a
sequence's events are ranked -- not the events, not their statuses, not their
scores -- so the deposited per-event output stays valid and the strain-level
summary can be recomputed from it exactly. This is a re-aggregation, not a
re-analysis, and it is the honest thing to do while OpenRDP cannot be
installed alongside CaliciTyper (numpy<2 against numpy 2.x; see
docs/INSTALL_TOOLS.md).

Rerun `analysis/slurm/recomb_table3.sh` once the Apptainer image builds, to
confirm end to end. The statuses should not move.

    python scripts/reaggregate_recomb_table3.py [--apply]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from calicityper.recomb import _status_rank, confidence_label   # noqa: E402

PANEL = ROOT / "analysis" / "benchmark_outputs" / "recomb_table3"
SUMMARY = PANEL / "table3_summary.csv"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--apply", action="store_true", help="rewrite the summary")
    a = ap.parse_args()

    if not SUMMARY.exists():
        sys.exit(f"missing {SUMMARY}")
    rows = list(csv.DictReader(open(SUMMARY)))
    changed = []
    for r in rows:
        f = PANEL / r["accession"] / "recomb" / "recombination.json"
        if not f.exists():
            print(f"  no per-event output for {r['accession']}; left as deposited")
            continue
        events = json.loads(f.read_text())["results"][0]["recombination"]["events"]
        if not events:
            continue
        top = sorted(events, key=lambda e: (-_status_rank(e["status"]), -e["score"]))[0]
        if top["status"] != r["status"]:
            changed.append((r["accession"], r["status"], top["status"]))
        r["status"], r["score"] = top["status"], f"{top['score']:g}"

    print(f"{len(rows)} strains; {len(changed)} change status under the new ranking")
    for acc, old, new in changed:
        print(f"  {acc:11s} {old:9s} -> {new}")

    # the property that matters: re-ranking may only raise a status, and must
    # not raise the negative control to a positive call
    for acc, old, new in changed:
        if _status_rank(new) < _status_rank(old):
            sys.exit(f"{acc} was lowered ({old} -> {new}); re-ranking can only raise")
    control = next((r for r in rows if r["accession"] == "JN420370"), None)
    if control and _status_rank(control["status"]) >= _status_rank("strong"):
        sys.exit(f"the negative control is now {control['status']}: precision is lost")

    if not a.apply:
        print("\nreport only; pass --apply to rewrite the summary")
        return 0
    with open(SUMMARY, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {SUMMARY}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
