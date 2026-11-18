"""The reported status must be the strongest evidence found, not the best score.

`recomb.py` used to sort a sequence's events by score alone and report the
first one. The score is a weighted blend of five components and is not a
measure of status, so a single-method event with a clean breakpoint and
distinct parents could outscore a four-method event and hide it. MZ488271 in
the Zhao et al. Table 3 panel was reported `weak` while containing a `strong`
event. Events are now ordered by (status rank, score).
"""
import json
from pathlib import Path

import pytest

from calicityper.recomb import STATUSES, _status_rank

ROOT = Path(__file__).resolve().parent.parent
PANEL = ROOT / "analysis" / "benchmark_outputs" / "recomb_table3"


def order(events):
    """The ordering recomb.py applies; events[0] becomes the reported status."""
    return sorted(events, key=lambda e: (-_status_rank(e["status"]), -e["score"]))


def test_status_rank_is_a_total_order_over_event_statuses():
    ranks = [_status_rank(s) for s in STATUSES[:5]]
    assert ranks == sorted(ranks) == [0, 1, 2, 3, 4]
    # the trailing three mean "not answered" and are never set on an event
    assert all(s not in STATUSES[:5] for s in STATUSES[5:])


def test_a_strong_event_is_never_hidden_by_a_better_scoring_weak_one():
    events = [{"status": "weak", "score": 0.90}, {"status": "strong", "score": 0.10}]
    assert order(events)[0]["status"] == "strong"


def test_score_still_breaks_ties_within_a_status():
    events = [{"status": "strong", "score": 0.40}, {"status": "strong", "score": 0.80}]
    top = order(events)[0]
    assert top["status"] == "strong" and top["score"] == 0.80


@pytest.mark.skipif(not (PANEL / "MZ488271" / "recomb" / "recombination.json").exists(),
                    reason="requires the deposited Table 3 panel")
def test_mz488271_is_strong_on_the_deposited_evidence():
    """The case that exposed the bug, checked against the real stored events.

    The evidence is unchanged -- only the ranking of it is -- so this reads
    the deposited per-event output rather than re-running the scanners.
    """
    f = PANEL / "MZ488271" / "recomb" / "recombination.json"
    events = json.loads(f.read_text())["results"][0]["recombination"]["events"]
    assert any(e["status"] == "strong" for e in events), \
        "fixture no longer contains a strong event; the test has lost its point"
    assert order(events)[0]["status"] == "strong"


@pytest.mark.skipif(not PANEL.exists(), reason="requires the deposited Table 3 panel")
def test_the_negative_control_is_not_promoted():
    """JN420370 is the one strain Zhao et al. report as not recombinant.

    Re-ranking may only raise a status, so the thing to check is that it does
    not raise this one to or past `strong`, which is where the positive call
    is made.
    """
    f = PANEL / "JN420370" / "recomb" / "recombination.json"
    if not f.exists():
        pytest.skip("negative control not deposited")
    events = json.loads(f.read_text())["results"][0]["recombination"]["events"]
    assert _status_rank(order(events)[0]["status"]) < _status_rank("strong")
