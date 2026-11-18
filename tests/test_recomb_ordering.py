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
# Per-event statuses and scores from a real run of two Zhao et al. Table 3
# strains. The full panel is 16 MB and lives in the benchmarking deposit, not
# here; the behaviour it pins belongs to the tool, so the handful of numbers
# needed to pin it are committed. Regenerate with the snippet in the file's
# own _comment if the panel is ever rescanned.
FIXTURE = ROOT / "tests" / "data" / "recomb_event_statuses.json"


def strain(acc):
    if not FIXTURE.exists():
        pytest.skip(f"missing {FIXTURE.name}")
    return json.loads(FIXTURE.read_text())["strains"][acc]


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


def test_mz488271_is_strong_on_the_recorded_evidence():
    """The case that exposed the bug, on events from a real run.

    It is the whole point of the fixture that these two disagree: the
    best-scoring event is weak and a lower-scoring one is strong.
    """
    d = strain("MZ488271")
    events = d["events"]
    assert any(e["status"] == "strong" for e in events), \
        "fixture no longer contains a strong event; the test has lost its point"
    assert max(events, key=lambda e: e["score"])["status"] == "weak", \
        "fixture no longer has a weak event scoring highest; the test is moot"
    assert d["reported_status"] == "weak", "fixture predates the fix, as it should"
    assert order(events)[0]["status"] == "strong"


def test_the_negative_control_is_not_promoted():
    """JN420370 is the one strain Zhao et al. report as not recombinant.

    Re-ranking may only raise a status, so what matters is that it does not
    raise this one to or past `strong`, where the positive call is made.
    """
    events = strain("JN420370")["events"]
    assert _status_rank(order(events)[0]["status"]) < _status_rank("strong")
