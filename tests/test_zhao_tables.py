"""The sapovirus tables must still agree with the paper they are built from.

scripts/check_against_zhao_tables.py cross-checks our generated VP1 and RdRp
taxonomy CSVs against the two other tables Zhao et al. 2024 publishes, which
assign the same strains independently. It was a script you had to remember to
run; here it is a test, so a RELABEL that quietly contradicts a published
table fails the suite instead of surviving to a release.

Deviations that have been adjudicated live in that script's ACCEPTED dict with
the reasoning for each. Anything else is a mismatch and fails.
"""
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "scripts" / "check_against_zhao_tables.py"


def _load():
    """Import the script by path; scripts/ is not a package."""
    spec = importlib.util.spec_from_file_location("check_against_zhao_tables", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


needs_tables = pytest.mark.skipif(
    not (ROOT / "validation" / "sapovirus_vp1_taxonomy.csv").exists(),
    reason="requires the generated sapovirus taxonomy CSVs",
)


@pytest.fixture(scope="module")
def check():
    if not SCRIPT.exists():
        pytest.skip("check_against_zhao_tables.py not present")
    return _load()


@needs_tables
def test_no_unadjudicated_mismatch(check):
    problems, _, _, _ = check.compare()
    assert problems == [], (
        "our tables disagree with a published Zhao et al. table. Resolve each "
        "against Supplemental Table 1 and, if the deviation is deliberate, add "
        "it to ACCEPTED in scripts/check_against_zhao_tables.py with the "
        "evidence:\n" + "\n".join(
            f"  {a} {label}: published {want}, ours {mine}"
            for a, label, want, mine in problems))


@needs_tables
def test_accepted_deviations_are_still_deviations(check):
    """ACCEPTED must not accumulate entries that no longer apply.

    An entry that stops firing means the underlying disagreement was fixed, and
    leaving it behind makes the record of what the paper got wrong untrue.
    """
    _, _, _, accepted = check.compare()
    fired = {(a, label) for a, label, _, _ in accepted}
    stale = sorted(set(check.ACCEPTED) - fired)
    assert not stale, (
        f"ACCEPTED entries that no longer correspond to a real disagreement: {stale}. "
        "Remove them, or explain why the check stopped seeing the strain.")


@needs_tables
def test_known_adjudications_are_held(check):
    """The three strains checked by hand against the published tables.

    These were each read off Supplemental Table 1 directly, so a regression in
    the parser or in RELABEL shows up here as a wrong label rather than as a
    silently different reference set.
    """
    vp1, rdrp, _ = check.ours()
    # Supplement says GNA4/GNA4.P; Table 2 and the pig/bat table say GXV/GXV.P,
    # confirmed by hand. The relabel is deliberate and GNA4 must end up empty.
    assert vp1["MH003837"] == "GXV"
    assert rdrp["MH003837"] == "GXV.P"
    assert "GNA4" not in set(vp1.values())
    assert "GNA4.P" not in set(rdrp.values())
    # Supplement row: VP1 GXXX / P-type GXIII.P, shared with JN387134.
    assert vp1["JN420370"] == "GXXX"
    assert rdrp["JN420370"] == "GXIII.P"
    # Supplement row: VP1 GXXXI.NA1 / P-type GV.PNA2.
    assert vp1["LC215885"] == "GXXXI.NA1"
    assert rdrp["LC215885"] == "GV.PNA2"
