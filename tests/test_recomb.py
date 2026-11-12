"""Tests for `calicityper recomb`.

Two layers:

* Unit tests over clustering, scoring, status assignment, coordinate mapping,
  p-value correction, parsing and QC. These are deterministic, need no external
  tool, and always run.
* End-to-end tests over the six seeded synthetic datasets in tests/synthetic.py,
  which need OpenRDP and are skipped without it (`CALICITYPER_OPENRDP`). These
  are the datasets BENCHMARKING.md section 4 quotes, so they are what keeps
  those numbers honest.

Real reference genomes cannot serve as negative controls here -- a published
norovirus strain frequently IS a recombinant -- so specificity is measured
against sequences that are non-recombinant by construction.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import synthetic                                            # noqa: E402

from calicityper import recomb as R                          # noqa: E402
from calicityper import recomb_config as RC                  # noqa: E402

ORF1_ORF2_JUNCTION = synthetic.ORF1_ORF2_JUNCTION


@pytest.fixture
def cfg():
    return RC.load(None)


def hit(method, start, end, p=0.001, rec="q", p1="A", p2="B"):
    return R.Hit(method=method, start=start, end=end, recombinant=rec,
                 parent1=p1, parent2=p2, pvalue=p)


# ----------------------------------------------------------------- coordinates

def test_coordinate_map_counts_ungapped_positions():
    cmap = R.coordinate_map("AC-GT")
    assert list(cmap) == [1, 2, 2, 3, 4]        # the gap holds the last real position


def test_coordinate_map_leading_gap_is_zero():
    """A breakpoint before the query's first base has no genome position."""
    assert list(R.coordinate_map("--AC")) == [0, 0, 1, 2]


# -------------------------------------------------------------------- parsing

def test_parse_openrdp_normalises_and_survives_bad_rows():
    text = (
        "Method,Start,End,Recombinant,Parent1,Parent2,Pvalue\n"
        "MaxChi,100,200,q,A,B,0.001\n"
        "Bootscan,150.0,250.0,q,-,B,1e-4\n"          # floats, '-' parent
        "Geneconv,oops,250,q,A,B,0.01\n"             # malformed: dropped
        ",1,2,q,A,B,0.01\n"                          # no method: dropped
    )
    hits = R.parse_openrdp(text)
    assert [h.method for h in hits] == ["maxchi", "bootscan"]
    assert hits[1].parent1 == ""                     # '-' becomes empty, not '-'
    assert hits[0].start == 100 and hits[1].end == 250


# ------------------------------------------------------------------ breakpoints

def test_narrow_hit_reports_its_centre_and_wide_hit_its_edges():
    """Window scanners report the window containing a breakpoint; GENECONV and
    3Seq report a converted fragment, whose edges are the breakpoints."""
    assert hit("maxchi", 1000, 1100).breakpoints(narrow_span=600) == [1050.0]
    assert hit("geneconv", 1000, 4000).breakpoints(narrow_span=600) == [1000.0, 4000.0]


# ------------------------------------------------------------------ clustering

def test_cluster_hits_groups_one_locus_as_one_event(cfg):
    """Five methods agreeing at one locus is one event, not five."""
    hits = [(hit(m, 5000, 5100), "query")
            for m in ("maxchi", "chimaera", "bootscan", "siscan", "rdp")]
    clusters = R.cluster_hits(hits, cfg)
    assert len(clusters) == 1
    assert len(clusters[0]["methods"]) == 5


def test_cluster_hits_separates_distant_loci(cfg):
    hits = [(hit("maxchi", 1000, 1100), "query"),
            (hit("chimaera", 6000, 6100), "query")]
    clusters = R.cluster_hits(hits, cfg)
    assert len(clusters) == 2


def test_cluster_hits_drops_non_significant_and_counts_a_hit_once(cfg):
    """A wide hit proposes two breakpoints; if both land in one cluster it must
    still count once toward that event's method support."""
    wide = hit("geneconv", 5000, 5000 + cfg.max_breakpoint_distance // 2)
    clusters = R.cluster_hits([(wide, "query"),
                               (hit("maxchi", 5000, 5100, p=0.9), "query")], cfg)
    assert len(clusters) == 1
    assert clusters[0]["methods"] == ["geneconv"]          # the p=0.9 hit is gone
    assert len(clusters[0]["hits"]) == 1                   # counted once, not twice


# -------------------------------------------------------------------- p-values

def test_bonferroni_applies_only_across_tests_performed(cfg):
    hits = [hit("maxchi", 1, 2, p=0.01)]
    assert R.correct_pvalues(hits, cfg, n_tests=1)[id(hits[0])] == 0.01
    corrected = R.correct_pvalues(hits, cfg, n_tests=5)[id(hits[0])]
    assert corrected == pytest.approx(0.05)


# ---------------------------------------------------------------------- status

def test_status_is_conservative_for_a_single_method(cfg):
    """One method plus a weak phylogeny hint produced 15 spurious 'possible'
    events on a non-recombinant control, so it must stay 'weak'."""
    assert R.assign_status(1, phylo_score=0.6, distinctness=0.9, cfg=cfg) == "weak"
    assert R.assign_status(1, phylo_score=0.7, distinctness=0.9, cfg=cfg) == "possible"
    assert R.assign_status(1, phylo_score=0.7, distinctness=0.1, cfg=cfg) == "weak"


def test_status_rises_with_independent_methods(cfg):
    weak = R.assign_status(1, 0.0, 0.0, cfg)
    many = R.assign_status(5, 0.9, 0.9, cfg)
    assert weak == "weak"
    assert many in ("likely", "strong")
    order = ["none", "weak", "possible", "likely", "strong"]
    assert order.index(many) > order.index(weak)


def test_score_is_bounded_and_weighted(cfg):
    assert R.score_event({}, cfg) == 0.0
    allmax = R.score_event({k: 1.0 for k in cfg.scoring.weights()}, cfg)
    assert 0.0 < allmax <= 1.0
    assert R.score_event({k: 99.0 for k in cfg.scoring.weights()}, cfg) == 1.0
    assert R.confidence_label(0.8) == "high"
    assert R.confidence_label(0.6) == "medium"
    assert R.confidence_label(0.2) == "low"


def test_biological_plausibility_prefers_the_junction(cfg):
    """Norovirus recombines at the ORF1/ORF2 junction, so a breakpoint there
    scores above one in the middle of ORF1."""
    at_junction = R.biological_plausibility_score(
        ORF1_ORF2_JUNCTION, ORF1_ORF2_JUNCTION, 7500, cfg)
    mid_orf1 = R.biological_plausibility_score(2000, 2000, 7500, cfg)
    assert at_junction > mid_orf1


# -------------------------------------------------------------------------- QC

def test_too_short_is_refused_rather_than_guessed(cfg, tmp_path):
    """Below the informative-site floor the scanners cannot say anything, and
    guessing would be worse than refusing."""
    p = tmp_path / "short.fasta"
    synthetic.too_short(str(p))
    aln = R.read_alignment(str(p))
    ok, reason, stats = R.qc_query(aln, aln.names[-1], cfg)
    assert not ok
    assert "min_sequence_length" in reason, reason
    assert stats["ungapped_length"] == 900


def test_informative_sites_ignores_gaps_and_singletons(tmp_path):
    p = tmp_path / "a.fasta"
    synthetic.write_fasta(str(p), [("s1", "ACGT"), ("s2", "ACGT"),
                                   ("s3", "ACGA"), ("s4", "ACG-")])
    aln = R.read_alignment(str(p))
    # column 4 has T,T,A,- : A appears once, so it is not informative
    assert R.informative_sites(aln) == 0


# ------------------------------------------------------- synthetic end-to-end

needs_openrdp = pytest.mark.skipif(
    not RC.load(None).openrdp(),
    reason="OpenRDP not installed; set CALICITYPER_OPENRDP")
# ~50 min per dataset (measured, job 1126058): deselected unless -m slow.
slow = pytest.mark.slow

# dataset -> (statuses allowed, expected breakpoint or None)
EXPECTED = {
    "non_recombinant": ({"none", "weak"}, None),
    "clear_recombinant": ({"possible", "likely", "strong"}, ORF1_ORF2_JUNCTION),
    "weak_signal": ({"none", "weak", "possible"}, None),
    "conflicting": ({"weak", "possible", "likely", "strong"}, None),
    # Near-identical donors leave too few informative sites, so the QC gate
    # refuses rather than guessing -- the correct answer, and the reason this
    # dataset exists. Measured: status insufficient_sequence, 0 events.
    "orf_mismatch_only": ({"none", "weak", "possible", "likely",
                           "insufficient_sequence"}, None),
}


def _run(name, tmp_path, cfg):
    """Returns (recombination block, provenance). A result is
    {sequence_id, qc, recombination{status, score, confidence, events, ...}}."""
    tmp_path = Path(tmp_path)
    tmp_path.mkdir(parents=True, exist_ok=True)
    p = tmp_path / f"{name}.fasta"
    getattr(synthetic, name)(str(p))
    aln = R.read_alignment(str(p))
    query = aln.names[-1]
    results, provenance = R.analyse(str(p), [query], cfg, workdir=str(tmp_path / "work"))
    return results[0]["recombination"], provenance


@needs_openrdp
@slow
@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_synthetic_dataset_status(name, tmp_path, cfg):
    allowed, _ = EXPECTED[name]
    res, _ = _run(name, tmp_path, cfg)
    detail = (f"{name}: status={res['status']!r} not in {sorted(allowed)}; "
              f"score={res['score']} confidence={res['confidence']} "
              f"events={len(res['events'])}; "
              f"methods per event="
              f"{[e.get('supporting_methods') for e in res['events'][:4]]}")
    assert res["status"] in allowed, detail


@needs_openrdp
@slow
def test_clear_recombinant_breakpoint_is_localised(tmp_path, cfg):
    """The planted breakpoint must be recovered, and within
    max_breakpoint_distance of the truth -- localisation, not just detection."""
    res, _ = _run("clear_recombinant", tmp_path, cfg)
    events = res.get("events") or []
    assert events, res
    # an event carries a LIST of breakpoints, each with query coordinates
    spans = [(b["query_start"], b["query_end"]) for e in events for b in e["breakpoints"]]
    assert spans, events
    assert any(lo - cfg.max_breakpoint_distance <= ORF1_ORF2_JUNCTION
               <= hi + cfg.max_breakpoint_distance for lo, hi in spans), spans


@needs_openrdp
@slow
def test_orf_mismatch_only_never_reaches_strong(tmp_path, cfg):
    """Distinct ORF1/ORF2 ancestry with near-identical donors has low parental
    distinctness, which must cap the status."""
    res, _ = _run("orf_mismatch_only", tmp_path, cfg)
    assert res["status"] != "strong", res


@needs_openrdp
@slow
def test_same_seed_gives_the_same_result(tmp_path, cfg):
    """Seeded end to end, so a failure is reproducible and a number quoted in
    the docs can be re-derived.

    If this fails, OpenRDP is not deterministic between runs (3Seq and BootScan
    may seed internally), and BENCHMARKING.md section 4's claim that the
    synthetic results are reproducible is wrong for the OpenRDP path. The
    message therefore reports exactly what differed.
    """
    a, _ = _run("clear_recombinant", tmp_path / "a", cfg)
    b, _ = _run("clear_recombinant", tmp_path / "b", cfg)
    summary = (f"run A: status={a['status']} score={a['score']} "
               f"events={len(a['events'])}\n"
               f"run B: status={b['status']} score={b['score']} "
               f"events={len(b['events'])}")
    assert a["status"] == b["status"], f"status differs between identical runs\n{summary}"
    assert a["score"] == b["score"], f"score differs between identical runs\n{summary}"
    assert len(a["events"]) == len(b["events"]), f"event count differs\n{summary}"


@needs_openrdp
@slow
def test_multi_query_reports_each_query(tmp_path, cfg):
    p = tmp_path / "panel.fasta"
    synthetic.clear_recombinant(str(p))
    aln = R.read_alignment(str(p))
    queries = [aln.names[-1], aln.names[0]]
    results, _ = R.analyse(str(p), queries, cfg, workdir=str(tmp_path / "w"))
    assert [r["sequence_id"] for r in results] == queries


def test_missing_openrdp_is_reported_not_fatal(tmp_path, cfg, monkeypatch):
    """Without OpenRDP the command still runs and says which methods it skipped,
    rather than failing."""
    monkeypatch.setenv("CALICITYPER_OPENRDP", str(tmp_path / "nope"))
    cfg2 = RC.load(None)
    p = tmp_path / "x.fasta"
    synthetic.clear_recombinant(str(p))
    aln = R.read_alignment(str(p))
    results, provenance = R.analyse(str(p), [aln.names[-1]], cfg2,
                                    workdir=str(tmp_path / "w"))
    assert results[0]["recombination"]["status"] in ("none", "weak")
    skipped = {t["tool"] for t in provenance["tools_skipped"]}
    assert "openrdp" in skipped, provenance["tools_skipped"]


def test_unknown_query_id_is_an_error(tmp_path, cfg):
    p = tmp_path / "x.fasta"
    synthetic.non_recombinant(str(p))
    with pytest.raises(R.AlignmentError):
        R.analyse(str(p), ["not_in_alignment"], cfg, workdir=str(tmp_path / "w"))


# ------------------------------------------------------------------ sapovirus

SAPO_JUNCTION = synthetic.SAPOVIRUS_RDRP_VP1_JUNCTION


@pytest.fixture
def sapo_cfg():
    cfg = RC.load(None)
    cfg.virus = "sapovirus"
    cfg.regions = {}                 # re-resolved from `virus` by validate()
    cfg.validate()
    return cfg


def test_sapovirus_model_is_selected_and_differs_from_norovirus(sapo_cfg, cfg):
    assert sapo_cfg.virus == "sapovirus"
    assert any("RdRp-VP1" in r for r in sapo_cfg.regions)
    assert not any("RdRp-VP1" in r for r in cfg.regions)
    assert sapo_cfg.regions != cfg.regions


def test_sapovirus_junction_is_labelled_as_such_not_as_orf1_orf2(sapo_cfg, cfg):
    """The point of a per-virus model. Both viruses score nt 5139 as a hotspot,
    but in sapovirus VP1 sits inside ORF1, so calling it an 'ORF1-ORF2
    junction' would be wrong even though the coordinate is similar."""
    s_score, s_where = R.biological_plausibility_score(
        SAPO_JUNCTION, SAPO_JUNCTION, synthetic.SAPOVIRUS_LENGTH, sapo_cfg)
    n_score, n_where = R.biological_plausibility_score(
        SAPO_JUNCTION, SAPO_JUNCTION, synthetic.SAPOVIRUS_LENGTH, cfg)
    assert s_score == 1.0 and "RdRp-VP1" in s_where
    assert n_score == 1.0 and "ORF1-ORF2" in n_where          # same score, wrong name
    assert s_where != n_where


def test_sapovirus_mid_ns_breakpoint_is_not_rewarded_as_a_junction(sapo_cfg):
    hot, _ = R.biological_plausibility_score(
        SAPO_JUNCTION, SAPO_JUNCTION, synthetic.SAPOVIRUS_LENGTH, sapo_cfg)
    mid, where = R.biological_plausibility_score(
        synthetic.SAPOVIRUS_MID_NS, synthetic.SAPOVIRUS_MID_NS,
        synthetic.SAPOVIRUS_LENGTH, sapo_cfg)
    assert mid < hot
    assert "nonstructural" in where.lower(), where


def test_sapovirus_vp1_region_lies_inside_orf1(sapo_cfg):
    """A coordinate in VP1 must report a VP1 region that is part of ORF1, not
    a separate ORF2 the way norovirus has."""
    _, where = R.biological_plausibility_score(
        6000, 6000, synthetic.SAPOVIRUS_LENGTH, sapo_cfg)
    assert "VP1" in where and "ORF1" in where, where


@needs_openrdp
@slow
@pytest.mark.parametrize("name,expect_detected", [
    ("sapovirus_non_recombinant", False),
    ("sapovirus_clear_recombinant", True),
    ("sapovirus_mid_ns_recombinant", True),
])
def test_sapovirus_synthetic_end_to_end(name, expect_detected, tmp_path, sapo_cfg):
    p = tmp_path / f"{name}.fasta"
    _, truth = getattr(synthetic, name)(str(p))
    aln = R.read_alignment(str(p))
    results, _ = R.analyse(str(p), [aln.names[-1]], sapo_cfg,
                           workdir=str(tmp_path / "work"))
    rec = results[0]["recombination"]
    if not expect_detected:
        assert rec["status"] in ("none", "weak"), rec
        return
    assert rec["status"] in ("possible", "likely", "strong"), rec
    spans = [(b["query_start"], b["query_end"]) for e in rec["events"] for b in e["breakpoints"]]
    assert spans, rec["events"]
    bp = truth["breakpoint"]
    assert any(lo - sapo_cfg.max_breakpoint_distance <= bp
               <= hi + sapo_cfg.max_breakpoint_distance for lo, hi in spans), (spans, bp)


@needs_openrdp
@slow
def test_sapovirus_junction_event_carries_the_sapovirus_region_name(tmp_path, sapo_cfg):
    p = tmp_path / "s.fasta"
    synthetic.sapovirus_clear_recombinant(str(p))
    aln = R.read_alignment(str(p))
    results, provenance = R.analyse(str(p), [aln.names[-1]], sapo_cfg,
                                    workdir=str(tmp_path / "w"))
    regions = [e.get("biological_region", "") for e in results[0]["recombination"]["events"]]
    assert any("RdRp-VP1" in r for r in regions), regions
    assert provenance["parameters"]["virus"] == "sapovirus"
