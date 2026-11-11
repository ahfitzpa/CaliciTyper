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

# dataset -> (statuses allowed, expected breakpoint or None)
EXPECTED = {
    "non_recombinant": ({"none", "weak"}, None),
    "clear_recombinant": ({"possible", "likely", "strong"}, ORF1_ORF2_JUNCTION),
    "weak_signal": ({"none", "weak", "possible"}, None),
    "conflicting": ({"weak", "possible", "likely", "strong"}, None),
    "orf_mismatch_only": ({"none", "weak", "possible", "likely"}, None),
}


def _run(name, tmp_path, cfg):
    """Returns (recombination block, provenance). A result is
    {sequence_id, qc, recombination{status, score, confidence, events, ...}}."""
    p = tmp_path / f"{name}.fasta"
    getattr(synthetic, name)(str(p))
    aln = R.read_alignment(str(p))
    query = aln.names[-1]
    results, provenance = R.analyse(str(p), [query], cfg, workdir=str(tmp_path / "work"))
    return results[0]["recombination"], provenance


@needs_openrdp
@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_synthetic_dataset_status(name, tmp_path, cfg):
    allowed, _ = EXPECTED[name]
    res, _ = _run(name, tmp_path, cfg)
    assert res["status"] in allowed, f"{name}: got {res['status']}, events={res.get('events')}"


@needs_openrdp
def test_clear_recombinant_breakpoint_is_localised(tmp_path, cfg):
    """The planted breakpoint must be recovered, and within
    max_breakpoint_distance of the truth -- localisation, not just detection."""
    res, _ = _run("clear_recombinant", tmp_path, cfg)
    events = res.get("events") or []
    assert events, res
    starts = [e["breakpoint"]["query_start"] for e in events
              if e.get("breakpoint", {}).get("query_start")]
    assert starts, events
    assert min(abs(s - ORF1_ORF2_JUNCTION) for s in starts) <= cfg.max_breakpoint_distance


@needs_openrdp
def test_orf_mismatch_only_never_reaches_strong(tmp_path, cfg):
    """Distinct ORF1/ORF2 ancestry with near-identical donors has low parental
    distinctness, which must cap the status."""
    res, _ = _run("orf_mismatch_only", tmp_path, cfg)
    assert res["status"] != "strong", res


@needs_openrdp
def test_same_seed_gives_the_same_result(tmp_path, cfg):
    """Seeded end to end, so a failure is reproducible and a number quoted in
    the docs can be re-derived."""
    a, _ = _run("clear_recombinant", tmp_path / "a", cfg)
    b, _ = _run("clear_recombinant", tmp_path / "b", cfg)
    assert a["status"] == b["status"]
    assert a["score"] == b["score"]


@needs_openrdp
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
