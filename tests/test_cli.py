"""End-to-end tests for the calicityper command line interface.

These complement tests/test_fast_paths.py (which checks the numeric kernels)
and tests/loo_validation.py (accuracy). The point here is behaviour a user
actually meets: does the CLI run, does it emit the documented columns, and
does it degrade sensibly on input that cannot be typed?

Each test needs BLAST and MAFFT on PATH and is skipped if they are absent.
"""
import csv
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO

from calicityper import config as cg

COLUMNS = ["query", "genogroup", "orf1_type", "orf1_support", "orf1_runners",
           "orf2_type", "orf2_support", "orf2_runners", "note",
           "orf1_putative_novel", "orf1_putative_novel_detail", "orf2_putative_novel", "orf2_putative_novel_detail"]

needs_tools = pytest.mark.skipif(
    not (shutil.which("blastn") and shutil.which("mafft")),
    reason="requires blastn and mafft on PATH",
)


def run_cli(fasta, out, *extra):
    """Invoke the CLI as a subprocess, the way a user does."""
    cmd = [sys.executable, "-m", "calicityper.cli", "type", str(fasta), "-o", str(out), *extra]
    p = subprocess.run(cmd, capture_output=True, text=True)
    assert p.returncode == 0, f"calicityper failed:\n{p.stdout}\n{p.stderr}"
    with open(out) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def reference_genomes(n):
    return list(SeqIO.parse(cg.BLAST_FASTA, "fasta"))[:n]


def test_reference_data_is_installed():
    """The data must ship with the package, not just exist in a source checkout."""
    cg.check_data()
    assert cg.BLAST_FASTA.exists()
    assert (cg.REF_DIR["orf1"] / "GI.aln").exists()
    assert (cg.REF_DIR["orf2"] / "GI.aln").exists()


@needs_tools
def test_types_known_genomes(tmp_path):
    """Two GI reference genomes should be typed as GI with high support."""
    q = tmp_path / "q.fa"
    SeqIO.write(reference_genomes(2), q, "fasta")
    rows = run_cli(q, tmp_path / "out.tsv", "--threads", "2")

    assert len(rows) == 2
    assert list(rows[0]) == COLUMNS, "output columns changed; README documents these"
    for r in rows:
        assert r["genogroup"] == "GI"
        assert r["orf2_type"].startswith("GI."), r
        assert float(r["orf2_support"]) >= 70.0, r


@needs_tools
def test_non_norovirus_input_is_reported_not_crashed(tmp_path):
    """Random sequence has no BLAST hit: one row, unassigned, explanatory note."""
    q = tmp_path / "junk.fa"
    q.write_text(">junk\n" + "ACGT" * 200 + "\n")
    rows = run_cli(q, tmp_path / "out.tsv")

    assert len(rows) == 1
    assert rows[0]["genogroup"] == "unassigned"
    assert rows[0]["orf1_type"] == "" and rows[0]["orf2_type"] == ""
    assert "no BLAST hit" in rows[0]["note"]


@needs_tools
def test_short_fragment_gets_overlap_note(tmp_path):
    """A 60 nt fragment is below the 100 nt ORF-overlap minimum for typing."""
    g = reference_genomes(1)[0]
    q = tmp_path / "frag.fa"
    q.write_text(f">frag\n{str(g.seq)[1000:1060]}\n")
    rows = run_cli(q, tmp_path / "out.tsv")

    assert len(rows) == 1
    r = rows[0]
    # Either no hit at all, or a hit with too little ORF overlap to type.
    assert r["orf1_type"] == "" or r["orf2_type"] == "" or "overlap" in r["note"]


@needs_tools
def test_every_input_record_yields_exactly_one_row(tmp_path):
    """Row count and order must track the input FASTA, typable or not."""
    recs = reference_genomes(2)
    q = tmp_path / "mixed.fa"
    with open(q, "w") as fh:
        SeqIO.write(recs[0], fh, "fasta")
        fh.write(">junk\n" + "ACGT" * 150 + "\n")
        SeqIO.write(recs[1], fh, "fasta")
    rows = run_cli(q, tmp_path / "out.tsv", "--threads", "2")

    assert len(rows) == 3
    assert [r["query"] for r in rows] == [recs[0].id, "junk", recs[1].id]
    assert rows[1]["genogroup"] == "unassigned"


@needs_tools
def test_results_are_deterministic_for_a_fixed_seed(tmp_path):
    """The bootstrap is seeded, so two identical runs must agree exactly."""
    q = tmp_path / "q.fa"
    SeqIO.write(reference_genomes(2), q, "fasta")
    a = run_cli(q, tmp_path / "a.tsv", "--seed", "7", "--bootstraps", "20")
    b = run_cli(q, tmp_path / "b.tsv", "--seed", "7", "--bootstraps", "20")
    assert a == b


@needs_tools
def test_runner_up_types_are_reported(tmp_path):
    """The ranked runner-up types must reach the TSV, best first."""
    q = tmp_path / "q.fa"
    SeqIO.write(reference_genomes(2), q, "fasta")
    rows = run_cli(q, tmp_path / "out.tsv", "--threads", "2", "--bootstraps", "30")

    for r in rows:
        for region in ("orf1", "orf2"):
            if r[f"{region}_type"] in ("", "unassigned"):
                continue
            parsed = [p.split(":") for p in r[f"{region}_runners"].split("|")]
            assert 1 <= len(parsed) <= 3, r[f"{region}_runners"]
            vals = [float(v) for _, v in parsed]
            assert vals == sorted(vals, reverse=True), "runners must be ranked"
            # the winning type is the best-supported one
            assert parsed[0][0] == r[f"{region}_type"]
            assert abs(vals[0] - float(r[f"{region}_support"])) < 0.05


@needs_tools
def test_report_writes_tsv_md_json_and_html(tmp_path):
    """--report emits all three extra formats without disturbing the TSV."""
    import json

    q = tmp_path / "q.fa"
    with open(q, "w") as fh:
        SeqIO.write(reference_genomes(2), fh, "fasta")
        fh.write(">junk\n" + "ACGT" * 150 + "\n")
    out = tmp_path / "out.tsv"
    prefix = tmp_path / "rep"
    rows = run_cli(q, out, "--threads", "2", "--bootstraps", "30", "--report", str(prefix))

    assert len(rows) == 3
    for ext in ("md", "json", "html"):
        f = prefix.with_suffix("." + ext)
        assert f.exists() and f.stat().st_size > 0, f"{ext} report missing"

    doc = json.loads(prefix.with_suffix(".json").read_text())
    assert doc["tool"] == "calicityper"
    assert doc["n_queries"] == 3
    assert len(doc["results"]) == 3
    assert doc["meta"]["bootstraps"] == 30
    # runner-ups are structured in JSON, not a packed string
    typed = [r for r in doc["results"] if r["orf2_type"] not in ("", "unassigned")]
    assert typed, "expected at least one typed query"
    assert isinstance(typed[0]["orf2_runners"], list)
    assert "support" in typed[0]["orf2_runners"][0]

    page = prefix.with_suffix(".html").read_text()
    assert "calicityper report" in page
    assert "validated for" in page, "HTML must carry the GI/GII validation caveat"


def test_report_module_works_without_bokeh(tmp_path, monkeypatch):
    """Markdown and JSON must not depend on bokeh, and HTML must degrade."""
    from calicityper import report

    monkeypatch.setattr(report, "HAVE_BOKEH", False)
    rows = [dict(query="q1", genogroup="GII", orf1_type="GII.P31", orf1_support=97.0,
                 orf1_runners="GII.P31:97.0|GII.P16:2.0", orf2_type="GII.4_Sydney",
                 orf2_support=99.0, orf2_runners="GII.4_Sydney:99.0", note="")]
    assert report.alignment_figure(["a"], ["ACGT"], "t") is None
    report.to_markdown(rows, tmp_path / "r.md")
    report.to_json(rows, tmp_path / "r.json")
    report.to_html(rows, tmp_path / "r.html", figures=[])
    assert (tmp_path / "r.md").stat().st_size > 0
    assert "GII.4_Sydney" in (tmp_path / "r.html").read_text()


def test_closest_reference_ignores_uncovered_columns():
    """A partial query is compared only where it has bases: the reference that is
    identical over the covered span wins even if it differs elsewhere."""
    from calicityper import report

    q = "----ACGTACGT----"
    refs = ["TTTTACGTACGTTTTT",      # identical over the covered 8 columns
            "ACGTACGAACGTACGT"]      # one difference inside the covered span
    k, p, n, diff = report.closest_reference(q, refs)
    assert (k, p, n, diff) == (0, 0.0, 8, 0)
    assert report.closest_reference("-" * 16, refs) is None


def test_difference_figure_draws_only_differences():
    from calicityper import report
    if not report.HAVE_BOKEH:
        return
    fig = report.difference_figure("q", "--ACGTTCGT-A", "r", "AAACGTACGTCA", "t")
    cells = [r for r in fig.renderers if "rect" in type(r.glyph).__name__.lower()][0]
    # one substitution (col 7) and one query deletion (col 11), each drawn in both rows
    assert sorted(set(cells.data_source.data["col"])) == [7, 11]


def test_complete_vp1_requires_start_stop_and_length():
    from calicityper import putative_novel as nv
    lens = [540, 560]
    orf = "ATG" + "GCT" * 549 + "TAA"                      # 550 aa, complete
    assert len(nv.complete_vp1("CC" + orf + "GG", lens)) == 550
    assert nv.complete_vp1("CC" + orf[:-3] + "GCT", lens) is None      # no stop: truncated
    assert nv.complete_vp1("ATG" + "GCT" * 100 + "TAA", lens) is None  # far too short


def test_classify_follows_per_cluster_bounds():
    from calicityper import putative_novel as nv
    t = {("orf2", "genotype", "GII.4"): 0.10, ("orf2", "genotype", "*pooled*"): 0.13,
         ("orf2", "genogroup", "GII"): 0.44, ("orf2", "genogroup", "*pooled*"): 0.44}
    assert nv.classify("orf2", {"GII.4": 0.05}, None, t)[0] == nv.KNOWN
    assert nv.classify("orf2", {"GII.4": 0.20}, {"GII": 0.30}, t)[0] == nv.NOVEL_TYPE
    assert nv.classify("orf2", {"GII.4": 0.50}, {"GII": 0.50}, t)[0] == nv.NOVEL_GROUP
    # an unlisted (small/tentative) cluster falls back to the pooled bound
    assert nv.classify("orf2", {"GII.NA1": 0.12}, None, t)[0] == nv.KNOWN


def test_putative_novel_columns_and_sapovirus_not_assessed(tmp_path):
    """Every row carries the putative novel columns; non-noroviruses are not assessed."""
    import pandas as pd
    from calicityper import putative_novel as nv
    ex = Path(__file__).resolve().parent.parent / "examples" / "example.fasta"
    out = tmp_path / "t.tsv"
    run_cli(ex, out)
    d = pd.read_csv(out, sep="\t", dtype=str).fillna("").set_index("query")
    for c in ("orf1_putative_novel", "orf2_putative_novel", "orf1_putative_novel_detail", "orf2_putative_novel_detail"):
        assert c in d.columns
    assert d.loc["MN509084.1", "orf2_putative_novel"] == nv.NA
    assert d.loc["MT451991.1", "orf2_putative_novel"] == nv.KNOWN
