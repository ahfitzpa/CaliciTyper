"""End-to-end tests for the norotyper command line interface.

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

import pytest
from Bio import SeqIO

from norotyper import config as cg

COLUMNS = ["query", "genogroup", "orf1_type", "orf1_support", "orf1_runners",
           "orf2_type", "orf2_support", "orf2_runners", "note"]

needs_tools = pytest.mark.skipif(
    not (shutil.which("blastn") and shutil.which("mafft")),
    reason="requires blastn and mafft on PATH",
)


def run_cli(fasta, out, *extra):
    """Invoke the CLI as a subprocess, the way a user does."""
    cmd = [sys.executable, "-m", "norotyper.cli", "type", str(fasta), "-o", str(out), *extra]
    p = subprocess.run(cmd, capture_output=True, text=True)
    assert p.returncode == 0, f"norotyper failed:\n{p.stdout}\n{p.stderr}"
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
    assert doc["tool"] == "norotyper"
    assert doc["n_queries"] == 3
    assert len(doc["results"]) == 3
    assert doc["meta"]["bootstraps"] == 30
    # runner-ups are structured in JSON, not a packed string
    typed = [r for r in doc["results"] if r["orf2_type"] not in ("", "unassigned")]
    assert typed, "expected at least one typed query"
    assert isinstance(typed[0]["orf2_runners"], list)
    assert "support" in typed[0]["orf2_runners"][0]

    page = prefix.with_suffix(".html").read_text()
    assert "norotyper report" in page
    assert "validated for" in page, "HTML must carry the GI/GII validation caveat"


def test_report_module_works_without_bokeh(tmp_path, monkeypatch):
    """Markdown and JSON must not depend on bokeh, and HTML must degrade."""
    from norotyper import report

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
