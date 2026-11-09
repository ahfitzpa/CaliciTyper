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

COLUMNS = ["query", "genogroup", "orf1_type", "orf1_support", "orf2_type", "orf2_support", "note"]

needs_tools = pytest.mark.skipif(
    not (shutil.which("blastn") and shutil.which("mafft")),
    reason="requires blastn and mafft on PATH",
)


def run_cli(fasta, out, *extra):
    """Invoke the CLI as a subprocess, the way a user does."""
    cmd = [sys.executable, "-m", "norotyper.cli", str(fasta), "-o", str(out), *extra]
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
