"""Sapovirus genotyping: distance cut-off calls (Zhao et al. 2024) and routing."""
import csv
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import AlignIO, SeqIO

from calicityper import config as cg
from calicityper import typing as ty

ROOT = Path(__file__).resolve().parent.parent
G = cg.SAPOVIRUS_PREFIX + "GII"

needs_tools = pytest.mark.skipif(
    not (shutil.which("blastn") and shutil.which("mafft")),
    reason="requires blastn and mafft on PATH",
)


def _row():
    return dict(query="q", genogroup=G, note="")


def _res(identity, nearest_type="GII.2", sites=1600, span=1650, best_type=None, supports=None):
    return dict(nearest="X", nearest_identity=identity, nearest_type=nearest_type, n_sites=sites,
                region_sites=span, best_type=best_type or nearest_type,
                supports=supports or {nearest_type: 90.0},
                all_supports=supports or {nearest_type: 90.0}, support=90.0)


def test_within_the_genotype_cutoff_is_called_and_accepted():
    row = _row()
    ty.apply_result(row, "orf2", _res(85.0), min_support=70.0)       # d = 0.150 < 0.161
    assert row["orf2_type"] == "GII.2"
    assert row["orf2_distance"] == 0.15
    assert row["orf2_accepted"] == "accepted"
    assert row["orf2_support"] == 90.0


def test_over_the_genotype_cutoff_is_untypeable_not_the_nearest_label():
    """The decided rule: outside every genotype cluster means untypeable."""
    row = _row()
    ty.apply_result(row, "orf2", _res(80.0), min_support=70.0)       # d = 0.200
    assert row["orf2_type"] == "untypeable"
    assert row["orf2_nearest_type"] == "GII.2"                       # still reported
    assert "possible novel genotype" in row["orf2_accepted"]


def test_over_the_genogroup_cutoff_says_the_genogroup_is_not_confirmed():
    row = _row()
    ty.apply_result(row, "orf2", _res(45.0), min_support=70.0)       # d = 0.550 > 0.503
    assert row["orf2_type"] == "untypeable"
    assert "genogroup is not confirmed" in row["orf2_accepted"]


def test_rdrp_uses_its_own_cutoff():
    row = _row()
    ty.apply_result(row, "orf1", _res(75.0, nearest_type="GII.P2"), min_support=70.0)  # 0.25 < 0.266
    assert row["orf1_type"] == "GII.P2" and row["orf1_accepted"] == "accepted"


def test_bootstrap_support_does_not_decide_the_call():
    """Distance decides; low bootstrap support is reported, not used to reject."""
    row = _row()
    ty.apply_result(row, "orf2", _res(90.0, supports={"GII.2": 20.0, "GII.1": 80.0},
                                      best_type="GII.1"), min_support=70.0)
    assert row["orf2_type"] == "GII.2"
    assert row["orf2_support"] == 20.0
    assert row["orf2_accepted"] == "accepted"
    assert "bootstrap favours GII.1" in row["note"]


def test_a_partial_region_is_typed_but_not_accepted():
    """The cut-offs were derived on complete regions (provisional coverage floor)."""
    row = _row()
    ty.apply_result(row, "orf2", _res(95.0, sites=600, span=1650), min_support=70.0)
    assert row["orf2_type"] == "GII.2"
    assert row["orf2_accepted"].startswith("rejected") and "complete regions" in row["orf2_accepted"]


def test_sapovirus_is_genogroup_only_unless_enabled(monkeypatch):
    hit = dict(genogroup=G, subject="x", aln_len=5000, sstart=1, send=5000, qstart=1, qend=5000,
               reverse=False)
    monkeypatch.setattr(cg, "TYPE_SAPOVIRUS", False)
    row, jobs = ty.plan_query("q", "ACGT" * 500, hit)
    assert jobs == [] and "genogroup only" in row["note"]


def test_a_genogroup_without_an_rdrp_set_says_so_instead_of_borrowing():
    assert cg.ref_paths("orf1", cg.SAPOVIRUS_PREFIX + "GIV") is None
    assert cg.ref_paths("orf2", cg.SAPOVIRUS_PREFIX + "GIV") is not None


@pytest.mark.parametrize("region", ["VP1", "RdRp"])
def test_alignment_names_match_the_labels(region):
    """MAFFT's CLUSTAL writer truncated names to 15 characters, so no aligned
    reference matched its label and every set loaded empty."""
    d = ROOT / "calicityper" / "data" / f"reference_sequences_sapovirus_{region}"
    if not d.exists():
        pytest.skip("sapovirus reference sets not installed")
    for aln in sorted(d.glob("*.aln")):
        names = {r.id.split("|")[-1] for r in AlignIO.read(aln, "clustal")}
        labels = {r.id.split("|")[2] for r in SeqIO.parse(aln.with_suffix(".orf"), "fasta")}
        assert names == labels, (aln.name, sorted(names ^ labels)[:5])


@needs_tools
def test_a_sapovirus_genome_is_typed_end_to_end(tmp_path):
    """MN509084 (RIVM: sapovirus GII.2), from examples/. It is itself a
    reference strain, so this checks the path, not accuracy (see the LOO)."""
    rec = next(r for r in SeqIO.parse(ROOT / "examples" / "example.fasta", "fasta")
               if r.id.startswith("MN509084"))
    q = tmp_path / "q.fa"
    SeqIO.write([rec], q, "fasta")
    out = tmp_path / "out.tsv"
    p = subprocess.run([sys.executable, "-m", "calicityper.cli", "type", str(q), "-o", str(out),
                        "--type-sapovirus"], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    (row,) = list(csv.DictReader(open(out), delimiter="\t"))
    assert row["genogroup"] == G
    assert (row["orf1_type"], row["orf2_type"]) == ("GII.P2", "GII.2")
    assert row["orf1_accepted"] == row["orf2_accepted"] == "accepted"
    # novelty is the distance cut-off's job; the norovirus 2x SD screen must not run
    from calicityper import putative_novel
    assert row["orf1_putative_novel"] == row["orf2_putative_novel"] == putative_novel.NA
