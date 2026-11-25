"""Regression tests for Milestone 1: numeric CLI validation, state isolation, URL safety, and atomic writes."""
from __future__ import annotations

from pathlib import Path
import pytest

from calicityper import antigenic as ag, cli, config as cg, putative_novel, setup_data


@pytest.mark.parametrize(
    "argv,expected",
    [
        (["type", "dummy.fa", "--threads", "0"], "--threads"),
        (["type", "dummy.fa", "--threads", "-2"], "--threads"),
        (["type", "dummy.fa", "--bootstraps", "0"], "--bootstraps"),
        (["type", "dummy.fa", "--min-support", "150"], "--min-support"),
        (["type", "dummy.fa", "--min-support", "-1"], "--min-support"),
        (["type", "dummy.fa", "--min-sites", "-5"], "--min-sites"),
        (["type", "dummy.fa", "--max-refs-per-type", "0"], "--max-refs-per-type"),
        (["recomb", "--alignment", "a.fa", "--openrdp-timeout", "0"], "--openrdp-timeout"),
        (["recomb", "--alignment", "a.fa", "--pvalue-threshold", "1.5"], "--pvalue-threshold"),
        (["recomb", "--alignment", "a.fa", "--min-length", "0"], "--min-length"),
        (["antigenic", "a.faa", "--threads", "0"], "--threads"),
    ],
)
def test_invalid_numeric_cli_args_exit_2(argv: list[str], expected: str, capsys):
    rc = cli.main(argv)
    err = capsys.readouterr().err
    assert rc == 2
    assert expected in err


def test_runtime_overrides_restores_module_globals(tmp_path: Path):
    fa = tmp_path / "empty.fa"
    fa.write_text("", encoding="utf-8")
    out = tmp_path / "out.tsv"

    orig_min_sites = cg.MIN_INFORMATIVE_SITES
    orig_orf1 = cg.MIN_INFORMATIVE_SITES_ORF1
    orig_sapo = cg.TYPE_SAPOVIRUS
    orig_dist = putative_novel.DISTANCE
    orig_fast_novel = putative_novel.FAST_NOVEL

    rc = cli.main([
        "type", str(fa), "-o", str(out),
        "--preset", "polcap",
        "--min-sites", "123",
        "--no-type-sapovirus",
        "--novel-distance", "p",
        "--fast-novel",
    ])
    assert rc == 0
    assert cg.MIN_INFORMATIVE_SITES == orig_min_sites
    assert cg.MIN_INFORMATIVE_SITES_ORF1 == orig_orf1
    assert cg.TYPE_SAPOVIRUS == orig_sapo
    assert putative_novel.DISTANCE == orig_dist
    assert putative_novel.FAST_NOVEL == orig_fast_novel


@pytest.mark.parametrize("bad_url", ["file:///etc/passwd", "ftp://example.com/ref.tar.gz"])
def test_setup_data_rejects_unsupported_url_schemes(bad_url: str, tmp_path: Path, capsys):
    rc = setup_data.install(data_dir=tmp_path / "ref", source=bad_url, force=True)
    err = capsys.readouterr().err
    assert rc == 1
    assert "unsupported URL scheme" in err


def test_write_pdimer_fasta_is_atomic(tmp_path: Path, monkeypatch):
    target = tmp_path / "pdimers.fa"
    target.write_text(">old\nAAA:AAA\n", encoding="utf-8")

    def fail_replace(_src, _dst):
        raise OSError("simulated failure")

    monkeypatch.setattr("calicityper.provenance.os.replace", fail_replace)
    with pytest.raises(OSError, match="simulated failure"):
        ag.write_pdimer_fasta(target, [("q1", "ACDEF")])

    assert target.read_text(encoding="utf-8") == ">old\nAAA:AAA\n"
    assert [p.name for p in tmp_path.iterdir()] == ["pdimers.fa"]
