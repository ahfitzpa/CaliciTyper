"""Regression tests for Milestone 3: structured logging, --quiet, antigenic stderr isolation, and report helpers."""
from __future__ import annotations

from pathlib import Path

from calicityper import cli, report


def test_type_quiet_suppresses_info_progress_logs(tmp_path: Path, capsys):
    fa = tmp_path / "empty.fa"
    fa.write_text("", encoding="utf-8")
    out = tmp_path / "out.tsv"

    rc = cli.main(["type", str(fa), "-o", str(out), "--quiet"])
    captured = capsys.readouterr()
    assert rc == 0
    assert captured.out == ""
    assert "[type] checking reference data" not in captured.err


def test_antigenic_diagnostics_go_to_stderr_not_stdout(tmp_path: Path, capsys):
    aln = tmp_path / "vp1.faa"
    ref_seq = "M" + "A" * 539
    qry_seq = "M" + "A" * 295 + "V" + "A" * 243
    aln.write_text(f">ref\n{ref_seq}\n>qry\n{qry_seq}\n", encoding="utf-8")
    out_tsv = tmp_path / "ag.tsv"

    rc = cli.main(["antigenic", str(aln), "-o", str(out_tsv), "--skip-input-check"])
    captured = capsys.readouterr()
    assert rc == 0
    assert captured.out == ""
    assert "annotated 1 sequence(s)" in captured.err

    # --show-sites intentionally still prints site definitions to stdout
    rc_sites = cli.main(["antigenic", "--show-sites"])
    captured_sites = capsys.readouterr()
    assert rc_sites == 0
    assert "A" in captured_sites.out


def test_report_builders_extracted_to_report_module():
    assert report.build_alignment_panels([], {}) == []
    assert report.build_tree_figures([], {}) == ([], [])
