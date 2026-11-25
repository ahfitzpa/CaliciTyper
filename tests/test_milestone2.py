"""Regression tests for Milestone 2: runtime provenance, antigenic provenance, and recomb validation/--dry-run."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from calicityper import cli, provenance


def test_runtime_environment_records_python_and_core_deps():
    rt = provenance.runtime_environment()
    assert isinstance(rt["python_version"], str) and "." in rt["python_version"]
    deps = rt["dependencies"]
    for pkg in ("numpy", "pandas", "biopython"):
        assert pkg in deps and deps[pkg]


def test_antigenic_provenance_records_input_sha256_and_runtime(tmp_path: Path):
    aln = tmp_path / "vp1.faa"
    ref_seq = "M" + "A" * 539
    qry_seq = "M" + "A" * 295 + "V" + "A" * 243
    payload = f">ref\n{ref_seq}\n>qry\n{qry_seq}\n"
    aln.write_text(payload, encoding="utf-8")
    out_tsv = tmp_path / "ag.tsv"

    rc = cli.main(["antigenic", str(aln), "-o", str(out_tsv), "--skip-input-check"])
    assert rc == 0

    prov_path = Path(str(out_tsv) + ".provenance.json")
    prov = json.loads(prov_path.read_text(encoding="utf-8"))
    assert prov["input_sha256"] == hashlib.sha256(payload.encode("utf-8")).hexdigest()
    assert "python_version" in prov["runtime"]
    assert "numpy" in prov["runtime"]["dependencies"]


def test_recomb_validation_errors_return_exit_2_without_creating_outdir(tmp_path: Path, capsys):
    aln = tmp_path / "bad.fa"
    aln.write_text(">a\nACGT\n>b\nACGTACGT\n", encoding="utf-8")
    outdir = tmp_path / "recomb_out"

    # Both --alignment and --input rejected
    rc = cli.main(["recomb", "--alignment", str(aln), "--input", str(aln), "-o", str(outdir)])
    assert rc == 2
    assert "not both" in capsys.readouterr().err

    # Neither --alignment nor --input rejected
    rc = cli.main(["recomb", "-o", str(outdir)])
    assert rc == 2
    assert "--alignment" in capsys.readouterr().err

    # Invalid method rejected cleanly (no ValueError traceback)
    rc = cli.main(["recomb", "--alignment", str(aln), "--methods", "bogus_method", "-o", str(outdir)])
    assert rc == 2
    assert "invalid recombination configuration" in capsys.readouterr().err

    # Unequal-length alignment rejected cleanly (no AlignmentError traceback)
    rc = cli.main(["recomb", "--alignment", str(aln), "--dry-run", "-o", str(outdir)])
    assert rc == 2
    assert not outdir.exists()
    assert "Error:" in capsys.readouterr().err


def test_recomb_dry_run_has_no_outdir_side_effects(tmp_path: Path, capsys):
    aln = tmp_path / "panel.fa"
    seq = "ACGT" * 50
    aln.write_text(f">r1\n{seq}\n>r2\n{seq}\n>q1\n{seq}\n", encoding="utf-8")
    outdir = tmp_path / "never_created"

    rc = cli.main(["recomb", "--alignment", str(aln), "--query", "q1", "--dry-run", "-o", str(outdir)])
    assert rc == 0
    assert not outdir.exists()
    assert "virus model" in capsys.readouterr().out


def test_recomb_malformed_config_exits_2(tmp_path: Path, capsys):
    bad_toml = tmp_path / "broken.toml"
    bad_toml.write_text("[recombination\nbroken = [", encoding="utf-8")
    rc = cli.main(["recomb", "--config", str(bad_toml), "--alignment", str(bad_toml)])
    assert rc == 2
    assert "malformed TOML" in capsys.readouterr().err

