"""Unit tests for atomic writes, provenance hashing/tool versions, and input validation."""
from __future__ import annotations

import gzip
import hashlib
from pathlib import Path
import pytest

from calicityper import cli, provenance, typing as ty


def test_atomic_write_text_replaces_atomically(tmp_path: Path):
    target = tmp_path / "out.tsv"
    target.write_text("old\n", encoding="utf-8")

    provenance.atomic_write_text(target, "new\n")
    assert target.read_text(encoding="utf-8") == "new\n"
    assert [p.name for p in tmp_path.iterdir()] == ["out.tsv"]


def test_atomic_write_text_preserves_target_on_failure(tmp_path: Path, monkeypatch):
    target = tmp_path / "out.tsv"
    target.write_text("original\n", encoding="utf-8")

    def boom(_src, _dst):
        raise OSError("disk full")

    monkeypatch.setattr(provenance.os, "replace", boom)
    with pytest.raises(OSError, match="disk full"):
        provenance.atomic_write_text(target, "partial\n")

    assert target.read_text(encoding="utf-8") == "original\n"
    assert [p.name for p in tmp_path.iterdir()] == ["out.tsv"]


def test_file_sha256_and_external_tool_versions(tmp_path: Path):
    fa = tmp_path / "q.fa"
    payload = b">q1\nACGTACGT\n"
    fa.write_bytes(payload)
    assert provenance.file_sha256(fa) == hashlib.sha256(payload).hexdigest()
    assert provenance.file_sha256(tmp_path / "missing.fa") == ""

    versions = provenance.external_tool_versions(("blastn", "mafft"))
    assert set(versions) == {"blastn", "mafft"}
    assert all(isinstance(v, str) and v for v in versions.values())


def test_read_fasta_rejects_gzip_with_clear_error(tmp_path: Path, capsys):
    gz = tmp_path / "reads.fa.gz"
    with gzip.open(gz, "wt", encoding="utf-8") as fh:
        fh.write(">s1\n" + "ACGT" * 60 + "\n")

    with pytest.raises(ValueError, match="gzip-compressed"):
        ty.read_fasta(gz)

    rc = cli.main(["type", str(gz), "-o", "-"])
    err = capsys.readouterr().err
    assert rc == 2
    assert "gzip-compressed" in err
