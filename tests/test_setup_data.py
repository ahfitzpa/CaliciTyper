import json
import tarfile

from calicityper import setup_data


def _sha256_16_bytes(data: bytes) -> str:
    import hashlib
    return hashlib.sha256(data).hexdigest()[:16]


def _write_fixture(tmp_path):
    # Minimal layout matching manifest paths under calicityper/data/...
    content = b">ref1\nACGT\n"
    rel = "calicityper/data/database/refseq_2.fasta"
    staged = tmp_path / "staged" / "database"
    staged.mkdir(parents=True)
    (staged / "refseq_2.fasta").write_bytes(content)

    archive = tmp_path / "calicityper-reference-fixture.tar.gz"
    # include a top-level wrapper dir (common for release tarballs)
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(tmp_path / "staged", arcname="reference-fixture")

    manifest = {
        "reference_set_version": "test-fixture",
        "doi": "10.5281/zenodo.example",
        "download_url": "https://example.invalid/fixture.tar.gz",
        "sha256_16": {rel: _sha256_16_bytes(content)},
    }
    mpath = tmp_path / "reference_manifest.json"
    mpath.write_text(json.dumps(manifest))
    return archive, manifest, mpath


def test_install_from_local_source_and_verify_only(tmp_path, capsys):
    archive, manifest, _ = _write_fixture(tmp_path)
    data_dir = tmp_path / "data"

    rc = setup_data.install(str(archive), data_dir, manifest=manifest, force=True)
    assert rc == 0
    assert (data_dir / "database" / "refseq_2.fasta").exists()

    ok, problems = setup_data.verify(data_dir, manifest=manifest)
    assert ok, problems
    assert problems == []

    # already present and verified -> no-op success
    rc2 = setup_data.install(str(archive), data_dir, manifest=manifest, force=False)
    assert rc2 == 0

    # verify-only path via main-like behavior
    ok2, problems2 = setup_data.verify(data_dir, manifest=manifest)
    assert ok2 and problems2 == []

    out = capsys.readouterr().out
    assert "installed reference set test-fixture" in out


def test_install_refuses_checksum_mismatch(tmp_path, capsys):
    archive, manifest, _ = _write_fixture(tmp_path)
    bad = dict(manifest)
    bad["sha256_16"] = dict(manifest["sha256_16"])
    rel = next(iter(bad["sha256_16"]))
    bad["sha256_16"][rel] = "0000000000000000"

    data_dir = tmp_path / "data-bad"
    rc = setup_data.install(str(archive), data_dir, manifest=bad, force=True)
    assert rc == 1
    assert not (data_dir / "database" / "refseq_2.fasta").exists()
    err = capsys.readouterr().err
    assert "does not match the manifest" in err
