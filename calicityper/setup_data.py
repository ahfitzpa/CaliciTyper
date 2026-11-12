"""Fetch and verify the reference data archive.

The reference set is distributed as a Zenodo dataset rather than inside the
wheel. Zenodo is the right home for it because it ships the BUILT artefacts --
the alignments, .orf files and BLAST database -- not just the accessions.
Typing compares a query against reference alignment *columns*, and all three
typing windows are column ranges, so a rebuilt alignment that differs from the
original silently changes calls. Downloading the exact built set avoids that,
and a Zenodo version DOI is immutable, so a run can state which reference set
produced it.

Every file is checked against `reference_manifest.json`, and a mismatch is an
error rather than a warning: carrying on with a reference set that is not the
one the manifest describes is the failure this module exists to prevent.

Note on licensing: putting data on Zenodo is still redistribution. This path is
appropriate for the 309 public GenBank records the current set is built from.
Sequences from a restricted source (CaliciNet, or any surveillance database
with data-sharing terms) must not be deposited either here or on GitHub --
reference such a set by accession and have users build it locally.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

from . import config as cg

MANIFEST = Path(__file__).resolve().parent / "data" / "reference_manifest.json"
# Members the typing path reads. Anything outside these is not extracted, so a
# tampered archive cannot drop files elsewhere on the filesystem.
ALLOWED_TOP = {"database", "reference_sequences_RdRp", "reference_sequences_VP1",
               "putative_novel_bounds.tsv", "reference_accessions.tsv",
               "reference_manifest.json"}


def load_manifest(path=None):
    p = Path(path) if path else MANIFEST
    if not p.exists():
        raise FileNotFoundError(f"reference manifest not found: {p}")
    return json.loads(p.read_text())


def sha256_16(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def verify(data_dir, manifest=None, quiet=False):
    """Check every artefact in the manifest. Returns (ok, [problems])."""
    m = manifest or load_manifest()
    sums = m.get("sha256_16", {})
    data_dir = Path(data_dir)
    problems = []
    for rel, want in sorted(sums.items()):
        # manifest paths are recorded relative to the package root
        p = data_dir / Path(rel).relative_to("calicityper/data") \
            if rel.startswith("calicityper/data") else data_dir / rel
        if not p.exists():
            problems.append(f"missing: {p}")
            continue
        got = sha256_16(p)
        if got != want:
            problems.append(f"checksum differs: {p} (expected {want}, got {got})")
    if not quiet:
        n = len(sums)
        print(f"verified {n - len(problems)}/{n} artefacts in {data_dir}")
    return (not problems), problems


def _safe_members(tar):
    for member in tar.getmembers():
        name = Path(member.name)
        parts = [p for p in name.parts if p not in (".", "")]
        if not parts:
            continue
        # strip a single leading directory, as release tarballs usually carry one
        if parts[0] not in ALLOWED_TOP and len(parts) > 1:
            parts = parts[1:]
            member.name = str(Path(*parts))
        if parts[0] not in ALLOWED_TOP:
            continue
        if member.islnk() or member.issym() or ".." in parts:
            continue                      # never follow links or escape the target
        yield member


def install(source, data_dir, manifest=None, force=False):
    """Download or copy `source` (URL or local tar) into `data_dir`, verified.

    MERGED into `data_dir`, not replacing it: the package data directory always
    already holds the manifest, the accession list and the bounds table, so
    "the directory is not empty" says nothing about whether the reference set
    is present. Treating it as a reason to refuse made `calicityper setup` exit
    1 on every fresh install, which is what broke CI.
    """
    data_dir = Path(data_dir)
    m = manifest or load_manifest()
    if data_dir.exists() and not force:
        ok, _ = verify(data_dir, m, quiet=True)
        if ok:
            print(f"reference data already present and verified in {data_dir}")
            return 0

    with tempfile.TemporaryDirectory(prefix="calicityper_db_") as td:
        archive = Path(td) / "reference.tar.gz"
        if "://" in str(source):
            print(f"downloading {source}")
            with urllib.request.urlopen(source, timeout=300) as r, open(archive, "wb") as fh:
                shutil.copyfileobj(r, fh)
        else:
            src = Path(source).expanduser()
            if not src.exists():
                print(f"no such archive: {src}", file=sys.stderr)
                return 1
            shutil.copy(src, archive)
        staged = Path(td) / "staged"
        staged.mkdir()
        with tarfile.open(archive, "r:*") as tar:
            tar.extractall(staged, members=_safe_members(tar))
        ok, problems = verify(staged, m, quiet=True)
        if not ok:
            print("downloaded reference data does not match the manifest:", file=sys.stderr)
            for p in problems[:10]:
                print(f"  {p}", file=sys.stderr)
            print("refusing to install: a different reference set would change calls "
                  "silently.", file=sys.stderr)
            return 1
        # merge, not replace: data_dir legitimately holds the manifest, the
        # accession list and the bounds table, none of which are in the archive
        data_dir.mkdir(parents=True, exist_ok=True)
        for item in staged.iterdir():
            target = data_dir / item.name
            if target.exists():
                shutil.rmtree(target) if target.is_dir() else target.unlink()
            shutil.move(str(item), str(target))
    print(f"installed reference set {m.get('reference_set_version', '?')} into {data_dir}")
    print(f"set CALICITYPER_DATA={data_dir} (or leave it, if that is the default path)")
    return 0


def default_data_dir(version=None):
    """Where a fetched set lives when the package data directory is not used.
    Versioned, so two reference sets can coexist and a run can name its own."""
    v = version or load_manifest().get("reference_set_version", "unknown")
    base = os.environ.get("XDG_CACHE_HOME") or (Path.home() / ".cache")
    return Path(base) / "calicityper" / f"reference-{v}"


def main(args):
    m = load_manifest(args.manifest)
    data_dir = Path(args.data_dir) if args.data_dir else (
        cg.DATA if args.in_package else default_data_dir())
    if args.verify_only:
        ok, problems = verify(data_dir, m)
        for p in problems[:20]:
            print(f"  {p}")
        return 0 if ok else 1
    source = args.source or m.get("download_url")
    if not source:
        print("no source given and the manifest has no download_url.\n"
              "Pass --source <url-or-tarball>, or add the Zenodo URL to "
              f"{args.manifest or MANIFEST}.", file=sys.stderr)
        return 1
    return install(source, data_dir, m, force=args.force)
