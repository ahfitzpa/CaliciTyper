#!/usr/bin/env python3
"""Build the reference-set archive that `calicityper setup` installs.

The 0.7.0 archive was tarred by hand, so there was no way to reproduce it or
to check that it held what the manifest said it held. This builds it from the
manifest itself, so the two cannot disagree.

    python scripts/build_reference_archive.py                 # dry run
    python scripts/build_reference_archive.py --out-dir dist --apply

Writes `calicityper-reference-<version>.tar.gz` plus a `.sha256` file. The
version comes from the manifest's `reference_set_version`, not from the
package version, because the two move independently.

The archive contains exactly the top-level entries the manifest checksums.
Every member is verified against the manifest before it goes in, and the
finished archive is read back and verified again, so a corrupt or stale file
cannot be published.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from calicityper import config as cg, setup_data as sd   # noqa: E402


def manifest_tops(manifest):
    """Top-level entries under calicityper/data/ that the manifest covers."""
    tops = set()
    for key in manifest["sha256_16"]:
        rel = key.split("calicityper/data/", 1)[1]
        tops.add(rel.split("/")[0])
    return sorted(tops)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out-dir", default="dist")
    ap.add_argument("--apply", action="store_true", help="write the archive")
    a = ap.parse_args()

    manifest = json.loads((cg.DATA / "reference_manifest.json").read_text())
    version = manifest["reference_set_version"]
    tops = manifest_tops(manifest)

    # Everything the manifest lists must be installable, or the archive is
    # useless on arrival -- this is the failure that hid the sapovirus sets.
    not_installable = sorted(set(tops) - sd.ALLOWED_TOP)
    if not_installable:
        raise SystemExit(
            "these are in the manifest but setup_data.ALLOWED_TOP would discard "
            f"them on install: {not_installable}")

    print(f"reference set version {version}")
    print(f"{len(manifest['sha256_16'])} artefacts over {len(tops)} top-level entries:")
    missing = []
    for t in tops:
        p = cg.DATA / t
        n = sum(1 for k in manifest["sha256_16"] if f"calicityper/data/{t}/" in k)
        print(f"  {t:42s} {n:4d} files" + ("" if p.exists() else "   MISSING"))
        if not p.exists():
            missing.append(t)
    if missing:
        raise SystemExit(f"not on disk: {missing}")

    # verify on disk before packing
    ok, problems = sd.verify(cg.DATA, manifest, quiet=True)
    if not ok:
        print(f"\n{len(problems)} artefact(s) do not match the manifest:")
        for p in problems[:10]:
            print(f"  {p}")
        raise SystemExit("refusing to build an archive from a set that does not "
                         "match its own manifest; rebuild the manifest first")
    print("\nall artefacts match the manifest")

    out_dir = Path(a.out_dir)
    archive = out_dir / f"calicityper-reference-{version}.tar.gz"
    if not a.apply:
        print(f"\ndry run; would write {archive}\npass --apply to build")
        return 0

    out_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "w:gz") as tar:
        for t in tops:
            tar.add(cg.DATA / t, arcname=t)
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    (archive.parent / (archive.name + ".sha256")).write_text(
        f"{digest}  {archive.name}\n")
    print(f"wrote {archive} ({archive.stat().st_size / 1e6:.1f} MB)")
    print(f"sha256 {digest}")

    # read it back and verify, so a bad archive cannot be published
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        with tarfile.open(archive) as tar:
            tar.extractall(td, members=sd._safe_members(tar))
        ok, problems = sd.verify(Path(td), manifest, quiet=True)
    if not ok:
        raise SystemExit(f"archive failed verification after writing: {problems[:5]}")
    print("archive re-verified after writing: every artefact matches the manifest")

    print("\nBefore release, set the manifest's download_url to the uploaded "
          "asset and rebuild the manifest if anything changed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
