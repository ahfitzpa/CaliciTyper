#!/usr/bin/env python3
"""Rebuild calicityper/data/reference_manifest.json from the files on disk.

`calicityper setup` refuses to install a reference set whose checksums do not
match this manifest, which is the point of it -- a reference set that is not
the one the manifest describes changes calls silently. So whenever the
reference set changes, this has to be rerun, and the version bumped: an old
cached copy must not be mistaken for the new one.

Metadata already in the manifest (download_url, doi, citation) is preserved;
only the checksums and, with --version, the version are rewritten.

    python scripts/build_reference_manifest.py [--version 0.7.0] [--apply]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from calicityper import config as cg, setup_data as sd          # noqa: E402

# Artefacts the typing path reads. The manifest covers what a user must have
# byte-identical, not every file in the directory: .dnd guide trees are legacy
# and never read, and the small tables ship inside the wheel instead.
SKIP_SUFFIX = {".dnd", ".bak", ".npy"}
SKIP_NAME = {"reference_manifest.json", "reference_accessions.tsv",
             "putative_novel_bounds.tsv", "sapovirus_excluded.tsv"}
SKIP_DIR = {".cache"}


def artefacts(data_dir):
    for p in sorted(data_dir.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(data_dir)
        if rel.parts[0] in SKIP_DIR or p.suffix in SKIP_SUFFIX or p.name in SKIP_NAME:
            continue
        yield rel, p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", help="new reference_set_version (e.g. 0.7.0)")
    ap.add_argument("--download-url", metavar="URL",
                    help="set download_url to the released asset (what "
                         "`calicityper setup` fetches); checked with a HEAD "
                         "request unless --no-check-url")
    ap.add_argument("--no-check-url", action="store_true",
                    help="skip the reachability check on --download-url")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()

    path = cg.DATA / "reference_manifest.json"
    m = json.loads(path.read_text())
    old = dict(m.get("sha256_16", {}))
    sums = {f"calicityper/data/{rel.as_posix()}": sd.sha256_16(p) for rel, p in artefacts(cg.DATA)}

    added = sorted(set(sums) - set(old))
    removed = sorted(set(old) - set(sums))
    changed = sorted(k for k in set(sums) & set(old) if sums[k] != old[k])
    print(f"{len(sums)} artefacts ({len(old)} before)")
    for label, keys in (("added", added), ("removed", removed), ("checksum changed", changed)):
        for k in keys:
            print(f"  {label:17} {k}")
    if not (added or removed or changed) and not a.version and not a.download_url:
        print("nothing to do")
        return 0

    m["sha256_16"] = dict(sorted(sums.items()))
    # n_artefacts is quoted in the docs and the setup output, so keep it in step
    # with the checksum table rather than letting it drift as it did at 0.7.0.
    if m.get("n_artefacts") != len(sums):
        print(f"  n_artefacts {m.get('n_artefacts')} -> {len(sums)}")
    m["n_artefacts"] = len(sums)
    if a.version:
        print(f"version {m.get('reference_set_version')} -> {a.version}")
        m["reference_set_version"] = a.version
    if a.download_url:
        if not a.no_check_url:
            import urllib.request, urllib.error
            req = urllib.request.Request(a.download_url, method="HEAD")
            try:
                with urllib.request.urlopen(req, timeout=30) as r:
                    size = r.headers.get("Content-Length")
                    print(f"  URL reachable: HTTP {r.status}"
                          + (f", {int(size)/1e6:.1f} MB" if size else ""))
            except Exception as e:
                print(f"  WARNING: could not reach the URL ({e}).")
                print("  Set it anyway with --no-check-url if you are sure, but a "
                      "download_url that 404s makes `calicityper setup` fail for "
                      "everyone who is not installing from a local archive.")
                return 1
        print(f"download_url:\n  old {m.get('download_url')}\n  new {a.download_url}")
        m["download_url"] = a.download_url
    elif m.get("download_url") and a.version:
        print(f"  NOTE: download_url still points at the previous deposit:\n"
              f"    {m['download_url']}\n  Pass --download-url once the asset is uploaded.")
    if not a.apply:
        print("\nreport only; pass --apply to write")
        return 0
    path.write_text(json.dumps(m, indent=2) + "\n")
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
