#!/usr/bin/env python3
"""Shared helpers for the F1-F4 development scripts. No tool behaviour changes."""
import os
import sys
from pathlib import Path

import numpy as np

# Resolve from this file rather than an assumed checkout layout: the earlier
# version assumed .../calicityper/dev/gii4_epitopes and broke in the HPC
# checkout, which sits at a different depth.
HERE = Path(__file__).resolve().parent.parent       # .../gii4_epitopes
ROOT = HERE.parent.parent                           # tool checkout (gii4_epitopes -> dev -> calicityper)


def _find_tool_checkout():
    """Locate a calicityper checkout so scripts run from a flat HPC folder.

    Order: already importable -> CALICITYPER_REPO env var -> parents of HERE
    containing calicityper/__init__.py. Returns None if not found; callers
    then get a clear error instead of ModuleNotFoundError.
    """
    try:
        import calicityper  # noqa: F401
        return Path(calicityper.__file__).resolve().parent.parent
    except ImportError:
        pass
    env = os.environ.get("CALICITYPER_REPO")
    if env and (Path(env) / "calicityper" / "__init__.py").exists():
        return Path(env).resolve()
    for parent in HERE.parents:
        if (parent / "calicityper" / "__init__.py").exists():
            return parent
    return None


_TOOL = _find_tool_checkout()
if _TOOL and str(_TOOL) not in sys.path:
    sys.path.insert(0, str(_TOOL))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
INPUTS = HERE / "inputs"
RESULTS = HERE / "results"
BASES = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T"}
BASES.update({k + "-": k for k in ("A", "C", "G", "T")})


def find_data_dir():
    """Find an installed reference set without importing calicityper.

    Order: CALICITYPER_DATA -> calicityper.config.DATA (if importable) -> common
    HPC cache locations. Returns None if nothing looks like a reference set.
    """
    import os
    env = os.environ.get("CALICITYPER_DATA")
    candidates = [Path(env)] if env else []
    if _TOOL:
        candidates.append(_TOOL / "calicityper" / "data")
    home = Path.home()
    cache = Path(os.environ.get("XDG_CACHE_HOME", home / ".cache")) / "calicityper"
    if cache.exists():
        candidates += sorted(cache.glob("reference-*"), reverse=True)
    for c in candidates:
        if c and (c / "database" / "refseq_2.fasta").exists():
            return c
    return None


def check_data_dir():
    """Return the verified reference-set directory, or raise with a setup hint.

    Deliberately does NOT import calicityper: the dev scripts must run from a
    flat HPC folder where the tool may not be on PYTHONPATH.
    """
    p = find_data_dir()
    if p is None:
        raise SystemExit(
            "Reference data not found. Run `calicityper setup`, set CALICITYPER_DATA "
            "to a verified set, or set CALICITYPER_REPO to the tool checkout. "
            "The repository does not ship the alignment files."
        )
    return p


def read_tsv(path, required=None):
    import csv
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if required:
        for r in rows:
            for key in required:
                if key not in r or r[key] in (None, ""):
                    raise SystemExit(f"{path}: missing required value for '{key}' in row {r}")
    return rows


def write_tsv(path, rows, columns):
    import csv
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {path}")


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.upper())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def write_fasta(path, seqs):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for n, s in seqs.items():
            fh.write(f">{n}\n{s}\n")
    print(f"wrote {path}")


def translate_nt(seq):
    """Translate the first complete codon-reading frame; return AA string."""
    from Bio.Seq import Seq
    return str(Seq(seq).translate(to_stop=False))


def p_distance(a, b):
    """Uncorrected p-distance over columns where both sequences carry a base.

    Truncates to the shared prefix length first: the mask is built over the
    zipped pairs, so indexing the original (possibly longer) strings would
    raise an IndexError.
    """
    n = min(len(a), len(b))
    x, y = a.upper()[:n], b.upper()[:n]
    ok = np.array([(p in BASES and q in BASES) for p, q in zip(x, y)])
    if not ok.any():
        return float("nan"), 0
    return float((np.array(list(x))[ok] != np.array(list(y))[ok]).mean()), int(ok.sum())


def provenance():
    """Manifest identity for generated tables, computed without the tool.

    Reads the small reference_manifest.json (from CALICITYPER_DATA or the tool
    checkout) and hashes its canonical JSON, matching calicityper.provenance.
    """
    import hashlib
    import json
    candidates = []
    d = find_data_dir()
    if d:
        candidates.append(d / "reference_manifest.json")
    if _TOOL:
        candidates.append(_TOOL / "calicityper" / "data" / "reference_manifest.json")
    for m in candidates:
        if m.exists():
            manifest = json.loads(m.read_text(encoding="utf-8"))
            canonical = json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode()
            return {
                "available": True,
                "manifest_path": str(m),
                "doi": manifest.get("doi", ""),
                "reference_set_version": manifest.get("reference_set_version", ""),
                "manifest_sha256": hashlib.sha256(canonical).hexdigest(),
            }
    return {"available": False, "manifest_path": "", "note": "no manifest found; set CALICITYPER_DATA or CALICITYPER_REPO"}
