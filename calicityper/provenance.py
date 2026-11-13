"""Provenance helpers for machine-readable outputs.

These functions collect audit-relevant metadata (tool version, timestamp,
parameters, and reference-set identity from the shipped manifest) so outputs
can record exactly how they were produced.
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from . import __version__


def _now_iso():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _manifest_path(manifest_path=None):
    if manifest_path:
        return Path(manifest_path)
    # shipped manifest (small metadata) lives in package data
    return Path(__file__).resolve().parent / "data" / "reference_manifest.json"


def reference_set_info(manifest_path=None):
    """Return reference-set identity/provenance from the manifest, if available.

    Includes DOI, reference_set_version, manifest hash, and artefact counts.
    The manifest hash is the SHA-256 of canonical JSON
    (sorted keys, compact separators) so it is stable and auditable.
    """
    path = _manifest_path(manifest_path)
    info = {"manifest_path": str(path), "available": False}
    if not path.exists():
        return info
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except Exception as e:
        info["error"] = f"{type(e).__name__}: {e}"
        return info
    canonical = json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode("utf-8")
    info.update({
        "available": True,
        "doi": manifest.get("doi", ""),
        "reference_set_version": manifest.get("reference_set_version", ""),
        "download_url": manifest.get("download_url", ""),
        "n_accessions": manifest.get("n_accessions"),
        "n_records": manifest.get("n_records"),
        "n_artefacts": manifest.get("n_artefacts"),
        "manifest_sha256": hashlib.sha256(canonical).hexdigest(),
    })
    return info


def base_provenance(tool="calicityper", command=None, parameters=None, include_reference_set=True,
                    manifest_path=None):
    prov = {
        "tool": tool,
        "tool_version": __version__,
        "generated": _now_iso(),
    }
    if command:
        prov["command"] = command
    if parameters is not None:
        prov["parameters"] = parameters
    if include_reference_set:
        prov["reference_set"] = reference_set_info(manifest_path=manifest_path)
    return prov
