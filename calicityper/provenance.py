"""Provenance and atomic output helpers for machine-readable results."""
from __future__ import annotations

from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
from datetime import datetime, timezone
from typing import Any

from . import __version__
from . import config as cg


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def atomic_write_text(path: str | Path, text: str, encoding: str = "utf-8") -> Path:
    """Write text to a sibling temp file and atomically commit via os.replace()."""
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding=encoding,
            newline="",
            dir=dest.parent,
            prefix=f".{dest.name}.",
            suffix=".tmp",
            delete=False,
        ) as fh:
            tmp_path = Path(fh.name)
            fh.write(text)
            fh.flush()
        os.replace(tmp_path, dest)
        tmp_path = None
        return dest
    finally:
        if tmp_path is not None and tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                pass


def file_sha256(path: str | Path | None) -> str:
    """Return hex SHA-256 of a file, or empty string if unreadable/absent."""
    if not path:
        return ""
    p = Path(path)
    if not p.is_file():
        return ""
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


@lru_cache(maxsize=None)
def _tool_version(name: str, exe_path: str) -> str:
    flag = "-version" if name in ("blastn", "makeblastdb") else "--version"
    try:
        proc = subprocess.run([exe_path, flag], capture_output=True, text=True, timeout=5)
        out = (proc.stdout or "") + "\n" + (proc.stderr or "")
        for line in out.splitlines():
            line = line.strip()
            if line:
                return line
    except Exception:
        pass
    return "unknown"


def external_tool_versions(names: tuple[str, ...] = ("blastn", "mafft")) -> dict[str, str]:
    """Return version strings for external CLI tools used by a run."""
    versions: dict[str, str] = {}
    for name in names:
        try:
            exe_path = cg.exe(name)
        except FileNotFoundError:
            versions[name] = "not_found"
            continue
        versions[name] = _tool_version(name, exe_path)
    return versions


@lru_cache(maxsize=1)
def runtime_environment() -> dict[str, Any]:
    """Return Python version and installed core/optional dependency versions."""
    import importlib.metadata
    import platform

    deps: dict[str, str] = {}
    for pkg in ("biopython", "numpy", "pandas", "numba", "bokeh"):
        try:
            deps[pkg] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            pass
    return {
        "python_version": platform.python_version(),
        "dependencies": deps,
    }


def _manifest_path(manifest_path: str | Path | None = None) -> Path:
    if manifest_path:
        return Path(manifest_path)
    return Path(__file__).resolve().parent / "data" / "reference_manifest.json"


def reference_set_info(manifest_path: str | Path | None = None) -> dict[str, Any]:
    """Return reference-set identity/provenance from the shipped manifest."""
    path = _manifest_path(manifest_path)
    info: dict[str, Any] = {"manifest_path": str(path), "available": False}
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


def base_provenance(
    tool: str = "calicityper",
    command: str | None = None,
    parameters: dict[str, Any] | None = None,
    include_reference_set: bool = True,
    manifest_path: str | Path | None = None,
    input_path: str | Path | None = None,
    tool_names: tuple[str, ...] | None = None,
) -> dict[str, Any]:
    prov: dict[str, Any] = {
        "tool": tool,
        "tool_version": __version__,
        "generated": _now_iso(),
        "runtime": runtime_environment(),
    }
    if command:
        prov["command"] = command
    if input_path is not None:
        prov["input_sha256"] = file_sha256(input_path)
    if tool_names is not None:
        prov["external_tools"] = external_tool_versions(tool_names)
    if parameters is not None:
        prov["parameters"] = parameters
    if include_reference_set:
        prov["reference_set"] = reference_set_info(manifest_path=manifest_path)
    return prov
