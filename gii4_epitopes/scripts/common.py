#!/usr/bin/env python3
"""Shared helpers for the F1-F4 development scripts. No tool behaviour changes."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]          # .../calicityper (scripts -> gii4_epitopes -> dev -> calicityper)
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
HERE = ROOT / "dev" / "gii4_epitopes"
INPUTS = HERE / "inputs"
RESULTS = HERE / "results"
BASES = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T"}
BASES.update({k + "-": k for k in ("A", "C", "G", "T")})


def check_data_dir():
    """Return the verified reference-set directory, or raise with setup hint."""
    import os
    from calicityper import config as cg
    p = cg.DATA
    missing = [str(x) for x in (cg.BLAST_FASTA, cg.ANNOTATION,
                                cg.REF_DIR["orf1"], cg.REF_DIR["orf2"]) if not x.exists()]
    if missing:
        raise SystemExit(
            "Reference data not found in " + str(p) + "\n"
            "Missing: " + ", ".join(missing) + "\n"
            "Run `calicityper setup` or set CALICITYPER_DATA to a verified set.\n"
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
    """Uncorrected p-distance over columns where both sequences carry a base."""
    ok = np.array([(x in BASES and y in BASES) for x, y in zip(a.upper(), b.upper())])
    if not ok.any():
        return float("nan"), 0
    aa = np.array(list(a.upper()))[ok]
    bb = np.array(list(b.upper()))[ok]
    return float((aa != bb).mean()), int(ok.sum())


def provenance():
    """Manifest identity for generated tables, so outputs can be audited."""
    from calicityper import provenance as prov
    return prov.reference_set_info()
