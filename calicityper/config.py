"""Locations of data and external programs. Nothing here is machine specific.

Data directory : $CALICITYPER_DATA, else the reference data bundled inside this
                 package (calicityper/data). Because the data lives *inside* the
                 package it is installed by pip along with the code, so a plain
                 `pip install calicityper` gives a working tool -- it does not
                 depend on the repository being present.
Programs       : $CALICITYPER_BLASTN / _MAKEBLASTDB / _CLUSTALW / _MAFFT, else whatever is on PATH

External programs are given real filesystem paths (BLAST and MAFFT are
subprocesses, so the data cannot be read from a zipped import).
"""
import os, shutil
from pathlib import Path

PKG = Path(__file__).resolve().parent
DATA = Path(os.environ.get("CALICITYPER_DATA", PKG / "data"))
BLAST_DB = DATA / "database" / "refseqdb_2"
BLAST_FASTA = DATA / "database" / "refseq_2.fasta"
ANNOTATION = DATA / "database" / "refseq_annotation.csv"
REF_DIR = {"orf1": DATA / "reference_sequences_RdRp", "orf2": DATA / "reference_sequences_VP1"}

NOROVIRUS_GROUPS = ["GI", "GII", "GIII", "GIV", "GV", "GVI", "GVII", "GVIII", "GIX", "GX"]


def check_data():
    """Raise a clear error if the reference data is missing or incomplete."""
    missing = [str(p) for p in (BLAST_FASTA, ANNOTATION, REF_DIR["orf1"], REF_DIR["orf2"]) if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "calicityper reference data not found: " + ", ".join(missing) +
            f"\nLooked in {DATA}. Set CALICITYPER_DATA to a copy of the data directory, "
            "or reinstall the package (the data ships inside it)."
        )


def exe(name):
    env = os.environ.get(f"CALICITYPER_{name.upper()}")
    candidates = {"clustalw": ["clustalw2", "clustalw"], "blastn": ["blastn"], "makeblastdb": ["makeblastdb"], "mafft": ["mafft"]}[name]
    path = env or next((shutil.which(c) for c in candidates if shutil.which(c)), None)
    if not path:
        raise FileNotFoundError(f"'{name}' not found on PATH; install it (conda install -c bioconda blast mafft clustalw) "
                                f"or set CALICITYPER_{name.upper()}")
    return path
