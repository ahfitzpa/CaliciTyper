"""Locations of data and external programs. Nothing here is machine specific.

Data directory : $NOROTYPER_DATA, else <repo>/new_data
Programs       : $NOROTYPER_BLASTN / _MAKEBLASTDB / _CLUSTALW / _MAFFT, else whatever is on PATH
"""
import os, shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = Path(os.environ.get("NOROTYPER_DATA", ROOT / "new_data"))
BLAST_DB = DATA / "database" / "refseqdb_2"
BLAST_FASTA = DATA / "database" / "refseq_2.fasta"
ANNOTATION = DATA / "database" / "refseq_annotation.csv"
REF_DIR = {"orf1": DATA / "reference_sequences_RdRp", "orf2": DATA / "reference_sequences_VP1"}

NOROVIRUS_GROUPS = ["GI", "GII", "GIII", "GIV", "GV", "GVI", "GVII", "GVIII", "GIX", "GX"]


def exe(name):
    env = os.environ.get(f"NOROTYPER_{name.upper()}")
    candidates = {"clustalw": ["clustalw2", "clustalw"], "blastn": ["blastn"], "makeblastdb": ["makeblastdb"], "mafft": ["mafft"]}[name]
    path = env or next((shutil.which(c) for c in candidates if shutil.which(c)), None)
    if not path:
        raise FileNotFoundError(f"'{name}' not found on PATH; install it (conda install -c bioconda blast mafft clustalw) "
                                f"or set NOROTYPER_{name.upper()}")
    return path
