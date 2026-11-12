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
# Genotyping is validated for GI and GII only, and the reference sets for the
# animal genogroups hold two or three sequences each -- enough to return a
# label, not enough for it to mean anything. Typing those spends the bootstrap
# (~95% of runtime) to produce a call the docs then tell users not to trust, so
# they are reported to genogroup and stop there. CALICITYPER_TYPE_ALL_GROUPS=1
# restores the old behaviour for reference-set development.
# RIVM/NoroNet operational typing windows: one per gene, not the A-E amplicon
# nomenclature. Region C is measured from the 5' end of VP1; the RdRp window is
# the 3' end of RdRp, which is what the ORF1 reference alignment already spans.
# Lengths are genogroup-dependent (GI.9 is a genotype, so its 267 nt only
# applies once a GI.9 call exists; GI's 264 is used until then).
RIVM_WINDOW_NT = {
    "orf2": {"GI": 264, "GII": 252, "GIX": 255, "_default": 252, "GI.9": 267},
    "orf1": {"GIV": 173, "_default": 172},
}
# GII.4 variant designation needs a third window. Region C is the VP1 shell
# domain, conserved BETWEEN GII.4 variants by construction -- measured on the
# bundled references, only 4 of 135 variant-diagnostic columns fall inside it,
# so a region C amplicon cannot support a variant call however long it is.
# The signal is 3' of there: P2 (aa 279-405) holds 41% of those columns, and P2
# through the C-terminus holds 70%, which is what this window spans. Strict P2
# was rejected for leaving more than half the signal out.
GII4_VARIANT_WINDOW_AA = (279, None)      # None = to the 3' end of VP1

# Fraction of a window a query must cover for a call to stand.
RIVM_WINDOW_MIN_COVERAGE = float(os.environ.get("CALICITYPER_MIN_WINDOW_COVERAGE", 0.80))
# Nothing shorter than the shortest window can cover any window, so it is
# rejected before BLAST rather than typed and then explained.
MIN_INPUT_NT = min(min(v for k, v in d.items() if k != "GI.9") for d in RIVM_WINDOW_NT.values())

TYPEABLE_GROUPS = ["GI", "GII"]
if os.environ.get("CALICITYPER_TYPE_ALL_GROUPS", "").lower() in ("1", "true", "yes"):
    TYPEABLE_GROUPS = list(NOROVIRUS_GROUPS)


def check_data():
    """Raise a clear error if the reference data is missing or incomplete."""
    missing = [str(p) for p in (BLAST_FASTA, ANNOTATION, REF_DIR["orf1"], REF_DIR["orf2"]) if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "calicityper reference data not found: " + ", ".join(missing) +
            f"\nLooked in {DATA}. Run `calicityper setup` to fetch and verify the "
            "reference set (Zenodo doi:10.5281/zenodo.24161236), or set "
            "CALICITYPER_DATA to an existing verified copy."
        )


def exe(name):
    env = os.environ.get(f"CALICITYPER_{name.upper()}")
    candidates = {"clustalw": ["clustalw2", "clustalw"], "blastn": ["blastn"], "makeblastdb": ["makeblastdb"], "mafft": ["mafft"]}[name]
    path = env or next((shutil.which(c) for c in candidates if shutil.which(c)), None)
    if not path:
        raise FileNotFoundError(f"'{name}' not found on PATH; install it (conda install -c bioconda blast mafft clustalw) "
                                f"or set CALICITYPER_{name.upper()}")
    return path
