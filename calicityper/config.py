"""Locations of data and external programs.

Data directory: $CALICITYPER_DATA, else /data inside a container, else calicityper/data.
Programs: $CALICITYPER_{BLASTN,MAKEBLASTDB,CLUSTALW,MAFFT,IQTREE}, else PATH.
"""
import os
import shutil
from pathlib import Path

PKG = Path(__file__).resolve().parent


def _default_data():
    """Return /data inside Apptainer/Singularity (read-only image), else package data/."""
    if os.environ.get("APPTAINER_CONTAINER") or os.environ.get("SINGULARITY_CONTAINER"):
        return Path("/data")
    return PKG / "data"


# Use `or` so an empty CALICITYPER_DATA="" falls back instead of resolving to Path("").
DATA = Path(os.environ.get("CALICITYPER_DATA") or _default_data())
BLAST_DB = DATA / "database" / "refseqdb_2"
BLAST_FASTA = DATA / "database" / "refseq_2.fasta"
# Fixed BLAST search-space size (v0.7.0 60-subject DB) so adding genomes cannot shift -max_hsps 1 calls.
BLAST_DBSIZE = 452925
# Minimum BLAST alignment length for a genogroup call (<100 nt yields false matches on non-caliciviruses).
MIN_GENOGROUP_HIT_NT = int(os.environ.get("CALICITYPER_MIN_GENOGROUP_HIT_NT", 100))
ANNOTATION = DATA / "database" / "refseq_annotation.csv"
REF_DIR = {"orf1": DATA / "reference_sequences_RdRp", "orf2": DATA / "reference_sequences_VP1"}

# Sapovirus: both typed regions lie in ORF1; output "orf1"=RdRp (P-type), "orf2"=VP1 (genotype).
SAPOVIRUS_PREFIX = "Caliciviridae_Sapovirus_"
SAPOVIRUS_REF_DIR = {
    "orf1": DATA / "reference_sequences_sapovirus_RdRp",
    "orf2": DATA / "reference_sequences_sapovirus_VP1",
}
# Zhao et al. 2024 (BMC Genomics 25:1131) nucleotide p-distance cut-offs: (same genotype, same genogroup).
SAPOVIRUS_CUTOFFS = {"orf1": (0.266, 0.531), "orf2": (0.161, 0.503)}
# Reported, not gated by default (truncation LOO showed 1.000 VP1 precision at floor 0.0).
SAPOVIRUS_MIN_REGION_COVERAGE = float(os.environ.get("CALICITYPER_SAPOVIRUS_MIN_COVERAGE", 0.0))
TYPE_SAPOVIRUS = os.environ.get("CALICITYPER_TYPE_SAPOVIRUS", "1").lower() in ("1", "true", "yes")


def is_sapovirus(group):
    return str(group).startswith(SAPOVIRUS_PREFIX)


def ref_paths(region, group):
    """(aln, orf) paths of the reference set typing `group` in `region`, or None.

    Norovirus falls back to the all-genogroup set; sapovirus does not, because
    its calls use within-genogroup distance cut-offs.
    """
    if is_sapovirus(group):
        d, stem = SAPOVIRUS_REF_DIR[region], group[len(SAPOVIRUS_PREFIX):]
        aln, orf = d / f"{stem}.aln", d / f"{stem}.orf"
        return (aln, orf) if aln.exists() and orf.exists() else None
    d = REF_DIR[region]
    if (d / f"{group}.aln").exists():
        return d / f"{group}.aln", d / f"{group}.orf"
    if (d / "all.aln").exists():
        return d / "all.aln", d / "all.orf"
    return None


NOROVIRUS_GROUPS = ["GI", "GII", "GIII", "GIV", "GV", "GVI", "GVII", "GVIII", "GIX", "GX"]
# RIVM/NoroNet operational typing windows (nt): VP1 5' Region C and RdRp 3' window.
RIVM_WINDOW_NT = {
    "orf2": {"GI": 264, "GII": 252, "GIX": 255, "_default": 252, "GI.9": 267},
    "orf1": {"GIV": 173, "_default": 172},
}
# GII.4 variant window: P2 through VP1 C-terminus (aa 279-end); Region C lacks variant-diagnostic signal.
GII4_VARIANT_WINDOW_AA = (279, None)

# Window coverage is reported only; RIVM_WINDOW_MODE restores the hard 0.80 gate for RIVM replication.
RIVM_WINDOW_MIN_COVERAGE = float(os.environ.get("CALICITYPER_MIN_WINDOW_COVERAGE", 0.80))
RIVM_WINDOW_MODE = os.environ.get("CALICITYPER_RIVM_WINDOW_MODE", "").lower() in ("1", "true", "yes")

# Acceptance floor in aligned informative sites (calls below are still reported with `_accepted` reason).
MIN_INFORMATIVE_SITES = int(os.environ.get("CALICITYPER_MIN_SITES", 300))
MIN_INFORMATIVE_SITES_ORF1 = (
    int(os.environ["CALICITYPER_MIN_SITES_ORF1"])
    if os.environ.get("CALICITYPER_MIN_SITES_ORF1") else None
)
MIN_INFORMATIVE_SITES_ORF2 = (
    int(os.environ["CALICITYPER_MIN_SITES_ORF2"])
    if os.environ.get("CALICITYPER_MIN_SITES_ORF2") else None
)
PRESETS = {
    "default": {},
    "polcap": {"min_sites_orf1": 170, "min_sites_orf2": 250},
}


def min_sites_for(region):
    """Return the minimum aligned informative sites floor for `region` ('orf1' or 'orf2')."""
    if region == "orf1" and MIN_INFORMATIVE_SITES_ORF1 is not None:
        return MIN_INFORMATIVE_SITES_ORF1
    if region == "orf2" and MIN_INFORMATIVE_SITES_ORF2 is not None:
        return MIN_INFORMATIVE_SITES_ORF2
    return MIN_INFORMATIVE_SITES


# Minimum input length (nt); matches the 100-site alignment floor.
MIN_INPUT_NT = int(os.environ.get("CALICITYPER_MIN_INPUT_NT", 100))

# Validated for GI/GII; set CALICITYPER_TYPE_ALL_GROUPS=1 to include sparse animal genogroups.
DEFAULT_TYPEABLE_GROUPS = ["GI", "GII"]
TYPEABLE_GROUPS = list(DEFAULT_TYPEABLE_GROUPS)
if os.environ.get("CALICITYPER_TYPE_ALL_GROUPS", "").lower() in ("1", "true", "yes"):
    TYPEABLE_GROUPS = list(NOROVIRUS_GROUPS)


def check_data():
    """Raise a clear error if the reference data is missing or incomplete."""
    missing = [str(p) for p in (BLAST_FASTA, ANNOTATION, REF_DIR["orf1"], REF_DIR["orf2"]) if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "calicityper reference data not found: " + ", ".join(missing) +
            f"\nLooked in {DATA}. Run `calicityper setup` to fetch and verify the "
            "reference set (Zenodo doi:10.5281/zenodo.24354774), or set "
            "CALICITYPER_DATA to an existing verified copy."
        )


def exe(name):
    env = os.environ.get(f"CALICITYPER_{name.upper()}")
    candidates = {
        "clustalw": ["clustalw2", "clustalw"],
        "blastn": ["blastn"],
        "makeblastdb": ["makeblastdb"],
        "mafft": ["mafft"],
        "iqtree": ["iqtree2", "iqtree"],
    }[name]
    path = env or next((shutil.which(c) for c in candidates if shutil.which(c)), None)
    if not path:
        raise FileNotFoundError(
            f"'{name}' not found on PATH; install it (conda install -c bioconda "
            f"blast mafft clustalw iqtree) or set CALICITYPER_{name.upper()}"
        )
    return path
