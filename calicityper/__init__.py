"""calicityper: RIVM-style norovirus genotyping."""
import os as _os

# Pin each process's BLAS/OpenMP pool to one thread before NumPy is imported.
# `--threads N` already runs N worker processes; left alone, each one's NumPy
# starts a pool sized to every core, and N x cores threads fight over the
# machine. A controlled A/B (BENCHMARKING.md, section 5) measured that as 62.9x
# slower with byte-identical calls. setdefault: a value the user exported wins.
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    _os.environ.setdefault(_var, "1")

__version__ = "0.7.0"
