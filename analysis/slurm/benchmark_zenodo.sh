#!/usr/bin/env bash
#SBATCH --job-name=zen_bench
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --time=180
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/zen_bench_%j.out
set -uo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd /home/people/fitzpatria/scratch/calicityper
echo "### host $(hostname) $(date) | commit $(git rev-parse --short HEAD)"
python scripts/benchmark_zenodo_amplicons.py \
  --data-dir /home/people/fitzpatria/scratch/calicityper_runs/zenodo_7271067 \
  --out-dir results/benchmark_zenodo --threads "${SLURM_CPUS_PER_TASK:-16}"
echo "### DONE $(date)"
