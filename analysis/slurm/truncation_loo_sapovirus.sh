#!/usr/bin/env bash
#SBATCH --job-name=trunc_sapo
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --time=120
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/trunc_sapo_%j.out
#
# Sapovirus truncation leave-one-out (scripts/truncation_loo_sapovirus.py).
# Every sapovirus reference, cut to 150-1200 nt at its 5' end, centre and 3'
# end (plus the whole region), typed against the other references of its
# genogroup with the Zhao et al. 2024 distance cut-offs, scored against its
# published type. Gives the coverage floor that replaces the provisional
# SAPOVIRUS_MIN_REGION_COVERAGE = 0.9 in calicityper/config.py.
#
#     sbatch truncation_loo_sapovirus.sh
#
# Outputs (calicityper/results/, gitignored):
#     truncation_loo_sapovirus.csv          one row per fragment
#     truncation_loo_sapovirus.floors.csv   precision/recall per coverage floor
set -euo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
cd /home/people/fitzpatria/scratch/calicityper
mkdir -p results
echo "### host $(hostname) $(date)"
echo "### commit $(git rev-parse --short HEAD)$(git diff --quiet || echo ' (+uncommitted changes)')"
python scripts/truncation_loo_sapovirus.py \
  --out results/truncation_loo_sapovirus.csv \
  --bootstraps 20 --threads "${SLURM_CPUS_PER_TASK:-16}"
echo "### DONE $(date)"
