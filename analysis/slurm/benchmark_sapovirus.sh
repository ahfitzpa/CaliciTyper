#!/usr/bin/env bash
#SBATCH --job-name=bench_sapo
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --time=120
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/bench_sapo_%j.out
#
# Sapovirus genotyping benchmark (scripts/benchmark_sapovirus.py):
#   1. within-genotype p-distances of the reference sets vs the ranges
#      published by Zhao et al. 2024, and separation between genotypes
#   2. all RIVM sapovirus sequences typed with --type-sapovirus, VP1 genotype
#      compared with RIVM's capsid type
# Runs independently of the truncation job: coverage and distance are recorded
# for every call, so acceptance is re-applied once the floor is chosen
# (python scripts/benchmark_sapovirus.py ... --floor <x>, or re-read the CSV).
# If the truncation table already exists, its full-length rows (the
# leave-one-out over the reference strains) are summarised too.
#
#     sbatch benchmark_sapovirus.sh
#
# Outputs (calicityper/results/benchmark_sapovirus/, gitignored):
#     distance_ranges.csv, rivm_comparison.csv, rivm_sapovirus_typed.tsv
set -euo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
cd /home/people/fitzpatria/scratch/calicityper
mkdir -p results/benchmark_sapovirus
echo "### host $(hostname) $(date)"
echo "### commit $(git rev-parse --short HEAD)$(git diff --quiet || echo ' (+uncommitted changes)')"
python scripts/benchmark_sapovirus.py \
  --rivm-fasta /home/people/fitzpatria/scratch/calicityper_runs/rivm_genogroup/rivm_all.fasta \
  --out-dir results/benchmark_sapovirus \
  --loo results/truncation_loo_sapovirus.csv \
  --threads "${SLURM_CPUS_PER_TASK:-16}"
echo "### DONE $(date)"
