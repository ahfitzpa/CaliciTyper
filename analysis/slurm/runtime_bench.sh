#!/usr/bin/env bash
#SBATCH --job-name=cal_runtime
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --time=240
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/runtime_%j.out
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --mail-type=END,FAIL
set -uo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
cd /home/people/fitzpatria/scratch/calicityper
R=/home/people/fitzpatria/scratch/calicityper_runs/rivm_genogroup/rivm_all.fasta
echo "### host $(hostname)  $(date)  cpus=${SLURM_CPUS_PER_TASK:-?}"
# Real, GII-dominated input: what a surveillance user actually submits.
python scripts/benchmark_runtime.py --source "$R" \
    --sizes 1 4 96 --threads 1 4 16 \
    --out results/runtime_rivm.tsv
# Bundled references (GI-heavy): reproducible by any user with no extra data.
python scripts/benchmark_runtime.py \
    --sizes 1 4 --threads 1 16 \
    --out results/runtime_bundled.tsv
echo "### DONE $(date)"
