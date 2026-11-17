#!/usr/bin/env bash
#SBATCH --job-name=trunc_loo
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --time=360
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/trunc_loo_%j.out
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
cd /home/people/fitzpatria/scratch/calicityper
mkdir -p results
python scripts/truncation_loo.py \
  --out results/truncation_loo.csv \
  --bootstraps 100 --min-support 70 --threads 16
echo "=== DONE ==="
