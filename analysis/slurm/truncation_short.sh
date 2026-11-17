#!/usr/bin/env bash
#SBATCH --job-name=trunc_short
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --time=360
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/trunc_short_%j.out
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
cd /home/people/fitzpatria/scratch/calicityper
mkdir -p results
python scripts/truncation_loo.py \
  --out results/truncation_loo_short.csv \
  --bootstraps 100 --min-support 70 --threads 16 --no-input-floor --lengths 80,100,120,140,160,172 --offsets 0
echo "=== DONE ==="
