#!/usr/bin/env bash
#SBATCH --job-name=recomb_retry
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --time=300
#SBATCH --array=1,4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/recomb_retry_%A_%a.out
#
# Re-run the two Table 3 strains that produced no scanner output at all:
#   index 1 GQ261222  openrdp timed out after 7200s
#   index 4 MZ488271  openrdp exit 1, IndexError on the last alignment column
# Both now benefit from the per-method retry: a shorter combined budget fails
# fast, then each method runs on its own, so one crashing scanner no longer
# discards the other six.
set -uo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
# OpenRDP lives in its own env (numpy<2); without this it is simply 'not found'
export CALICITYPER_OPENRDP=/home/people/fitzpatria/scratch/openrdp_env/bin/openrdp
cd /home/people/fitzpatria/scratch/calicityper
echo "### host $(hostname) $(date) | task ${SLURM_ARRAY_TASK_ID} | commit $(git rev-parse --short HEAD)"
python scripts/validate_recombination_table3.py --index "${SLURM_ARRAY_TASK_ID}" \
    --out-dir results/recomb_table3 --threads "${SLURM_CPUS_PER_TASK:-8}" \
    --openrdp-timeout 1800
echo "### exit $? $(date)"
