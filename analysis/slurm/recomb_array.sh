#!/usr/bin/env bash
#SBATCH --job-name=cal_recomb
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --array=0-12
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/recomb_%A_%a.out
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --mail-type=END,FAIL
# One array task per slow test. The previous single job hit its 8 h limit and
# was killed before pytest printed any failure detail, so two failures went
# unexplained; one test per task means a slow dataset cannot starve the others
# and each failure is written to its own log as it happens.
set -uo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export CALICITYPER_OPENRDP=/home/people/fitzpatria/scratch/openrdp_env/bin/openrdp
export PYTHONPATH=/home/people/fitzpatria/scratch/calicityper
export PYTHONUNBUFFERED=1
cd /home/people/fitzpatria/scratch/calicityper
mapfile -t TESTS < <(python -m pytest tests/test_recomb.py -m slow -q -p no:cacheprovider --co 2>/dev/null | grep '::')
T="${TESTS[$SLURM_ARRAY_TASK_ID]:-}"
if [ -z "$T" ]; then echo "no test at index $SLURM_ARRAY_TASK_ID"; exit 0; fi
echo "### $(hostname) $(date)"
echo "### task $SLURM_ARRAY_TASK_ID: $T"
python -m pytest "$T" -m slow -v -p no:cacheprovider --tb=long -rA --durations=0
echo "### exit $? $(date)"
