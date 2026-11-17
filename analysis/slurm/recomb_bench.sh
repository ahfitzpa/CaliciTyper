#!/usr/bin/env bash
#SBATCH --job-name=cal_recomb
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --time=480
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/recomb_bench_%j.out
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --mail-type=END,FAIL
set -uo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export CALICITYPER_OPENRDP=/home/people/fitzpatria/scratch/openrdp_env/bin/openrdp
export PYTHONPATH=/home/people/fitzpatria/scratch/calicityper
cd /home/people/fitzpatria/scratch/calicityper
echo "### host $(hostname)  $(date)  openrdp=$CALICITYPER_OPENRDP"
echo "### sapovirus recombination, end to end over the three sapovirus datasets"
python -m pytest tests/test_recomb.py -v -p no:cacheprovider \
    -k "sapovirus and (end_to_end or carries)" --durations=0
rc_sapo=$?
echo "### norovirus recombination, end to end over the six norovirus datasets"
python -m pytest tests/test_recomb.py -v -p no:cacheprovider \
    -k "synthetic_dataset or clear_recombinant or orf_mismatch or same_seed or multi_query" \
    --durations=0
rc_noro=$?
echo "### sapovirus exit $rc_sapo | norovirus exit $rc_noro"
echo "### DONE $(date)"
