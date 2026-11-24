#!/bin/bash
#SBATCH --exclude=sonic59
#SBATCH --job-name=cal_build
#SBATCH --output=logs/build_%j.log
#SBATCH --error=logs/build_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#
# Build calicityper-tools.sif then calicityper.sif from the repo root.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
mkdir -p logs containers/cache

module purge
module load apptainer/1.3.4-gcc-11.5.0-ojp6nts

SCRATCH_TMP="${TMPDIR:-/scratch/${USER}/tmp}"
export APPTAINER_TMPDIR="${SCRATCH_TMP}/apptainer-${SLURM_JOB_ID:-$$}"
export APPTAINER_CACHEDIR="/scratch/${USER}/tmp/apptainer-cache"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"
trap 'rm -rf "$APPTAINER_TMPDIR"' EXIT

echo "[build] $(date): $(apptainer --version)"
status=0
for def in calicityper-tools calicityper; do
    echo "[build] $(date): $def"
    apptainer build --force "$def.sif" "containers/$def.def"
    rc=$?
    echo "[build] $(date): $def exit $rc"
    [ $rc -eq 0 ] || { status=1; break; }
done

if [ $status -eq 0 ]; then
    for def in calicityper-tools calicityper; do
        echo "[build] test $def.sif"
        apptainer test "$def.sif" || status=1
    done
    ls -lh calicityper-tools.sif calicityper.sif
fi
echo "[build] $(date): done, status $status"
exit $status
