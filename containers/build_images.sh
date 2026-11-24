#!/bin/bash
#SBATCH --exclude=sonic59
#SBATCH --job-name=cal_build
#SBATCH --output=logs/build_%j.log
#SBATCH --error=logs/build_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#
# 1. Build both images from the repository root, tools first: calicityper.def
#    bootstraps from calicityper-tools.sif and %files paths are resolved
#    against the current directory.
# 2. Apptainer's scratch goes to $TMPDIR, not home: a conda solve plus two
#    image squashes will not fit in a home quota.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
mkdir -p logs
module purge
module load apptainer/1.3.4-gcc-11.5.0-ojp6nts
export APPTAINER_TMPDIR="${TMPDIR:-/scratch/fitzpatria/tmp}/apptainer-$SLURM_JOB_ID"
export APPTAINER_CACHEDIR="/scratch/fitzpatria/tmp/apptainer-cache"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"
trap 'rm -rf "$APPTAINER_TMPDIR"' EXIT

echo "[build] $(date): apptainer $(apptainer --version)"
status=0
for def in calicityper-tools calicityper; do
    echo "[build] $(date): $def"
    apptainer build --force "$def.sif" "containers/$def.def"
    rc=$?
    echo "[build] $(date): $def exit $rc"
    [ $rc -eq 0 ] || { status=1; break; }
done

# 3. %test runs at build time, but run it again on the finished images so the
#    log records what the artifacts actually do.
if [ $status -eq 0 ]; then
    for def in calicityper-tools calicityper; do
        echo "[build] test $def.sif"
        apptainer test "$def.sif" || status=1
    done
    ls -la calicityper-tools.sif calicityper.sif
fi
echo "[build] $(date): done, status $status"
exit $status
