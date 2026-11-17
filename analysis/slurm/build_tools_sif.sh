#!/usr/bin/env bash
#SBATCH --job-name=tools_sif
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --time=120
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/tools_sif_%j.out
#
# Build the pinned external-tools image and run its %test block.
#
# The first attempt (1159831) failed because apt-get cannot drop privileges in
# a root-mapped namespace. The definition no longer uses apt at all. A second
# failure on the login node was mksquashfs "Failed to create thread", a
# resource limit there, so this runs on a compute node.
set -uo pipefail
source /etc/profile.d/modules.sh
module load apptainer/1.3.4-gcc-11.5.0-ojp6nts
cd /home/people/fitzpatria/scratch/calicityper
export APPTAINER_CACHEDIR=/scratch/fitzpatria/.apptainer_cache
export APPTAINER_TMPDIR=/scratch/fitzpatria/.apptainer_tmp
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
SIF=/scratch/fitzpatria/calicityper_runs/calicityper-tools.sif
echo "### host $(hostname) $(date) | $(apptainer --version) | cpus ${SLURM_CPUS_PER_TASK:-?}"
echo "### ulimit -u: $(ulimit -u)"

# UCD's Apptainer guidance (KB0011720) gives:
#   apptainer build --fakeroot --sandbox --bind="$TMPDIR:/tmp" out.sif in.def
# --fakeroot is omitted: fakeroot is not installed and this account has no
# /etc/subuid range, so it cannot work here. --bind of TMPDIR is adopted.
export TMPDIR="${TMPDIR:-$APPTAINER_TMPDIR}"
apptainer build --force --bind="$TMPDIR:/tmp" "$SIF" containers/calicityper-tools.def
rc=$?
echo "### build exit $rc"

# --sandbox writes a directory tree instead of a SIF, so it never calls
# mksquashfs -- which is what failed on the login node with "Failed to create
# thread". Fall back to it, then convert, so a squashfs limit is not fatal.
if [ $rc -ne 0 ]; then
    echo "### SIF build failed; retrying as a sandbox (KB0011720 --sandbox)"
    SBX=/scratch/fitzpatria/calicityper_runs/calicityper-tools.sandbox
    rm -rf "$SBX"
    apptainer build --force --sandbox --bind="$TMPDIR:/tmp" "$SBX" \
        containers/calicityper-tools.def
    src=$?
    echo "### sandbox build exit $src"
    if [ $src -eq 0 ]; then
        echo "### sandbox built; the tools are usable directly from it:"
        apptainer exec "$SBX" blastn -version 2>&1 | head -1
        echo "### converting the sandbox to a SIF"
        apptainer build --force "$SIF" "$SBX" && rc=0
        echo "### convert exit $?"
    fi
fi
if [ $rc -eq 0 ]; then
    echo "### sha256:"; sha256sum "$SIF"
    echo "### size: $(du -h "$SIF" | cut -f1)"
    echo "### versions from inside the image:"
    apptainer exec "$SIF" blastn -version 2>&1 | head -1
    apptainer exec "$SIF" makeblastdb -version 2>&1 | head -1
    apptainer exec "$SIF" mafft --version 2>&1 | head -1
    apptainer exec "$SIF" openrdp --help >/dev/null 2>&1 && echo "openrdp: runs"
    apptainer exec "$SIF" /opt/openrdp-venv/bin/python -c "import numpy;print('openrdp numpy',numpy.__version__)" 2>&1 | tail -1
else
    echo "### FAILED -- the definition stays marked unverified"
fi
echo "### done $(date)"
