#!/usr/bin/env bash
#SBATCH --job-name=cal_figs
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --time=90
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/cal_figs_%j.out
#
# Regenerate all manuscript figures from the current benchmark outputs. No
# LaTeX: the manuscript is a Word document, built by manuscript/scripts/.
#
# figenv is treated as a cache and recreated here if absent, so it can be
# deleted whenever disk is needed (it is ~1.4 GB).
set -uo pipefail
M=/home/people/fitzpatria/scratch/calicityper_runs/manuscript
ENV=$M/figenv
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export MPLCONFIGDIR=$XDG_CACHE_HOME/matplotlib
mkdir -p "$MPLCONFIGDIR" "$M/figures"
echo "### host $(hostname) $(date)"
# figenv is a cache, not source: it is deleted when disk is reclaimed and
# rebuilt here on demand. Keep this in step with manuscript_figures.sh.
if [ ! -x "$ENV/bin/python" ]; then
    echo "### creating the plotting environment (figenv absent)"
    export CONDA_PKGS_DIRS=$XDG_CACHE_HOME/conda_pkgs
    source /opt/software/el9/anaconda3/2025.12-1/etc/profile.d/conda.sh
    rm -rf "$ENV"
    mamba create -y -p "$ENV" -c conda-forge python=3.12 matplotlib pandas numpy biopython \
        || conda create -y -p "$ENV" -c conda-forge python=3.12 matplotlib pandas numpy biopython \
        || { echo "### could not create figenv"; exit 1; }
fi
"$ENV/bin/python" -c "import matplotlib,pandas,Bio;print('matplotlib',matplotlib.__version__,'| pandas',pandas.__version__)" || exit 1
cd "$M"
"$ENV/bin/python" scripts/make_figures.py --out-dir "$M/figures"
echo "### figures now:"
ls -la "$M/figures"/*.png | awk '{print $NF, $5"B", $6, $7, $8}'
echo "### DONE $(date)"
