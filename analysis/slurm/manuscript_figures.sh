#!/usr/bin/env bash
#SBATCH --job-name=cal_manuscript
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --time=120
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/manuscript_%j.out
#
# Build the benchmarking manuscript: figures from the benchmark outputs, then
# the PDF. Everything runs here rather than on the login node.
#
# Creates ../manuscript/figenv on first run (matplotlib/pandas; the rebuilt
# calicityper_env has neither, and bokeh is for the tool's own HTML reports).
#
#     sbatch manuscript_figures.sh
set -uo pipefail
M=/home/people/fitzpatria/scratch/calicityper_runs/manuscript
ENV=$M/figenv
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
export CONDA_PKGS_DIRS=$XDG_CACHE_HOME/conda_pkgs
export MPLCONFIGDIR=$XDG_CACHE_HOME/matplotlib
mkdir -p "$MPLCONFIGDIR" "$M/figures"
echo "### host $(hostname) $(date)"

if [ ! -x "$ENV/bin/python" ]; then
    echo "### creating the plotting environment"
    source /opt/software/el9/anaconda3/2025.12-1/etc/profile.d/conda.sh
    rm -rf "$ENV"
    mamba create -y -p "$ENV" -c conda-forge python=3.12 matplotlib pandas numpy biopython \
        || conda create -y -p "$ENV" -c conda-forge python=3.12 matplotlib pandas numpy biopython
fi
if ! "$ENV/bin/python" -c "import Bio" 2>/dev/null; then
    echo "### adding biopython (needed for the distance-distribution figure)"
    source /opt/software/el9/anaconda3/2025.12-1/etc/profile.d/conda.sh
    mamba install -y -p "$ENV" -c conda-forge biopython || "$ENV/bin/pip" install biopython
fi
"$ENV/bin/python" -c "import matplotlib, pandas, Bio; print('matplotlib', matplotlib.__version__, '| pandas', pandas.__version__, '| biopython', Bio.__version__)" || exit 1

echo "### figures  $(date)"
cd "$M"
"$ENV/bin/python" scripts/make_figures.py --out-dir "$M/figures"
ls -la "$M/figures"

echo "### LaTeX  $(date)"
export PATH=/home/people/fitzpatria/scratch/tools/texlive2024/bin/x86_64-linux:$PATH
pdflatex -interaction=nonstopmode main.tex >/dev/null
bibtex main >/dev/null 2>&1 || true
pdflatex -interaction=nonstopmode main.tex >/dev/null
pdflatex -interaction=nonstopmode main.tex | tail -5
if [ -s main.pdf ]; then
    echo "### built main.pdf ($(stat -c%s main.pdf) bytes, $(pdfinfo main.pdf 2>/dev/null | awk '/^Pages/{print $2}') pages)"
else
    echo "### LaTeX FAILED -- see main.log"; grep -m10 "^!" main.log
fi
grep -c "Citation.*undefined" main.log 2>/dev/null | sed 's/^/undefined citations: /'
echo "### DONE $(date)"
