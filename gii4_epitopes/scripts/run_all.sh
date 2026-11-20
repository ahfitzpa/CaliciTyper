#!/usr/bin/env bash
# F1-F4 scaffolding driver. Works from ANY checkout layout: the local repo
# (.../calicityper/dev/gii4_epitopes) or a flat HPC copy
# (~/scratch/calicityper_runs/gii4_epitopes).
#
# On UCD Sonic, activate the calicityper conda environment first so that
# calicityper, biopython, numpy and MAFFT are available, then run:
#     bash scripts/run_all.sh
set -euo pipefail
cd "$(dirname "$0")/.."
HERE="$(pwd)"

# HPC modules already loaded in a conda env? Do not override; just make TMPDIR safe.
mkdir -p "${TMPDIR:-/tmp}"
export TMPDIR="${TMPDIR:-/tmp}"

# Tell every script where the tool checkout is, whether or not calicityper is
# importable in this shell. common.py also auto-detects it.
if [ -n "${CALICITYPER_REPO:-}" ]; then
    export PYTHONPATH="${CALICITYPER_REPO}:${PYTHONPATH:-}"
fi

echo "== 0. validate inputs =="
python scripts/check_inputs.py

echo "== F1. build GII.4 translated VP1 AA alignment (needs reference set + MAFFT) =="
python scripts/build_gii4_aa_reference.py || echo "SKIP F1 build: requires installed reference set and MAFFT"

echo "== F1/F2. LOO threshold evidence (review results/f1_variant_loo.csv) =="
python scripts/validate_variant_thresholds.py || echo "SKIP: needs AA alignment from F1"

echo "== F2. sub-lineage review table =="
python scripts/build_sublineage_table.py

echo "== F3. representative selection + LOO comparison =="
python scripts/select_representatives.py
python scripts/validate_representatives.py || echo "SKIP F3 LOO: needs reference set"

echo "== F4. immunotype/regime lookup =="
python scripts/build_immunotype_lookup.py

echo "== extras: aa x year heatmap from the fetched S3 proteins =="
if [ -f inputs/parra_s3/vp1_proteins.faa ]; then
    python scripts/aa_year_heatmap.py \
        inputs/parra_s3/accessions.tsv inputs/parra_s3/vp1_proteins.faa \
        results/heatmaps/ || echo "SKIP heatmap"
fi

echo "== cross-check Parra S4 accessions against bundled reference set =="
python scripts/crosscheck_parra_accessions.py || echo "SKIP cross-check: needs reference set"

echo "DONE. Outputs in results/ (untracked). Review docs/DESIGN.md promotion gate."
