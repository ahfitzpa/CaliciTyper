#!/usr/bin/env bash
# F1-F4 scaffolding driver, mirroring ~/scratch/calicityper_runs/scripts conventions.
# Run from the repo root on the UCD Sonic login/compute node with the verified
# reference set installed (CALICITYPER_DATA set by `calicityper setup`).
set -euo pipefail
cd "$(dirname "$0")/.."
module load python mafft 2>/dev/null || true   # adapt to Sonic module names
export TMPDIR="${TMPDIR:-/scratch/tmp}"

echo "== 0. validate inputs =="
python scripts/check_inputs.py

echo "== F1. build GII.4 translated VP1 AA alignment =="
python scripts/build_gii4_aa_reference.py

echo "== F1/F2. LOO threshold evidence (review results/f1_variant_loo.csv) =="
python scripts/validate_variant_thresholds.py

echo "== F2. sub-lineage review table =="
python scripts/build_sublineage_table.py

echo "== F3. representative selection + LOO comparison =="
python scripts/select_representatives.py
python scripts/validate_representatives.py

echo "== F4. immunotype/regime lookup (draft JSON) =="
python scripts/build_immunotype_lookup.py

echo "DONE. All outputs in results/ (untracked). Review docs/DESIGN.md promotion gate."
