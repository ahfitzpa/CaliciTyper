#!/usr/bin/env bash
#SBATCH --job-name=cal_build
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --time=90
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=amy.fitzpatrick@ucd.ie
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/cal_build_%j.out
#
# Build the 0.8.0 release artefacts and prove they work from a clean install:
#   1. reference archive, built from the manifest and re-verified
#   2. sdist + wheel
#   3. wheel audit: no reference data may be inside the wheel
#   4. fresh venv, install the wheel, setup from the archive, type a sequence
set -uo pipefail
export PATH=/home/people/fitzpatria/scratch/calicityper_env/bin:$PATH
export XDG_CACHE_HOME=/scratch/fitzpatria/.cache
cd /home/people/fitzpatria/scratch/calicityper
echo "### host $(hostname) $(date) | commit $(git rev-parse --short HEAD)"

echo; echo "### 1. reference archive"
python scripts/build_reference_archive.py --out-dir dist --apply || exit 1
ARCHIVE=$(ls -t dist/calicityper-reference-*.tar.gz | head -1)

echo; echo "### 2. sdist + wheel"
python -m pip install --quiet --upgrade build 2>&1 | tail -2
rm -rf dist/*.whl dist/*.tar.gz.tmp build/ *.egg-info
python -m build --outdir dist 2>&1 | tail -5
WHEEL=$(ls -t dist/*.whl | head -1)
echo "wheel: $WHEEL"

echo; echo "### 3. wheel audit (no reference data inside)"
python - "$WHEEL" <<'PY'
import sys, zipfile
names = zipfile.ZipFile(sys.argv[1]).namelist()
bad = [n for n in names if "/data/reference_sequences" in n or "/data/database/" in n]
big = [(n, z) for n, z in
       ((i.filename, i.file_size) for i in zipfile.ZipFile(sys.argv[1]).infolist())
       if z > 1_000_000]
print(f"{len(names)} files in the wheel")
if bad:
    print(f"FAIL: reference data inside the wheel ({len(bad)}): {bad[:5]}")
    sys.exit(1)
for n, z in big:
    print(f"  large member: {n} {z/1e6:.1f} MB")
print("PASS: no reference sequence sets or BLAST database in the wheel")
PY
[ $? -ne 0 ] && { echo "### wheel audit FAILED"; exit 1; }

echo; echo "### 4. fresh-environment smoke test"
SMOKE=/scratch/fitzpatria/.cache/smoke_$$
rm -rf "$SMOKE"; python -m venv "$SMOKE" || exit 1
"$SMOKE/bin/pip" install --quiet "$WHEEL" 2>&1 | tail -3
echo "-- version:"; "$SMOKE/bin/calicityper" --version 2>&1 | head -2
echo "-- setup from the archive into a clean data dir:"
DATA="$SMOKE/refdata"
"$SMOKE/bin/calicityper" setup --source "$(readlink -f "$ARCHIVE")" --data-dir "$DATA" 2>&1 | tail -4
echo "-- verify:"
CALICITYPER_DATA="$DATA" "$SMOKE/bin/calicityper" setup --verify-only --data-dir "$DATA" 2>&1 | tail -3
echo "-- sapovirus sets present after install:"
for d in reference_sequences_sapovirus_VP1 reference_sequences_sapovirus_RdRp; do
    n=$(ls "$DATA/$d" 2>/dev/null | wc -l); echo "   $d: $n files"
done
echo "-- type the bundled example (needs blastn/mafft on PATH):"
CALICITYPER_DATA="$DATA" "$SMOKE/bin/calicityper" type examples/example.fasta \
    -o "$SMOKE/out.tsv" --threads 4 2>&1 | tail -5
[ -s "$SMOKE/out.tsv" ] && { echo "   rows: $(($(wc -l < "$SMOKE/out.tsv") - 1))"; head -2 "$SMOKE/out.tsv" | cut -c1-140; }
rm -rf "$SMOKE"
echo; echo "### DONE $(date)"
