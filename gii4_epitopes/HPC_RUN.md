# Running F1–F4 on UCD Sonic

Works from either layout: the repo (`.../calicityper/dev/gii4_epitopes`) or a
flat HPC copy (`~/scratch/calicityper_runs/gii4_epitopes`).

## One-time setup on the login node

```bash
# Activate the environment that has calicityper + biopython + numpy + MAFFT
source ~/scratch/calicityper_env/bin/activate      # adapt to your env path

# Tell the scripts where the tool checkout and verified reference set are
export CALICITYPER_REPO=~/scratch/calicityper/calicityper
export CALICITYPER_DATA=~/.cache/calicityper/reference-0.8.1   # from `calicityper setup --verify-only`
mkdir -p /tmp/calicityper_tmp && export TMPDIR=/tmp/calicityper_tmp
```

If the reference set lives elsewhere, point `CALICITYPER_DATA` at the directory
containing `database/refseq_2.fasta`, `reference_accessions.tsv` and
`reference_sequences_VP1/`.

## Run

```bash
cd ~/scratch/calicityper_runs/gii4_epitopes     # or your copy
bash scripts/run_all.sh
```

`run_all.sh` is tolerant: steps needing MAFFT or the reference set print
`SKIP ...` instead of aborting, so the input-only steps always complete.

## Run one step directly

```bash
python scripts/check_inputs.py                  # no external deps
python scripts/build_sublineage_table.py        # no external deps
python scripts/select_representatives.py        # no external deps
python scripts/build_immunotype_lookup.py       # no external deps (manifest if present)
python scripts/aa_year_heatmap.py inputs/parra_s3/accessions.tsv \
       inputs/parra_s3/vp1_proteins.faa results/heatmaps/
python scripts/validate_variant_thresholds.py   # needs gii4_vp1_aa_279.aln (F1)
python scripts/validate_representatives.py      # needs reference set
python scripts/crosscheck_parra_accessions.py   # needs reference set
python scripts/check_immunotype_consistency.py  # needs tool + MAFFT + reference set
```

The optional checks that need MAFFT or the tool emit a clear one-line message
and exit non-zero if their dependency is missing — set the env vars above or
run them inside the conda environment.

## Verified locally (2026-11-19)

From a **flat copy** with no parent repo and no `PYTHONPATH`:
- `check_inputs.py` — 0 placeholder warnings
- `build_sublineage_table.py`, `select_representatives.py`, `build_immunotype_lookup.py` — pass
- with `CALICITYPER_REPO` set: `validate_representatives.py` — 0/3 distance increases
- with `CALICITYPER_DATA` set: `crosscheck_parra_accessions.py` — writes TSV
- all 16 scripts compile; 7 antigenic tests pass in the repo
