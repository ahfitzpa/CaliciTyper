# Analysis provenance

Everything needed to trace a number in the manuscript back to the code that
produced it and the file it came out of. Each row below names the claim, the
script that generated it, and the output deposited here.

All analyses ran on the UCD "Sonic" HPC cluster through SLURM. The batch
scripts in `slurm/` are the ones that were actually submitted, unedited.

```
analysis/
  slurm/                     the SLURM batch scripts, as submitted
  benchmark_outputs/         the scored output tables the figures and text use
  make_figures.py            generates every manuscript figure from those tables
  calicityper.bib            the manuscript bibliography
```

The tool's own `results/` directory is working output and is not tracked; the
curated tables the manuscript depends on are copied here instead.

## Manuscript claim to output file

| Claim in the manuscript | Script | Output |
|---|---|---|
| 98.1% genogroup / 99.6% ORF2 / 99.8% ORF1 agreement with the Norovirus Typing Tool | `scripts/validate_vs_rivm.py` via `slurm/rivm_benchmark_070.sh` | `benchmark_outputs/norovirus_0.7.0/rivm_comparison_all_0.7.0.csv` (+ the four `*_discordant.csv`) |
| Rebuilt-environment comparison, 68 genogroup changes | `slurm/rivm_benchmark_dev.sh` | `benchmark_outputs/rivm_comparison_all_dev.csv` (+ `*_discordant.csv`) |
| Norovirus leave-one-out, 267/270 (98.9%) | `scripts/truncation_loo.py` via `slurm/truncation_loo.sh` | `benchmark_outputs/norovirus_0.7.0/loo_0.7.0.csv` — the 270 rows with `truncated=0`; `sum(call==truth)` is 267. Run record: `benchmark_outputs/loo_default.log` |
| Per-type concordance (Fig. 3) | `scripts/validate_vs_rivm.py` | the `*_discordant.csv` sets above and `benchmark_outputs/benchmark_sapovirus/rivm_comparison.csv` |
| Pairwise distance distributions (Fig. 4) | `make_figures.py`, computed from the reference alignments in `calicityper/data/` | figure only; no intermediate table |
| Published sapovirus ranges: VP1 31/32, RdRp 19/20 | `scripts/benchmark_sapovirus.py` via `slurm/rebuild_sapovirus.sh` | `benchmark_outputs/benchmark_sapovirus/distance_ranges.csv` |
| Sapovirus vs the Norovirus Typing Tool, 265/286 and 198/207 | `scripts/benchmark_sapovirus.py` | `benchmark_outputs/benchmark_sapovirus/rivm_comparison.csv`, `rivm_sapovirus_typed.tsv` |
| Truncation: VP1 82.9–87.2%, RdRp 91.2–93.1%; precision 1.000 at every floor | `scripts/truncation_loo_sapovirus.py` via `slurm/truncation_loo_sapovirus.sh` | `benchmark_outputs/truncation_loo_sapovirus.csv`, `truncation_loo_sapovirus.floors.csv` |
| Amplicon benchmark: 98.6% accuracy, 917/918 accepted correct | `scripts/benchmark_zenodo_amplicons.py` via `slurm/benchmark_zenodo.sh` | `benchmark_outputs/benchmark_zenodo/zenodo_comparison.csv`, `zenodo_per_genotype.csv`, `zenodo_typed.tsv` |
| Recombination, Zhao Table 3: 7 strong / 4 possible / 1 weak, control possible | `scripts/validate_recombination_table3.py` via `slurm/recomb_table3.sh`, `slurm/rerun_recomb_failures.sh` | `benchmark_outputs/recomb_table3/table3_summary.csv`, plus per-strain `result.json`, `panel_hits.tsv` and `recomb/recombination.{tsv,json}` |
| Acceptance threshold kept at 70% | `scripts/calibrate_thresholds.py` via `slurm/calibrate.sh` | `benchmark_outputs/calibration/calibration_sweep.csv` |
| GI.P5 is two P-types; GI.P6 restored | `scripts/analyse_gi5_ptype.py` via `slurm/analyse_gi5.sh` | `benchmark_outputs/gi5_ptype/gi5_clusters.json`, `gi5_sequences.csv` |
| Cross-check against both published Zhao tables | `scripts/check_against_zhao_tables.py` via `slurm/check_zhao_tables.sh` | printed report; see `docs/SAPOVIRUS_STATUS.md` for the adjudicated deviations |
| Temporal hold-out of the novelty screen | `scripts/temporal_holdout_novelty.py` via `slurm/temporal_holdout.sh` | `benchmark_outputs/temporal_holdout*/temporal_holdout.csv` (four cut-offs) |
| Leave-one-genotype-out novelty screen: ORF2 sensitivity 0.950, ORF1 0.639 | `scripts/logo_novelty.py` via `slurm/logo_novelty.sh` | `benchmark_outputs/logo_novelty.csv` (17,238 rows); run record `benchmark_outputs/logo_novelty.log` |
| Our parse of Zhao Supplemental Table 1 against an independent one | `scripts/third_party/convert_sapovirus_table.py` | `validation/zhao_supplement_independent_parse.csv`; asserted by `tests/test_zhao_independent_parse.py` |
| Synthetic recombination datasets, 9/9 expected status | `tests/test_recomb.py` (`needs_openrdp`, `slow`) | test assertions; worked examples in `examples/recombination/` |

## Reproducing

Order matters: the reference sets must exist before anything is typed.

```bash
calicityper setup                                   # reference data
sbatch analysis/slurm/fetch_zenodo_benchmark.sh     # amplicon benchmark data
sbatch analysis/slurm/rebuild_sapovirus.sh "reason" # tables, refsets, benchmarks
sbatch analysis/slurm/benchmark_zenodo.sh
sbatch analysis/slurm/recomb_table3.sh              # needs OpenRDP
```

Paths in the batch scripts are absolute and site-specific
(`/home/people/fitzpatria/scratch/...`). They are deposited unedited, as run,
rather than generalised after the fact; adapt the paths for another system.

`make_figures.py` reads from `~/scratch/calicityper` and
`~/scratch/calicityper_runs` for the same reason. Point `REPO` and `BENCH070`
at `analysis/benchmark_outputs/` to regenerate the figures from this deposit.

## Known reproducibility caveats

Stated plainly, because they affect whether these numbers can be reproduced
exactly rather than approximately.

- **OpenRDP is not on PyPI.** `pip install openrdp` fails outright. It is
  installed from `github.com/PoonLab/OpenRDP`, and the environment used for the
  results here was built from a local checkout that no longer exists, so
  byte-identical reproduction cannot be guaranteed. Its recorded metadata
  (`numpy<2.0.0,>=1.17.4`) matches commit `b42d9e5` on master, not the `v0.1.0`
  tag, which does not cap numpy and therefore resolves numpy 2.x and breaks
  OpenRDP's isolation. `containers/calicityper-tools.def` pins `b42d9e5` and
  asserts at build time that numpy resolves below 2. Anyone reproducing the
  recombination analysis should expect the same statuses but should verify
  rather than assume.
- **OpenRDP must be present or the recombination tests are skipped.** Without
  it, 13 tests are deselected and a result is reported as `not_assessed`, which
  is not the same as "no recombination found". See `docs/INSTALL_TOOLS.md`.
- **The reference set version is part of the result.** Every output carries a
  `*.provenance.json` recording the tool version, parameters and the reference
  manifest hash. A different reference set changes calls, which is why
  `calicityper setup` refuses to install one that does not match its manifest.
- **The 0.7.0 norovirus figures were produced before the sapovirus corrections.**
  They are deposited under `benchmark_outputs/norovirus_0.7.0/` and are
  unaffected by those corrections, which touch only sapovirus tables. The
  manuscript says where each number comes from.
