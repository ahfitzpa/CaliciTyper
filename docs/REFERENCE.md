# CaliciTyper detailed reference

This document contains detailed usage, methods, outputs, configuration, and validation/reference information. For installation and a minimal first run, see [`../README.md`](../README.md).

## Contents
- [Supported vs unsupported](#supported-vs-unsupported)
- [Input and output behavior](#input-and-output-behavior)
- [Commands and options](#commands-and-options)
- [Output fields and report files](#output-fields-and-report-files)
- [Method details](#method-details)
- [Genogroup/genotype coverage and acceptance](#genogroupgenotype-coverage-and-acceptance)
- [Recombination analysis details](#recombination-analysis-details)
- [Advanced configuration and reproducibility](#advanced-configuration-and-reproducibility)
- [Reference data and offline/cluster installs](#reference-data-and-offlinecluster-installs)
- [Performance and scalability](#performance-and-scalability)
- [Troubleshooting](#troubleshooting)
- [How to cite and how results were validated](#how-to-cite-and-how-results-were-validated)

---

## Supported vs unsupported

| Capability | Scope | Status | Notes |
|---|---|---|---|
| Genogroup detection | Caliciviridae (including norovirus and non-norovirus genera in bundled DB) | **Working** | BLAST-based detection/assignment against bundled reference database. |
| Norovirus dual genotyping (ORF1 P-type + ORF2 capsid type) | **GI, GII** | **Validated/supported** | Primary supported typing use case. |
| Norovirus genotyping attempts | GIII–GX | **Detected / not validated** | Calls may be returned but are not validated (sparse reference coverage); use with caution. |
| GII.4 variant resolution | 11 variants (bundled set) | **Validated (within bundled references/method limits)** | Variant carried in `orf2_type` (e.g. `GII.4_Sydney`). |
| Putative novel ORF1/ORF2 type screening | Norovirus only (`type` screening fields, `novel`) | **Screening / tentative** | 2×SD screen; flags are not genotype designations. `not_assessable` means required region/sequence basis is insufficient. |
| Recombination detection/scoring | Norovirus, sapovirus (alignment-based) | **Working with caveats / experimental aspects** | Evidence-weighted results; depends on alignment quality, sequence length/informative sites, and optional external tools (OpenRDP). |
| Sapovirus genotyping | GI–GXIX (planned) | **Not implemented** | Detected to genus/genogroup only; taxonomy tables drafted, typing/reference pipeline not implemented. |
| FASTQ / raw-read typing | — | **Unsupported** | Assemble/cluster first; input must be FASTA contigs/genomes. |

**Capability/status terms used in outputs/docs:** validated/supported, detected/not validated, screened/tentative, not implemented, experimental/working with caveats.

---

## Input and output behavior

### Input
- **FASTA input is required** for `type` and `novel` (`fasta` positional argument). CaliciTyper types contigs/genomes, **not raw sequencing reads (FASTQ)**.
- `recomb` requires either `--alignment` (aligned FASTA) or `--input` (unaligned FASTA; CaliciTyper aligns with MAFFT first into the output directory).
- **Stdin is not supported** as a FASTA source for `type`/`novel`/`recomb --input`; provide a file path.

### Output streams and files
- `type` and `novel` write the primary TSV table to **stdout by default** (`-o -`) or to a file with `-o/--out <path>`.
- Scope/status/progress summaries are written to **stderr** (so TSV on stdout remains machine-readable).
- When `type`/`novel` write a TSV **file**, CaliciTyper also writes a machine-readable provenance sidecar: `<output>.provenance.json`.
- `recomb` writes outputs into `-o/--output <dir>` (default `recomb_out`), including `recombination.json`, `recombination.tsv`, `recombination.tsv.provenance.json`, and `work/` artifacts (and `aligned.fasta` if `--input` was used).

### Predictable `--report` naming (`type`)
`--report PREFIX` produces:
- `PREFIX.md` (Markdown report)
- `PREFIX.json` (JSON report with structured runners and provenance metadata)
- `PREFIX.html` (self-contained HTML report; figures require `bokeh`)

You may also request formats individually with `--json`, `--md`, `--html`.

### Auditability
Machine-readable outputs record provenance for reproducibility/auditability (tool version, effective parameters, reference-set DOI/version/manifest hash). See [`OUTPUT_SCHEMA.md`](OUTPUT_SCHEMA.md).

---

## Commands and options

### `calicityper type`
Assign genogroup and norovirus GI/GII ORF1/ORF2 types (validated scope). Other norovirus genogroups/non-norovirus Caliciviridae may be detected but genotyping is not validated; sapovirus genotyping is not implemented.

**Usage:** `calicityper type FASTA [-o OUT] [options]`

Common options:
- `-o, --out PATH|-` TSV output (default stdout)
- `--threads N`, `--seed N`, `--bootstraps N`, `--min-support PCT`, `--max-refs-per-type N`
- Reports: `--report PREFIX`, `--json PATH`, `--md PATH`, `--html PATH`, `--report-max-panels N`, `--no-tree`, `--html-cdn`
- Screening: `--screen-all` (screen every call, not only unaccepted calls)
- Advanced (correctness/reproducibility; prefer CLI flags over env vars): `--min-sites N`, `--novel-distance {ml,nj,p}`, `--type-all-groups/--no-type-all-groups`, `--rivm-window-mode/--no-rivm-window-mode`

**Examples**
```bash
calicityper type contigs.fasta -o typed.tsv --threads 8
calicityper type contigs.fasta -o typed.tsv --report run01 --threads 8
calicityper type contigs.fasta -o typed.tsv --report run01 --report-max-panels 200 --no-tree
```

### `calicityper novel`
Chhabra et al. 2019 2×SD putative-novel screen (norovirus only, tentative). Requires complete VP1 (ORF2, amino-acid basis) and/or ≥95% of the 762 nt 3′ RdRp typing region (ORF1, nucleotide basis); otherwise `not_assessable`.

**Usage:** `calicityper novel FASTA [-o OUT] [--region {orf1,orf2,both}] [options]`

**Examples**
```bash
calicityper novel sequences.fasta -o novel.tsv
calicityper novel sequences.fasta --region orf2 -o novel_vp1.tsv --threads 4
```

### `calicityper recomb`
Detect and score recombination in norovirus or sapovirus from an alignment/panel.

**Usage:** `calicityper recomb (--alignment ALN | --input FASTA) [options]`

Common options:
- `--alignment PATH` (aligned FASTA) or `--input PATH` (unaligned FASTA; MAFFT alignment step performed)
- `--query ID` (repeatable; omit to analyze all sequences in alignment)
- `-o, --output DIR` (default `recomb_out`)
- `--virus {norovirus,sapovirus}`, `--reference-metadata PATH`, `--config PATH`
- Method/threshold controls: `--methods ...`, `--min-length`, `--pvalue-threshold`, `--max-breakpoint-distance`
- Runtime/reporting: `--threads N` (MAFFT step when `--input`), `--skip-rdp`, `--skip-gard`, `--no-segment-trees`, `--dry-run`, `--verbose`

**Examples**
```bash
calicityper recomb --alignment panel.fasta --query SAMPLE_01 -o recomb_out/
calicityper recomb --input samples.fasta -o recomb_out/ --threads 8 --virus norovirus
calicityper recomb --alignment panel.fasta --query SAMPLE_01 --dry-run
```

### `calicityper setup`
Fetch and verify the reference data set required for typing (`type`/`novel`). Refuses to install data that does not match the shipped manifest.

**Usage:** `calicityper setup [--source URL|TARBALL] [--data-dir PATH] [--verify-only] [--force] [--in-package] [--manifest PATH]`

**Examples**
```bash
calicityper setup
calicityper setup --verify-only
calicityper setup --source /path/to/reference.tar.gz --data-dir /shared/calicityper-ref
CALICITYPER_DATA=/shared/calicityper-ref calicityper type contigs.fasta -o typed.tsv
```

---

## Output fields and report files

- **Schemas, column semantics, JSON structures, provenance fields/sidecars, and versioning:** see [`OUTPUT_SCHEMA.md`](OUTPUT_SCHEMA.md).
- **`type` TSV:** one row per query with genogroup, ORF1/ORF2 calls/support/runners/sites/nearest/acceptance/window coverage, putative-novel fields, and notes.
- **`novel` TSV:** one row per screened region per query with flag/closest/bound/distance/basis/detail and context fields.
- **`recomb` outputs:** `recombination.json` (full results + provenance), `recombination.tsv` (event/summary rows), and `recombination.tsv.provenance.json`.

---

## Method details

### Typing pipeline (`type`)
1. **Genogroup detection:** BLAST against bundled whole-genome/genogroup reference database.
2. **Region selection/planning:** regions typed only when BLAST hit overlaps annotated ORF1/ORF2 sufficiently (subject to implementation thresholds and QC).
3. **Alignment:** query cropped/oriented and added to pre-aligned references with MAFFT (`--addfragments --6merpair --keeplength`) so reference columns remain fixed.
4. **Tree and bootstrap:** distance/NJ-based assignment with bootstrap replicates; support and acceptance criteria determine reported/accepted calls.
5. **Putative-novel screening:** 2×SD criterion (Chhabra et al. 2019) applied as a tentative screen (by default focused on unaccepted calls in `type`; configurable).

### Putative-novel screen
- ORF2 (VP1 genotype/genogroup): complete VP1 amino-acid basis as required by the criterion.
- ORF1 (RdRp P-type): partial RdRp nucleotide typing region basis (≥95% coverage requirement in `novel` command).
- Distance measures/bounds are measure-specific (`ml` default, `nj`, `p`); see `--novel-distance` and bundled bounds data.
- Flags are **putative/screening** results, not genotype designations.

### Recombination
- Sanitizes/validates alignment, runs available scanners (OpenRDP when configured/available), clusters breakpoint evidence, infers candidate parents, evaluates segment/phylogenetic support where enabled, scores evidence components, and assigns status/confidence.
- Uses a virus-specific genome model (`--virus norovirus|sapovirus`) for biological plausibility/junction interpretation.
- Can return `insufficient_sequence`/non-recombinant/borderline/recombinant outcomes depending on evidence and QC.

---

## Genogroup/genotype coverage and acceptance

- **Genogroup detection:** broad Caliciviridae detection via BLAST database.
- **Validated typing:** norovirus GI/GII. GIII–GX may return calls but are **not validated** (restricted by default unless `--type-all-groups`).
- **GII.4 variants:** resolved at variant level in bundled set (variant suffix in `orf2_type`).
- **Acceptance vs reported call:** genotype/type may be reported while acceptance is recorded separately (support/sites/window-mode conditions). Check `*_accepted`, `*_support`, `*_sites`, runner-ups, and notes.
- **Typing windows:** window coverages are reported (`orf1_window_cov`, `orf2_window_cov`, `gii4_window_cov`). RIVM-style gating is off by default and can be enabled with `--rivm-window-mode` for comparison/replication.
- **Sapovirus typing:** not implemented (detected to genus/genogroup only).

---

## Recombination analysis details

### Status labels
Typical labels include `none`/`non_recombinant`, `weak`/`possible`/`likely`/`strong` (evidence-weighted), and `insufficient_sequence` or `not_assessed` depending on QC/tool availability and scoring. One method alone is insufficient for strong support; short/ambiguous/low-information inputs should not be forced into confident calls.

### Scoring components (configurable)
Default evidence weighting combines method support, phylogenetic/segment support, breakpoint consistency, parental distinctness, and biological plausibility. See [`../calicityper/recomb_config.py`](../calicityper/recomb_config.py) and `--config`.

### OpenRDP dependency
OpenRDP scanners are external and may be unavailable; use a separate compatible environment and set `CALICITYPER_OPENRDP` (see README/quickstart). `--dry-run` reports tool availability and whether queries would be assessed.

---

## Advanced configuration and reproducibility

Prefer CLI flags over environment variables for correctness-affecting settings (CLI flags take precedence). Effective settings are recorded in machine-readable outputs/provenance.

| Setting | CLI flag | Env fallback | Default | Effect |
|---|---|---|---|---|
| Minimum aligned informative sites | `--min-sites N` | `CALICITYPER_MIN_SITES` | 300 | Below threshold, call reported but not accepted (per acceptance logic). |
| Novel-screen distance measure | `--novel-distance {ml,nj,p}` | `CALICITYPER_NOVEL_DISTANCE` | `ml` | Changes putative-novel distance basis; bounds are measure-specific. |
| Type all norovirus genogroups | `--type-all-groups` / `--no-type-all-groups` | `CALICITYPER_TYPE_ALL_GROUPS=1` | off | Enables GIII–GX typing attempts (not validated). |
| RIVM-style window gating | `--rivm-window-mode` / `--no-rivm-window-mode` | `CALICITYPER_RIVM_WINDOW_MODE=1` | off | For RIVM comparison/replication rather than maximizing correct calls. |

Other key flags: `--bootstraps`, `--min-support`, `--max-refs-per-type`, `--seed`, `--threads`.
External paths/env vars: `CALICITYPER_DATA`, `CALICITYPER_BLASTN`, `CALICITYPER_MAFFT`, `CALICITYPER_CLUSTALW`, `CALICITYPER_OPENRDP`.

---

## Reference data and offline/cluster installs

Typing requires a verified reference set matching the manifest shipped with your installed CaliciTyper version.

```bash
calicityper setup
calicityper setup --verify-only
calicityper setup --source /path/to/reference.tar.gz --data-dir /shared/calicityper-ref
export CALICITYPER_DATA=/shared/calicityper-ref
```

- `calicityper setup` downloads/verifies (or installs from local `--source`) and refuses mismatched data (checksums/manifest mismatch protection).
- **Offline/air-gapped/cluster:** obtain the correct manifest-matching archive out-of-band, install once to shared storage, `--verify-only`, and set `CALICITYPER_DATA` in jobs.
- Reference manifest: `calicityper/data/reference_manifest.json` (version/DOI/source/checksums metadata). Installed set identity/provenance should be recorded with run outputs for auditability.

---

## Performance and scalability

### Threading model
- `--threads N` controls **worker processes** for typing/alignment work (multiprocessing pool over chunks/jobs), not BLAS-level parallelism inside each worker.
- On import, CaliciTyper pins common BLAS/OpenMP thread environment variables to `1` by default (`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `NUMEXPR_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`) via `setdefault`, so a user-exported value takes precedence. This avoids oversubscription when running multiple worker processes.
- `recomb --threads` is used for the MAFFT alignment step when `--input` is provided (unaligned FASTA).

### Memory and disk expectations
- Memory/disk use scales with number and length of input sequences, number of region/genogroup alignment jobs, bootstrap/tree work, and report generation.
- Large HTML reports (embedded Bokeh JS/resources and many alignment panels/trees) can be significantly larger and slower to generate than TSV/JSON-only runs.
- `recomb` writes intermediate/working artifacts under `<output>/work/` and may write an aligned FASTA when starting from `--input`.
- For very large batches, prefer writing TSV/JSON outputs, limit or disable heavy report artifacts, and run on a machine/filesystem with sufficient RAM and scratch space.

### Recommended options for large batches
- Use `--report-max-panels N` to cap the number of alignment panels included in reports.
- Use `--no-tree` to skip per-region trees in HTML reports.
- Omit HTML (or use `--html-cdn` if appropriate) when standalone visual reports are not needed.
- Run with an appropriate `--threads` value for available CPU cores and I/O; validate on a subset before full batch runs.

---

## Troubleshooting

### Reference data missing or verification/install fails
- **Symptom:** `calicityper reference data not found ... Run \`calicityper setup\`` or checksum/manifest mismatch errors.
- **Actions:**
  - Run `calicityper setup` (with network access) or install from a manifest-matching local archive via `--source`.
  - Verify with `calicityper setup --verify-only [--data-dir ...]`.
  - Ensure `CALICITYPER_DATA` points to the intended verified directory in all shells/jobs.
  - Do not mix archives/manifests from different reference-set/tool versions.

### External tool not found (BLAST/MAFFT/CLUSTALW/OpenRDP)
- **Symptom:** errors indicating `blastn`/`mafft`/other executable not found, or OpenRDP reported as unavailable.
- **Actions:**
  - Install required tools (conda/mamba recommended) and ensure they are on `PATH`, or set `CALICITYPER_BLASTN`, `CALICITYPER_MAFFT`, `CALICITYPER_CLUSTALW`, `CALICITYPER_OPENRDP`.
  - For recombination, keep OpenRDP in a separate compatible environment (`numpy<2`) and use `calicityper recomb --dry-run` to verify wiring.

### `recomb` returns `insufficient_sequence` / `not_assessed` or no events
- Check input is an appropriate **alignment/panel** with sufficient length and informative sites (defaults are conservative).
- Confirm query IDs exist in the alignment (`--query` values must match sequence headers).
- Review `recombination.json`/TSV notes and `tools_skipped` in provenance for missing tools or skipped analyses.

### Typing output is unassigned, rejected, or detected-but-not-validated
- Inspect `note`, `*_accepted`, `*_support`, `*_sites`, `*_window_cov`, and runner-up/nearest-reference fields.
- Calls outside GI/GII may be reported but are not validated by default; `--type-all-groups` enables attempts for development/reference work only.
- Putative-novel flags are tentative screens; `not_assessable` indicates required region/basis/coverage was insufficient.

### Large run/report generation is slow or too large
- Use `--report-max-panels N` and/or `--no-tree`, avoid unnecessary HTML, and prefer TSV/JSON/provenance sidecars for batch processing.

---

## How to cite and how results were validated

### How to cite
- Cite CaliciTyper software (repository/release/tag/commit as used) and include the tool version and run parameters/provenance from machine-readable outputs.
- Cite the typing method and norovirus classification standards underpinning the approach:
  - Kroneman A, Vennema H, Deforche K, *et al.* An automated genotyping tool for enteroviruses and noroviruses. *J Clin Virol* 2011;51(2):121-125. doi:10.1016/j.jcv.2011.03.006
  - Chhabra P, de Graaf M, Parra GI, *et al.* Updated classification of norovirus genogroups and genotypes. *J Gen Virol* 2019;100(10):1393-1406. doi:10.1099/jgv.0.001318
- For putative-novel screening, cite Chhabra *et al.* 2019 (2×SD criterion context) and interpret flags as tentative screening results.
- For recombination analysis, cite OpenRDP and the underlying scanner/method publications as appropriate.

### How results were validated
- **Validated/supported scope:** norovirus **GI/GII** ORF1/ORF2 genotyping and Caliciviridae genogroup detection within documented method/reference limits.
- **In-repo automated checks:** `pytest tests/ -v` covers core behavior in this repository (environment-dependent; see CI configuration for tested baselines).
- **Benchmarking/extended validation materials:** detailed benchmarking artifacts, figures, and extended validation workflows are maintained outside the public repository during active development and are intended to be shared as curated records when available. (Historical/backup materials referenced elsewhere may be out of date; rely on released/curated records and recorded run provenance.)
- **Result provenance:** for any reported result, retain machine-readable provenance (tool version, effective parameters, and reference-set DOI/version/manifest hash) and the installed reference-set/manifest context (`calicityper setup --verify-only`, `CALICITYPER_DATA`) to support reproducibility and audit.

See also: [`../README.md`](../README.md) (quickstart/summary), [`OUTPUT_SCHEMA.md`](OUTPUT_SCHEMA.md), [`../RELEASING.md`](../RELEASING.md), and [`../CHANGELOG.md`](../CHANGELOG.md).