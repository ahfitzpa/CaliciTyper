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
| Putative novel ORF1/ORF2 type screening | Norovirus ORF1/ORF2 + sapovirus VP1 (`type` screening fields, `novel`) | **Screening / tentative** | Norovirus 2×SD patristic screen and sapovirus VP1 cut-off novelty; flags are not genotype designations. **Sapovirus RdRp novelty screening is not implemented** (`not_assessable`). |
| Recombination detection/scoring | Norovirus, sapovirus (alignment-based) | **Working with caveats** | Evidence-weighted results across 7 OpenRDP scanners; depends on alignment quality, sequence length/informative sites, and external OpenRDP (`numpy<2` env or Apptainer container). |
| Sapovirus dual genotyping (`VP1` + `RdRp`) | Multi-genogroup (Zhao et al. 2024 scheme) | **On by default, benchmarked** | Nucleotide p-distance cut-offs (`VP1 < 0.161/0.503`; `RdRp < 0.266/0.531`) with 95% Wilson CI and second-type distance margin. RdRp novelty screening is not implemented, so P-types are provisional. |
| GII.4 antigenic-site profiling (`antigenic`) | Norovirus GII VP1 | **Annotation only** | Reports site substitutions, BLOSUM62 `site_radicality` (substitution-score summary, not antigenic distance), and `sequence_antigenic_alert`. Must NOT be used to infer immune escape or vaccine relevance without experimental evidence. |
| FASTQ / raw-read typing | — | **Unsupported** | Assemble/cluster first; input must be FASTA contigs/genomes. |

**Capability/status terms used in outputs/docs:** validated/supported, detected/not validated, screened/tentative, annotation only, not implemented, working with caveats.

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
Assign genogroup and norovirus GI/GII ORF1/ORF2 types (validated scope). Other norovirus genogroups/non-norovirus Caliciviridae may be detected but genotyping is not validated; sapovirus VP1 genotyping and RdRp P-typing are implemented and enabled by default, reported as a dual type; RdRp cannot reliably flag a novel P-type, so P-types are provisional (see docs/SAPOVIRUS_STATUS.md).

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
- **Sapovirus typing:** VP1 genotyping and RdRp P-typing (Zhao et al. 2024 distance cut-offs) are enabled by default (`--no-type-sapovirus` reverts to genogroup only). Each sapovirus call reports observed p-distance, its 95% Wilson CI (`orf{1,2}_distance_ci`), and the nearest competing type and margin (`orf{1,2}_second_type`, `orf{1,2}_second_distance`, `orf{1,2}_distance_margin`). Sapovirus RdRp novelty screening is **not implemented** (`not_assessable`), so RdRp P-types should be treated as provisional.

## Expected behaviour on edge cases

| Edge case | Behaviour in CaliciTyper | Recommended action |
|---|---|---|
| **Mixed infections / coinfections** | CaliciTyper types assembled consensus FASTA contigs, not raw reads. A chimeric or consensus-collapsed mixed infection may show elevated IUPAC ambiguity codes (masked pairwise, lowering comparable sites and bootstrap support) or conflicting `ORF1`/`ORF2` calls (`dual_type_discordant`). | Deconvolve haplotypes or run read-level classifiers (e.g. `NorotyperII`, `MGtree`) upstream before passing consensus contigs to CaliciTyper. |
| **Recombinant genomes** | `type` classifies `ORF1` (RdRp) and `ORF2` (VP1) independently and emits dual notation (`Genotype[P-type]`), flagging inter-genogroup incongruence in `dual_type_discordant`. `recomb` tests aligned whole/near-complete genomes (`>=5,000 nt`) across 7 OpenRDP scanners and genus-specific hotspot models. | Verify discordant dual types or suspected intra-genogroup recombinants with `calicityper recomb` and independent flanking-region phylogenies. |
| **Low-coverage / short amplicons** | Inputs `<100` unambiguous nt are rejected before BLAST. For norovirus, alignments with `<100` informative columns return `unassigned`; `100–299` columns emit the best candidate marked `rejected` under default `--min-sites 300` (use `--preset polcap` for `>=170`/`250` nt dual PolCap or Region A–D amplicons). For sapovirus, no default coverage floor is enforced (`>=100` aligned sites and `>=50` comparable bases still required), and `orf{1,2}_distance_ci` widens automatically to flag cut-off uncertainty. | Inspect `*_sites`, `*_accepted`, `*_distance_ci`, and `gii4_window_cov`; do not call GII.4 variants from Region C alone when `gii4_window_cov` is low. |
| **Sequences containing ambiguous bases (IUPAC `N`, `R`, `Y`, etc.)** | Non-`ACGT` IUPAC codes are accepted by input validation (non-IUPAC characters are rejected cleanly without crashing other batch entries), mapped to a masked state, and excluded pairwise from TN93 distances, p-distances, and identity denominators. In `recomb`, IUPAC codes are masked to `N` and recorded in provenance; queries exceeding `5%` ambiguous bases return `insufficient_sequence`. | Check `*_sites` (effective unambiguous overlap) and `recomb` `qc.ambiguous_fraction`; high ambiguity reduces statistical support rather than forcing erroneous calls. |

---

## Recombination analysis details

### Status labels
Typical labels include `none`, `weak`/`possible`/`likely`/`strong` (evidence-weighted status-first monotonic ladder), and `insufficient_sequence`, `not_assessed`, or `failed` depending on QC/tool availability. One method alone never reaches `likely` or `strong`; short (`<5,000 nt`), low-information (`<500` parsimony-informative sites), or ambiguous (`>5%` N) inputs return `insufficient_sequence` rather than forced calls.

### Scoring components (configurable)
Default evidence weighting combines method support, phylogenetic/segment support, breakpoint consistency, parental distinctness, and biological plausibility. See [`../calicityper/recomb_config.py`](../calicityper/recomb_config.py) and `--config`.

### OpenRDP dependency
OpenRDP scanners are external and require a separate `numpy<2` environment (or the Apptainer image `containers/calicityper-tools.def`); set `CALICITYPER_OPENRDP` to the executable. `--dry-run` reports tool availability and whether queries would be assessed.

---

## Advanced configuration and reproducibility

Prefer CLI flags over environment variables for correctness-affecting settings (CLI flags take precedence). Effective settings are recorded in machine-readable outputs/provenance.

| Setting | CLI flag | Env fallback | Default | Effect |
|---|---|---|---|---|
| Minimum aligned informative sites | `--min-sites N` | `CALICITYPER_MIN_SITES` | 300 | Below threshold, call reported but not accepted (per acceptance logic). |
| Per-region site floor / preset | `--preset polcap`, `--min-sites-orf1`, `--min-sites-orf2` | `CALICITYPER_MIN_SITES_ORF1/2` | unset | Sets region-specific floors (`polcap` = `170`/`250`). |
| Novel-screen distance measure | `--novel-distance {ml,nj,p}` | `CALICITYPER_NOVEL_DISTANCE` | `ml` | Changes putative-novel distance basis; bounds are measure-specific. |
| Type all norovirus genogroups | `--type-all-groups` / `--no-type-all-groups` | `CALICITYPER_TYPE_ALL_GROUPS=1` | off | Enables GIII–GX typing attempts (not validated). |
| Sapovirus dual typing | `--type-sapovirus` / `--no-type-sapovirus` | `CALICITYPER_TYPE_SAPOVIRUS=1` | on | Enables Zhao et al. 2024 sapovirus VP1/RdRp dual typing. |
| RIVM-style window gating | `--rivm-window-mode` / `--no-rivm-window-mode` | `CALICITYPER_RIVM_WINDOW_MODE=1` | off | For RIVM comparison/replication rather than maximizing correct calls. |

Other key flags: `--bootstraps`, `--min-support`, `--max-refs-per-type`, `--seed`, `--threads`.
External paths/env vars: `CALICITYPER_DATA`, `CALICITYPER_BLASTN`, `CALICITYPER_MAFFT`, `CALICITYPER_CLUSTALW`, `CALICITYPER_OPENRDP`.

### Determinism definition
"Byte-identical determinism" means identical SHA-256 digests of TSV output files (and identical JSON result fields excluding UTC `generated` timestamps) under fixed seeds (`--seed 1` for `type`/`novel`, `seed=3` for `recomb`), pinned reference bundle (`v0.8.1`), across thread counts (`--threads 1` vs `8`), and across independently rebuilt Linux x86_64 environments (`Python 3.12`/`3.15`, `NumPy 2.x`). Cross-architecture BLAS implementations or differing upstream OpenRDP builds may exhibit minor floating-point or scanner-version differences.

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

## Performance, hardware requirements and scalability

### Typical hardware requirements and effect of optional dependencies
- **Baseline (`type` / `antigenic` TSV/JSON):** runs comfortably on a standard laptop or workstation (`>=2` CPU cores, `>=1 GB` RAM, `~150 MB` disk for the verified reference bundle; measured peak RSS `<430 MB` across `5,986` sequences at `16` threads).
- **`numba` (optional, recommended):** accelerates NJ bootstrap clade evaluation `~3–5x` with minor one-time JIT disk-cache warmup (`~30–50 MB` RAM).
- **`iqtree` (optional, for `novel` ML distances):** increases per-query CPU time for ML patristic tree inference (`~50–150 MB` additional RAM); without IQ-TREE, `novel` falls back to NJ patristic distances (`--novel-distance nj`).
- **`bokeh` (optional, for `--html` reports):** embedding BokehJS and hundreds of alignment panels/trees increases RAM (`>=2–4 GB` recommended for large HTML batches) and HTML file size; use `--report-max-panels` and `--no-tree` on large runs.
- **`openrdp` (optional, for `recomb`):** requires an isolated `numpy<2` environment or the Apptainer container (`~335 MB` image) and dominates runtime (`~20–50 min` CPU per 9-genome `7.2 kb` alignment).

### Threading model
- `--threads N` controls **worker processes** for typing/alignment work (multiprocessing pool over chunks/jobs), not BLAS-level parallelism inside each worker.
- On import, CaliciTyper pins common BLAS/OpenMP thread environment variables to `1` by default (`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `NUMEXPR_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`) via `setdefault`, so a user-exported value takes precedence. This avoids oversubscription when running multiple worker processes.
- `recomb --threads` is used for the MAFFT alignment step when `--input` is provided (unaligned FASTA).

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