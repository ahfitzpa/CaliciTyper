<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/calicityper-lockup-horizontal-dark.png">
    <img src="assets/calicityper-lockup-horizontal.png" width="420" alt="CaliciTyper">
  </picture>
</p>

<p align="center">Command-line typing for <em>Caliciviridae</em></p>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)
[![CI tested](https://img.shields.io/badge/CI-3.10%20%7C%203.12-blue.svg)](.github/workflows/ci.yml)
[![Status](https://img.shields.io/badge/status-research_code-orange.svg)](#how-to-cite--how-results-were-validated)

# CaliciTyper

CaliciTyper (`calicityper`) assigns norovirus genogroup, ORF1 polymerase (P-) type, and ORF2 capsid type to assembled contigs/genomes (FASTA), and provides evidence-weighted recombination analysis for norovirus/sapovirus alignments.

- **Validated/supported genotyping scope:** norovirus **GI/GII** (ORF1 + ORF2) and broad Caliciviridae genogroup detection.
- **Not validated / not implemented:** non-GI/GII norovirus genotyping calls may be reported but are not validated; sapovirus genotyping is not implemented.
- **Research software:** independent reimplementation (not produced by or affiliated with RIVM).

Detailed reference: [`docs/REFERENCE.md`](docs/REFERENCE.md) • Output schemas/provenance: [`docs/OUTPUT_SCHEMA.md`](docs/OUTPUT_SCHEMA.md)

## Quickstart

### 1) Install (from source, recommended conda/mamba)

```bash
git clone https://github.com/ahfitzpa/calicityper.git
cd calicityper
conda env create -f environment.yml
conda activate calicityper
calicityper setup   # required before first typing use (installs/verifies reference set)
```

Pip/source install is also supported when Python and external tools are already available (`blastn`, `mafft` on `PATH`): `pip install .` (or from Git) then `calicityper setup`.

> Distribution status: source install only at present (not yet published on PyPI/bioconda). `environment.yml` is a convenience environment (not a lock file); for reproducibility pin/export resolved environments.

### 2) Type example sequences

```bash
calicityper type examples/example.fasta -o example_typed.tsv --threads 4 --report example
```

- Primary table: `example_typed.tsv` (TSV). Compare with [`examples/example_typed.tsv`](examples/example_typed.tsv).
- Reports from `--report example`: `example.md`, `example.json`, `example.html` (figures require `bokeh`).
- File TSV outputs also write a provenance sidecar (`example_typed.tsv.provenance.json`).

### 3) Common commands

```bash
# Typing (FASTA contigs/genomes; stdin not supported)
calicityper type contigs.fasta -o typed.tsv --threads 8

# Reports (predictable naming: PREFIX.md/.json/.html)
calicityper type contigs.fasta -o typed.tsv --report run01 --threads 8

# Putative-novel screen (norovirus only, tentative)
calicityper novel sequences.fasta -o novel.tsv

# Recombination (aligned FASTA, or --input to align with MAFFT)
calicityper recomb --alignment panel.fasta --query SAMPLE_01 -o recomb_out/
calicityper recomb --alignment panel.fasta --query SAMPLE_01 --dry-run

# Reference data install/verify (offline/local archive supported)
calicityper setup
calicityper setup --verify-only
calicityper setup --source /path/to/reference.tar.gz --data-dir /shared/calicityper-ref
```

Run `calicityper <command> --help` for full options and examples.

## Input/output notes

- **Input:** FASTA contigs/genomes for `type`/`novel` (not raw FASTQ). `recomb` requires `--alignment` or `--input`. **Stdin is not supported** for sequence input.
- **Output:** `type`/`novel` write TSV to stdout by default (`-o -`) or to a file; status/progress messages go to stderr.
- **`--report PREFIX`:** writes `PREFIX.md`, `PREFIX.json`, `PREFIX.html`.
- **Provenance/auditability:** file-based TSV outputs get `*.provenance.json`; JSON outputs include `tool_version`, effective parameters, and reference-set DOI/version/manifest hash. See [`docs/OUTPUT_SCHEMA.md`](docs/OUTPUT_SCHEMA.md).

## Supported vs unsupported

| Capability | Scope | Status | Notes |
|---|---|---|---|
| Genogroup detection | Caliciviridae | Working | BLAST-based detection against bundled references. |
| Norovirus dual genotyping (ORF1 + ORF2) | GI, GII | **Validated/supported** | Primary supported typing use case. |
| Norovirus genotyping attempts | GIII–GX | **Not validated** | May be reported; restricted by default. |
| GII.4 variant resolution | 11 variants | Supported (within bundled set) | Variant suffix in `orf2_type`. |
| Putative novel ORF1/ORF2 screening | Norovirus | **Tentative screening** | Not a genotype designation; `not_assessable` possible. |
| Recombination analysis | Norovirus, sapovirus | Working with caveats | Alignment-dependent; optional OpenRDP/external tools. |
| Sapovirus genotyping | GI–GXIX | **Not implemented** | Genogroup detection only. |
| Raw reads (FASTQ) | — | Unsupported | Assemble first. |

More detail: [`docs/REFERENCE.md`](docs/REFERENCE.md#supported-vs-unsupported).

## Troubleshooting

- **Reference data missing/mismatch:** run `calicityper setup` (or `--source <manifest-matching tarball>` + `--verify-only`) and set `CALICITYPER_DATA` if using a non-default location.
- **`blastn`/`mafft` not found:** install tools (conda-forge/bioconda recommended) or set `CALICITYPER_BLASTN`/`CALICITYPER_MAFFT`.
- **OpenRDP unavailable (recomb):** install in a separate `numpy<2` environment and set `CALICITYPER_OPENRDP`; use `--dry-run` to verify. Recomb can still run with skipped scanners.
- **`insufficient_sequence` / not assessed (recomb):** input may be too short/ambiguous or have too few informative sites; check dry-run/QC and provenance notes.
- **Slow/large reports on big batches:** use `--report-max-panels N` and/or `--no-tree`; prefer TSV/JSON/provenance outputs for batch processing.
- **Calls outside GI/GII or novel flags:** treat as not validated / tentative screening; inspect `*_accepted`, support/sites, coverage, notes, and provenance.

Full troubleshooting and performance guidance: [`docs/REFERENCE.md`](docs/REFERENCE.md#troubleshooting), [`docs/REFERENCE.md`](docs/REFERENCE.md#performance-and-scalability).

## How to cite & how results were validated

### How to cite
- Cite the CaliciTyper software version/commit/release used and retain run provenance (tool version, parameters, reference-set DOI/version/manifest hash) from machine-readable outputs.
- Cite the method/classification context:
  - Kroneman A, Vennema H, Deforche K, *et al.* An automated genotyping tool for enteroviruses and noroviruses. *J Clin Virol* 2011;51(2):121-125. [doi:10.1016/j.jcv.2011.03.006](https://doi.org/10.1016/j.jcv.2011.03.006)
  - Chhabra P, de Graaf M, Parra GI, *et al.* Updated classification of norovirus genogroups and genotypes. *J Gen Virol* 2019;100(10):1393-1406. [doi:10.1099/jgv.0.001318](https://doi.org/10.1099/jgv.0.001318)
- If using `calicityper recomb`, cite [OpenRDP](https://github.com/PoonLab/OpenRDP) and relevant underlying methods.

```bibtex
@software{fitzpatrick_calicityper,
  author = {Amy Heather Fitzpatrick},
  title  = {CaliciTyper: command-line calicivirus genotyping and recombination analysis},
  url    = {https://github.com/ahfitzpa/calicityper},
  year   = {2026}
}
```

### How results were validated
- **Validated/supported scope:** norovirus GI/GII ORF1/ORF2 genotyping and Caliciviridae genogroup detection within documented method/reference limits.
- **Automated checks:** `pytest tests/ -v` covers core in-repo behavior; CI tests a baseline on Ubuntu with Python 3.10/3.12 (see `.github/workflows/ci.yml`).
- **Extended benchmarking/validation artifacts:** detailed benchmarking outputs, figures, and extended validation materials are maintained separately during active development and will be shared as curated records when available (public repository HEAD no longer tracks those detailed artifacts). Use released/curated records plus recorded run provenance for reproducibility.
- **Audit trail:** record `calicityper setup` reference-set/manifest context and per-run machine-readable provenance (`*.json` / `*.provenance.json`: tool version, parameters, reference-set DOI/version/manifest hash). See [`docs/OUTPUT_SCHEMA.md`](docs/OUTPUT_SCHEMA.md).

## Performance and scalability (summary)

- `--threads N` uses **worker processes** for typing/alignment chunks; CaliciTyper pins BLAS/OpenMP threads per process by default to avoid oversubscription (user-exported values take precedence). `recomb --threads` applies to MAFFT when aligning `--input`.
- For large batches, prefer TSV/JSON outputs, use `--report-max-panels N` and/or `--no-tree`, and ensure adequate RAM/disk for alignments, trees, and report artifacts. See [`docs/REFERENCE.md`](docs/REFERENCE.md#performance-and-scalability).

## Dependencies and reproducible environments

- **Mandatory:** Python `>=3.10` (CI: 3.10/3.12), `biopython`, `numpy`, `pandas`; external `blastn` (BLAST+) and `mafft`.
- **Recommended:** `numba` (bootstrap performance; version must be compatible with Python/NumPy).
- **Optional:** `bokeh` (HTML figures), `iqtree`/`clustalw` (specific workflows), OpenRDP (recomb scanners, separate `numpy<2` env).
- Pin/export resolved environments (e.g. `conda env export`) or use containers for strong reproducibility. OpenRDP workflow and version guidance: [`docs/REFERENCE.md`](docs/REFERENCE.md) and installation notes in this repo.

## Documentation

- [`docs/REFERENCE.md`](docs/REFERENCE.md) — detailed usage, methods, options, coverage/acceptance, recombination, troubleshooting, performance
- [`docs/OUTPUT_SCHEMA.md`](docs/OUTPUT_SCHEMA.md) — TSV/JSON schemas, column semantics, provenance/sidecars, versioning
- [`RELEASING.md`](RELEASING.md), [`CHANGELOG.md`](CHANGELOG.md) — release/versioning
- [`ATTRIBUTION.md`](ATTRIBUTION.md), [`LICENSE`](LICENSE)

## Acknowledgements and licence

Fork and substantial rewrite of [zhuzhanji/NorovirusGenotyping](https://github.com/zhuzhanji/NorovirusGenotyping) (Yang Ji). See [`ATTRIBUTION.md`](ATTRIBUTION.md).

MIT License ([`LICENSE`](LICENSE)). CaliciTyper name/logo notices: [`assets/NOTICE.md`](assets/NOTICE.md).

Author: Amy Fitzpatrick, University College Dublin.