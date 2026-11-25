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
- **Sapovirus:** VP1 genotyping and RdRp P-typing per Zhao *et al.* 2024, on by default, reported as a dual type (`GI.4[P1]`) with 95% Wilson CIs on observed p-distances and second-nearest type distance margins. VP1 is reliable; RdRp assigns accurately among known P-types, but **sapovirus RdRp novelty screening is not implemented** (P-types are provisional) — see [`docs/SAPOVIRUS_STATUS.md`](docs/SAPOVIRUS_STATUS.md). Disable with `--no-type-sapovirus` or `CALICITYPER_TYPE_SAPOVIRUS=0`.
- **GII.4 antigenic sites:** `calicityper antigenic` annotates published antigenic-site differences, BLOSUM62 `site_radicality` (a substitution-score summary, not antigenic distance), and emits a sequence-based antigenic-change alert (`sequence_antigenic_alert`). Annotation only — feeds no typing, novelty or recombination call, refuses non-GII input, and must **not** be used to infer immune escape, vaccine relevance or phenotypic antigenicity without neutralisation, blockade or other experimental evidence.
- **Not validated:** non-GI/GII norovirus genotyping calls may be reported (`--type-all-groups`) but are not validated.
- **Research software:** independent reimplementation (not produced by or affiliated with RIVM).

Detailed reference: [`docs/REFERENCE.md`](docs/REFERENCE.md) • Output schemas/provenance: [`docs/OUTPUT_SCHEMA.md`](docs/OUTPUT_SCHEMA.md)

## Quickstart

### 1) Install

**From the GitHub release (recommended).** Each release carries a wheel, an
sdist and the matching reference-set archive.

```bash
pip install https://github.com/ahfitzpa/CaliciTyper/releases/download/v0.8.1/calicityper-0.8.1-py3-none-any.whl
calicityper setup --source calicityper-reference-0.8.1.tar.gz
```

Or straight from Git:

```bash
pip install git+https://github.com/ahfitzpa/CaliciTyper@v0.8.1
calicityper setup
```

**From source, with conda/mamba:**

```bash
git clone https://github.com/ahfitzpa/CaliciTyper.git
cd CaliciTyper
conda env create -f environment.yml
conda activate calicityper
calicityper setup
```

`environment.yml` is a convenience environment, not a lock file. Pin or export
a resolved environment if you need reproducibility, or use the Apptainer image
described in [`docs/INSTALL_TOOLS.md`](docs/INSTALL_TOOLS.md), which fixes
BLAST+, MAFFT and OpenRDP at known versions in a single file you can checksum.

> Not yet on PyPI or bioconda.

### `calicityper setup`

**Required before first use.** The wheel deliberately contains no reference
data; `setup` fetches the reference archive, verifies every artefact against
the manifest's SHA-256 list, and refuses to install a set that does not match,
because a different reference set changes calls silently.

```bash
calicityper setup                       # download, verify, install
calicityper setup --verify-only         # check an existing install
calicityper setup --source ref.tar.gz   # install from a local archive (offline)
calicityper setup --data-dir /shared/calicityper-ref   # shared location
```

Set `CALICITYPER_DATA` to use a non-default location afterwards. For
air-gapped systems, copy the archive across and use `--source`; nothing in
`setup` requires network access when a local archive is given.

External tools (`blastn`, `makeblastdb`, `mafft`, and `openrdp` for
recombination) are found on `PATH` or via `CALICITYPER_*` variables — see
[`docs/INSTALL_TOOLS.md`](docs/INSTALL_TOOLS.md). Note that **OpenRDP cannot
share an environment with CaliciTyper**: it requires numpy<2.

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

# Typing with sapovirus dual typing off (genogroup only)
calicityper type contigs.fasta -o typed.tsv --no-type-sapovirus

# Short dual PolCap amplicons (570/579 bp): per-region aligned-site floors
calicityper type amplicons.fasta -o typed.tsv --preset polcap

# Putative-novel screen (tentative). Norovirus: mean + 2 SD patristic bounds.
# Sapovirus: the VP1 distance cut-off's outcome; RdRp is not screened.
calicityper novel sequences.fasta -o novel.tsv

# Recombination (aligned FASTA, or --input to align with MAFFT).
# --virus picks the genome model and DEFAULTS TO norovirus: pass
# --virus sapovirus for sapovirus, otherwise junctions are reported against
# the norovirus ORF1/ORF2 boundary rather than the sapovirus RdRp/VP1 one.
calicityper recomb --alignment panel.fasta --query SAMPLE_01 -o recomb_out/
calicityper recomb --alignment panel.fasta --query SAMPLE_01 --virus sapovirus -o recomb_out/
calicityper recomb --alignment panel.fasta --query SAMPLE_01 --dry-run

# GII.4 antigenic-site annotation. Norovirus GII VP1 only (complete ORF2);
# --unaligned translates and aligns against the bundled GII.4 variants.
calicityper antigenic vp1.faa -o antigenic.tsv
calicityper antigenic orf2_cds.fasta --unaligned -o antigenic.tsv --threads 8
calicityper antigenic --show-sites

# Reference data install/verify (offline/local archive supported)
calicityper setup
calicityper setup --verify-only
calicityper setup --source /path/to/reference.tar.gz --data-dir /shared/calicityper-ref
```

Run `calicityper <command> --help` for full options and examples.

## How it works

A query goes through one pipeline per genus, and both end at a call or an
explicit refusal to call.

1. **Genogroup.** `blastn` against the bundled genogroup database decides
   whether the query is norovirus, sapovirus or another *Caliciviridae*, and
   which genogroup. Everything after this is per region: ORF1 (polymerase) and
   ORF2 (capsid).
2. **Norovirus — phylogenetic typing.** The query is aligned with MAFFT onto
   the genogroup's pre-aligned reference set for that region, then placed by
   neighbour joining with 100 bootstrap replicates. The nearest reference's
   type is the call; the bootstrap value is reported beside it and gates
   acceptance at **70%**. Queries overlapping fewer than **300** alignment
   columns are not called. GII.4 variants need sequence from VP1 residue 279
   to the C-terminus, so a short amplicon gets a genotype but no variant.
3. **Sapovirus — distance typing.** Uncorrected nucleotide p-distance (Zhao
   *et al.* 2024) over shared unambiguous columns, with gaps and ambiguities
   deleted pairwise. Cut-offs: VP1 `0.161` genotype and `0.503` genogroup,
   RdRp `0.266` P-type and `0.531` P-group. Beyond the cut-off the query is
   **`untypeable`**, never assigned to its nearest neighbour. Bootstrap support
   is reported but does not decide the call.
4. **Putative-novel screening.** For norovirus, a query's mean patristic
   distance to each reference cluster is compared against that cluster's
   `mean + 2 SD` bound (`calicityper/data/putative_novel_bounds.tsv`, ML
   patristic by default, IQ-TREE with a pinned seed). For sapovirus it is the
   cut-off outcome from step 3 — beyond `0.161` but within `0.503` is a
   possible novel genotype. Sapovirus RdRp is not screened: the `0.266`
   cut-off overlaps the between-type distances.
5. **Recombination.** `calicityper recomb` runs the OpenRDP scanners over an
   alignment and ranks events by categorical evidence rules — agreement across
   `>= 3` scanners plus distinct candidate parents — rather than by a composite
   score. `--virus` selects the junction model.
6. **Provenance.** Every file output carries the tool version, parameters,
   seeds and the SHA-256 digest of the reference manifest, which
   `calicityper setup` verifies on installation.

Design rationale and the measurements behind each threshold:
[`docs/REFERENCE.md`](docs/REFERENCE.md).

## Input/output notes

- **Input:** FASTA contigs/genomes for `type`/`novel` (not raw FASTQ). `recomb` requires `--alignment` or `--input`. `antigenic` takes an aligned VP1 amino-acid FASTA, or, with `--unaligned`, unaligned ORF2 nucleotide or VP1 amino-acid sequence; either way it must be a **norovirus GII complete VP1** and anything else is refused, because the site coordinates are GII VP1 numbering resolved through a GII.4 anchor. **Stdin is not supported** for sequence input.
- **Output:** `type`/`novel` write TSV to stdout by default (`-o -`) or to a file; status/progress messages go to stderr.
- **`--report PREFIX`:** writes `PREFIX.md`, `PREFIX.json`, `PREFIX.html`.
- **Provenance/auditability:** file-based TSV outputs get `*.provenance.json` written atomically; JSON outputs record `tool_version`, `runtime` (Python + dependency versions), `input_sha256`, `external_tools` versions, effective parameters, and reference-set DOI/version/manifest hash. See [`docs/OUTPUT_SCHEMA.md`](docs/OUTPUT_SCHEMA.md).
- **Exit codes:** `0` success; `1` runtime/external-tool/download failure; `2` CLI/config/input validation or missing-file error; `130` interrupted (`SIGINT`); `143` terminated (`SIGTERM`).

## Supported vs unsupported

| Capability | Scope | Status | Notes |
|---|---|---|---|
| Genogroup detection | Caliciviridae | Working | BLAST-based detection against bundled references. |
| Norovirus dual genotyping (ORF1 + ORF2) | GI, GII | **Validated/supported** | Primary supported typing use case. |
| Norovirus genotyping attempts | GIII–GX | **Not validated** | May be reported; restricted by default. |
| GII.4 variant resolution | 11 variants | Supported (within bundled set) | Variant suffix in `orf2_type`. |
| Putative novel ORF1/ORF2 screening | Norovirus only | **Tentative screening** | Not a genotype designation; `not_assessable` possible. Sapovirus is declined with a reason: the published RdRp cut-off overlaps between-type distances, so only 1 of 20 testable P-types is separated from every other P-type in its genogroup. |
| Recombination analysis | Norovirus, sapovirus | Working with caveats | Alignment-dependent; optional OpenRDP/external tools. **`--virus` defaults to `norovirus`** and sets the junction model, so pass `--virus sapovirus` for sapovirus. |
| Sapovirus dual typing (VP1 genotype + RdRp P-type) | Multi-genogroup, Zhao et al. 2024 scheme | **On by default**, benchmarked | Nucleotide p-distance cut-offs (VP1 `0.161`/`0.503`, RdRp `0.266`/`0.531`); beyond the cut-off a query is `untypeable`, never assigned to its nearest neighbour. `--no-type-sapovirus` reports genogroup only. RdRp P-types are provisional: the cut-off cannot flag a novel P-type. |
| GII.4 antigenic-site annotation (`antigenic`) | Norovirus GII VP1 | Annotation only | Reports which published antigenic-site residues differ from the nearest GII.4 variant. Needs a complete VP1; non-GII input is refused (`--skip-input-check` overrides). Feeds no genotype, P-type, novelty or recombination call. |
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
- Measured on one 16-core node (UCD Sonic, Xeon Gold 6542Y): **5,986 sequences in 20.8 min at `--threads 16`, 416 MB peak RSS**; 100 sequences in 23.5 s at 8 threads; a single sequence in 2.4 s, dominated by process start-up. Determinism is verified at every release: `--threads 1` and `--threads 8` give byte-identical output.
- Speed-ups available, all opt-in and off by default so results stay bit-reproducible:
  - `--fast-boot` stops bootstrapping at 25 of 100 replicates when support is already unanimous or hopeless. Measured gain is only ~1% on a full panel, and it changed a handful of calls, so it is not recommended for reportable runs.
  - `--fast-novel` reuses cached ML reference backbone trees during putative-novel screening.
  - `--preset polcap` (`--min-sites-orf1 170 --min-sites-orf2 250`) accepts dual `Genotype[P-type]` calls on 570/579 bp PolCap amplicons, where the polymerase half is shorter than the default 300-column floor. `--min-sites-orf1/--min-sites-orf2` set the two floors directly.
- For large batches, prefer TSV/JSON outputs, use `--report-max-panels N` and/or `--no-tree`, and ensure adequate RAM/disk for alignments, trees, and report artifacts. See [`docs/REFERENCE.md`](docs/REFERENCE.md#performance-and-scalability).

## Dependencies and reproducible environments

- **Mandatory:** Python `>=3.10` (CI: 3.10/3.12), `biopython>=1.81,<2`, `numpy>=2.0,<3`, `pandas>=2.0,<4`; external `blastn` (BLAST+) and `mafft`.
- **Recommended:** `numba>=0.59,<1` (bootstrap performance).
- **Optional:** `bokeh>=3.1,<4` (HTML figures), `iqtree`/`clustalw` (specific workflows), OpenRDP (recomb scanners, separate `numpy<2` env).
- Dependencies carry an upper bound at the next major, so a new upstream release cannot change behaviour without a deliberate pin change. Raise a bound only after CI has run against that major.
- Three levels of reproducibility, strongest last:
  - [`environment.yml`](environment.yml) — bounded spec, for a normal install. Resolves fresh, so versions differ between machines.
  - [`environment.lock.yml`](environment.lock.yml) — one solved, **test-verified** set (linux-64, Python 3.12): every version pinned, written only after `pytest tests/` passed in that exact environment. `conda env create -n calicityper-locked -f environment.lock.yml`.
  - [`containers/calicityper-tools.def`](containers/calicityper-tools.def) and [`containers/calicityper.def`](containers/calicityper.def) — Apptainer images (tools-only and all-in-one), isolating OpenRDP's `numpy<2` requirement from CaliciTyper's `numpy>=2`.
- CI deliberately resolves 3.10 and 3.12 fresh rather than installing the lock, so upstream breakage still surfaces; the upper bounds are what keep that from being noise. OpenRDP workflow and version guidance: [`docs/REFERENCE.md`](docs/REFERENCE.md).

## Codebase map

| Path | Role |
|---|---|
| `calicityper/cli.py` | CLI entry point (`type`, `novel`, `recomb`, `antigenic`, `setup`) and multiprocessing orchestration |
| `calicityper/config.py` | Reference paths, container `/data` default, thresholds, and external-tool resolution |
| `calicityper/align.py` | BLAST genogrouping and MAFFT profile alignment onto fixed reference columns |
| `calicityper/typing.py`, `tree.py`, `patristic.py` | Norovirus NJ/bootstrap typing, sapovirus p-distance typing, and patristic distance |
| `calicityper/putative_novel.py` | Norovirus `mean + 2 SD` patristic novelty screening |
| `calicityper/recomb.py`, `recomb_config.py` | OpenRDP wrapper, clustering, and categorical evidence-level rules |
| `calicityper/antigenic.py` | Norovirus GII.4 VP1 antigenic-site annotation and BLOSUM62 radicality |
| `calicityper/report.py`, `provenance.py`, `setup_data.py` | Markdown/JSON/HTML reporting, provenance sidecars, and verified Zenodo data setup |
| `containers/`, `scripts/`, `tests/`, `examples/` | Apptainer recipes, reference-build utilities, `pytest` suite, and worked examples |

## Documentation

- [`docs/REFERENCE.md`](docs/REFERENCE.md) — detailed usage, methods, options, coverage/acceptance, recombination, troubleshooting, performance
- [`docs/OUTPUT_SCHEMA.md`](docs/OUTPUT_SCHEMA.md) — TSV/JSON schemas, column semantics, provenance/sidecars, versioning
- [`RELEASING.md`](RELEASING.md), [`CHANGELOG.md`](CHANGELOG.md) — release/versioning
- [`ATTRIBUTION.md`](ATTRIBUTION.md), [`LICENSE`](LICENSE)

## Acknowledgements and licence

Fork and substantial rewrite of [zhuzhanji/NorovirusGenotyping](https://github.com/zhuzhanji/NorovirusGenotyping) (Yang Ji). See [`ATTRIBUTION.md`](ATTRIBUTION.md).

MIT License ([`LICENSE`](LICENSE)). CaliciTyper name/logo notices: [`assets/NOTICE.md`](assets/NOTICE.md).

Author: Amy Fitzpatrick, University College Dublin.