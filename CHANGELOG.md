# Changelog

All notable changes to this project should be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project aims to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
for released versions.

## [Unreleased]

## [0.8.0] - 2026-11-17

### Added
- Sapovirus genotyping, on by default (disable with `CALICITYPER_TYPE_SAPOVIRUS=0`). VP1 genotype in the `orf2_*` columns and RdRp P-type in `orf1_*` (both regions lie in sapovirus ORF1). Calls follow the Zhao et al. 2024 nucleotide p-distance cut-offs: within the genotype cut-off of the nearest reference the call is its type; otherwise `untypeable` (with "possible novel genotype" inside the genogroup cut-off). Bootstrap support is reported beside the call, not used to decide it. New `orf1_distance`/`orf2_distance` columns. Coverage of the typed region is reported but not gated on: the truncation benchmark showed VP1 precision of 1.000 at every coverage floor from 0.0 to 0.95, so a floor only withholds correct calls. A genogroup with no reference set for a region (e.g. GIV RdRp) says so instead of borrowing another genogroup's.
- Sapovirus genogroup database: 39 subjects over 35 genogroups (was 12/12).
- First-class documentation for offline, air-gapped, and cluster reference-data installs (`--source`, `--data-dir`, `--verify-only`, `CALICITYPER_DATA`).
- `tests/test_setup_data.py`: local-source/fixture validation of reference-set install, verification, and checksum-mismatch rejection.
- Release/versioning and packaging guidance (`RELEASING.md`, changelog, README/distribution notes).
- Formalized output schema/provenance documentation (`docs/OUTPUT_SCHEMA.md`) covering `type`/`novel`/`recomb` TSV and JSON structures, column semantics, provenance fields, and schema versioning.
- Machine-readable provenance for auditability: `type` JSON/provenance sidecars, `novel`/`type` TSV sidecars (`*.provenance.json`), and enriched `recomb` JSON/provenance sidecar records including tool version, effective parameters, and shipped reference-set DOI/version/manifest hash.
- Detailed user/reference documentation (`docs/REFERENCE.md`) with supported-vs-unsupported matrix, troubleshooting, performance/scalability guidance, and citation/validation guidance.

- Sapovirus dual VP1/RdRp typing is on by default, reported as `GI.4[P1]`, and as `GIV.1[GII.P10]` where the P-group differs from the capsid genogroup (`dual_type`, `dual_type_discordant`).
- `--openrdp-timeout` for the `recomb` command. A failed combined OpenRDP run is retried method by method, so one crashing scanner no longer discards the evidence from the others.
- `scripts/check_against_zhao_tables.py`, cross-checking the generated sapovirus tables against Table 2 and the pig/bat table of Zhao et al. 2024, with an explicit list of adjudicated deviations.
### Changed
- Benchmarking artefacts are tracked again, under `analysis/`: the 32 SLURM batch scripts as submitted, the scored output tables behind every figure and reported number, the figure generator, and a README mapping each claim in the paper to the script and file that produced it. `BENCHMARKING.md` is no longer gitignored.
- README restructured into a quickstart-focused entry point with links to detailed reference and schema docs.
- CLI/UX polish: consistent subcommand help/descriptions and examples, clearer input/output behavior (stdin not supported; stdout TSV and stderr status/progress), predictable `--report` naming (`PREFIX.md/.json/.html`), improved progress messages, and clearer actionable errors/hints for missing reference data or external tools.

- Sapovirus minimum region coverage is 0 by default: coverage is reported, not gated on.
- The sapovirus reference sets are now included in the reference manifest and checksummed (previously skipped).
- `recomb` reports `not_assessed` when no scanner ran and `failed` when a query raised, instead of `none`. A tool failure is no longer indistinguishable from "no recombination found".
### Fixed
- Genogroup BLAST uses a fixed search-space size, so growing the database no longer moves calls (MN567957.1 GII.11 flipped to GIX when sapovirus subjects were added).
- A BLAST hit shorter than 100 nt no longer names a genogroup; on the RIVM set every such call was wrong (non-calicivirus records called norovirus GVI/GX).
- Sapovirus alignments were written with names truncated to 15 characters by MAFFT's CLUSTAL writer, so they could not be joined to their labels; now named by accession.
- The six sapovirus strains listed under two genogroups are resolved in the taxonomy generator, so regenerating the CSV no longer undoes it.
- Stopped tracking generated reference-data cache artifacts and ignored `calicityper/data/.cache/`.

- MH003837 is GXV/GXV.P, not GNA4/GNA4.P, confirmed against the published tables.
- GNA2.P and GNA3.P were crossed between OP963623 and OQ818319 in Supplemental Table 1; uncrossed per the paper's pig/bat table.
- AJ606693 is GI.4[GI.P1], not GI.P2 as the supplement lists.
- GI.P5 was two P-types: Supplemental Table 1 has no GI.P6 row and lists its three strains under GI.P5. Restoring GI.P6 returns GI.P5's within-type maximum from 0.255 to 0.050, inside its published 0-0.052 range, and adds three complete RdRp references.
- The BLAST database search space is pinned, so a call cannot change as the reference database grows.
## [0.7.0] - 2026 (baseline)

Baseline packaged version currently described by `pyproject.toml`. See git history/tags for detailed changes leading to this version.

> Note: this changelog was introduced after 0.7.0; earlier history may be incomplete here. When cutting releases, move relevant `Unreleased` entries under a dated version heading and tag the release commit.
