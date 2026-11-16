# Changelog

All notable changes to this project should be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project aims to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
for released versions.

## [Unreleased]

### Added
- Experimental sapovirus genotyping (`--type-sapovirus`, or `CALICITYPER_TYPE_SAPOVIRUS=1`; off by default until benchmarked). VP1 genotype in the `orf2_*` columns and RdRp P-type in `orf1_*` (both regions lie in sapovirus ORF1). Calls follow the Zhao et al. 2024 nucleotide p-distance cut-offs: within the genotype cut-off of the nearest reference the call is its type; otherwise `untypeable` (with "possible novel genotype" inside the genogroup cut-off). Bootstrap support is reported beside the call, not used to decide it. New `orf1_distance`/`orf2_distance` columns. Acceptance provisionally requires >=90% coverage of the region until a sapovirus truncation curve exists. A genogroup with no reference set for a region (e.g. GIV RdRp) says so instead of borrowing another genogroup's.
- Sapovirus genogroup database: 39 subjects over 35 genogroups (was 12/12).
- First-class documentation for offline, air-gapped, and cluster reference-data installs (`--source`, `--data-dir`, `--verify-only`, `CALICITYPER_DATA`).
- `tests/test_setup_data.py`: local-source/fixture validation of reference-set install, verification, and checksum-mismatch rejection.
- Release/versioning and packaging guidance (`RELEASING.md`, changelog, README/distribution notes).
- Formalized output schema/provenance documentation (`docs/OUTPUT_SCHEMA.md`) covering `type`/`novel`/`recomb` TSV and JSON structures, column semantics, provenance fields, and schema versioning.
- Machine-readable provenance for auditability: `type` JSON/provenance sidecars, `novel`/`type` TSV sidecars (`*.provenance.json`), and enriched `recomb` JSON/provenance sidecar records including tool version, effective parameters, and shipped reference-set DOI/version/manifest hash.
- Detailed user/reference documentation (`docs/REFERENCE.md`) with supported-vs-unsupported matrix, troubleshooting, performance/scalability guidance, and citation/validation guidance.

### Changed
- Public `main` no longer tracks detailed benchmarking artifacts/docs in-repo (moved out of HEAD; historical blobs remain in git history until a future history-rewrite, if performed).
- README restructured into a quickstart-focused entry point with links to detailed reference and schema docs.
- CLI/UX polish: consistent subcommand help/descriptions and examples, clearer input/output behavior (stdin not supported; stdout TSV and stderr status/progress), predictable `--report` naming (`PREFIX.md/.json/.html`), improved progress messages, and clearer actionable errors/hints for missing reference data or external tools.

### Fixed
- Genogroup BLAST uses a fixed search-space size, so growing the database no longer moves calls (MN567957.1 GII.11 flipped to GIX when sapovirus subjects were added).
- A BLAST hit shorter than 100 nt no longer names a genogroup; on the RIVM set every such call was wrong (non-calicivirus records called norovirus GVI/GX).
- Sapovirus alignments were written with names truncated to 15 characters by MAFFT's CLUSTAL writer, so they could not be joined to their labels; now named by accession.
- The six sapovirus strains listed under two genogroups are resolved in the taxonomy generator, so regenerating the CSV no longer undoes it.
- Stopped tracking generated reference-data cache artifacts and ignored `calicityper/data/.cache/`.

## [0.7.0] - 2026 (baseline)

Baseline packaged version currently described by `pyproject.toml`. See git history/tags for detailed changes leading to this version.

> Note: this changelog was introduced after 0.7.0; earlier history may be incomplete here. When cutting releases, move relevant `Unreleased` entries under a dated version heading and tag the release commit.
