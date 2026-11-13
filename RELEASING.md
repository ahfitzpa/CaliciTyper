# Releasing CaliciTyper

This document defines the release/versioning conventions and a reproducible checklist for maintainers.

## Versioning policy

- **Scheme:** [Semantic Versioning](https://semver.org/) (`MAJOR.MINOR.PATCH`).
  - **MAJOR:** incompatible CLI/output/schema or behavior changes that require user action.
  - **MINOR:** backward-compatible features, new capabilities, or meaningful non-breaking improvements.
  - **PATCH:** backward-compatible bug fixes, docs/metadata fixes, and low-risk maintenance updates.
- **Pre-1.0 note:** while research-grade and evolving, breaking changes should still be communicated clearly in release notes/changelog and, when practical, grouped into MINOR releases with explicit migration notes.
- **Source of truth for version:** `pyproject.toml` `[project].version` and `calicityper.__version__` must stay in sync.

## Release artifacts and channels

- **Git tags:** every public release is tagged from a clean commit on `main` as `vX.Y.Z` (annotated tag preferred).
- **Release notes:** publish notes for each tag (GitHub Release) summarizing user-facing changes, compatibility notes, and provenance (including reference-set/manifest information when relevant). Keep `CHANGELOG.md` updated.
- **Python packages:** build `sdist` and `wheel` from the tagged commit. Intended distribution channels: PyPI and/or bioconda (when ready). Until published, document and use source/conda-from-repo install paths.
- **Reference data:** distributed separately as a versioned, checksummed set described by `calicityper/data/reference_manifest.json`; releases must not bundle heavy reference sequences/alignments/BLAST DB in the wheel/repo.

## Reproducible release checklist

Run from a clean checkout of the release commit/tag candidate (no uncommitted changes, no local generated artifacts).

### 1) Prepare and version
- [ ] Confirm target version (`X.Y.Z`) and scope; ensure `CHANGELOG.md` has curated entries under `Unreleased`.
- [ ] Update `pyproject.toml` version (and `calicityper.__version__` if maintained separately) to `X.Y.Z`.
- [ ] Move `Unreleased` entries in `CHANGELOG.md` to `## [X.Y.Z] - YYYY-MM-DD`; leave a fresh empty `Unreleased` section.
- [ ] Verify user-facing docs/metadata are consistent (README install/distribution status, `CITATION.cff` if applicable, project URLs).

### 2) Reference data / manifest provenance
- [ ] Confirm `calicityper/data/reference_manifest.json` is correct for the release (`reference_set_version`, DOI/record, download source, checksums/artefact counts as applicable).
- [ ] If the reference set changes: update manifest/checksums and publish/update the versioned dataset/DOI record; document provenance. (Do not ship restricted-source data publicly.)
- [ ] Verify setup flows: `calicityper setup --source <verified-tarball> --data-dir <tmpdir>` and `--verify-only`; ensure offline/cluster instructions remain accurate.

### 3) Build
- [ ] Build artifacts in a clean environment (e.g. `python -m build`) producing `sdist` and `wheel`.
- [ ] Record tool versions/environment used for the build if needed for reproducibility.

### 4) Test matrix
- [ ] Run unit/CLI/regression tests: `pytest tests/ -v` (default, non-slow) on supported Python versions (at minimum CI-covered versions, e.g. 3.10 and 3.12).
- [ ] Run/confirm CI is green for the release commit (OS/Python matrix as configured).
- [ ] Smoke-test external dependencies where relevant (BLAST/MAFFT present) and basic commands (`calicityper --help`/subcommands).
- [ ] If applicable, run separately tracked slow/integration validations (e.g. OpenRDP-dependent recombination tests) and record outcomes.

### 5) Wheel-data audit and install smoke test
- [ ] Inspect built wheel contents: ensure **no** heavy reference data is bundled (no `calicityper/data/database/` or `calicityper/data/reference_sequences_*`) and required small metadata files are present (`reference_manifest.json`, `reference_accessions.tsv`, `putative_novel_bounds.tsv`).
- [ ] Inspect sdist for unintended generated/cache/build artifacts.
- [ ] In a fresh environment, install from the built wheel/sdist, run `calicityper setup --source <verified-tarball-or-fixture> --data-dir <tmpdir>` (or documented path), and smoke-test e.g. `calicityper type examples/example.fasta -o out.tsv`.

### 6) Tag, publish, and announce
- [ ] Create annotated tag `vX.Y.Z` on the release commit and push tag/branch as appropriate.
- [ ] Create GitHub Release with release notes (include highlights, breaking/compat notes, install instructions, and reference-set/DOI provenance).
- [ ] Publish to PyPI and/or bioconda if/when approved (credentials/access required); verify install from the published channel.
- [ ] Post-release: bump to next development version if desired and reopen `Unreleased` in `CHANGELOG.md`.

## Notes

- Do not commit generated artifacts (`build/`, `dist/`, `.pytest_cache/`, `__pycache__/`, local data caches, benchmark/run outputs) except intentionally curated, documented assets.
- Prefer reproducible inputs: pinned/recorded environments where practical, versioned manifests/DOIs for data, and explicit test/build commands in CI and this checklist.
