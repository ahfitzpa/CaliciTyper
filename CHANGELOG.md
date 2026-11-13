# Changelog

All notable changes to this project should be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project aims to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
for released versions.

## [Unreleased]

### Added
- First-class documentation for offline, air-gapped, and cluster reference-data installs (`--source`, `--data-dir`, `--verify-only`, `CALICITYPER_DATA`).
- `tests/test_setup_data.py`: local-source/fixture validation of reference-set install, verification, and checksum-mismatch rejection.
- Release/versioning and packaging guidance (`RELEASING.md`, changelog, README/distribution notes).

### Changed
- Public `main` no longer tracks detailed benchmarking artifacts/docs in-repo (moved out of HEAD; historical blobs remain in git history until a future history-rewrite, if performed).
- README validation/install/distribution guidance clarified for current release status.

### Fixed
- Stopped tracking generated reference-data cache artifacts and ignored `calicityper/data/.cache/`.

## [0.7.0] - 2026 (baseline)

Baseline packaged version currently described by `pyproject.toml`. See git history/tags for detailed changes leading to this version.

> Note: this changelog was introduced after 0.7.0; earlier history may be incomplete here. When cutting releases, move relevant `Unreleased` entries under a dated version heading and tag the release commit.
