# Output schemas and provenance

This document formalizes CaliciTyper's machine-readable outputs: column semantics for TSV files, JSON structures, provenance fields, and schema/versioning policy.

## Scope and stability

- **Primary machine-readable outputs:** TSV files written by `calicityper type` (`-o/--out`), `calicityper novel` (`-o/--out`), and `calicityper recomb` (`recombination.tsv` in `-o/--output` directory).
- **Structured machine-readable outputs:** JSON reports/documents produced by `calicityper type --json/--report` and `calicityper recomb` (`recombination.json`).
- **Provenance sidecars:** machine-readable `*.provenance.json` files written alongside file-based TSV outputs to ensure auditability even when the primary table is TSV.
- **Human-readable outputs** (Markdown/HTML) are not governed as stable machine schemas, but JSON/TSV schemas and provenance are.

Schema objects include a `schema` field with `name` and `version` where applicable. Within a schema major version, fields documented here are intended to remain stable; additive fields may be introduced in minor updates. Breaking changes to meaning, type, or required fields should increment the schema major version and be recorded in `CHANGELOG.md`.

## Common provenance fields (auditability)

Machine-readable provenance (JSON reports and `*.provenance.json` sidecars) records, at minimum:

| Field | Meaning |
|---|---|
| `tool` | Tool identifier (`calicityper`). |
| `tool_version` | CaliciTyper package version used for the run. |
| `command` | Subcommand that produced the output (`type`, `novel`, `recomb`). |
| `generated` | UTC timestamp (ISO 8601). |
| `parameters` | Effective run parameters/settings (inputs, thresholds, flags, seeds, threads, modes, etc.) relevant to reproducibility/interpretation. CLI flags take precedence over env-var fallbacks for correctness-affecting settings where exposed. |
| `reference_set` | Shipped reference-manifest identity: `doi`, `reference_set_version`, `manifest_sha256` (SHA-256 of canonical manifest JSON), `download_url`, counts (`n_accessions`, `n_records`, `n_artefacts`), `manifest_path`, and availability/error info. |

Notes:
- `manifest_sha256` identifies the **manifest metadata** shipped with the tool build (canonicalized JSON hash), supporting auditable linkage to a declared reference set/DOI/version. It is not, by itself, a substitute for per-file data verification performed by `calicityper setup --verify-only`.
- For `recomb`, the bundled genotyping reference set is not used for scoring; the shipped reference-manifest identity is still recorded for auditability/scope clarity, alongside alignment-specific provenance.

---

## `calicityper type` outputs

### TSV (`-o/--out`, tab-separated; stdout when `-o -`)

One row per input query sequence (FASTA record), in input order. Columns are emitted by the typing pipeline; the set below reflects the current implementation (additional columns may be appended in future minor schema updates).

| Column | Type / values | Meaning |
|---|---|---|
| `query` | string | Query sequence identifier (FASTA record id). |
| `genogroup` | string | Assigned/detected genogroup (e.g. `GI`, `GII`, other norovirus genogroups, or other Caliciviridae labels) or `unassigned`. |
| `note` | string | Semicolon-separated notes/reasons (detection limits, unsupported/not-validated scope, rejection reasons, coverage notes, QC issues, etc.). |
| `orf1_type` | string | Best ORF1 (RdRp/P-type) call, or empty/`unassigned` when not called. |
| `orf1_support` | number/string | Bootstrap support (%) for the ORF1 call (rounded as reported). |
| `orf1_runners` | string | Runner-up ORF1 candidates/support values when available (implementation-formatted string; see JSON for structured runners). |
| `orf1_sites` | integer/string | Number of aligned informative sites used for ORF1 typing. |
| `orf1_nearest` | string | Nearest reference sequence/accession label (ORF1). |
| `orf1_nearest_identity` | number/string | Identity/similarity metric to nearest ORF1 reference (as computed/reported). |
| `orf1_nearest_type` | string | Type/genotype of nearest ORF1 reference. |
| `orf1_accepted` | string | `accepted` or `rejected: ...` describing why the call is not accepted (e.g. low support, insufficient sites, window-mode rejection when enabled). |
| `orf1_window_cov` | string | Coverage of the RIVM ORF1 typing window as `covered/width` (reported; gating depends on `--rivm-window-mode`). |
| `orf1_putative_novel` | string | Putative-novel screen status for ORF1 (e.g. `known`, `putative_novel_type`, `putative_novel_group`, `not_assessable`/NA marker, or empty depending on context). Tentative screening only. |
| `orf1_putative_novel_detail` | string | Explanation/detail for the ORF1 putative-novel screen result. |
| `orf2_type` | string | Best ORF2 (VP1/capsid type, including GII.4 variant suffix when applicable) or empty/`unassigned`. |
| `orf2_support` | number/string | Bootstrap support (%) for the ORF2 call. |
| `orf2_runners` | string | Runner-up ORF2 candidates/support values when available. |
| `orf2_sites` | integer/string | Number of aligned informative sites used for ORF2 typing. |
| `orf2_nearest` | string | Nearest reference sequence/accession label (ORF2). |
| `orf2_nearest_identity` | number/string | Identity/similarity metric to nearest ORF2 reference. |
| `orf2_nearest_type` | string | Type/genotype of nearest ORF2 reference. |
| `orf2_accepted` | string | `accepted` or `rejected: ...` for ORF2. |
| `orf2_window_cov` | string | Coverage of the RIVM ORF2 typing window as `covered/width`. |
| `orf2_putative_novel` | string | Putative-novel screen status for ORF2 (tentative). |
| `orf2_putative_novel_detail` | string | Explanation/detail for the ORF2 putative-novel screen result. |
| `gii4_window_cov` | string | Coverage of the GII.4 variant window as `covered/width` when a GII.4 variant call is considered; otherwise empty. |

Interpretation caveats (see README/capability terms): norovirus GI/GII genotyping is the validated/supported scope; other genogroup/genotype calls may be detected but not validated; putative-novel flags are screening signals, not genotype designations.

### JSON report (`--json` or `--report PREFIX` → `PREFIX.json`)

Schema: `calicityper.type.report` v1.0

```json
{
  "schema": {"name": "calicityper.type.report", "version": "1.0"},
  "tool": "calicityper",
  "tool_version": "<version>",
  "command": "type",
  "generated": "<ISO8601 UTC>",
  "n_queries": <int>,
  "reference_set": { "doi": "...", "reference_set_version": "...", "manifest_sha256": "...", ... },
  "parameters": { "input": "...", "bootstraps": 100, "min_support": 70.0, "min_sites": 300, "novel_distance": "ml", "typeable_groups": ["GI","GII"], "rivm_window_mode": false, "threads": 1, "seed": 1, ... },
  "results": [ { /* per-query record; same fields as TSV plus structured runners */ } ]
}
```

- `results[].orf1_runners` / `results[].orf2_runners`: structured lists (parsed from TSV runner strings) with entries like `{ "type": "...", "support": <float> }`.
- `results[]` contains the per-query typing fields described for TSV (including acceptance, coverage, nearest-reference, and putative-novel fields).

### Provenance sidecar for TSV (`<out>.provenance.json` when `-o/--out` is a file path)

Schema: `calicityper.type.provenance` v1.0

Contains `schema`, `tool`, `tool_version`, `command` (`type`), `generated`, `n_queries`, `reference_set`, and `parameters` (effective settings/inputs) as defined above.

### Markdown/HTML reports

`--md`, `--html`, and `--report` also produce human-readable reports. These are not versioned as strict machine schemas; use JSON/provenance sidecars for auditability.

---

## `calicityper novel` outputs

`calicityper novel` performs the Chhabra et al. 2019 2×SD putative-novel screen (norovirus only). Results are **tentative screening outputs**, not genotype designations.

### TSV (`-o/--out`, tab-separated; stdout when `-o -`)

One row per screened region per query (regions selected by `--region`: `orf1`, `orf2`, or `both`).

| Column | Type / values | Meaning |
|---|---|---|
| `query` | string | Query sequence identifier. |
| `genogroup` | string | Detected/assigned genogroup context for the query. |
| `region` | `orf1` \| `orf2` | Screened region. |
| `basis` | string | Sequence basis required/used for the screen (e.g. nucleotides for ORF1 typing region, amino acids for complete VP1/ORF2). |
| `distance` | string | Effective distance measure used (`ml`, `nj`, or `p`) as configured. |
| `type` | string | Current/best type call context for the region (if available). |
| `support` | number/string | Bootstrap/support value for the type context (if available). |
| `accepted` | string | Acceptance status for the type context (if available). |
| `nearest` | string | Nearest reference (if available). |
| `nearest_identity` | number/string | Identity/similarity to nearest reference (if available). |
| `flag` | string | Screen outcome: e.g. known/not-flagged, `putative_novel_type`, `putative_novel_group`, or `not_assessable` (with reason). |
| `closest_type` | string | Closest known cluster/type considered (when applicable). |
| `distance` (value column) | number/string | Observed distance to closest/assessed cluster as reported by the screen (column name is `distance` in output; interpretation depends on `distance` measure/basis). |
| `bound` | number/string | Applicable 2×SD bound for the assessed cluster/measure (when applicable). |
| `detail` | string | Human-readable explanation (including why not assessable, or evidence for tentative flag). |

> Note: the output includes a column describing the selected distance **measure** and the numeric distance/bound values used by the screen; consult the header row of produced TSVs and `docs/OUTPUT_SCHEMA.md` version notes if field naming evolves. Current implementation emits `basis`, `distance` (measure), and numeric `distance`/`bound` result fields as shown by the command.

### Provenance sidecar for TSV (`<out>.provenance.json` when `-o/--out` is a file path)

Schema: `calicityper.novel.provenance` v1.0

Contains `schema`, `tool`, `tool_version`, `command` (`novel`), `generated`, `n_rows`, `reference_set`, and `parameters` (including `input`, selected `region`, `screen_regions`, effective typing/screening settings such as `min_sites`, `novel_distance`, `typeable_groups`, `rivm_window_mode`, `bootstraps`, `min_support`, `max_refs_per_type`, `seed`, `threads`).

---

## `calicityper recomb` outputs

`calicityper recomb` analyzes a user-supplied alignment (or aligns `--input`) and reports recombination evidence. It does **not** perform genotyping against CaliciTyper's bundled typing reference set; however, shipped reference-manifest identity is recorded for auditability.

### JSON (`recombination.json`)

Schema: `calicityper.recomb` v1.0

```json
{
  "schema": {"name": "calicityper.recomb", "version": "1.0"},
  "tool": "calicityper",
  "tool_version": "<version>",
  "command": "recomb",
  "generated": "<ISO8601 UTC>",
  "provenance": { /* detailed provenance */ },
  "results": [ /* per-query recombination results */ ]
}
```

`provenance` includes: `tool`, `tool_version`, `command`, `generated`, `alignment_file`, `alignment_sequences`, `alignment_length`, `ambiguity_chars_masked`, `reference_metadata`, `tools_used`, `tools_skipped`, `parameters` (virus/regions, thresholds, correction, requested methods, scoring weights, segment-tree settings, seeds/timeouts, requested queries, etc.), and `reference_set` (shipped manifest DOI/version/hash identity, scope-noted).

`results[]` entries include `sequence_id`, `input_qc`/QC details (as applicable), and a `recombination` object with overall status/score/confidence/notes and `events[]` containing breakpoint coordinates (alignment/query), supporting methods/p-values, parental candidates, phylogenetic/biological/plausibility details, and per-event scores/components.

### TSV (`recombination.tsv`)

Tab-separated event-oriented table (one row per event; queries with no events emit a summary row). Columns:

| Column | Meaning |
|---|---|
| `sequence_id` | Query sequence id. |
| `status` | Overall or event status (e.g. non_recombinant/recombinant/borderline/not_assessed/insufficient_sequence depending on context). |
| `score` | Overall or event recombination score. |
| `confidence` | Confidence label derived from score/status. |
| `event_id` | Event identifier (empty/absent semantics for no-event rows). |
| `breakpoint_alignment_start` | Breakpoint start in alignment coordinates. |
| `breakpoint_alignment_end` | Breakpoint end in alignment coordinates. |
| `breakpoint_query_start` | Breakpoint start in query coordinates. |
| `breakpoint_query_end` | Breakpoint end in query coordinates. |
| `supporting_methods` | Comma-separated supporting methods/scanners. |
| `candidate_major_parent` | Candidate major parent label/id (lineage label when metadata available, else sequence id). |
| `candidate_minor_parent` | Candidate minor parent label/id. |
| `phylogenetic_support` | Phylogenetic/segment support summary for the event. |
| `biological_region` | Genome region/context of breakpoint plausibility. |
| `notes` | Additional notes/warnings (including QC/skip reasons). |

### Provenance sidecar (`recombination.tsv.provenance.json`)

Schema: `calicityper.recomb.provenance` v1.0

Contains `schema`, `tool`, `tool_version`, `command` (`recomb`), `generated`, `n_results`, `n_events`, and the full `provenance` object described above.

## Versioning policy

- Schema `name` identifies the artifact family; `version` follows `MAJOR.MINOR` for JSON/provenance documents.
- **MAJOR**: breaking semantic or structural changes.
- **MINOR**: backward-compatible additions (new optional fields/columns, added detail).
- `tool_version` records the CaliciTyper software version; schema versions may evolve independently but changes must be documented in `CHANGELOG.md` and reflected here.
- Consumers should treat unknown/additional fields as forward-compatible extensions and rely on documented required/core semantics for interpretation.
