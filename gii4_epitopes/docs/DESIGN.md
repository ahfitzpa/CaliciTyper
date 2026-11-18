# F1–F4 design and evidence register

Companion to `../README.md` and the manuscript note
`GII4_VARIANTS_AND_EPITOPES_NOTE.md`. This register is the single place where a
number is either a published fact or a proposal — never both.

## Published facts (verified 2026-11-18)

| Value | Source | Meaning | Safe to encode as a threshold? |
|---|---|---|---|
| <5% within, >=5% between, nt or aa | Parra 2017 Results + Table 1 footnote | Definition of an intra-genotypic variant (>=2 clustered strains) | No — decide nt vs aa and mean vs max first |
| 15% to next-nearest genotype (14.1% within) | Kroneman 2013; Chhabra 2019 | Genotype designation cut-off | No — different taxonomic level |
| ~>=20% aa between groups | Parra 2017 Fig 5B | Immunotype grouping of genotypes | As a static lookup only, not a query call |
| GII.17 Kawasaki C vs D: 5.4 +/- 1.1% aa | Parra 2017 Results | Candidate divergent lineage | As a case for validation, not a rule |
| GII.6 A/B/C up to ~57 aa apart | Parra 2017 Fig 2A | Static genotype can hold divergent co-circulating variants | Motivates F3 representative choice |
| GII.4 >10 variants, 3–8 y (mean 5.3) | Parra 2017 Table 1 | Evolving genotype | Supports current GII.4-specific window |
| P2 dN/dS >=2x whole ORF2 in GI.3/GI.4/GII.4/GII.13/GII.17 | Parra 2017 Results | Surface-exposed antigenic variation | Supports P2-inclusive window |
| Antigenic sites A (E368A), E (S412N) | Parra 2017 Fig 3 | GII.4 examples | Not a substitution set for a caller |
| GII.17 "Static" in Table 1 | Parra 2017 | Regime label | Annotate `static_emerging_candidate`, not `evolving` |

## Input files required (all in `../inputs/`)

| File | Purpose | Status |
|---|---|---|
| `gii4_variant_references.tsv` | 20 GII.4 VP1 accessions + 11 variant labels for F1/F2 | Derived from bundled v0.8.1 reference_accessions.tsv — complete |
| `variant_thresholds.tsv` | Candidate thresholds with `status=PROPOSED/REFERENCE` and evidence | Complete; PROPOSED values NOT approved |
| `sublineage_definitions.tsv` | Candidate divergent sub-lineages for F2 | Partly complete — 5 rows need real accessions from Fig 5C / Table 1 |
| `intragenotypic_clusters.tsv` | Cluster assignments for F3 | Partly complete — GII.17/GII.6 rows need real accessions |
| `immunotype_regime.tsv` | Genotype→immunotype/regime for F4 | 8 supported, 12 to verify from Fig 5C |

Required external inputs NOT in the repo: the verified reference set itself
(`CALICITYPER_DATA`; alignments are installed by `calicityper setup`, not
shipped in the repo/wheel). MAFFT is needed to produce the AA alignment.

## Outputs (generated, `../results/`, untracked)

- `gii4_vp1_aa_279.fasta`, `gii4_vp1_aa_279.aln`, `gii4_variant_index.tsv` (F1)
- `f1_variant_loo.csv` (F1/F2 threshold evidence)
- `f2_sublineage_candidates.tsv` (F2 review table)
- `f3_current_representatives.tsv`, `f3_proposed_representatives.tsv`, `f3_representative_loo.csv` (F3)
- `immunotype_regime.json` (F4, draft with provenance)

## Promotion gate (before any tool change)

1. `check_inputs.py` passes with zero unresolved placeholders.
2. F1: AA LOO accuracy is reported; a threshold is chosen by the maintainer, not by the script.
3. F3: `f3_representative_loo.csv` shows no increase in nearest-reference distance; full re-typing confirms no call moves.
4. F4: every genotype in `immunotype_regime.tsv` is `SUPPORTED`; JSON is added to the reference set and manifest with a new version/hash.
5. F2: a divergent-sublineage rule is defined, validated on the GII.17 case, and documented before any flag is added to output schemas.

Only then change `align.py`/`tree.py`/`cli.py`/`OUTPUT_SCHEMA.md`, with a changelog entry, version bump, and a fresh leave-one-out + benchmark run in `~/scratch/calicityper_runs`.
