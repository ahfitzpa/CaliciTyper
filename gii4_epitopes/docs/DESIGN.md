# F1–F4 design and evidence register

> **Checkout merge (2026-11-19):** files from the parallel HPC checkout at
> `C:\Users\Amy\Downloads\gii4_epitopes` were merged in: verbatim-quoted antigenic
> site TSV, Kendra coevolving-pair TSV, fetched S3/S4 protein alignments, and
> reconstruction scripts/results. The two directories were then reconciled
> file-by-file: 33 files identical; HPC-unique files all copied here; where
> files differed, the **locally verified** version was kept (immunotype letters
> from Fig S5 + paper text, resolved sub-lineage accessions from Figs 2A/S2,
> expanded threshold and literature notes, F3 results covering all three
> genotypes). `scripts/common.py` was merged: HPC's portable path resolution
> (works in either checkout) plus local's `p_distance` IndexError fix.
> **Authority rule:** where a reconstruction
> disagrees with the published Fig S5 / paper text, the published figure wins.
> Specifically, the HPC reconstruction puts GI.6 as a singleton and GII.22 in a
> cluster; the paper's Fig S5 shows GI.6 in B ("B (GI.2, GI.4, GI.5, GI.6)") and
> GII.22 as singleton L, and the paper text names GII.6 (not GI.6) as a
> singleton. The shipped `antigenic.IMMUNOTYPE` therefore keeps the figure-based
> letters. The reconstructions are retained in `results/` as evidence of the
> independent check, not as the source of truth.

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

## Parra 2017 supplementary tables (added 2026-11-19)

Downloaded publisher supplementary files (not committed) and parsed into
`inputs/` by `scripts/extract_parra_tables.py`:

| File | Content | Parsed TSV |
|---|---|---|
| S3 Table (`ppat.1006136.s009.docx`) | 131 ORF2 accessions, 31 genotypes, with host/year/country — the sequences behind the immunotype tree (Fig 5) | `parra2017_s3_immunotype_sequences.tsv` |
| S4 Table (`ppat.1006136.s010.docx`) | 47 GII.6, 26 GII.12, 144 GII.17 accessions + strain names | `parra2017_s4_gii6_gii12_gii17.tsv` |
| S6 Fig (`ppat.1006136.s006.tiff`) | re-infection matrix (immunotypes A–L). Viewed: it encodes same-immunotype pairs, **not** genotype→letter assignments | not parsed |
| S1/S5 Fig | GII.12 diversity; ML version of Fig 5 | not parsed |

**Resolved 2026-11-19 from additional supplied figures:**
- **Fig S2** (GII.17 trees + diversity heatmap): the four clades are labelled
  A–D, and the heatmap's dashed line sits at **~27% aa difference**. Read from
  the tree: Variant_C contains AB983218 (Kawasaki323/2013/JPN) and LC043139/
  LC043305; Variant_D contains LC037415 (Kawasaki308/2015/JPN) and KT589391;
  Variant_B contains KJ156329, DQ438972, GQ266696/7; Variant_A contains
  AY502009, KJ196286. GII.17 C/D placeholders are now filled with AB983218
  and LC037415. The ~27% cut-off is recorded in `variant_thresholds.tsv` as
  `REFERENCE_FROM_FIGURE` — it is NOT the 5% within-variant definition and is
  not hard-coded into any caller.
- **Fig S5** (ML version of Fig 5, same topology) and **Table 1**: the complete
  12-immunotype mapping is now verified — A GI.1; B GI.2/4/5/6; C GI.3/7/8/9;
  D GII.1/12/16; E GII.2/5/10; F GII.3; G GII.4/20; H GII.6; I GII.7/8/9/14;
  J GII.13/17/21; K GII.15; L GII.22. Table 1 confirms evolving=GII.4 only;
  all other pattern-assigned genotypes are static. Eight earlier guessed
  letters in `antigenic.py` were wrong (GI.2/4/6, GII.1/3/7/12/13) and were
  corrected. The lookup now covers 31 genotypes with zero unverified rows.

**Immunotype B — text confirmation (2026-11-19, `document.pdf` = the paper):**
A parallel checkout disputed whether GI.6 belongs to immunotype B. The paper
text resolves it: "others appeared as **single** genotypes (e.g. **GI.1, GII.3,
GII.6**; Fig 5A, S3 Table)" — so GII.6 (immunotype H) is a singleton, while
GI.6 is grouped. The B label in Fig S5 reads "B (GI.2, GI.4, GI.5, GI.6)", and
the text separately confirms GI.3/GI.7/GI.8/GI.9 cluster as immunotype C and
GII.4+GII.20 as G. The shipped lookup (B = GI.2/4/5/6, GI.6 → B) is therefore
correct; the reconstruction that moved GI.6 to a singleton was wrong. The
paper also states the variant definition verbatim: "a group of strains (>=2)
that clustered together ... <5% difference in their nt **or** aa sequences,
but >=5% difference compared to other strains" — confirming F2's unresolved
aa-vs-nt and mean-vs-max caveat.

**Resolved 2026-11-19 (Fig 2A supplied):** GII.6 variants A/B/C assigned from
the Figure 2A tree: A = AB039777 (SaitamaU4/1997/JPN), B = AF414410
(Miami/1994/USA), C = HM633213 (8913 Shizuoka/2000/JPN). The figure shows
inter-variant amino-acid differences of 5.2%, 7.2% and 9.2% and within-variant
values of 0.7–1.2%, consistent with the ~5% variant definition. All
`FILL_FROM_PARRA` placeholders are now resolved: `check_inputs.py` reports
**0 placeholder warnings**. F3 representative selection and LOO ran on the
bundled reference set: 6 proposed representatives across 3 genotypes, 0
differences from the current first-two rule, and 0/3 increases in nearest
reference distance. Also fixed a p_distance IndexError in the dev helper.

`scripts/crosscheck_parra_accessions.py` reports overlap with the installed
reference set: GII.6 8/47, GII.12 3/26, GII.17 6/144. The GII.17 overlap
includes the Kawasaki strains `KT589391`, `AB983218` and `LC037415`, but the
paper does not state which is clade C vs D, so the C/D labels remain pending.

## Promotion gate (before any tool change)

1. `check_inputs.py` passes with zero unresolved placeholders.
2. F1: AA LOO accuracy is reported; a threshold is chosen by the maintainer, not by the script.
3. F3: `f3_representative_loo.csv` shows no increase in nearest-reference distance; full re-typing confirms no call moves.
4. F4: every genotype in `immunotype_regime.tsv` is `SUPPORTED`; JSON is added to the reference set and manifest with a new version/hash.
5. F2: a divergent-sublineage rule is defined, validated on the GII.17 case, and documented before any flag is added to output schemas.

Only then change `align.py`/`tree.py`/`cli.py`/`OUTPUT_SCHEMA.md`, with a changelog entry, version bump, and a fresh leave-one-out + benchmark run in `~/scratch/calicityper_runs`.
