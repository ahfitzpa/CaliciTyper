# F1–F4 development scaffolding: GII.4 variants, epitopes and immunotypes

**This directory is research/development scaffolding, not shipped tool behaviour.** Nothing here changes a call, threshold or output column. Scripts build candidate inputs and validation tables; a maintainer must review the results and decide whether to promote any of F1–F4 into the tool.

Follows the benchmarking-folder convention: `scripts/` (runnable steps), `inputs/` (TSV/CSV required), `results/` (generated, do not commit), `docs/` (design + evidence).

## Status and dependencies

| Item | Status | Needs before it can be implemented |
|---|---|---|
| F1 GII.4 AA variant caller | scaffolding only | GII.4-only AA alignment; validated variant threshold; LOO re-run |
| F2 divergent-sublineage flag | scaffolding only | sublineage table; validated distance rule; **5%/15% not yet approved as flags** |
| F3 cluster-aware representatives | scaffolding only | intra-genotypic cluster table; LOO proof calls do not move |
| F4 immunotype/regime lookup | **table recovered** | letters for 9 of 12 clades; manifest version bump |

External tools: Python 3.10+, `biopython` (translation/alignment parsing), `MAFFT` (F1 alignment), `iqtree` (optional, ML distances). All analyses were designed for the UCD Sonic HPC (`~/scratch/calicityper_runs`), matching the manuscript benchmarking layout.

## Order of operations

1. `python scripts/check_inputs.py` — validate every input TSV before use.
2. `python scripts/build_gii4_aa_reference.py` (F1) — build the GII.4-only translated VP1 alignment from the installed reference set. **Requires `CALICITYPER_DATA`** pointing at a verified reference set (the 0.8.x set), because the bundled repo does not carry the alignments.
3. `python scripts/validate_variant_thresholds.py` (F1/F2) — LOO recovery across candidate distance thresholds; writes `results/`. **Review before choosing any threshold.**
4. `python scripts/build_sublineage_table.py` (F2) — derive candidate divergent sub-lineages from the AA alignment; produces a reviewable TSV.
5. `python scripts/select_representatives.py` (F3) — pick two representatives per genotype from distinct clusters; compare against current first-two rule.
6. `python scripts/validate_representatives.py` (F3) — LOO with current vs proposed representatives; **must show no existing call moves** before adopting.
7. `python scripts/build_immunotype_lookup.py` (F4) — curate genotype→immunotype/regime table and emit JSON with provenance; add to the reference set and manifest.
8. `python scripts/audit_reference_tables.py` — optional cross-check that every genotype in the lookup exists in the reference set.

## Critical evidence notes (do not bypass)

- Parra et al. 2017 (doi:10.1371/journal.ppat.1006136; PMC5283768; CC0 1.0) define an intra-genotypic **variant** as a cluster of >=2 strains differing by **<5% within, >=5% between**, on **nt or aa** sequence. The abstract's "<=5" is a *count of variants per genotype*, not a percentage. The AA/nt basis and mean-vs-maximum question must be settled before any threshold enters code.
- The **15%** VP1 aa distance is the *genotype* cut-off to the next-nearest genotype (Kroneman 2013 / Chhabra 2019; within-genotype divergence 14.1%). It is not a sub-lineage flag.
- Immunotypes A–L are grouped at approximately **>=20% VP1 aa difference between groups**; GI.3/GI.7/GI.8/GI.9 form one immunotype, GII.4+GII.20 another. These are **published genotype groupings**, not properties of a query sequence.
- Parra's Table 1 labels **GII.17 "Static"**; the Kawasaki_2014 C/D 5.4 +/- 1.1% aa difference is discussed as a possible emerging lineage. Do not tag GII.17 as "evolving" without a decision.

## Provenance

Generated tables carry the reference-set manifest SHA-256 and the source DOI where applicable. Generated results are written to `results/` and are intentionally untracked.
