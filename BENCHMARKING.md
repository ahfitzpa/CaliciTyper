# Benchmarking

Updated 2026-11-17, against calicityper 0.8.0.

Four independent evaluations, because no single one is sufficient. Agreement
with another tool measures agreement, not correctness. Leave-one-out measures
the algorithm but only over types it already knows. The published distance
ranges test the reference sets rather than the code. Only the amplicon
benchmark has a ground truth known by construction.

All runs are SLURM batch jobs; scripts are in `scripts/`, job logs in
`calicityper_runs/logs/`.

## 1. Against the RIVM Norovirus Typing Tool

Panel: 5,986 Caliciviridae sequences retrieved from GenBank by **versioned**
accession in February 2024 and submitted to the Norovirus Typing Tool. 4,295
norovirus, 515 sapovirus, 1,176 other or unassigned.

| Level | Agreement | n typed by both |
|---|---|---|
| Genogroup, whole panel | 98.1% | 5,887 |
| Genogroup, norovirus only | 98.4% | — |
| ORF2 capsid genotype | 99.6% | 4,042 |
| ORF1 P-type | 99.8% | 3,787 |
| GII.4 variant | 96.2% | — |

Most residual difference is disagreement about which sequences each tool
*declines*, not conflicting assignments. GVI, GVIII, GIX and GX are absent
from the comparator's reference set, so those sequences count as
disagreements however we type them.

This is agreement, not accuracy. Where the classification has moved on, the
comparator is scored as right and we are scored as wrong: 62 sapovirus
sequences reassigned under Zhao et al. 2024 pull whole-panel genogroup
agreement down from 98.06% to 97.59%, while norovirus agreement rises.

Script: `scripts/validate_vs_rivm.py`. Any new comparison code must reuse its
`genotype()` and `norm()` helpers, which strip the legacy alias RIVM carries
in parentheses (`GII.P31 (GII.Pe)`) and our GII.4 variant suffix
(`GII.4_Sydney`). Comparing raw strings invents disagreements that are purely
nomenclature.

## 2. Leave-one-out recovery

Each reference strain is removed from both the alignment and the tree, then
re-typed against what remains, and scored against its published type. This
isolates the typing algorithm from the composition of the reference set.

**Norovirus:** 267 of 270 region/strain combinations recovered (98.9%), with
**no incorrect assignment**. Complete for ORF1 GI (36/36), ORF2 GI (34/34) and
ORF2 GII (91/91); 97.2% for ORF1 GII (106/109). All three failures were
declines, all GII.P4 (X76716.1, AF145896.1, GQ845368.2), and all three had
bootstrap support of 45-59%, below the 70% acceptance threshold. Median
support across successful recoveries was 100%.

**Sapovirus, complete length:** RdRp 68 of 72 (94.4%), VP1 106 of 125 (84.8%),
no discordant call in either region.

What it cannot measure: recovery of a *known* type is not detection of a novel
one. A removed genotype evolved inside the surviving radiation and sits on a
shorter branch than genuinely novel lineage would, so these are lower bounds.
The singleton and temporal hold-out tests below address that.

## 3. Against the published distance ranges

Recomputes every pairwise within-genotype distance from the aligned reference
sets and compares with the ranges published by Zhao et al. 2024.

| | Within published range | Separated from every other type in its genogroup |
|---|---|---|
| VP1 | 31 of 32 | **22 of 22** |
| RdRp | 19 of 20 | **1 of 20** |

The VP1 exception is GXVI (0.330 against a published 0.273), represented by
only two strains. The RdRp exception is GI.P1 at 0.144 against 0-0.085, which
is AJ606693 sitting atypically within its own type.

The 1-of-20 RdRp separation figure is the headline limitation; see
`docs/SAPOVIRUS_STATUS.md` for why it does not mean the 0.266 cut-off is
wrong, and what it does mean in practice.

### Singleton test

A genotype with one reference cannot be recovered once that reference is
removed, so `untypeable` is the only correct answer. This measures how often a
genuinely novel type would be handed an existing label.

- VP1: **6 of 6** correctly untypeable.
- RdRp: **0 of 5**. All five singleton P-types were given another known
  P-type's label.

## 4. Surveillance-length amplicons

The only benchmark with ground truth known by construction: the simulated
capsid amplicon datasets of Fitzpatrick et al. (doi:10.5281/zenodo.7271067,
CC BY 4.0). Ten simulations, 1,069 amplicons of 330-344 nt from the GI
(GISKF/GISKR) and GII (G2SKF/G2SKR) primer pairs, 1,057 carrying an expected
genotype. Norovirus capsid only.

| | Result |
|---|---|
| Genogroup | 1,057 of 1,067 (99.1%) |
| Genotype, all amplicons | 1,042 of 1,057 (98.6%) |
| Genotype, among typed | 98.8% |
| **Precision among accepted calls** | **917 of 918 (99.9%)** |
| Correct calls withheld | 125 |

Ten of the 13 errors are GI.8 amplicons typed as GI.3, from ten distinct
source accessions, with one reciprocal GI.3 typed as GI.8. Across the 314
aligned sites these amplicons provide, GI.3 and GI.8 are not reliably
separable. That is a limit of the amplicon, not of the implementation: no
method can separate types that differ at too few sites within the amplified
region.

Twelve of the 13 errors had support of 0-66% and were declined. The one that
got through had 71% support, one point over the threshold.

### Why the threshold stays at 70%

Calibration across three datasets put every error at support <= 71, making 75
the lowest threshold with no errors anywhere. The entire margin is that single
observation at 71%, and raising the threshold would withhold correct calls to
catch it. Decision taken 2026-11-17: keep 70, and report support beside every
call so a borderline result is visible. Revisit with more than one observation
in the gap.

## 5. Truncation behaviour (sapovirus)

Each of 208 references cut to 150-1,200 nt at its 5' end, centre and 3' end,
then typed leave-one-out: 5,824 assignments.

- VP1: 82.9-87.2% correct across the whole range, against 84.8% at complete
  length. No fragment at any length produced a discordant genotype, so
  precision is 1.000 at every coverage floor from 0.0 to 0.95 while recall
  falls from 1.000 to 0.030.
- RdRp: 91.2-93.1%, against 93.1% at complete length.

This is why no minimum coverage is applied to sapovirus calls. Coverage is
recorded for every call; a floor would only withhold correct ones.

## 6. Recombination

Nine synthetic datasets with planted breakpoints, three on a sapovirus genome
model: all nine returned the expected status. A tenth reported the same
planted event as an RdRp-VP1 junction under the sapovirus model and an
ORF1-ORF2 junction under the norovirus model, which is the point of having
genus-specific models.

Worked examples, with output, in `examples/recombination/`.

**Determinism is verified, not assumed:** two runs of the same sapovirus
dataset with the same seed return identical status, score and event counts.

Against the 13 strains of Zhao et al. Table 3, scored on that paper's own
expectation. Note that the validation script counts only `likely` and `strong`
as detections, so `possible` is scored as a non-detection.

A tool failure is never reported as a negative. If no scanner ran, the status
is `not_assessed`, not `none`; a query that raised is `failed`. Both sit in
`UNASSESSED` alongside `insufficient_sequence`, and the reason is recorded in
`provenance.tools_skipped`.

## Reproducing

OpenRDP is not in `calicityper_env`; it needs its own environment (numpy<2)
and `CALICITYPER_OPENRDP` pointing at the binary. Without it, 13 tests --
every real-scanner test -- are silently deselected. See `docs/INSTALL_TOOLS.md`.
