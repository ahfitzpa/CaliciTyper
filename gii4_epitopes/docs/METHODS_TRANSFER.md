# What can be transferred from the four papers' methods

Sources: pasted-7/8 (Tohma 2017 mSphere GII.2), pasted-9 (Kendra 2021 PNAS
GII.4 antigenic cartography), pasted-10 (Parra 2017 PLoS Pathog genomics).

## Tier A — adoptable in CaliciTyper now (no new wet-lab data)

1. **Aa-vs-nt divergence heatmap (pasted-8, pasted-10).** Parra's Python script
   builds a matrix of pairwise aa differences × isolation-year differences and
   plots counts. We already have the S3 Table with year/country for 131
   accessions — the exact input this needs. adoptable as a benchmark figure:
   `dev/gii4_epitopes/scripts/aa_year_heatmap.py`, consuming
   `parra2017_s3_immunotype_sequences.tsv`. Reinforces the static/evolving
   distinction and the synonymous-nt-vs-aa argument (F1 rationale).
2. **Shannon entropy per VP1 residue (pasted-9).** Kendra used Shannon Entropy-One
   on 1,572 GII.4 VP1 sequences, with entropy ≤0.3 = conserved. Pure computation
   over an alignment; adoptable as a `calicityper antigenic` annotation
   (`entropy` per site) and as a benchmark column. Uses only alignment columns.
3. **dN/dS per genotype and per P2 (pasted-10).** Parra reported dN/dS over VP1
   and separately P2. Computable from codon-aligned nt sequences; adoptable as
   a benchmark metric for the GII.4 AA-window decision (P2 dN/dS ≥2× VP1 in
   five genotypes). Needs codon-aware alignment (existing pipeline).
4. **Coevolving-residue consensus (pasted-9).** Kendra required agreement of two
   methods (BGM edge support >0.95 in HyPhy AND CAPS2 bootstrap >0.80) before
   calling a pair coevolving. Our `COEVOLVING` list is from Kendra Table 1 —
   already the *consensus* set. We can add the two-support criterion as a
   documented provenance field; no new analysis needed.
5. **Micro-averaged concordance + Wilson CIs (already in benchmark).** Matches
   the papers' count-based comparison style; keep.

## Tier B — adoptable in benchmarking only (needs data we have or can fetch)

6. **Root-to-tip / temporal signal (pasted-7).** Tohma used TempEst on ML trees
   with isolation years (R²=0.96 for GII.2 VP1). Adoptable as a benchmark
   diagnostic: our S3 Table gives years; a root-to-tip regression on the NJ/ML
   tree would show whether the reference sets carry clock-like signal. Cheap,
   no external lab work.
7. **Best-model selection by AICc (pasted-7).** GTR+G+I chosen for GII.2. Our
   tool uses Tamura–Nei by design (reproducing RIVM). Adoptable as a
   *sensitivity analysis* in the benchmark: compare TN vs GTR+G+I NJ trees and
   show calls unchanged — a reviewer-friendly robustness check.
8. **ALRT / bootstrap node support (pasted-7).** PhyML aLRT vs our 100-replicate
   bootstrap; report both on the same trees as a robustness comparison.
9. **Selection scans iFEL/MEME (pasted-8).** Needs HyPhy and a codon alignment;
   adoptable in a *future study*, not in the tool (runtime/dependencies).
10. **BEAST MCMC evolutionary rates (pasted-7, pasted-10).** Strict/relaxed
    clock, constant/skyline priors, ESS>200, 10% burn-in. Heavy; benchmark-only
    if at all; not for the CLI.

## Tier C — NOT adoptable from sequence alone (experimental)

11. **Antigenic cartography (pasted-9).** EC50 blockade matrices → ACMACS 2D/3D
    maps, 1 unit = twofold difference, stress 0.1, Mantel correlation of
    antigenic vs aa distance. Requires VLPs + hyperimmune sera (24 VLPs, 9
    mouse sera, 4 mice each); impossible from sequence. Our ≥18-aa rule is the
    sequence-only shadow of this — keep the caveat explicit.
12. **HBGA blockade and HIE neutralisation (pasted-9).** Wet-lab; defines the
    ground truth that our `candidate_antigenic_shift` screen approximates.
13. **Coevolution via HyPhy BGM + CAPS2 (pasted-9).** Could be reproduced
    computationally ( HyPhy available) but needs 308-sequence subsample and
    two tools; future study only.
14. **Structural mapping to PDB 4RPB / 4QPX / 7K6V (pasted-8, pasted-9).**
    Needs Chimera and structures; future study.

## Concrete next actions (smallest first)

- A1: `aa_year_heatmap.py` from S3 Table with years (Tier A1) — uses data already parsed.
- A2: add Shannon-entropy annotation to `antigenic.py` (Tier A2) with a synthetic test.
- B1: TN vs GTR+G+I sensitivity row in benchmark (Tier B7) — HPC job.
- C1: write a "future study" paragraph in the manuscript: antigenic cartography
  on CaliciTyper-called GII.4 variants would require VLP/sera data; our screen
  is not a substitute.

Keep every adopted metric provenance-tagged with its source DOI and whether it
is sequence-derived or assay-derived.
