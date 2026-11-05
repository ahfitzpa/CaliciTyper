# norotyper: command-line norovirus genotyping (RIVM method)

Fork of `zhuzhanji/NorovirusGenotyping` (UCD diploma project). The original scripts are kept in `scripts/`; the
maintained tool is the `norotyper` package.

## What changed
- **No machine-specific paths.** Data: `$NOROTYPER_DATA` or `<repo>/new_data`. Programs: found on `PATH`
  (override with `NOROTYPER_BLASTN`, `NOROTYPER_MAKEBLASTDB`, `NOROTYPER_CLUSTALW`). Legacy scripts no longer
  contain personal cluster paths or depend on the directory they are run from; temporary files are unique per run.
- **Follows the published RIVM method**: BLAST against whole genomes (genogroup) -> >=100 nt overlap with ORF1 or
  ORF2 -> ClustalW profile alignment to pre-aligned references -> neighbour-joining tree -> **100 bootstrap
  replicates** -> a type is assigned only with bootstrap **>= 70%** (RIVM: "Only with bootstrap > 70").
  The old `phylo.R` built a maximum-likelihood tree with **5** bootstraps; it is now 100 and only used by the legacy workflow.
- **Finishes the job**: reports genogroup, ORF1 (P-type) and ORF2 (capsid) type with bootstrap support for any
  number of sequences, in parallel (`--threads`).

## Use
```
conda install -c conda-forge -c bioconda blast clustalw biopython numpy pandas
pip install -e .
norotyper queries.fasta -o typed.tsv --threads 8        # --bootstraps 100 --min-support 70 are the defaults
```
Output columns: query, genogroup, orf1_type, orf1_support, orf2_type, orf2_support, note.

## Validation
- `python tests/loo_validation.py`: leave-one-out on the bundled reference sets (no internet).
- `scripts/fetch_rivm_test_set.sh` then `norotyper ... ` then `python scripts/validate_vs_rivm.py typed.tsv`: compares
  with the 5,986 RIVM calls in `new_data/results_rivm.csv` (needs NCBI access; not yet run).

## Known differences from RIVM (to check with the RIVM comparison)
- Distance: Tamura-Nei (analytic) instead of PAUP* HKY85; tree: NJ, as RIVM describes, but tie-breaking and
  distance estimation differ slightly.
- "Assigned" is defined as the query forming a clade with all (up to 2) reference sequences of one type in >=70%
  of replicates. RIVM's exact clade rule and its choice of representative references are not published; the
  subset of 2 per type is inferred from their description and from leave-one-out accuracy.
- GII.4 variant assignment (RIVM step 3) is not implemented. Reference sets are the 2019 (Chhabra) sets.
- The ORF region windows are the full annotated ORFs, not RIVM's fixed typing windows.
- Licence: none yet; to be set by UCD.

## Leave-one-out results (bundled references, 100 bootstraps, >=70% support, 2 references per type)
| Region | Queries | Assigned | Correct | Wrong |
|---|---|---|---|---|
| ORF2 capsid (GI, GII) | 125 | 100% | 100% | 0 |
| ORF1 P-type (GI, GII) | 147 | 93.9% | 93.2% | 1 (GII.PNA4 called GII.P11 at 99%) |
Speed: about 5 s per sequence per core (ClustalW dominates; 100 bootstrap NJ trees take about 1 s); runs in parallel with `--threads`.
