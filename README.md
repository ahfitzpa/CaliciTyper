# norotyper: command-line norovirus genotyping (RIVM method)

Fork of `zhuzhanji/NorovirusGenotyping` (UCD diploma project). The original scripts are kept in `scripts/`; the
maintained tool is the `norotyper` package.

## What changed
- **No machine-specific paths.** Data: `$NOROTYPER_DATA` or `<repo>/new_data`. Programs: found on `PATH`
  (override with `NOROTYPER_BLASTN`, `NOROTYPER_MAKEBLASTDB`, `NOROTYPER_CLUSTALW`). Legacy scripts no longer
  contain personal cluster paths or depend on the directory they are run from; temporary files are unique per run.
- **Follows the published RIVM method**: BLAST against whole genomes (genogroup) -> >=100 nt overlap with ORF1 or
  ORF2 -> MAFFT alignment to pre-aligned references -> neighbour-joining tree -> **100 bootstrap
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
## Speed (v0.3)
| Version | Aligner | Time per whole genome (1 core) |
|---|---|---|
| 0.2 | ClustalW profile alignment, one call per query and region | about 10 s (ORF1 + ORF2) |
| 0.3 | MAFFT `--addfragments --6merpair --keeplength`, queries batched per region/genogroup, cropped to the ORF window | about 0.5 s |

What changed in 0.3, with no change to the calls:
- **MAFFT replaces ClustalW** as default. The reference alignment is never recomputed and its columns stay fixed. `--6merpair` matters: MAFFT's default pairwise mode took 52 s for 20 ORF queries, `--6merpair` took 0.23 s.
- **Queries are cropped** to the ORF region (BLAST hit coordinates + 300 nt) before alignment, so whole genomes are no longer aligned in full.
- **Bootstrap:** identical alignment columns are collapsed to weighted patterns, base indicators are precomputed once (float32; counts are integers, so exact), and neighbour-joining plus the clade test run in numba (optional; pure numpy fallback).
- **Parallel:** BLAST uses `--threads`, and alignment and bootstrapping run in the same worker pool, in chunks.

Measured on a 2-core machine (100 bootstraps): 100 genomes in 26 s, 400 genomes in 99 s. Throughput scales with cores (each query is independent), so 16 cores should do 400 genomes in roughly 15 s.

Regression checks against the v0.2 ClustalW results: leave-one-out, 272 sequences, 0 differing calls (ORF1 93.9% assigned / 93.2% correct, ORF2 100% / 100%); the 400-genome CLI output was byte-identical before and after float32 indicators, cropping and chunked parallelism (compared with the MAFFT run before those changes). Use `--backend clustalw` in `tests/loo_validation.py` to reproduce the old path.

**Illumina runs:** the tool types contigs or genomes, not reads. Map or BLAST reads to norovirus references, assemble or cluster the hits (e.g. vsearch at 99%), then give the representatives to `norotyper`.
