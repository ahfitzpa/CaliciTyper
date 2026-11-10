# Validation data (not tracked in git)

The large validation datasets are deliberately **not** committed: they are
third-party outputs and intermediate test results, not part of the tool, and
they would make every clone of a genotyping tool 1.2 MB heavier for no benefit.
`.gitignore` excludes them.

To reproduce the RIVM comparison you need one file here:

### `results_rivm.csv`

Calls produced by the [RIVM Norovirus Typing
Tool](https://www.rivm.nl/mpf/typingtool/norovirus/) for 4,392 GenBank
accessions (1,868 of them GII.4), used as the reference answer set.

Columns used by `scripts/validate_vs_rivm.py`: `name`, `polymerase type`,
`capsid type`, and the two matching `* support` columns.

It is **not redistributed here.** Obtain it either from the upstream project
(<https://github.com/zhuzhanji/NorovirusGenotyping>, as
`new_data/results_rivm.csv`), or regenerate it by submitting the accession set
to the RIVM typing tool and exporting the results as CSV.

Once it is in place:

```bash
python scripts/fetch_rivm_test_set.py rivm_test_set.fasta   # ~4,392 seqs from GenBank
calicityper rivm_test_set.fasta -o rivm_typed.tsv --threads 16
python scripts/validate_vs_rivm.py rivm_typed.tsv           # -> results/rivm_comparison.csv
```

Only the small summary written to `results/` is committed, not the inputs.

### Also ignored

- `test_result_3.csv` — an intermediate test output from the original project;
  not needed to run or validate the tool.
- `rivm_test_set.fasta` — the downloaded sequences; regenerate with the fetch
  script rather than storing them.

### Kept in git

`rivm_accessions.tsv` lists every accession in the RIVM benchmark, all 5,986
of them, **with the GenBank version** (`MW460210.1`, not `MW460210`). GenBank
records get revised, and a fetch by bare accession returns the newest version,
which is not necessarily the one RIVM typed. `genotype_set` is `yes` for the
4,392 records that have a capsid or polymerase call (the genotype comparison)
and `no` for the rest (non-norovirus caliciviruses, used in the genogroup
comparison). Accession IDs are public GenBank identifiers; RIVM's calls are
not included.

`rdrp_taxonomy.csv` and `vp1_taxonomy.csv` (about 7 kB each) stay tracked.
They are not validation data: the legacy `scripts/legacy/msa_reference_data.py` reads
them to rebuild the reference sequence sets, so they are inputs to the
reference data rather than test fixtures.
