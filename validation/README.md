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

## Known reference-strain caveats

### JF781268.1 — genuinely discordant, kept and documented

`JF781268.1` (CU081210E) is capsid genotype **GIV.2** but polymerase group
**GVI.P1**, so it sits in the ORF1 GIV reference set while its polymerase
belongs elsewhere. That inflates the apparent within-cluster spread of ORF1
GIV in any statistic computed per genogroup.

It is kept rather than excluded or special-cased. The discordance is real
biology, not a curation error, and the published classification lists it this
way; removing it would make our reference set disagree with the scheme the
tool implements, and special-casing it would hide a case users will meet in
their own data. Read ORF1 GIV cluster statistics with this in mind.

### Sapovirus P-types cannot yet flag a novel type

The RdRp P-type cut-off of 0.266 overlaps the body of the between-type
distance distribution, and most P-types have a neighbouring P-type closer than
the cut-off. The cut-off was derived as the mean minus two standard deviations
of a pooled between-type peak, so a tail below it is expected by construction
and this is not evidence that 0.266 is wrong. It is, however, a practical
limit: where two known P-types are closer together than the threshold, a query
from an unrepresented P-type can fall inside the threshold of one or both, and
the threshold cannot then establish that it belongs to neither.

Measured directly: of five singleton RdRp P-types, **none** was reported as
untypeable once its only reference was removed; all five were given the label
of a different known P-type. VP1 behaves correctly under the same test, with
6 of 6 singleton genotypes reported untypeable.

So treat sapovirus P-type assignments as provisional, particularly for a query
distant from a dense cluster of references, and do not rely on the RdRp
threshold alone to identify a novel P-type. Assignment *among known* P-types
is accurate: 68 of 72 recovered under leave-one-out, with no discordant call.
