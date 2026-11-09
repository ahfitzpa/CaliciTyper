# norotyper

Command-line norovirus genotyping, replicating the published method of the
[RIVM Norovirus Typing Tool](https://www.rivm.nl/mpf/typingtool/norovirus/).

`norotyper` assigns, for each input sequence:

- **genogroup** (BLAST against whole norovirus genomes)
- **ORF1 P-type** (polymerase) with bootstrap support
- **ORF2 capsid type** with bootstrap support

A type is reported only when it reaches **>= 70% bootstrap support** over 100
replicates, following RIVM's rule.

> **Status: research code, validated on reference sequences only.** See
> [Validation](#validation) before relying on the calls. The comparison against
> RIVM's own calls has not yet been run.

## Install

```bash
conda create -n norotyper -c conda-forge -c bioconda \
    python=3.12 blast mafft numpy pandas biopython numba
conda activate norotyper
pip install .
```

BLAST and MAFFT must be on `PATH`. The reference data ships inside the package,
so no extra download or configuration is needed. `numba` is optional and makes
the bootstrap roughly ten times faster; without it a pure-NumPy fallback is used.

For development, `conda env create -f environment.yml` installs the same stack
plus `pytest` and an editable install.

## Use

```bash
norotyper contigs.fasta -o typed.tsv --threads 8
```

Output is a TSV with one row per input sequence, in input order:

| column | meaning |
|---|---|
| `query` | FASTA record id |
| `genogroup` | `GI`, `GII`, … or `unassigned` |
| `orf1_type` | P-type, `unassigned`, or empty if not typeable |
| `orf1_support` | bootstrap % for the ORF1 call |
| `orf2_type` | capsid type, `unassigned`, or empty if not typeable |
| `orf2_support` | bootstrap % for the ORF2 call |
| `note` | why a region was not typed (no BLAST hit, `<100 nt` overlap, support below threshold) |

Useful options: `--bootstraps` (default 100), `--min-support` (default 70),
`--threads`, `--seed` (the bootstrap is seeded, so runs are reproducible).

### Input must be assembled

`norotyper` types **contigs or whole genomes, not sequencing reads.** For an
Illumina run: map or BLAST reads against norovirus references, assemble or
cluster the hits (for example `vsearch` at 99% identity), then pass the
representative sequences to `norotyper`.

### Speed

About **0.5 s per whole genome on one core**; 100 genomes in 26 s and 400 in
99 s on two cores. Queries are independent, so throughput scales with
`--threads`.

## Method

1. **Genogroup** — `blastn` against whole genomes (`norotyper/data/database`).
2. **Region selection** — a region is typed only if the BLAST hit overlaps the
   annotated ORF1 or ORF2 by at least 100 nt.
3. **Alignment** — the query is cropped to the ORF window (+/- 300 nt) and added
   to the pre-aligned reference set with MAFFT
   (`--addfragments --6merpair --keeplength`), so reference columns never move.
4. **Tree** — Tamura-Nei distances, neighbour-joining.
5. **Bootstrap** — 100 replicates. A type is assigned when the query forms a
   clade with all references of exactly one type in >= 70% of replicates.

Reference sets follow the 2019 classification (Chhabra *et al.*).

## Validation

**Leave-one-out on the bundled reference set** (`python tests/loo_validation.py`,
no internet needed). Each reference ORF is typed against all other references of
its genogroup; its own label is the answer.

| Region | Queries | Genogroups | Distinct types | Assigned | Correct | Wrong |
|---|---|---|---|---|---|---|
| ORF2 capsid | 125 | GI, GII | 42 | 100% | 100% | 0 |
| ORF1 P-type | 147 | GI, GII | 49 | 93.9% | 93.2% | 1 |

Two caveats that matter:

- **The single error is a confident one.** `AB074893.1` (true `GII.PNA4`) was
  called `GII.P11` at **99% support**. The 70% threshold did not catch it, so
  bootstrap support must not be read as a guarantee of correctness.
- **Leave-one-out flatters the tool.** Queries are curated, full-length
  reference sequences. It shows the method is wired up correctly; it does not
  predict accuracy on partial, divergent or chimeric environmental contigs.

**Not yet validated:** against RIVM's own calls; on any genogroup outside
GI/GII; on assembled environmental or wastewater contigs.

### Running the RIVM comparison

`validation/results_rivm.csv` holds RIVM Typing Tool calls for **4,392
GenBank accessions** (1,868 of them GII.4). You do **not** need to run anything
through the RIVM web tool — the calls are already here. You need
`entrez-direct` and outbound NCBI access:

```bash
scripts/fetch_rivm_test_set.sh rivm_test_set.fasta    # ~4,392 sequences from GenBank
norotyper rivm_test_set.fasta -o rivm_typed.tsv --threads 16
python scripts/validate_vs_rivm.py rivm_typed.tsv     # -> results/rivm_comparison.csv
```

The comparison reports, per region, how often the two tools agree where both
assign a type, plus the cases where only one of them does.

## Known differences from RIVM

- **Distance and tree**: Tamura-Nei computed analytically with neighbour-joining,
  rather than PAUP\* HKY85. Tie-breaking and distance estimation differ slightly.
- **Clade rule**: "assigned" means the query groups with all (up to 2) references
  of one type in >= 70% of replicates. RIVM's exact clade rule and its choice of
  representative references are not published; the subset of two per type was
  inferred from their description and from leave-one-out accuracy.
- **GII.4 variant assignment (RIVM step 3) is not implemented.** This matters:
  GII.4 is the most common capsid type in surveillance data.
- **ORF windows** are the full annotated ORFs, not RIVM's fixed typing windows.

## Repository layout

```
norotyper/            the maintained package
  data/               reference sequences and BLAST database (installed with the code)
scripts/              the original upstream pipeline, kept for provenance
tests/                unit tests, CLI tests, leave-one-out validation
validation/           RIVM comparison calls and taxonomy tables
results/              leave-one-out output
```

## Attribution and licence

MIT licensed. This is a fork and substantial rewrite of
[zhuzhanji/NorovirusGenotyping](https://github.com/zhuzhanji/NorovirusGenotyping)
by Yang Ji, and it reimplements the method of the RIVM Norovirus Typing Tool.
See **[ATTRIBUTION.md](ATTRIBUTION.md)** for what came from where, the citations
for RIVM and the 2019 classification, and an important note on upstream
licensing. `norotyper` is not affiliated with or endorsed by RIVM.
