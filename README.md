# norotyper

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python Version](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![Classification Standard](https://img.shields.io/badge/Classification_Standard-Chhabra_2019-green.svg)](https://doi.org/10.1099/jgv.0.001318)
[![Status](https://img.shields.io/badge/status-research_code-orange.svg)](#validation)

`norotyper` is a command-line tool for dual-typing norovirus **assembled
contigs and whole genomes**, replicating the published method of the
[RIVM Norovirus Typing Tool](https://www.rivm.nl/mpf/typingtool/norovirus/).

It assigns a **genogroup**, an **ORF2 capsid genotype** and an **ORF1
polymerase (P-) type**, each with bootstrap support, following the updated
classification of [Chhabra et al. (2019)](https://doi.org/10.1099/jgv.0.001318).

> **Status: research code.** Genotyping is validated for GI and GII only, and
> the comparison against RIVM's own calls has not yet completed. Read
> [Genogroup and genotype coverage](#genogroup-and-genotype-coverage) and
> [Validation](#validation) before relying on a call.

As far as we are aware, every other norovirus typing tool (the RIVM tool,
CDC's HuCaT) is web-only. `norotyper` is intended to be the scriptable
equivalent: no uploading of surveillance sequences to a third-party form, and
it drops into a pipeline.

---

## Key Features

* **Dual-typing.** ORF1 (RdRp / P-type) and ORF2 (VP1 / capsid) are typed
  independently, so recombinants are visible as a mismatched pair.
* **GII.4 variant resolution.** The ORF2 reference set is labelled at variant
  level, so GII.4 queries resolve to `GII.4_Sydney`, `GII.4_Den_Haag` and so on
  rather than stopping at `GII.4`.
* **RIVM's assignment rule.** 100 bootstrap replicates; a type is reported only
  at **>= 70% support**.
* **Reproducible.** The bootstrap is seeded (`--seed`), so repeated runs on the
  same input give identical output.
* **Self-contained.** The curated reference sets and BLAST database are
  installed with the package; there is no separate database build step.
* **Parallel.** Queries are independent and distributed across `--threads`.

---

## Installation

### Option 1: Conda / Mamba (recommended)

```bash
git clone https://github.com/ahfitzpa/norotyper.git
cd norotyper

conda env create -f environment.yml
conda activate norotyper
```

### Option 2: pip, with dependencies already present

Ensure `blast` and `mafft` are on `PATH`, then:

```bash
pip install git+https://github.com/ahfitzpa/norotyper.git
```

The reference data ships inside the package, so this gives a working install
with no further setup.

---

## Requirements

* **Python** >= 3.9
* **External programs** (must be on `PATH`, or set `NOROTYPER_BLASTN` /
  `NOROTYPER_MAFFT`):
  * [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) - genogroup assignment
  * [MAFFT](https://mafft.cbrc.jp/alignment/software/) - aligns queries onto the
    fixed reference alignment
* **Python packages:** biopython, numpy, pandas
* **Strongly recommended:** `numba`. It is nominally optional, but without it
  the bootstrap runs roughly 8x slower. Install it.
* **Optional:** `clustalw`, only for the legacy `--backend clustalw` path;
  `entrez-direct` is *not* required (the RIVM fetch script uses plain HTTP).

---

## Quick Start

```bash
norotyper contigs.fasta -o typed.tsv --threads 8
```

That is the whole interface: one input FASTA, one output TSV. There is no
database build step and no subcommand.

### Input must be assembled

`norotyper` types **contigs or genomes, not sequencing reads.** There is no
FASTQ mode. For an Illumina run, assemble first: map or BLAST reads against
norovirus references, assemble or cluster the hits (for example `vsearch` at
99% identity), and pass the representative sequences in.

---

## Output

A single TSV, one row per input record, in input order:

| Column | Meaning |
|---|---|
| `query` | FASTA record id |
| `genogroup` | `GI`, `GII`, ... `GX`, or `unassigned` |
| `orf1_type` | P-type, `unassigned`, or empty if the region was not typeable |
| `orf1_support` | bootstrap % for the ORF1 call |
| `orf2_type` | capsid type (variant-level for GII.4), `unassigned`, or empty |
| `orf2_support` | bootstrap % for the ORF2 call |
| `note` | why a region was not typed: no BLAST hit, `<100 nt` ORF overlap, or support below threshold |

### Example

| query | genogroup | orf1_type | orf1_support | orf2_type | orf2_support | note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| NC_001959.2 | GI | GI.P1 | 100.0 | GI.1 | 100.0 | |
| JX459908.1 | GII | GII.P31 | 97.0 | GII.4_Sydney | 99.0 | |
| frag_07 | GII | | | GII.3 | 88.0 | orf1: <100 nt overlap |
| junk_01 | unassigned | | | | | no BLAST hit (not a norovirus / Caliciviridae sequence in reference set) |

Note that for GII.4 the variant is carried **inside** `orf2_type`. Strip the
suffix (`GII.4_Sydney` -> `GII.4`) if you need genotype-level labels, as RIVM
reports them in a separate column.

### Options

```text
usage: norotyper [-h] [-o OUT] [--bootstraps N] [--min-support PCT]
                 [--max-refs-per-type N] [--threads N] [--seed N] [--version]
                 fasta

  -o, --out               output TSV, or "-" for stdout [default: -]
  --bootstraps            bootstrap replicates [default: 100, RIVM's value]
  --min-support           minimum bootstrap % to assign a type [default: 70]
  --max-refs-per-type     reference sequences per type used in the tree [default: 2]
  --threads               worker processes [default: 1]
  --seed                  bootstrap seed, for reproducibility [default: 1]
```

---

## Genogroup and genotype coverage

These are two different steps with very different reliability, and the
distinction matters.

**Genogroup detection works for all genogroups.** It is a BLAST search against
whole genomes, independent of the typing alignments, and resolves any of
GI-GX. Non-norovirus caliciviruses are reported with the genogroup only, plus a
note.

**Genotype assignment is only validated for GI and GII.** Typing requires a
curated reference alignment, and only GI and GII have one. Every other
genogroup falls back to `all.aln`, a single alignment of all genogroups mixed
together, in which they are represented by a handful of sequences:

| Genogroup | ORF1 refs | ORF2 refs | Genotyping |
|---|---|---|---|
| **GI** | 41 | 34 | validated (leave-one-out) |
| **GII** | 121 | 96 | validated (leave-one-out) |
| GIII | 2 | 3 | returns a call, not validated |
| GIV | 4 | 12 | returns a call, not validated |
| GV | 2 | 2 | returns a call, not validated |
| GVI | 8 | 8 | returns a call, not validated |
| GVII | 2 | 2 | returns a call, not validated |
| GVIII | - | 2 | ORF2 only, not validated |
| GIX | - | 2 | ORF2 only, not validated |
| GX | 2 | 2 | returns a call, not validated |

A query outside GI/GII **will** return a genotype and a bootstrap percentage,
and that percentage may look convincing. Treat it as unverified: two or three
references cannot support a reliable assignment, and leave-one-out covered GI
and GII only. GVIII and GIX cannot be ORF1-typed at all.

In practice this is adequate for human and wastewater surveillance, which is
almost entirely GI, GII and GIX. It is **not** adequate for the animal
genogroups: GIII (bovine), GV (murine), GVI/GVII (canine), GX (bat).

### GII.4 variants

GII.4 is resolved to variant level, from 11 variants in the ORF2 reference set
(20 sequences): Bristol, US95_96, Camberwell, Farmington_Hills, Hunter, Asia,
Yerseke, Den_Haag, Osaka, New_Orleans, Sydney.

Leave-one-out assigns all 17 evaluable GII.4 references correctly, at 81-100%
support. This is narrower than RIVM's variant set (which includes
Apeldoorn_2007, absent here) and RIVM appends the emergence year
(`Sydney_2012`), which `norotyper` does not.

---

## Method

1. **Genogroup** - `blastn` against whole genomes.
2. **Region selection** - a region is typed only if the BLAST hit overlaps the
   annotated ORF1 or ORF2 by at least 100 nt.
3. **Alignment** - the query is cropped to the ORF window (+/- 300 nt) and added
   to the pre-aligned reference set with MAFFT
   (`--addfragments --6merpair --keeplength`), so reference columns never move.
4. **Tree** - Tamura-Nei distances, neighbour-joining.
5. **Bootstrap** - 100 replicates. A type is assigned when the query forms a
   clade with all references of exactly one type in >= 70% of replicates.

---

## Validation

### Leave-one-out on the bundled reference set

`python tests/loo_validation.py` (no internet needed). Each reference ORF is
typed against all other references of its genogroup; its own label is the
answer.

| Region | Queries | Genogroups | Distinct types | Assigned | Correct | Wrong |
|---|---|---|---|---|---|---|
| ORF2 capsid | 125 | GI, GII | 42 | 100% | 100% | 0 |
| ORF1 P-type | 147 | GI, GII | 49 | 93.9% | 93.2% | 1 |

Two caveats that matter:

* **The single error is a confident one.** `AB074893.1` (true `GII.PNA4`) was
  called `GII.P11` at **99% support**. The 70% threshold did not catch it, so
  bootstrap support must not be read as a guarantee of correctness.
* **Leave-one-out flatters the tool.** Queries are curated, full-length
  reference sequences. It shows the method is wired up correctly; it does not
  predict accuracy on partial, divergent or chimeric environmental contigs.

### RIVM comparison (not yet complete)

The answer set is RIVM Typing Tool calls for 4,392 GenBank accessions, 1,868 of
them GII.4. It is **not redistributed here** - see
[validation/README.md](validation/README.md) for how to obtain it.

```bash
python scripts/fetch_rivm_test_set.py rivm_test_set.fasta   # ~4,392 seqs from GenBank
norotyper rivm_test_set.fasta -o rivm_typed.tsv --threads 16
python scripts/validate_vs_rivm.py rivm_typed.tsv           # -> results/rivm_comparison.csv
```

The comparison reconciles the two labelling conventions before scoring: RIVM
puts plain `GII.4` in `capsid type` with the variant in `capsid subtype`, while
`norotyper` carries the variant inside the type. Genotype agreement is measured
on the genotype stem, and the GII.4 variant is scored separately.

### Not yet validated

Against RIVM's calls; on any genogroup outside GI/GII; on assembled
environmental or wastewater contigs.

---

## Performance

Cost is dominated by the bootstrap, not by alignment. Measured per query
(GII, 16 threads, one core's share):

| Step | Per query | Share |
|---|---|---|
| BLAST | ~42 ms | ~2% |
| MAFFT | ~40 ms | ~4% |
| **Bootstrap, 100 replicates** | **~900 ms** | **~95%** |

Within the bootstrap, the pairwise distance matrix is ~98% of the time and
neighbour-joining ~2%. Cost scales with the *square* of the genogroup's
reference count, so genogroup drives runtime: GII reference sets are about
three times the size of GI.

**Budget for GII, not GI.** An earlier figure of 0.5 s per genome was measured
on GI and does not generalise. Real surveillance data is overwhelmingly GII.
Install `numba`; without it the bootstrap is about 8x slower.

---

## Known differences from RIVM

* **Distance and tree**: Tamura-Nei computed analytically with
  neighbour-joining, rather than PAUP\* HKY85. Tie-breaking and distance
  estimation differ slightly.
* **Clade rule**: "assigned" means the query groups with all (up to 2)
  references of one type in >= 70% of replicates. RIVM's exact clade rule and
  its choice of representative references are not published; the subset of two
  per type was inferred from their description and from leave-one-out accuracy.
* **ORF windows** are the full annotated ORFs, not RIVM's fixed typing windows.
* **GII.4 variants** are a narrower set than RIVM's, and without the emergence
  year suffix.
* **Reference sets** are the 2019 (Chhabra) sets.

---

## Roadmap

Not implemented. Listed so the gap between this tool and a surveillance-ready
one is explicit.

* **Visual alignment report.** Per assignment, plot the query against the
  reference sequences of the called genotype, using a Bokeh sequence viewer
  along the lines of
  [Farrell's Bokeh sequence aligner](https://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner).
  One panel per typed region, so a call can be inspected rather than trusted.
* **Structured output.** JSON alongside the TSV, carrying the alignment and the
  ranked runner-up types that are currently computed and discarded.
* **Explicit recombinant reporting.** ORF1 and ORF2 are already typed
  independently; a `GII.4_Sydney[P16]`-style dual-type label and a
  recombinant flag would make the mismatch explicit.
* **Reference curation for the animal genogroups**, to lift GIII-GX above two
  or three sequences.
* **FASTQ entry point**, wrapping read QC and assembly.
* **Distribution** via PyPI and a bioconda recipe. Currently install is from
  this repository only.

---

## Repository layout

```
norotyper/            the maintained package
  data/               reference sequences and BLAST database (installed with the code)
scripts/              the original upstream pipeline, kept for provenance
tests/                unit tests, CLI tests, leave-one-out validation
validation/           taxonomy tables; RIVM comparison inputs go here (untracked)
results/              leave-one-out output
```

---

## Citation

### `norotyper`

```bibtex
@software{fitzpatrick_norotyper,
  author = {Amy Heather Fitzpatrick},
  title  = {norotyper: command-line norovirus dual-typing following the RIVM method},
  url    = {https://github.com/ahfitzpa/norotyper},
  year   = {2026}
}
```

### Method and classification standard

```bibtex
@article{kroneman2011automated,
  author  = {Kroneman, Annelies and Vennema, Harry and Deforche, Koen and
             van der Avoort, Harrie and Pe{\~n}aranda, Silvia and Oberste, M Steven
             and Vinj{\'e}, Jan and Koopmans, Marion},
  title   = {An automated genotyping tool for enteroviruses and noroviruses},
  journal = {Journal of Clinical Virology},
  volume  = {51},
  number  = {2},
  pages   = {121--125},
  year    = {2011},
  doi     = {10.1016/j.jcv.2011.03.006}
}

@article{chhabra2019updated,
  author  = {Chhabra, Preeti and de Graaf, Miranda and Parra, Gabriel I and
             Chan, Martin Chi-Wai and Green, Kim and Martella, Vito and Wang, Qiuhong
             and White, Patricia A and Katayama, Kazuhiko and Vinj{\'e}, Jan},
  title   = {Updated classification of norovirus genogroups and genotypes},
  journal = {Journal of General Virology},
  volume  = {100},
  number  = {10},
  pages   = {1393--1406},
  year    = {2019},
  doi     = {10.1099/jgv.0.001318}
}
```

---

## Attribution and licence

MIT licensed. This is a fork and substantial rewrite of
[zhuzhanji/NorovirusGenotyping](https://github.com/zhuzhanji/NorovirusGenotyping)
by Yang Ji, and it reimplements the method of the RIVM Norovirus Typing Tool.
See **[ATTRIBUTION.md](ATTRIBUTION.md)** for what came from where and an
important note on upstream licensing. `norotyper` is not produced by, endorsed
by, or affiliated with RIVM.
