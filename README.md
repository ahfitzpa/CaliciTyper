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
* **For reports:** `bokeh`, needed only for the figures in the HTML report.
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
| `orf1_runners` | up to 3 best-supported types, best first: `GII.P31:97.0\|GII.P16:2.0` |
| `orf2_type` | capsid type (variant-level for GII.4), `unassigned`, or empty |
| `orf2_support` | bootstrap % for the ORF2 call |
| `orf2_runners` | as `orf1_runners`, for ORF2 |
| `note` | why a region was not typed: no BLAST hit, `<100 nt` ORF overlap, or support below threshold |

### Example

| query | genogroup | orf1_type | orf1_support | orf1_runners | orf2_type | orf2_support | note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| NC_001959.2 | GI | GI.P1 | 100.0 | GI.P1:100.0\|GI.P10:0.0 | GI.1 | 100.0 | |
| NC_039476.1 | GII | GII.P16 | 72.0 | GII.P16:72.0\|GII.P1:0.0 | GII.2 | 100.0 | |
| frag_07 | GII | | | | GII.3 | 88.0 | orf1: <100 nt overlap |
| junk_01 | unassigned | | | | | | no BLAST hit (not a norovirus / Caliciviridae sequence in reference set) |

The second row is the case the runner-up column exists for: `GII.P16` cleared the
70% threshold by two points, and nothing else drew meaningful support.

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

reports (the TSV stays the primary output):
  --report PREFIX         write PREFIX.md, PREFIX.json and PREFIX.html too
  --json PATH             results as JSON, with runner-ups as structured lists
  --md PATH               Markdown report
  --html PATH             self-contained HTML with alignment panels and trees
  --report-max-panels N   most alignment panels to draw [default: 20]
  --no-tree               skip the per-region trees in the HTML report
```

### Reports

The TSV is the machine-readable output and never changes shape. For looking at
a run:

```bash
norotyper contigs.fasta -o typed.tsv --threads 8 --report run01
# -> typed.tsv, run01.md, run01.json, run01.html
```

* **`.json`** - every field from the TSV, with runner-ups as structured lists
  (`[{"type": "GII.P31", "support": 97.0}, ...]`) plus the run parameters.
* **`.md`** - the table and a genogroup tally, for pasting into a notebook or PR.
* **`.html`** - a self-contained page carrying the table plus two kinds of figure:
  * **Alignment panels.** Per call, the query against the references of the type
    it was assigned *and* the references of the runner-up type, drawn with
    [Bokeh](https://bokeh.org) following
    [Farrell's sequence aligner](https://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner).
    Pan and zoom; letters appear once legible. A bootstrap percentage is not
    inspectable, but the alignment it came from is.
  * **Trees.** Per region and genogroup, a neighbour-joining tree placing the
    queries (highlighted) among the references, so you can see *where* a run
    landed rather than reading it row by row.

Panels are capped at `--report-max-panels` (highest support first) and trees at
120 taxa, with a note in the report when either truncates. `bokeh` is needed
for the figures; without it the Markdown and JSON are unaffected and the HTML is
written without them.

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

## Datasets

Everything the tool ships with, and everything it has been tested on.

### Bundled reference data (`norotyper/data/`)

Installed with the package; no download step.

**BLAST database for genogroup detection** - 60 whole genomes, deliberately
including other genera so that non-noroviruses are *rejected* rather than
forced into a genogroup:

| | Sequences |
|---|---|
| Norovirus GI-GX | 30 (GI 5, GII 8, GIX 3, GVI 4, GIV 2, GV 2, GVII 2, GX 2, GIII 1, GVIII 1) |
| Other Caliciviridae | 30 (Sapovirus GI-GXII, Lagovirus, Vesivirus, Nebovirus, Nacovirus, Recovirus, Valovirus, and unclassified) |

**Typing reference alignments** - norovirus only, 2019 (Chhabra) sets. GI and
GII have their own curated alignments; all other genogroups fall back to a
mixed `all.aln`:

| Alignment | ORF1 (RdRp) | ORF2 (VP1) |
|---|---|---|
| `GI.aln` | 41 seqs, 806 cols | 34 seqs, 1,672 cols |
| `GII.aln` | 121 seqs, 873 cols | 96 seqs, 1,719 cols |
| `all.aln` (fallback) | 185 seqs, 874 cols | 166 seqs, 1,973 cols |

GII.4 is labelled at variant level: 11 variants, 20 sequences (Bristol,
US95_96, Camberwell, Farmington_Hills, Hunter, Asia, Yerseke, Den_Haag, Osaka,
New_Orleans, Sydney).

### Test set 1: leave-one-out on the bundled references

Source: the bundled alignments themselves. Each reference ORF is typed against
all other references of its genogroup; its own label is the answer. Only types
with >= 2 references are evaluable. No internet needed.

| | ORF1 | ORF2 |
|---|---|---|
| Queries | 147 | 125 |
| Genogroups | GI (36), GII (111) | GI (34), GII (91) |
| Distinct types | 49 | 42 |

**This set is not independent.** Queries are drawn from the same curated,
full-length sequences the tool types against. It demonstrates the method is
wired up correctly; it does not estimate field accuracy.

### Test set 2: RIVM Norovirus Typing Tool calls

Source: 5,986 GenBank accessions with calls from the RIVM tool, in
`validation/results_rivm.csv`. **Not redistributed here** - see
[validation/README.md](validation/README.md). Sequences are fetched from NCBI
with `scripts/fetch_rivm_test_set.py`.

This is the independent set, and it is used two ways:

| Subset | n | Used for |
|---|---|---|
| Records with a capsid or polymerase call | 4,392 | genotype agreement (default fetch) |
| **All records** | **5,986** | genogroup detection (`--all`) |

Composition of the full 5,986:

| | n |
|---|---|
| **Norovirus** | **4,295** |
| &nbsp;&nbsp;GII | 3,771 |
| &nbsp;&nbsp;GI | 339 |
| &nbsp;&nbsp;GV | 100 |
| &nbsp;&nbsp;GIV | 51 |
| &nbsp;&nbsp;GIII | 29 |
| &nbsp;&nbsp;GVII | 5 |
| **Non-norovirus Caliciviridae** | **1,691** |
| &nbsp;&nbsp;Lagovirus (RHDV + unclassified) | 765 |
| &nbsp;&nbsp;Sapovirus (GI-GIV and others) | ~420 |
| &nbsp;&nbsp;Vesivirus (feline calicivirus and others) | ~190 |
| &nbsp;&nbsp;Unassigned / other genera | ~320 |

Of the capsid calls, 1,868 are GII.4 - 42% of the typed subset - which is why
the GII.4 variant handling matters. Note the genogroup evaluation does reach
the animal genogroups (GV 100, GIV 51, GIII 29), even though *genotyping* for
those is unvalidated.

### Not tested on anything else

In particular: no assembled environmental or wastewater contigs, no partial or
chimeric sequences from a real run, and no sequencing reads. Every number
below comes from curated GenBank records.

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

The comparison reports four things:

| Row | What it measures |
|---|---|
| **Genogroup** | every matched record, including the non-norovirus caliciviruses |
| **Genogroup (norovirus only)** | records RIVM calls norovirus at all |
| **ORF2 capsid** / **ORF1 P-type** | genotype agreement, variant suffixes stripped |
| **GII.4 variant** | against RIVM's separate `capsid subtype` column |

Genogroup is scored over a wider set than genotyping: the BLAST database holds
30 norovirus references (GI-GX) and 30 other Caliciviridae (Sapovirus,
Lagovirus, Vesivirus, Nebovirus and others), so detection is doing real
discrimination and must also *reject* non-noroviruses correctly. Use
`--all` on the fetch script to pull all 5,986 RIVM records rather than only the
4,392 carrying a type call, which is what that evaluation needs.

Label conventions are reconciled before scoring. RIVM puts plain `GII.4` in
`capsid type` with the variant in `capsid subtype`, while `norotyper` carries
the variant inside the type; and RIVM records genogroup as a lineage string
(`Caliciviridae Norovirus GII`) where `norotyper` reports a bare `GII`.

### Not yet validated

Against RIVM's calls (genogroup or genotype -- the run has not completed); on
genotyping outside GI/GII; on assembled environmental or wastewater contigs.

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
