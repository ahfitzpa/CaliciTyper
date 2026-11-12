<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/calicityper-lockup-horizontal-dark.png">
    <img src="assets/calicityper-lockup-horizontal.png" width="420" alt="CaliciTyper">
  </picture>
</p>

<p align="center">Command-line typing for <em>Caliciviridae</em></p>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)
[![CI tested](https://img.shields.io/badge/CI-3.10%20%7C%203.12-blue.svg)](.github/workflows/ci.yml)
[![Classification](https://img.shields.io/badge/classification-Chhabra_2019-green.svg)](https://doi.org/10.1099/jgv.0.001318)
[![Status](https://img.shields.io/badge/status-research_code-orange.svg)](BENCHMARKING.md)

# CaliciTyper: command-line calicivirus genotyping and recombination analysis

CaliciTyper (command: `calicityper`) assigns norovirus genogroup, ORF1 polymerase (P-) type and ORF2
capsid type to assembled contigs and whole genomes, replicating the published
method of the [RIVM Norovirus Typing
Tool](https://www.rivm.nl/mpf/typingtool/norovirus/). It also scans norovirus
genomes for recombination and reports evidence-weighted events rather than a
yes/no answer.

Every other norovirus typing tool we are aware of -- the RIVM tool, CDC's
HuCaT -- is web-only. 

### What works today

| Capability | Scope | State |
|---|---|---|
| Genogroup detection | all Caliciviridae genera | working |
| Norovirus genotyping (ORF1 + ORF2) | GI, GII | validated |
| Norovirus genotyping | GIII-GX | returns a call, **not validated** |
| GII.4 variant resolution | 11 variants | validated |
| Putative novel ORF1/ORF2 types (RdRp P-type, VP1 genotype / genogroup) | norovirus | working, see [caveats](#putative-novel-orf1orf2-types) |
| Recombination analysis | norovirus, sapovirus | working, see caveats |
| **Sapovirus genotyping** | GI-GXIX | **not implemented** -- tables drafted, reference sets not built |

---

## Installation

### Conda / Mamba (recommended)

```bash
git clone https://github.com/ahfitzpa/calicityper.git
cd calicityper
conda env create -f environment.yml
conda activate calicityper
```

### pip, with dependencies already present

Needs `blast` and `mafft` on `PATH`:

```bash
pip install git+https://github.com/ahfitzpa/calicityper.git
```

The reference data currently ships inside the package, so this gives a working
install with no database build step. That is changing: see
[Reference data](#reference-data).

### OpenRDP (only for `calicityper recomb`)

OpenRDP pins `numpy<2` while `calicityper` needs `numpy>=2`. Install it separately. 

```bash
conda create -n openrdp -c conda-forge python=3.11 "numpy<2" "scipy>=1.5" "h5py>=3.8,<3.11"
conda activate openrdp
git clone https://github.com/PoonLab/OpenRDP.git && pip install ./OpenRDP
conda deactivate

export CALICITYPER_OPENRDP=$(conda info --base)/envs/openrdp/bin/openrdp
```

OpenRDP is unpackaged (no PyPI or bioconda release). 

## Quick start

Type the four bundled example genomes (three noroviruses, one of them flagged
as a putative novel RdRp P-type, and one sapovirus):

```bash
calicityper type examples/example.fasta -o example_typed.tsv --threads 4 --report example
```

This takes under a minute. Compare your `example_typed.tsv` with
[examples/example_typed.tsv](examples/example_typed.tsv); they should match
exactly. 

## Test

```bash
pytest tests/ -v                      # unit, CLI and regression tests
python tests/loo_validation.py        # leave-one-out, no internet needed
python tests/putative_novel_validation.py   # putative novel ORF1/ORF2 type screen, leave-one-out
```

---

## Running CaliciTyper

### Simple taxonomy assignment

```bash
% calicityper type contigs.fasta -o typed.tsv --threads 8
```

### Taxonomy assignment with visual report

```bash
# also writes run01.md, run01.json and run01.html
% calicityper type contigs.fasta -o typed.tsv --threads 8 --report run01
```

For very large runs cap alignments with
`--report-max-panels`, and use `--html-cdn` to link Bokeh instead of embedding. 

### Recombination, from an alignment

```bash
% calicityper recomb --alignment panel.fasta --query SAMPLE_01 -o recomb_out/
```

### Recombination, from unaligned sequences

```bash
# aligned with MAFFT first, then scanned
% calicityper recomb --input samples.fasta -o recomb_out/ --threads 8
```

### Overview to CaliciTyper commands

```bash
% calicityper recomb --alignment panel.fasta --query SAMPLE_01 --dry-run
alignment        : panel.fasta
  sequences      : 9
  columns        : 7500
  informative    : 2166
queries          : 1
methods requested: rdp geneconv bootscan maxchi chimaera siscan threeseq
OpenRDP          : /path/to/openrdp
  SAMPLE_01: OK
```

---

## Commands

### `calicityper type`

Assigns genus for Caliciviridae and for performs dual-genotyping for norovirus and sapovirus. 

#### Output files

| File | Description |
|---|---|
| `typed.tsv` | primary machine-readable output, one row per input sequence |
| `PREFIX.json` | the same fields, with runner-ups as structured lists, plus run parameters |
| `PREFIX.md` | table and genogroup tally, for a notebook or a PR |
| `PREFIX.html` | self-contained page with alignment panels and trees |

#### Output table fields

| Column | Description |
|---|---|
| `query` | FASTA record id |
| `genogroup` | `GI`-`GX`, a non-norovirus lineage, or `unassigned` |
| `orf1_type` | P-type, `unassigned`, or empty if the region was not typeable |
| `orf1_support` | bootstrap % for the ORF1 call |
| `orf1_runners` | up to 3 best-supported types, best first: `GII.P31:97.0\|GII.P16:2.0` |
| `orf2_type` | capsid type (variant-level for GII.4) |
| `orf2_support` | bootstrap % for the ORF2 call |
| `orf2_runners` | as `orf1_runners` |
| `note` | why a region was not typed: no BLAST hit, `<100 nt` ORF overlap, support below threshold |
| `orf1_putative_novel` | putative novel ORF1 (RdRp) type: `known`, `putative_novel_type` (P-type), `putative_novel_group` (P-group), `not_assessable` (see [Putative novel ORF1/ORF2 types](#putative-novel-orf1orf2-types)) |
| `orf1_putative_novel_detail` | the distance, the nearest cluster and its bound, or why it was not assessable |
| `orf2_putative_novel` | as `orf1_putative_novel`, for ORF2 (VP1): genotype and genogroup |
| `orf2_putative_novel_detail` | as `orf1_putative_novel_detail` |
| `orf1_window_cov` / `orf2_window_cov` | nucleotides of the RIVM typing window covered, as `covered/window` |
| `gii4_window_cov` | as above for the GII.4 variant window; empty unless a variant was in play |

Example, the real output for [examples/example.fasta](examples/example.fasta)
(runner-up columns shortened to the top two; detail columns omitted):

| query | genogroup | orf1_type | orf1_support | orf1_runners | orf1_putative_novel | orf2_type | orf2_support | orf2_putative_novel | note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| MW305639.1 | GII | GII.P16 | 100.0 | `GII.P16:100.0\|GII.P1:0.0` | known | GII.4_Sydney | 100.0 | known | |
| MT451991.1 | GI | GI.P13 | 100.0 | `GI.P13:100.0\|GI.P1:0.0` | known | GI.3 | 100.0 | known | |
| MN509084.1 | Caliciviridae_Sapovirus_GII | | | | not_assessable | | | not_assessable | other Caliciviridae: genogroup only |
| AB242258.2 | GII | GII.P21 | 100.0 | `GII.P21:100.0\|GII.P1:0.0` | **putative_novel_type** | GII.3 | 100.0 | known | |

#### Borderline flag: AB242258.2

**AB242258.2** types as **GII.P21** at 100%, matching RIVM, but its RdRp is 0.047 nt from the GII.P21 references—just outside the cluster bound of 0.044.

GII.P21 references are very close together, so the cluster bound is tight. Typing asks which type the tree supports; the screen asks whether the sequence falls outside every known cluster. A clean-typing sequence this close to the bound is the screen’s main false-positive mode.

#### Options

| Argument | Type | Description |
|---|---|---|
| `fasta` | path | input FASTA (positional, required) |
| `-o, --out` | path | output TSV, or `-` for stdout [`-`] |
| `--bootstraps` | int | bootstrap replicates [100, RIVM's value] |
| `--min-support` | float | minimum bootstrap % to assign a type [70] |
| `--max-refs-per-type` | int | references per type used in the tree [2] |
| `--threads` | int | worker processes [1] |
| `--seed` | int | bootstrap seed, for reproducibility [1] |
| `--report` | prefix | write `.md`, `.json` and `.html` as well |
| `--json` / `--md` / `--html` | path | write one format individually |
| `--report-max-panels` | int | cap on alignment panels, in input order [no cap: every call] |
| `--no-tree` | flag | skip the per-region trees in the HTML |

### `calicityper recomb`

Detects and scores recombination in **norovirus or sapovirus** (`--virus`).

This command works on an alignment you supply, so it is independent of whether a virus is
genotypeable.

`--virus` selects the genome model used to judge whether a breakpoint is in a
plausible place. 

| | Genome organisation | Recombination hotspot |
|---|---|---|
| norovirus | ORF1 nonstructural, ORF2 VP1, ORF3 VP2 | ORF1/ORF2 junction |
| sapovirus | **ORF1 carries NS proteins *and* VP1**; ORF2 is VP2 | **RdRp(NS7)/VP1 junction, inside ORF1** |

#### Output files

| File | Description |
|---|---|
| `recombination.json` | full result: events, breakpoints, score components, provenance |
| `recombination.tsv` | one row per candidate event |
| `aligned.fasta` | written only when `--input` was used |
| `work/` | the sanitised alignment and raw OpenRDP CSV |

#### Status labels

| Status | Rule |
|---|---|
| `none` | no method passes the significance threshold |
| `weak` | one method passes, evidence limited |
| `possible` | two or more methods agree, **or** one method with strong phylogenetic support and distinct parents |
| `likely` | three or more methods agree **and** segment phylogenies are incongruent |
| `strong` | as `likely`, plus clearly distinct parental lineages |
| `insufficient_sequence` | too short, too ambiguous, or too few informative sites |

`strong` is never reached on one method, and a short or ambiguous sequence
returns `insufficient_sequence` rather than a forced call.

#### Scoring

```
score = 0.35 * method_support
      + 0.20 * phylogenetic_support
      + 0.20 * breakpoint_consistency
      + 0.15 * parental_distinctness
      + 0.10 * biological_plausibility
```

All weights and thresholds are configurable; see
[`calicityper/recomb_config.py`](calicityper/recomb_config.py) for the defaults
and pass your own with `--config`.

#### Options

| Argument | Type | Description |
|---|---|---|
| `--alignment` | path | aligned FASTA with the query and a reference panel |
| `--input` | path | unaligned FASTA; aligned with MAFFT first |
| `--query` | id | query id; repeat for several, omit for all |
| `-o, --output` | dir | output directory [`recomb_out`] |
| `--reference-metadata` | path | TSV/CSV/JSON mapping reference ids to lineages |
| `--config` | path | YAML or TOML configuration |
| `--methods` | list | subset of rdp geneconv bootscan maxchi chimaera siscan threeseq |
| `--min-length` | int | override `min_sequence_length` [5000] |
| `--pvalue-threshold` | float | override `pvalue_threshold` [0.05] |
| `--max-breakpoint-distance` | int | clustering distance [250] |
| `--skip-rdp` / `--skip-gard` | flag | skip the OpenRDP scanners / GARD |
| `--no-segment-trees` | flag | skip the phylogeny check (faster, loses a component) |
| `--threads` | int | threads for the MAFFT step [1] |
| `--dry-run` | flag | validate inputs, report what would run, stop |
| `--verbose` | flag | progress to stderr |

---

## Input must be a fasta file

`calicityper` types **contigs and genomes, not raw sequencing reads.**. 


Recombination analysis additionally needs **near-complete genomes**: the
default `min_sequence_length` is 5,000 nt and `min_informative_sites` 500.
Most wastewater contigs will correctly return `insufficient_sequence`.

---

## Genogroup and genotype coverage

Two different steps, with very different reliability.

**Genogroup detection works across Caliciviridae.** It is a BLAST search
against 60 whole genomes -- 30 norovirus GI-GX and 30 other genera (Sapovirus
GI-GXII, Lagovirus, Vesivirus, Nebovirus, Nacovirus, Recovirus, Valovirus). 

**Genotype assignment is validated for norovirus GI and GII only.** 

| Genogroup | ORF1 refs | ORF2 refs | Genotyping |
|---|---|---|---|
| **GI** | 41 | 34 | validated |
| **GII** | 121 | 96 | validated |
| GIII | 2 | 3 | returns a call, not validated |
| GIV | 4 | 12 | returns a call, not validated |
| GV | 2 | 2 | returns a call, not validated |
| GVI | 8 | 8 | returns a call, not validated |
| GVII | 2 | 2 | returns a call, not validated |
| GVIII | - | 2 | ORF2 only, not validated |
| GIX | - | 2 | ORF2 only, not validated |
| GX | 2 | 2 | returns a call, not validated |

### GII.4 variants

GII.4 resolves to variant level from 11 variants (20 sequences): Bristol,
US95_96, Camberwell, Farmington_Hills, Hunter, Asia, Yerseke, Den_Haag, Osaka,
New_Orleans, Sydney. The variant is carried **inside** `orf2_type`
(`GII.4_Sydney`); strip the suffix for genotype-level labels, as RIVM reports
the variant in a separate column.

### Typing windows and when a call is withheld

Three windows decide whether the evidence supports a call. Each is a set of
columns in the reference alignment, so an insertion in a few references cannot
shift a window the way genomic coordinates would.

| Window | Region | Length | Used for |
|---|---|---|---|
| **VP1 region C** | 5' end of VP1 | 252 nt (GII), 264 (GI), 267 (GI.9), 255 (GIX) | capsid genotype |
| **RdRp 3' window** | 3' end of RdRp | 172 nt (173 for GIV) | P-type |
| **GII.4 variant window** | P2 to the VP1 C-terminus | 789 nt | GII.4 variant |

The first two are the RIVM/NoroNet operational windows. The third exists
because region C is the VP1 **shell** domain, which is conserved *between*
GII.4 variants: measured on the bundled references, only **4 of 135
variant-diagnostic columns** fall inside region C, so a region C amplicon
cannot support a variant call however long it is. P2 alone (aa 279-405) holds
41% of those columns; P2 through the C-terminus holds 70%, which is the window
used.

Three rules follow:

- **Input shorter than 172 nt is rejected** before BLAST, with the reason. No
  window can be covered, so nothing can be typed.
- **A call is withheld when its window is under 80% covered.** A 3' VP1
  fragment that misses region C entirely gets no capsid type, even if it
  matches a reference closely — RIVM could not type that input either.
- **A GII.4 variant is reported only when the variant window is covered.**
  Otherwise the genotype stem is returned: `GII.4`, not `GII.4_Sydney`.

Coverage is reported per call in `orf1_window_cov`, `orf2_window_cov` and
`gii4_window_cov` as `covered/window`, so a withheld call can be understood
without rerunning. The 80% floor is set by `CALICITYPER_MIN_WINDOW_COVERAGE`
and is provisional until the truncation leave-one-out curve fixes it.

Worked example, one GII.4 Sydney genome cut three ways:

| Input | region C | variant window | Call |
|---|---|---|---|
| region C only (360 nt) | 241/252 | — | `unassigned` (bootstrap 69% < 70%) |
| region C + P2 (1,340 nt) | 241/252 | 508/789 (64%) | **`GII.4`** — variant withheld |
| complete VP1 (1,690 nt) | 241/252 | 775/789 (98%) | `GII.4_Sydney` |

### Genotyping is restricted to GI and GII

Sequences from the animal genogroups (GIII-GX) and from other *Caliciviridae*
are reported to genogroup and stop there. Their reference sets hold two or
three sequences each -- enough to return a label, not enough for it to mean
anything -- and typing them spends the bootstrap, ~95% of runtime, on a call
this README would tell you not to trust. `CALICITYPER_TYPE_ALL_GROUPS=1`
restores the old behaviour for reference-set development.

### Sapovirus typing: in preparation, not yet available

Reference tables are built directly from Supplemental Table 1 of Zhao *et al.*
2024 (BMC Genomics 25:1131) by
[`scripts/build_sapovirus_taxonomy.py`](scripts/build_sapovirus_taxonomy.py):
53 VP1 genotypes over 156 strains and 45 RdRp P-types over 139 strains, in
[`validation/sapovirus_vp1_taxonomy.csv`](validation/sapovirus_vp1_taxonomy.csv)
and
[`validation/sapovirus_rdrp_taxonomy.csv`](validation/sapovirus_rdrp_taxonomy.csv).

Generating them from the spreadsheet rather than transcribing matters: an
earlier hand-transcription put six accessions in the wrong genotype or changed
a character, and a test now fails if the CSVs drift from the published table.
Every accession was also checked against NCBI, which found four errors in the
published table itself (`L504335`, `KC215898`, `KC469052` and `LC215696` are
not sapovirus records; the `note` column gives each correction) and one that
cannot be resolved: `KC309412`, listed for VP1 genogroup GXXVI, is a 366 nt
*Pseudomonas aeruginosa* 16S record. It stays flagged and is excluded from any
reference set.

**Sapovirus genotyping is not implemented.** The tables are step one of five;
the sequences have not had their ORFs extracted, reference sets are not built,
and nothing in the typing path is sapovirus-aware. Sapovirus input is currently
detected to genus and genogroup only.
---

## Method

1. **Genogroup** -- `blastn` against whole genomes.
2. **Region selection** -- a region is typed only if the BLAST hit overlaps the
   annotated ORF1 or ORF2 by at least 100 nt.
3. **Alignment** -- the query is cropped to the ORF window (+/- 300 nt) and
   added to the pre-aligned references with MAFFT
   (`--addfragments --6merpair --keeplength`), so reference columns never move.
4. **Tree** -- Tamura-Nei distances, neighbour-joining.
5. **Bootstrap** -- 100 replicates; a type is assigned when the query forms a
   clade with all references of exactly one type in >= 70% of them.

### Putative novel ORF1/ORF2 types

Every typed region is also checked for whether it is a **putative novel ORF1
(RdRp) or ORF2 (VP1) type**: whether it falls outside every known genotype or
genogroup, using the classification rule of Chhabra *et al.*
2019. That paper assigns capsid genotypes and genogroups on uncorrected
p-distances of the **complete VP1 amino acid** sequence and P-types on
**partial RdRp nucleotide** sequence, and recognises a new cluster when "the
average distance between all sequences within a new genogroup or genotype and
its nearest cluster(s) should not overlap within two standard deviations".

Applied to one query:

1. Each reference genotype and genogroup has a bound: the mean distance between
   its members + 2 SD.
2. The query is placed by its **average distance to the members** of each
   cluster.
3. It is `known` if that average is inside some genotype's bound;
   `putative_novel_type` (a putative novel RdRp P-type or VP1 genotype) if it
   is outside every genotype's; `putative_novel_group` (P-group or VP1
   genogroup) if it is also outside every genogroup's.
4. It is `not_assessable` if it lacks a complete
   VP1 open reading frame (ATG to stop, no internal stop, within 10% of the
   reference lengths), or at least 95% of the RdRp typing region.

Clusters with fewer than three references, and the tentative (NA) clusters,
have no SD of their own; they use the pooled distribution of all clusters at
that level, with far outliers dropped. Genogroups always use the pooled bound,
because their internal diversity varies too much between so few groups (GVI
alone would otherwise absorb most queries). The per-cluster table ships as
`calicityper/data/putative_novel_bounds.tsv` and is rebuilt by
`scripts/build_putative_novel_bounds.py`.

**Read a flag as "putative"** The paper names
a cluster from a single sequence, or from one location, as not assigned (NA).
A flag can also sit just outside a very tight cluster's bound, as in the
AB242258.2 example above, or point at a problem in the reference set itself:
an earlier flag on Chhabra's own GII.P8 reference strain exposed two corrupted
RdRp references, since fixed (see [BENCHMARKING.md](BENCHMARKING.md#reference-set-correction-061)).

How well it works, from leaving each reference out in turn
([BENCHMARKING.md](BENCHMARKING.md#6-putative-novel-orf1orf2-types)):

| | Known type correctly `known` | Removed type correctly flagged |
|---|---|---|
| VP1 genotype | 94% | 95% |
| RdRp P-type | 96% | 50% |
| VP1 genogroup | | 18% |


### Recombination

For recombination: sanitise and validate the alignment, run the OpenRDP
scanners, cluster hits into candidate events by breakpoint proximity, infer
parents, test segment phylogenies either side of each breakpoint, score, and
assign a status.

---

## Validation

Full detail, including what each test does *not* show, is in
**[BENCHMARKING.md](BENCHMARKING.md)**.

| Region | Queries | Assigned | Correct | Wrong |
|---|---|---|---|---|
| ORF2 capsid | 125 | 100% | 100% | 0 |
| ORF1 P-type | 147 | 97.3% | 96.6% | 1 |

Two caveats that matter:

- **Leave-one-out flatters the tool.** Queries are the same curated,
  full-length references the tool types against.
- **The one error is a confident error.** `AB074893.1` (true `GII.PNA4`) was
  called `GII.P11` at **100% support**. Bootstrap support is not a guarantee of
  correctness.

The comparison against RIVM's own calls for 4,392 accessions has **not yet
completed**, and the genogroup comparison over all 5,986 is **queued**.

---

## Known differences from RIVM

- **Distance and tree**: Tamura-Nei computed analytically with
  neighbour-joining, not PAUP\* HKY85. Tie-breaking differs slightly.
- **Clade rule**: the query must group with all (up to 2) references of one
  type in >= 70% of replicates. RIVM's exact rule is not published; the subset
  of two per type was inferred from their description and from leave-one-out.
- **ORF windows** are the full annotated ORFs, not RIVM's fixed typing windows.
- **GII.4 variants** are a narrower set, without the emergence year suffix.

`calicityper` is an independent reimplementation. It is **not** produced by,
endorsed by, or affiliated with RIVM, and its calls are not guaranteed to match
the RIVM web tool.

---

## FAQ

**Can it type sapovirus?** No. It detects sapovirus to genogroup only. The
genotype tables are transcribed and verified, but the reference sets have not
been built and the typing path is not sapovirus-aware.

**But `recomb --virus sapovirus` works?** Yes, and the two are independent:
`recomb` scores an user provided alignment. 

**Can it take FASTQ?** No. Assemble or cluster first.

**Why did recombination return `insufficient_sequence`?** The sequence is under
5,000 nt or has fewer than 500 informative sites. The classical scanners cannot
say anything useful below that.

**Why is a GIII call untrustworthy when it shows 98% support?** Support
measures how often the bootstrap recovers the same grouping, not whether the
reference set contains the right answer.

**Do I need OpenRDP?** Only for `recomb`. Without it the command still runs and
reports every method it skipped.

---

## Dependencies

**Mandatory**

- **Python** >= 3.9, with biopython, numpy, pandas
- **[BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)** -- genogroup assignment
- **[MAFFT](https://mafft.cbrc.jp/alignment/software/)** -- aligns queries onto the fixed reference alignment

**Strongly recommended**

- **[numba](https://numba.pydata.org/)** -- nominally optional, but the bootstrap is ~8x slower without it

**Optional**

- **[bokeh](https://bokeh.org)** -- figures in the HTML report
- **[OpenRDP](https://github.com/PoonLab/OpenRDP)** -- the seven recombination scanners
- **clustalw** -- only for the legacy `--backend clustalw` path

External programs are found on `PATH`, or via `CALICITYPER_BLASTN`,
`CALICITYPER_MAFFT`, `CALICITYPER_CLUSTALW`, `CALICITYPER_OPENRDP`. Reference
data location is overridable with `CALICITYPER_DATA`.

---

## Roadmap

- Build the sapovirus reference sets from the drafted taxonomy tables.
- HyPhy GARD as independent phylogenetic confirmation of breakpoints.
- Explicit recombinant labelling (`GII.4_Sydney[P16]`) and a recombinant flag
  in `calicityper type`.
- Reference curation for the animal genogroups, to lift GIII-GX above two or
  three sequences.
- Distribution via PyPI and bioconda.

---

## Reference data

The reference set is **309 GenBank records** across 10 genogroups, held as
pre-computed alignments for the VP1 and RdRp typing regions plus the BLAST
database used for genogroup assignment. The accessions are listed in
[`calicityper/data/reference_accessions.tsv`](calicityper/data/reference_accessions.tsv).

It is distributed as a **Zenodo dataset**, not inside this repository:

```bash
calicityper setup                      # fetch and verify, into ~/.cache/calicityper
calicityper setup --verify-only        # check an existing copy
calicityper setup --source set.tar.gz --data-dir /shared/calicityper-ref
```

`CALICITYPER_DATA` points the tool at a particular copy, so one download can be
shared across a cluster. Fetch on a login node if compute nodes have no
outbound network.

**Why the built alignments are distributed, not just the accessions.** Typing
compares a query against reference alignment *columns*, and all three typing
windows are column ranges: region C is "the first 252 covered columns of VP1".
An alignment rebuilt with a different aligner version can shift those windows
and change calls with no other symptom. Shipping the built set removes that
variable, and
[`reference_manifest.json`](calicityper/data/reference_manifest.json) holds a
SHA-256 for each of the 84 artefacts. `calicityper setup` verifies every one and
**refuses to install a set that does not match** — a reference set that is not
the one the manifest describes is an error, not a warning.

Each Zenodo version has its own DOI, recorded in the manifest, so a result can
name the exact reference set that produced it.

**A note on restricted data.** Zenodo is still public redistribution. This path
suits the public GenBank records above. A set built from a restricted
surveillance source such as CaliciNet must not be deposited on Zenodo or
GitHub: reference it by accession and have users build it locally, against the
same manifest discipline.

## Resources

### Papers

- **Typing method** this tool reimplements: Kroneman *et al.* 2011, *J Clin
  Virol* 51:121-125. [doi:10.1016/j.jcv.2011.03.006](https://doi.org/10.1016/j.jcv.2011.03.006)
- **Norovirus classification** (genogroups, P-types, capsid types; the
  reference strains bundled here): Chhabra *et al.* 2019, *J Gen Virol*
  100:1393-1406. [doi:10.1099/jgv.0.001318](https://doi.org/10.1099/jgv.0.001318)
- **Putative novel type criterion** (2x SD rule, NA clusters): as above, Chhabra *et al.*
  2019, Methods and Results.
- **Sapovirus dual typing** (VP1 and RdRp genogroups, genotypes and their
  cut-offs): Zhao *et al.* 2024, *BMC Genomics* 25:1131.
  [doi:10.1186/s12864-024-11048-w](https://doi.org/10.1186/s12864-024-11048-w)

### Related tools

- [RIVM Norovirus Typing Tool](https://www.rivm.nl/mpf/typingtool/norovirus/):
  the web tool whose method and calls `calicityper` is benchmarked against.
- [OpenRDP](https://github.com/PoonLab/OpenRDP): the recombination scanners
  behind `calicityper recomb`.
- [zhuzhanji/NorovirusGenotyping](https://github.com/zhuzhanji/NorovirusGenotyping):
  the upstream project this repository was forked from.

### Datasets

| Data | Where | Notes |
|---|---|---|
| Typing reference sequences | `calicityper/data/` (installed with the package) | from Chhabra *et al.* 2019 |
| Example input and outputs | [examples/](examples/) | three public GenBank genomes |
| RIVM benchmark accession list | [validation/rivm_accessions.tsv](validation/rivm_accessions.tsv) | 5,986 versioned GenBank accessions; `genotype_set` marks the 4,392 used for the genotype comparison |
| RIVM's calls for those accessions | not redistributed; see [validation/README.md](validation/README.md) | required to score a benchmark run |

The benchmark sequences themselves are not stored here. They are public
GenBank records, so the accession list is enough to download them;
[BENCHMARKING.md](BENCHMARKING.md) has the commands and the results.

## Citation

```bibtex
@software{fitzpatrick_calicityper,
  author = {Amy Heather Fitzpatrick},
  title  = {calicityper: command-line calicivirus genotyping and recombination analysis},
  url    = {https://github.com/ahfitzpa/calicityper},
  year   = {2026}
}
```

Please also cite the method and the classification standard:

> Kroneman A, Vennema H, Deforche K, *et al.* An automated genotyping tool for
> enteroviruses and noroviruses. *J Clin Virol* 2011;51(2):121-125.
> [doi:10.1016/j.jcv.2011.03.006](https://doi.org/10.1016/j.jcv.2011.03.006)

> Chhabra P, de Graaf M, Parra GI, *et al.* Updated classification of norovirus
> genogroups and genotypes. *J Gen Virol* 2019;100(10):1393-1406.
> [doi:10.1099/jgv.0.001318](https://doi.org/10.1099/jgv.0.001318)

If you use `recomb`, cite [OpenRDP](https://github.com/PoonLab/OpenRDP) and the
underlying methods.

## Licence

MIT. See [LICENSE](LICENSE).

> The CaliciTyper name and logo are © A. H. Fitzpatrick and are not covered by
> the software licence. See [assets/NOTICE.md](assets/NOTICE.md).

## Attribution

A fork and substantial rewrite of
[zhuzhanji/NorovirusGenotyping](https://github.com/zhuzhanji/NorovirusGenotyping)
by Yang Ji. See **[ATTRIBUTION.md](ATTRIBUTION.md)** for what came from where
and an important note on upstream licensing.

## Author

Amy Fitzpatrick, University College Dublin.
