[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![Classification](https://img.shields.io/badge/classification-Chhabra_2019-green.svg)](https://doi.org/10.1099/jgv.0.001318)
[![Status](https://img.shields.io/badge/status-research_code-orange.svg)](BENCHMARKING.md)

# calicityper: command-line calicivirus genotyping and norovirus recombination analysis

`calicityper` assigns norovirus genogroup, ORF1 polymerase (P-) type and ORF2
capsid type to assembled contigs and whole genomes, replicating the published
method of the [RIVM Norovirus Typing
Tool](https://www.rivm.nl/mpf/typingtool/norovirus/). It also scans norovirus
genomes for recombination and reports evidence-weighted events rather than a
yes/no answer.

Every other norovirus typing tool we are aware of -- the RIVM tool, CDC's
HuCaT -- is web-only. `calicityper` is the scriptable equivalent: no uploading
surveillance sequences to a third-party form, reproducible with a fixed seed,
and it drops into a pipeline.

> **Status: research code.** Genotyping is validated for norovirus GI and GII
> only. Read [Genogroup and genotype coverage](#genogroup-and-genotype-coverage)
> and [BENCHMARKING.md](BENCHMARKING.md) before trusting a call.

### What works today

| Capability | Scope | State |
|---|---|---|
| Genogroup detection | all Caliciviridae genera | working |
| Norovirus genotyping (ORF1 + ORF2) | GI, GII | validated |
| Norovirus genotyping | GIII-GX | returns a call, **not validated** |
| GII.4 variant resolution | 11 variants | validated |
| Recombination analysis | **norovirus only** | working, see caveats |
| **Sapovirus genotyping** | GI-GXIX | **reference tables drafted, not yet built** |

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

The reference data ships inside the package, so this gives a working install
with no database build step.

### OpenRDP (only for `calicityper recomb`)

OpenRDP pins `numpy<2` while `calicityper` needs `numpy>=2`, so **it cannot
share the environment**. Install it separately and point at the CLI, the same
way BLAST and MAFFT are found:

```bash
conda create -n openrdp -c conda-forge python=3.11 "numpy<2" "scipy>=1.5" "h5py>=3.8,<3.11"
conda activate openrdp
git clone https://github.com/PoonLab/OpenRDP.git && pip install ./OpenRDP
conda deactivate

export CALICITYPER_OPENRDP=$(conda info --base)/envs/openrdp/bin/openrdp
```

OpenRDP is unpackaged (no PyPI or bioconda release), so pin the git SHA for a
reproducible deployment. Without it, `recomb` still runs and reports which
methods were skipped and why.

## Test

```bash
pytest tests/ -v                      # unit, CLI and regression tests
python tests/loo_validation.py        # leave-one-out, no internet needed
```

---

## Invoking calicityper

### Typing: the common case

```bash
% calicityper type contigs.fasta -o typed.tsv --threads 8
```

### Typing with a visual report

```bash
# also writes run01.md, run01.json and run01.html
% calicityper type contigs.fasta -o typed.tsv --threads 8 --report run01
```

The HTML carries Bokeh alignment panels (query against the references of its
assigned type *and* of the runner-up) and a neighbour-joining tree per region.
A bootstrap percentage cannot be inspected; the alignment behind it can.

### Recombination, from an alignment

```bash
% calicityper recomb --alignment panel.fasta --query SAMPLE_01 -o recomb_out/
```

### Recombination, from unaligned sequences

```bash
# aligned with MAFFT first, then scanned
% calicityper recomb --input samples.fasta -o recomb_out/ --threads 8
```

### Check what would run, without running it

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

Assigns genogroup and, for norovirus, ORF1 and ORF2 types.

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

Example:

| query | genogroup | orf1_type | orf1_support | orf1_runners | orf2_type | note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| NC_001959.2 | GI | GI.P1 | 100.0 | `GI.P1:100.0\|GI.P10:0.0` | GI.1 | |
| NC_039476.1 | GII | GII.P16 | 72.0 | `GII.P16:72.0\|GII.P1:0.0` | GII.2 | |
| frag_07 | GII | | | | GII.3 | orf1: <100 nt overlap |
| junk_01 | unassigned | | | | | no BLAST hit |

The second row is what the runner-up column is for: `GII.P16` cleared the 70%
threshold by two points and nothing else drew support.

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
| `--report-max-panels` | int | most alignment panels to draw [20] |
| `--no-tree` | flag | skip the per-region trees in the HTML |

### `calicityper recomb`

Detects and scores recombination. **Norovirus only** -- the reference panels
and the ORF1/ORF2 junction model are norovirus-specific.

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

## Input must be assembled

`calicityper` types **contigs and genomes, not sequencing reads.** There is no
FASTQ mode. For an Illumina run: map or BLAST reads against norovirus
references, assemble or cluster the hits (for example `vsearch` at 99%
identity), and pass the representatives in.

Recombination analysis additionally needs **near-complete genomes**: the
default `min_sequence_length` is 5,000 nt and `min_informative_sites` 500.
Most wastewater contigs will correctly return `insufficient_sequence`.

---

## Genogroup and genotype coverage

Two different steps, with very different reliability.

**Genogroup detection works across Caliciviridae.** It is a BLAST search
against 60 whole genomes -- 30 norovirus GI-GX and 30 other genera (Sapovirus
GI-GXII, Lagovirus, Vesivirus, Nebovirus, Nacovirus, Recovirus, Valovirus) --
so non-noroviruses are *rejected* rather than forced into a norovirus
genogroup.

**Genotype assignment is validated for norovirus GI and GII only.** Typing
needs a curated reference alignment and only GI and GII have one. Everything
else falls back to `all.aln`, a mixed alignment in which they are represented
by a handful of sequences:

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

A query outside GI/GII **will** return a genotype and a confident-looking
percentage. Treat it as unverified: two or three references cannot support a
reliable assignment. This is adequate for human and wastewater surveillance
(GI, GII, GIX) and **not** for the animal genogroups: GIII (bovine), GV
(murine), GVI/GVII (canine), GX (bat).

### GII.4 variants

GII.4 resolves to variant level from 11 variants (20 sequences): Bristol,
US95_96, Camberwell, Farmington_Hills, Hunter, Asia, Yerseke, Den_Haag, Osaka,
New_Orleans, Sydney. The variant is carried **inside** `orf2_type`
(`GII.4_Sydney`); strip the suffix for genotype-level labels, as RIVM reports
the variant in a separate column.

### Sapovirus typing: in preparation, not yet available

Reference tables have been transcribed from the BMC Genomics dual VP1/RdRp
typing proposal -- 54 VP1 genotypes and 46 RdRp P-types over 198 accessions, in
[`validation/sapovirus_vp1_taxonomy.csv`](validation/sapovirus_vp1_taxonomy.csv)
and
[`validation/sapovirus_rdrp_taxonomy.csv`](validation/sapovirus_rdrp_taxonomy.csv).
Every accession was verified against NCBI, which caught four transcription
errors.

**The reference alignments have not been built from them yet**, so
`calicityper` does not currently assign sapovirus genotypes. Sapovirus
sequences are detected to genus and genogroup only. The tables are marked
DRAFT and need a human read-through before the sets are built.

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
| ORF1 P-type | 147 | 93.9% | 93.2% | 1 |

Two caveats that matter:

- **Leave-one-out flatters the tool.** Queries are the same curated,
  full-length references the tool types against.
- **The one error is a confident error.** `AB074893.1` (true `GII.PNA4`) was
  called `GII.P11` at **99% support**. Bootstrap support is not a guarantee of
  correctness.

The comparison against RIVM's own calls for 4,392 accessions is **running**;
the genogroup comparison over all 5,986 is **queued**. Neither number is
available yet.

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

**Can it type sapovirus?** Not yet. It detects sapovirus to genogroup; the
genotype reference sets are drafted but not built.

**Can it take FASTQ?** No. Assemble or cluster first.

**Why did recombination return `insufficient_sequence`?** The sequence is under
5,000 nt or has fewer than 500 informative sites. The classical scanners cannot
say anything useful below that, and guessing would be worse.

**Why is a GIII call untrustworthy when it shows 98% support?** Support
measures how often the bootstrap recovers the same grouping, not whether the
reference set contains the right answer. GIII has two or three references.

**Do I need OpenRDP?** Only for `recomb`. Without it the command still runs and
reports every method it skipped.

**Why is it slow?** The bootstrap is ~95% of runtime, and cost scales with the
square of the genogroup's reference count, so GII is much slower than GI.
Install `numba`; without it the bootstrap is roughly 8x slower.

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

## Citation

```bibtex
@software{fitzpatrick_calicityper,
  author = {Amy Heather Fitzpatrick},
  title  = {calicityper: command-line calicivirus genotyping and norovirus recombination analysis},
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

## Attribution

A fork and substantial rewrite of
[zhuzhanji/NorovirusGenotyping](https://github.com/zhuzhanji/NorovirusGenotyping)
by Yang Ji. See **[ATTRIBUTION.md](ATTRIBUTION.md)** for what came from where
and an important note on upstream licensing.

## Author

Amy Fitzpatrick, University College Dublin.
