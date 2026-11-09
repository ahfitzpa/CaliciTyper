# Benchmarking and validation

How `calicityper` has been tested, what each test does and does not show, and how
to reproduce every number quoted in the [README](README.md).

Read the caveats. Several of these tests are easy to misread as stronger
evidence than they are.

---

## Summary

| Test | What it measures | Scope | Status |
|---|---|---|---|
| [Leave-one-out](#1-leave-one-out-on-the-bundled-references) | internal consistency of the typing method | 272 queries, GI + GII | complete |
| [RIVM genotype comparison](#2-rivm-genotype-comparison) | agreement with the reference implementation | 4,392 GenBank accessions | running |
| [RIVM genogroup comparison](#3-rivm-genogroup-comparison) | genogroup detection, including rejecting non-noroviruses | 5,986 accessions | queued |
| [Synthetic recombination](#4-synthetic-recombination-datasets) | recombination sensitivity and specificity | 6 generated datasets | complete |
| [Performance](#5-performance) | runtime and where it goes | profiling + end-to-end | partial |

---

## 1. Leave-one-out on the bundled references

Each reference ORF is typed against all other references of its genogroup; its
own label is the answer. Only types with at least two references can be
evaluated, since a type with one reference cannot be recovered once it is
removed.

```bash
python tests/loo_validation.py --threads 8 --out results/loo.csv
```

No internet required.

| Region | Queries | Genogroups | Distinct types | Assigned | Correct | Wrong |
|---|---|---|---|---|---|---|
| ORF2 capsid | 125 | GI (34), GII (91) | 42 | 100% | 100% | 0 |
| ORF1 P-type | 147 | GI (36), GII (111) | 49 | 93.9% | 93.2% | 1 |

GII.4 is resolved to variant level; all 17 evaluable GII.4 references are
called correctly at 81-100% support.

### What this does not show

**It is not independent.** The queries are the same curated, full-length
sequences the tool types against. It demonstrates the method is wired up
correctly. It does not estimate accuracy on real data, and it cannot, by
construction.

**The one error is a confident error.** `AB074893.1` (true `GII.PNA4`) was
called `GII.P11` at **99% support**. The 70% threshold did not catch it.
Bootstrap support is not a guarantee of correctness, and no threshold makes it
one. Check the runner-up columns on any call that matters.

---

## 2. RIVM genotype comparison

The independent test: compare against the published calls of the
[RIVM Norovirus Typing Tool](https://www.rivm.nl/mpf/typingtool/norovirus/) for
4,392 GenBank accessions, 1,868 of them GII.4.

```bash
python scripts/fetch_rivm_test_set.py rivm_test_set.fasta   # ~4,392 seqs from NCBI
calicityper type rivm_test_set.fasta -o rivm_typed.tsv --threads 16
python scripts/validate_vs_rivm.py rivm_typed.tsv
```

`validation/results_rivm.csv` is **not redistributed** here; see
[validation/README.md](validation/README.md).

### Reconciling the labels first

The two tools disagree about where a GII.4 variant lives, and comparing the
fields naively scores every GII.4 record as a mismatch -- 42% of the test set.

| | genotype field | variant field |
|---|---|---|
| RIVM | `capsid type` = `GII.4` | `capsid subtype` = `Sydney_2012` |
| calicityper | `orf2_type` = `GII.4_Sydney` | (inside the type) |

The comparison therefore scores genotype agreement on the genotype stem
(`GII.4_Sydney` -> `GII.4`) and reports the GII.4 variant separately against
RIVM's subtype column, with the emergence year stripped.

**Result: not yet available.** The run has not completed.

---

## 3. RIVM genogroup comparison

A wider test than genotyping, because genogroup detection is a BLAST search
that does not need a curated alignment, and because it must also *reject*
non-noroviruses.

```bash
python scripts/fetch_rivm_test_set.py rivm_all.fasta --all   # all 5,986
calicityper type rivm_all.fasta -o all_typed.tsv --threads 16
python scripts/validate_vs_rivm.py all_typed.tsv --out results/rivm_comparison_all.csv
```

Composition of the 5,986 records:

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
| &nbsp;&nbsp;Lagovirus | 765 |
| &nbsp;&nbsp;Sapovirus | ~420 |
| &nbsp;&nbsp;Vesivirus | ~190 |
| &nbsp;&nbsp;Other / unassigned | ~320 |

This reaches the animal genogroups (GV, GIV, GIII) that genotyping cannot
validate, and the 1,691 non-noroviruses exercise rejection.

**Result: not yet available.** The run is queued.

---

## 4. Synthetic recombination datasets

### Why synthetic

Real norovirus reference genomes **cannot** serve as negative controls.
Norovirus recombines constantly around the ORF1-ORF2 junction, so a published
reference strain frequently *is* a recombinant, and a `possible` call on one may
be correct rather than a false positive. Measuring specificity requires
sequences that are non-recombinant by construction.

`tests/synthetic.py` generates six seeded datasets, each a reference panel of
four lineages (two members each) plus one query:

| Dataset | Construction | Expected |
|---|---|---|
| `non_recombinant` | query descends from one lineage | `none` or `weak` |
| `clear_recombinant` | lineage0 left of 5100, lineage2 right | `possible` or better, breakpoint near 5100 |
| `weak_signal` | 150 nt imported tract | not above `possible` |
| `conflicting` | two tracts from different donors | multiple events |
| `orf_mismatch_only` | distinct ORF1/ORF2 ancestry, near-identical donors | low parental distinctness, never `strong` |
| `too_short` | 900 nt | `insufficient_sequence` |

### Results

| Dataset | Status | Score | Top event |
|---|---|---|---|
| `non_recombinant` | `weak` | 0.38 | 1 event, single method, no phylogenetic support |
| `clear_recombinant` | `possible` | 0.84 | **query 5061-5351** (planted: 5100) |
| `too_short` | `insufficient_sequence` | - | correctly refused |

The true breakpoint is recovered to within 39 nt, and the non-recombinant does
not reach `possible`.

### Component validation

Each scoring component was checked against ground truth separately.

**Coordinate mapping.** Alignment column 3551 maps to query position 3201 for a
breakpoint planted at 3200.

**Segment phylogenies.** The first implementation compared everything left of a
breakpoint against everything right of it. That confirms recombination but
cannot localise it: on a real recombinant almost any split shows donor A against
donor B, so a decoy breakpoint at 1000 scored 1.00, identical to the true one.
Flanking windows (`segment_window`, default 1000) fixed it:

| Breakpoint tested | Score | |
|---|---|---|
| 3551 | 1.00 | true |
| 3400 | 1.00 | 151 nt from true, inside `max_breakpoint_distance` |
| 1000 | 0.00 | decoy |
| 2000 | 0.00 | decoy |
| 6000 | 0.00 | decoy |

**Breakpoint semantics.** Window scanners (MaxChi, Chimaera, BootScan, SiScan)
report the window *containing* a breakpoint, so the estimate is its centre.
GENECONV and 3Seq report a converted *fragment* spanning thousands of columns,
where the breakpoints are the **edges**. Clustering on midpoints alone hid the
true event entirely; span width now decides which rule applies.

**Status rules.** One method plus a phylogeny hint previously reached
`possible` and produced 15 spurious events on a control. A single method now
qualifies only with phylogenetic support >= 0.7 *and* distinct parents.

### What this does not show

The datasets are generated from random sequence with a simple substitution
model. They test that the machinery localises a known breakpoint and does not
invent one. They do **not** reproduce real norovirus sequence evolution, and
they say nothing about performance on partial or chimeric wastewater contigs.

---

## 5. Performance

### Where the time goes

Measured per query on GII (the common case), 100 bootstrap replicates:

| Step | Per query | Share |
|---|---|---|
| BLAST | ~42 ms | ~2% |
| MAFFT | ~40 ms | ~4% |
| **Bootstrap** | **~900 ms** | **~95%** |

Within the bootstrap, the pairwise distance matrix is **98%** of the time and
neighbour-joining 2%. Cost scales with the *square* of the genogroup's
reference count, so genogroup drives runtime:

| Genogroup | ORF1 refs | ORF2 refs |
|---|---|---|
| GI | 41 | 34 |
| GII | 121 | 96 |

**Replacing MAFFT is not worth doing** -- a perfect aligner saves under 5%.

### Unresolved

An earlier figure of 0.5 s per genome was measured on GI and does not
generalise; real data is overwhelmingly GII. A full run of 4,391 mostly-GII
sequences on 16 cores took far longer than the per-component measurements
predict, and **that discrepancy is not yet explained**. A controlled A/B test of
BLAS thread pinning is pending. Until it resolves, treat the per-step table as
reliable and any end-to-end projection as optimistic.

Install `numba`; without it the bootstrap is roughly 8x slower.

---

## Reproducing

```bash
conda env create -f environment.yml && conda activate calicityper
pytest tests/ -v                                   # unit, CLI and regression tests
python tests/loo_validation.py --threads 8         # section 1
python scripts/fetch_rivm_test_set.py rivm.fasta   # sections 2-3 (needs NCBI)
```

The synthetic datasets are regenerated by the test suite and are seeded, so a
failure is reproducible.
