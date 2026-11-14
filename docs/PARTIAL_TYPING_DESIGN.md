# Partial-segment typing: design

Status: proposal, not implemented. Written against calicityper 0.6.1.

## The problem

Typing currently assumes a query that covers a whole typing region. A partial
sequence — a surveillance amplicon, a low-coverage consensus, a fragmentary
contig — either fails to reach the 70% bootstrap threshold or is dropped with
`<100 nt overlap`. The tool then says nothing useful, when in practice the
fragment often carries enough signal for the same call the RIVM tool would make
from an amplicon.

Two reference standards are in play, and the design has to serve both:

| Standard | Capsid | P-type |
|---|---|---|
| Chhabra 2019 (formal) | complete VP1, **amino acids** | **762 nt** from the 3′ end of ORF1/RdRp |
| RIVM/NoroNet (operational) | VP1 region C: **252 nt** (GII), 264 (GI), 267 (GI.9), 255 (GIX), from the 5′ end of VP1 | final **172 nt** of RdRp (173 for GIV) |

Sapovirus is different again: Zhao 2024 defines both genotype and P-type on
**nucleotide** distances (amino acids explicitly could not separate sapovirus
genotypes statistically), on complete VP1 and the complete RdRp.

## What already exists and should not be rebuilt

The architecture is most of the way there, which constrains the design usefully:

1. `blastn` against whole genomes picks the genogroup and the best reference,
   and already reports orientation, so reverse-complement handling is solved.
2. The query is added to a **pre-computed reference alignment** with
   `mafft --addfragments --keeplength`. Reference columns therefore never move,
   and a query's coordinates are reference coordinates for free.
3. Comparison already happens only over reference columns. The ORF1 reference
   alignment spans **762 nt** — it is already the Chhabra P-type window, so the
   "full ORF vs fixed window" difference in the README is smaller in practice
   than it reads.
4. `putative_novel.rdrp_region_columns()` already selects a window **as a set of
   alignment columns** by reference occupancy, and `msa_distances()` already
   does pairwise-deletion distances over a column subset.

So window extraction should be expressed as **named column ranges over the
existing reference alignments**, not as genomic coordinates and not as a new
alignment step.

## Pipeline

### 1. Input handling

| Input | Treatment |
|---|---|
| Contigs / consensus | the normal path |
| Raw reads | reject with a clear message; assemble or cluster first (current behaviour, keep it) |
| Ambiguous bases | counted per window, not globally: a query can be clean in region C and unusable in RdRp |
| Reverse complement | from the BLAST HSP orientation, as now |
| Low complexity | reject if the window's distinct-3-mer count is below a floor; a homopolymer run can align anywhere |
| Mixed infection | see below — flag, never silently type |

**Mixed-infection detection.** Two cheap signals, both already available:
ambiguity codes concentrated at columns that are polymorphic between two
candidate types, and two or more references of *different* types within a small
distance margin of the query. Either triggers `possible_mixture`, which reports
both candidates and withholds a definitive call. This matters more for short
windows, where a mixture looks like a near-tie.

### 2. Reference database

Keep **separate VP1 and RdRp sets** — they must be, because a recombinant has a
different type in each, which is the entire point of dual typing.

Each reference record should carry:

- `accession` with version (GenBank records are revised; the 0.6.1 benchmark
  list already pins versions for exactly this reason)
- `virus` (norovirus / sapovirus), `genogroup`, `type`, and for GII.4 the
  `variant`
- `region` (VP1 / RdRp) and the alignment column span it occupies
- `role`: `prototype` (the designated reference strain from the publication) or
  `representative` (an additional curated member)
- `source`: the publication the assignment comes from, plus whether the
  accession was corrected and why — the `note` column added for the sapovirus
  tables is the pattern
- `complete`: whether the record covers the whole region, so partial references
  are never used to define a window

**Versioning.** A single `reference_manifest.json` holding a set version, the
publication each label derives from, per-alignment checksums, and the window
column ranges. Every result should record the reference set version, so a call
can be reproduced after the set changes. Today a changed reference silently
changes every call — the two corrupted GII references showed what that costs.

**Prototypes vs representatives.** Store both, mark the difference, and let the
window definitions depend only on `prototype` + `complete` records. Use
representatives for typing but never for defining coordinates.

### 3. Alignment strategy

Keep the two-stage approach rather than replacing it:

- **BLAST** for genogroup and reference selection. Add explicit gates: a hit
  needs ≥ 80% identity over ≥ 100 nt to be considered at all, and the
  best-scoring genogroup must beat the next by a bit-score margin (10%
  suggested) or the result is `genogroup_ambiguous`.
- **MAFFT `--addfragments --keeplength`** to place the query in the reference
  alignment. This is the right tool precisely because it is built for adding
  short fragments to a fixed alignment, and it is what makes coordinates
  reference-anchored.
- **minimap2** is the wrong fit here: it is built for long-read mapping against
  a genome, not for placing a fragment into a curated protein-coding alignment
  where column identity is the whole point.

**Tie-breaking** should be explicit and reported, not silent: if two references
of different types are within the distance margin, return both and downgrade
the call. Current behaviour ranks by support and reports runner-ups, which is
the right instinct — make it a defined state rather than something the user has
to read off a column.

### 4. Window extraction

Every window is a set of columns in a reference alignment, selected once from
`prototype` + `complete` references and stored in the manifest:

| Window | Definition |
|---|---|
| Complete VP1 (Chhabra) | all columns of the VP1 alignment; translate per query for the amino-acid distance |
| 762 nt RdRp (Chhabra) | the 762 columns of the ORF1 alignment, i.e. what it already spans |
| RIVM VP1 region C | the first 252/264/267/255 columns of VP1, **per genogroup**, measured from the VP1 start column |
| RIVM RdRp short | the final 172 columns (173 for GIV) of the RdRp window |

Extraction for a query is then: take its aligned row, slice those columns, and
measure coverage as the fraction of them carrying an unambiguous base. No
genomic arithmetic, so an insertion or deletion in the query cannot shift a
window — the failure mode that makes fixed genomic coordinates unsafe.

### 5. Scoring and assignment

Run the window ladder from most to least demanding and stop at the first that
passes, recording which was used:

1. **Definitive** — the Chhabra window is covered. Current method: NJ tree,
   100 bootstrap replicates, type assigned when the query clades with all
   references of exactly one type in ≥ 70% (both values confirmed from RIVM's
   own output: its support values are integers 70–100 and nothing is assigned
   below 70).
2. **Partial (RIVM-style)** — only a short window is covered. **Do not use
   bootstrap support here**: on 252 nt a bootstrap percentage is unstable and
   will read as confident. Use nearest-reference distance with a margin: assign
   the type of the closest reference if it is closer than the nearest reference
   of any *other* type by Δp ≥ 0.02, and if its own distance is within that
   type's 2× SD bound. Otherwise report the closest reference without a call.
3. **Closest reference only** — coverage below the partial minimum, or a tie
   within the margin. Report the closest reference, its distance, and the
   candidates. No type.
4. **Untypeable** — no window reaches its minimum, no BLAST hit, or the window
   is below the ambiguity/complexity gates.

**Ambiguous bases and gaps.** Exclude by pairwise deletion, as
`msa_distances()` already does: a position counts only where both sequences
have an unambiguous base. Report the number compared alongside every distance —
a p-distance over 40 positions is not the same evidence as one over 700, and
the number is what tells a reader which they have. Reject a window whose
ambiguity exceeds 5% of covered positions.

**Thresholds to validate, not assume.** Δp ≥ 0.02, the 10% bit-score margin
and the 5% ambiguity ceiling are starting points. Each should be set by the
leave-one-out harness: truncate every reference to each RIVM window, re-type it,
and choose the margin that maximises correct calls while keeping miscalls at
zero. That harness already exists for whole regions
(`tests/loo_validation.py`); extending it to truncated queries is the main
piece of validation work this design needs.

## Minimum-length rules

Lengths are of the **target window**, not the submitted sequence; a 7 kb genome
that misses region C still fails region C.

| Target region | Standard | Definitive | Partial / RIVM-style | Below the partial minimum | Label |
|---|---|---|---|---|---|
| VP1 capsid, norovirus | Chhabra: complete VP1 aa | ≥ 95% of VP1 columns (~1,520 of 1,600 nt) | — | falls through to region C | `definitive` |
| VP1 region C, GII | RIVM | — | ≥ 90% of 252 nt (227 nt) | `untypeable` for capsid | `partial` |
| VP1 region C, GI | RIVM | — | ≥ 90% of 264 nt (238 nt) | `untypeable` | `partial` |
| VP1 region C, GI.9 / GIX | RIVM | — | ≥ 90% of 267 / 255 nt | `untypeable` | `partial` |
| RdRp P-type, norovirus | Chhabra: 762 nt | ≥ 90% of 762 nt (686 nt) | — | falls through to the short window | `definitive` |
| RdRp short | RIVM: final 172 nt (173 GIV) | — | ≥ 90% (155 nt) | `untypeable` for P-type | `partial` |
| VP1, sapovirus | Zhao 2024: complete VP1 **nt** | ≥ 95% of VP1 columns | ≥ 50% (provisional only) | `untypeable` | `definitive` / `provisional` |
| RdRp, sapovirus | Zhao 2024: complete RdRp **nt** (512 aa / 1,536 nt, derived — see below) | ≥ 90% | ≥ 50% (provisional only) | `untypeable` | `definitive` / `provisional` |

`provisional` exists for sapovirus because there is no published short-window
convention to fall back on: a partial sapovirus call is an extrapolation of the
authors' cutoffs to a region they did not validate, and should be labelled as
such rather than presented as equivalent.

Every result should carry, in all output formats: the window used, columns
covered out of columns in the window, positions compared, the distance to the
closest reference and to the nearest other type, and the label. A user who sees
`partial` must be able to see *why* without rerunning anything.

## Sapovirus RdRp: how the 1,536 nt comes about

The paper says "complete RdRp region" without a coordinate, and no sapovirus
record in GenBank annotates `NS7` (checked: zero hits for the feature key).
Two routes were tried:

- **Sapovirus RefSeq NC_006269** — no mature-peptide annotation at all, only
  three CDS. Dead end.
- **Derive the NS6/NS7 boundary inside an annotated NS6-7.** MH541033 annotates
  NS6-7 as 665 aa (polyprotein aa 1052–1717), with VP1 starting immediately
  after. Three candidate E|[GAST] cleavage sites sit in the expected region;
  the one at offset 153 gives NS7 = **512 aa**, with the catalytic YGDD motif
  344 aa into it — the canonical position for a calicivirus polymerase, and
  NS6 = 153 aa, the expected protease size.

Independent check: the "RNA dependent RNA polymerase" annotations on U95643,
U95644 and U95645 all map onto MH541033 ending at **exactly** aa 1717, 0 aa
discrepancy — so the 3′ end of RdRp is confirmed by a second source. Those
three are 5′-truncated (293 aa, starting at their records' own first codon),
which is why they cannot define the start.

The 3′ end is therefore solid; the 5′ end rests on cleavage-site consensus plus
motif spacing, not on annotation. It should be labelled as derived, and is
worth confirming with the paper's authors before the sapovirus reference sets
are published.
