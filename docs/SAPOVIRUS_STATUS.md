# Sapovirus: where this stands

Updated 2026-11-12, against calicityper 0.6.1.

Short version: **recombination works and is tested. Genotyping does not exist
yet** — the reference sets are extracted but not built, and nothing in the
typing path is sapovirus-aware. One decision is needed before the RdRp set can
be built; see "The RdRp coverage problem".

## What the paper specifies

Zhao *et al.* 2024, BMC Genomics 25:1131, dual VP1/RdRp typing.

- **Nucleotide distances, not amino acid.** Explicitly: "the pairwise distance
  of the amino acid sequences for sapovirus could not be used to define the
  genotype range statistically". Trees were built on amino acids; the cutoffs
  are nucleotide. This is the opposite of norovirus, where the capsid genotype
  is defined on complete VP1 **amino acids**.
- **Cutoffs** (nucleotide p-distance): VP1 < 0.161 same genotype, < 0.503 same
  genogroup; RdRp < 0.266 same P-type, < 0.531 same P-group.
- **Regions**: complete VP1, and the "complete RdRp region".
- **Genome shape**: ORF1 encodes the nonstructural proteins *and* VP1 as one
  polyprotein, so both regions sit inside ORF1. The norovirus "take ORF1, take
  ORF2" approach does not transfer.

## Done

### Taxonomy tables, verified against the published spreadsheet

`scripts/build_sapovirus_taxonomy.py` generates both CSVs directly from
Supplemental Table 1 (`validation/12864_2024_11048_MOESM1_ESM.xlsx`), replacing
a transcription made from an image of that table. The transcription had six
accessions in the wrong row or with a changed character:

| | Accession | Was | Corrected to |
|---|---|---|---|
| VP1 | MN509083 | GI.7 | GI.6 |
| VP1 | LC215899 | GX | GXI |
| VP1 | MF766258 | GXXIII | GXXIV |
| VP1 | GNA2 | OP963622 (BtSY1) | OP963623 (BtSY2) |
| RdRp | GII.P10 | MG012401 (a GI.3 record) | MG012421 (the GIV.1 strain) |
| RdRp | MW450824 | GXIV.P | GXVI.P |

Counts now match the paper: 53 VP1 genotypes / 156 strains, 45 P-types / 139
strains. Two tests fail if the CSVs drift from the spreadsheet.

Four accessions are corrected against NCBI because the **published** value is
not a sapovirus record — `L504335` resolves to nothing, `KC215898` is a
scorpion 18S, `KC469052` an uncultured bacterium, `LC215696` a whitefly COI.
Each correction is recorded in the CSV `note` column.

One is unresolved: **`KC309412`**, listed for VP1 genogroup GXXVI, is a 366 nt
*Pseudomonas aeruginosa* 16S record. It is flagged and excluded. **Worth
raising with the authors**, along with the four typos.

### Region extraction

`scripts/extract_sapovirus_regions.py` takes the 196 GenBank flatfiles in
`~/scratch/calicityper_runs/sapovirus_build/sapovirus.gb` and extracts VP1 and
RdRp from each.

Only 18 of the 196 records annotate any boundary (`mat_peptide`), and most of
those annotate only VP1, so the method is: define the regions **once** on a
canonical record, project that definition onto the other annotated records so
all anchors describe the same region, then map onto each query through a
pairwise protein alignment, using whichever anchor aligns best.

Canonical record is **MH541033**, which annotates NS1–NS6-7 and VP1:

- **VP1** = polyprotein aa 1717–2277 (560 aa)
- **RdRp** = aa 1205–1717 (512 aa, 1,536 nt)

The RdRp 5′ boundary is **derived, not annotated**. No sapovirus record in
GenBank carries an `NS7` feature (checked: zero), and the sapovirus RefSeq
NC_006269 has no mature-peptide annotation at all. MH541033 annotates the
combined NS6-7 (665 aa); the NS6/NS7 cleavage site at offset 153 — an E|N
dipeptide — gives NS7 = 512 aa with the catalytic YGDD motif 344 aa in, the
canonical position for a calicivirus polymerase, and NS6 = 153 aa, the expected
protease size. The two alternative cleavage sites put YGDD at 358 and 306 aa,
both wrong.

The **3′ boundary is independently confirmed**: the annotated RdRp of U95643,
U95644 and U95645 all map onto MH541033 ending at exactly aa 1717, 0 aa
discrepancy.

Accuracy is measurable because 18 records annotate their own boundaries:
**exact agreement (0 aa) in 23 of 25 checks**. The exceptions are OQ818319, a
divergent bat sapovirus at 0.445 identity, off by 13 aa. Mean anchor identity
across all records is 0.74.

Results: **192 VP1** sequences and **110 RdRp**, from 193 of 196 records.
Not extracted: KC309412 and MN510437 have no CDS feature; OM480531 (the bamboo
rat GX-F1 strain, which the paper itself notes lacks definitive typing) fell
below the 0.35 identity floor.

### Recombination

Working and tested. `--virus sapovirus` selects a genome model whose hotspot is
the NS7(RdRp)/VP1 junction *inside* ORF1, not a junction between two ORFs.
Three synthetic datasets at the sapovirus hotspot (nt 5139) and away from it
(nt 2500) are in `tests/synthetic.py`.

All four sapovirus recombination tests pass end to end (job 1139087,
2026-11-12), roughly 20 minutes each:

| Dataset | Planted | Recovered | Methods |
|---|---|---|---|
| `sapovirus_clear_recombinant` | 5139 | 4949–5415 | bootscan, geneconv, rdp, threeseq |
| `sapovirus_mid_ns_recombinant` | 2500 | 2301–2541 | bootscan, geneconv, threeseq |
| `sapovirus_non_recombinant` | none | none | — |

A fourth test confirms the event is labelled an **RdRp-VP1 junction** rather
than an ORF1-ORF2 junction — the point of having a per-virus model, since both
viruses score nt 5139 as a hotspot but only one names it correctly.

## Not done

Building the per-group reference sets, aligning them, wiring
`plan_query` routing, and any genotyping benchmark. **No sapovirus genotype
call exists.** Sapovirus input is detected to genogroup and stops there.

### The RdRp coverage problem — decision needed

Joining the extracted sequences to the verified tables:

| | Extracted | In the table | Usable | Labels covered | Labels with ≥2 members |
|---|---|---|---|---|---|
| VP1 | 192 | 147 | **145** | 51 of 52 | 32 |
| RdRp | 110 | 138 | **81** | **29 of 45** | 20 |

VP1 is in good shape. **RdRp covers only 29 of 45 P-types**, because requiring
a *complete* RdRp (the paper's wording) excludes 75 records that are 5′
truncated — their mapped RdRp is 188–259 aa rather than ~512. The kept 110 are
tight and consistent: 496–517 aa, median 512.

Three options, in order of preference:

1. **Build the P-type set from complete RdRp only** (29 P-types). Consistent
   with the paper and internally comparable, but a query whose P-type is one of
   the missing 16 can never be typed correctly — it will be assigned the
   nearest available P-type, or flagged as putative novel. Honest if the
   coverage is documented.
2. **Admit partial RdRp references and compare with pairwise deletion**, which
   is what `msa_distances()` already does. Recovers coverage, but a distance
   computed over 200 aa is not comparable with the paper's cutoffs, which were
   derived on complete sequences.
3. **Fetch more sequences.** The 196 flatfiles came from one query; the P-types
   missing a complete RdRp may have complete records elsewhere in GenBank. This
   is the only option that fixes the coverage without compromising either
   consistency or the cutoffs, and should be tried first.

My recommendation is 3, then 1 for whatever is still missing, with the covered
P-types listed in the README so a user knows what the set can and cannot type.

## Next steps

1. Try to fill the 16 missing P-types from GenBank (option 3 above).
2. Build per-genogroup FASTAs for VP1 and RdRp, in the norovirus layout
   (`{group}.fasta`, `{group}.orf`, `all.orf`), align each with MAFFT into
   `{group}.aln`.
3. Wire the typing path: a sapovirus group list, per-virus `REF_DIR`, and
   genogroup→virus routing in `plan_query`.
4. Benchmark: leave-one-out over the reference strains, plus a check that
   computed within-cluster distances fall inside the paper's published ranges —
   an independent validation norovirus cannot offer, because the paper
   publishes a distance range per genotype.
