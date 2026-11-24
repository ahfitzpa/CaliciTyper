# Sapovirus: where this stands

Updated 2026-11-17, against calicityper 0.8.0.

Short version: **genotyping is implemented, benchmarked and on by default.**
VP1 is reliable. RdRp assigns accurately among the P-types it holds, but
cannot establish that a query belongs to none of them, so it cannot flag a
novel P-type. Recombination works and is tested.

Disable sapovirus typing with `CALICITYPER_TYPE_SAPOVIRUS=0`.

## What the paper specifies

Zhao *et al.* 2024, BMC Genomics 25:1131, dual VP1/RdRp typing. Both typed
regions lie inside sapovirus ORF1, which encodes the non-structural proteins
and VP1 as one polyprotein, so VP1 appears in the `orf2_*` columns and RdRp in
the `orf1_*` columns.

Nucleotide p-distance to the nearest reference decides the call:

| Region | Genotype / P-type | Genogroup / P-group |
|---|---|---|
| VP1 | < 0.161 | < 0.503 |
| RdRp | < 0.266 | < 0.531 |

A query outside the genotype cut-off of every reference is `untypeable`, not
given its nearest neighbour's label. Bootstrap support, a 95% Wilson score
confidence interval on the observed p-distance (`orf{1,2}_distance_ci`), and
the distance/margin to the nearest distinct alternative type
(`orf{1,2}_second_type`, `orf{1,2}_second_distance`, `orf{1,2}_distance_margin`)
are reported beside each call to quantify fragment-length and inter-type
uncertainty. Warning notes are appended when the 95% CI straddles the cut-off
or when a second distinct type also falls below the cut-off.

Results carry a dual type: `GI.4[P1]`, and `GIV.1[GII.P10]` where the P-group
differs from the capsid genogroup. Columns `dual_type`, `dual_type_discordant`.

## Reference sets

Generated from Supplemental Table 1 by `scripts/build_sapovirus_taxonomy.py`,
then built and aligned by `scripts/build_sapovirus_refsets.py`.

| | Records | Types | Genogroups |
|---|---|---|---|
| VP1 | 144 | 50 | 33 |
| RdRp | 82 | 31 | 14 |

The taxonomy tables hold more: 156 VP1 accessions over 52 genotypes, and 139
RdRp accessions over 45 P-types. The reference sets are smaller because only
records carrying the **complete** typed region are kept. The published
cut-offs were derived from complete sequences, so a reference truncated to 290
of 512 aa would make distances incomparable. Excluded records are listed in
`calicityper/data/sapovirus_excluded.tsv`.

Region boundaries were defined once on MH541033, the only GenBank record
annotated with both NS1-NS6-7 and VP1, and projected onto every reference by
pairwise protein alignment. VP1 is polyprotein aa 1717-2277, RdRp aa
1205-1717. The 5' RdRp boundary is inferred, because no sapovirus GenBank
record carries an NS7 feature.

### Corrections to the published tables

Four, each recorded as a `RELABEL` entry in the generator with its evidence,
so regenerating the tables cannot silently undo them.
`scripts/check_against_zhao_tables.py` cross-checks the result against Table 2
and the pig/bat table and reports no mismatches.

| Strain | Supplement says | We use | Why |
|---|---|---|---|
| AJ606693 | GI.P2 | **GI.P1** | Table 2 designates it GI.4[P1]; distances agree (0.144 to GI.P1, 0.272-0.281 to GI.P2) |
| OP963623 / OQ818319 | GNA2.P / GNA3.P, crossed | **uncrossed** | The pig/bat table gives the uncrossed pairing; these are the only two discordances in the supplement naming a P-group coined in that same paper |
| MH003837 | GNA4 / GNA4.P | **GXV / GXV.P** | The pig/bat table, confirmed manually against the published tables |
| LC380411, MG012443, MN509083 | GI.P5 | **GI.P6** | The supplement has no GI.P6 row at all, jumping GI.P5 to GI.P7 |

The GI.P6 correction is the substantive one. The five strains the supplement
lists under GI.P5 are two clusters -- DQ366345 + LC504396 at 0.050, and the
three above at 0.019-0.057 -- lying 0.249-0.255 apart. The published GI.P5
range of 0-0.052 fits the first cluster exactly and is impossible for all
five, whose spread is 0.255. The VP1 genotypes split on the same line (GI.5
against GI.6), and Table 2 names DQ366345 GI.5[P5] and MG012443 GI.6[P6], one
in each cluster. Restoring GI.P6 returns GI.P5 to its published range and adds
three complete RdRp references.

With GNA4.P emptied and GI.P6 restored, the P-type count is the paper's own
45. VP1 is 52, one short of the paper's 53, because GNA4 now has no members.
That is an open question for the authors.

## How well it works

Full detail in the benchmarking repository's `benchmarking/README.md`.

**Against the published distance ranges.** VP1 reproduces the published range
for 31 of 32 genotypes, and all 22 testable genotypes are separated from every
other genotype in their genogroup by at least the 0.161 cut-off. RdRp
reproduces 19 of 20; the one exception is GI.P1 at 0.144 against a published
0-0.085, which is AJ606693 sitting atypically within its own type.

**Leave-one-out, complete length.** RdRp recovers 68 of 72 references (94.4%),
VP1 106 of 125 (84.8%). Neither region produced a single discordant call:
every failure is a decline, not a mistake.

**Under truncation.** VP1 is 82.9-87.2% correct across 150-1,200 nt against
84.8% at complete length, so essentially flat. RdRp is 91.2-93.1%. No VP1
fragment at any length produced a discordant genotype, so VP1 precision is
1.000 at every coverage floor from 0.0 to 0.95. That is why no minimum
coverage is applied: a floor only withholds correct calls.

**Against the Norovirus Typing Tool.** Of 515 sapovirus sequences it assigned
a VP1 genotype to 294; we typed 286 and agreed on 265. Restricted to the 213
that are not our own reference strains, we typed 207 and agreed on 198. All
nine differences are reassignments under the newer classification.

## The RdRp limitation

**Only 1 of 20 testable P-types is separated from every other P-type in its
genogroup by the 0.266 cut-off**, and the cut-off overlaps the body of the
between-type distance distribution.

This is not evidence that 0.266 is wrong. It was derived as the mean minus two
standard deviations of a pooled between-type peak, so a tail below it is
expected by construction. It is a practical limit on novel-type detection:
where two known P-types are closer together than the cut-off, a query from an
unrepresented P-type can fall inside the cut-off of one or both, and the
cut-off cannot then establish that it belongs to neither.

Measured directly, that is exactly what happens. Of five singleton RdRp
P-types, **none** was reported untypeable once its only reference was removed;
all five were given another known P-type's label. VP1 passes the same test
6 of 6.

Two alternatives were tested and rejected. Tamura-Nei correction raises
separation from 1 to 4 of 19 but worsens agreement with the published ranges.
Assigning to the lowest-mean-distance cluster lowers accuracy for both regions
(RdRp 92%, VP1 78%) and introduces two discordant calls that nearest-neighbour
assignment does not make.

**So:** assignment among known P-types is sound; novel P-type detection from
RdRp alone is not. Treat P-types as provisional for any query distant from a
dense cluster of references. VP1 is the reliable region for novelty.

## Coverage gaps

14 of 45 P-types have no reference carrying a complete RdRp, so a query of
those types cannot be assigned. This is a coverage gap, not an accuracy
problem: the tool reports that no complete reference exists rather than
borrowing a neighbouring type's label. Complete-RdRp accessions for those
types are on the list of questions for the authors.

GIV has no RdRp reference set at all. The paper states that GIV RdRp clusters
within GII.P, so reporting "no RdRp set for this genogroup" is correct
behaviour rather than a gap -- pending confirmation from the authors.

## Recombination

Works and is tested; see `examples/recombination/`. The genome model matters:
the recombination hotspot upstream of the capsid is an ORF1-ORF2 junction in
norovirus but an RdRp-VP1 junction *within* ORF1 in sapovirus, and `--virus
sapovirus` reports it correctly.

Validated on nine synthetic datasets, all returning the expected status, and
against the 13 strains of Zhao Table 3. VP1/RdRp discordance alone is not
evidence of recombination, and the tool does not treat it as such.
