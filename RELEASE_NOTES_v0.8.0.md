# v0.8.0 — "Dual Typing"

Sapovirus VP1 genotyping and RdRp P-typing, reported together as a dual type,
and a fix that anyone running 0.7.0 needs.

## Upgrade now if you installed 0.7.0

`calicityper setup` silently discarded the sapovirus reference sets. They were
filtered out during extraction, so the archive and the manifest could both be
correct and the data still never arrived. Sapovirus typing could not have
worked from a released 0.7.0 install. A test now asserts that every top-level
entry the manifest checksums is actually installable, so this cannot recur.

## Sapovirus dual typing

On by default; disable with `CALICITYPER_TYPE_SAPOVIRUS=0`.

VP1 genotype and RdRp P-type are assigned by nucleotide p-distance to the
nearest reference, using the cut-offs of Zhao *et al.* 2024, and reported
together as `GI.4[P1]` — or `GIV.1[GII.P10]` where the polymerase group differs
from the capsid genogroup. A query outside the genotype cut-off of every
reference is `untypeable` rather than given its nearest neighbour's label.
Coverage of the typed region is reported but not gated on.

**What is validated, and what is not.** VP1 genotyping is supported by
leave-one-out recovery (106/127, no discordant call), by agreement with the
Norovirus Typing Tool on sequences that are not reference strains (198/207),
and by reproduction of the published within-genotype ranges (31/32). RdRp
P-typing is supported for assignment among the P-types the reference set
represents (68/72, no discordant call) — but **not** for recognising a P-type
that is absent: none of five singleton P-types was reported untypeable once its
only reference was removed, and only 1 of 20 testable P-types is separated from
every other in its genogroup by the 0.266 cut-off.

Read a dual type as two calls of unequal standing. A polymerase from an
unrepresented lineage will be confidently wrong rather than declined. Coverage
is the binding constraint in practice: 51 of 144 VP1 references are also in the
RdRp set, and on a 515-sequence panel only 224 (43.5%) received both halves.
See `docs/SAPOVIRUS_STATUS.md`.

## Corrections to the published reference tables

Rebuilding the Zhao *et al.* 2024 classification programmatically surfaced four
strain-level discrepancies between that paper's own tables. Each is recorded
with its evidence in the generator, so regenerating cannot silently reverse it,
and `scripts/check_against_zhao_tables.py` cross-checks the result against both
published tables.

The substantive one: **Supplemental Table 1 has no GI.P6 row at all**, jumping
from GI.P5 to GI.P7 and listing GI.P6's three strains under GI.P5. Those five
strains are two clusters 0.249–0.255 apart; the published GI.P5 range of
0–0.052 fits one of them exactly and is impossible for all five. Restoring
GI.P6 returns GI.P5's within-type maximum from 0.255 to 0.050 and adds three
complete RdRp references. The others: AJ606693 is GI.P1, not GI.P2; GNA2.P and
GNA3.P are crossed between OP963623 and OQ818319; MH003837 is GXV/GXV.P.

## Recombination

- A failure is no longer reported as a negative. When no scanner runs the
  status is `not_assessed`, and a query that raised is `failed`; both are
  distinct from `none`. A tool that times out previously looked like a clean
  sequence.
- A failed combined OpenRDP run is retried method by method, so one crashing
  scanner no longer discards the evidence from the other six.
- `--openrdp-timeout` makes the scanner budget configurable.

## Reproducibility

- `analysis/` holds the 32 SLURM batch scripts as submitted, the scored output
  tables behind every reported number, the figure generator, and a README
  mapping each claim to the script and file that produced it.
- The sapovirus reference sets are in the manifest and checksummed; the set is
  229 artefacts, up from 84.
- `scripts/build_reference_archive.py` builds the reference archive from the
  manifest and re-verifies it afterwards, replacing a hand-made tarball.
- Genogroup BLAST uses a fixed search-space size, so growing the database no
  longer moves calls, and a hit shorter than 100 nt no longer names a genogroup.
- `containers/calicityper-tools.def` pins BLAST+ 2.17.0, MAFFT 7.526 and
  OpenRDP 0.1.0 in one Apptainer image. See `docs/INSTALL_TOOLS.md`.

## Install

```bash
pip install calicityper-0.8.0-py3-none-any.whl
calicityper setup --source calicityper-reference-0.8.1.tar.gz
```

The wheel deliberately contains no reference data. `setup` verifies every
artefact against the manifest's SHA-256 list and refuses a set that does not
match, because a different reference set changes calls silently.

`calicityper-reference-0.8.1.tar.gz`
sha256 `47a79b5e24622cdd311b701d33553145f1528a54ed37178f2804520926e517f6`

**Use 0.8.1, not 0.8.0.** The 0.8.0 archive shipped three stale files,
`reference_sequences_sapovirus_VP1/GNA4.{orf,fasta,aln}`, left behind when
MH003837 was relabelled GNA4 to GXV: the builder wrote a file per genogroup
but never removed one whose genogroup had lost its last member. They put
MH003837 under a genotype this release retracts, while the pooled `all.orf`
did not, so the set disagreed with itself. Nothing could reach them, because
the genogroup BLAST database has no GNA4 entry and the reference set is
selected by genogroup name, so no call made with 0.8.0 is affected. The
builder now prunes such files and three tests assert the set is internally
consistent.

Requires BLAST+ and MAFFT on `PATH`. OpenRDP, needed only for `recomb`, cannot
share an environment with CaliciTyper: it requires numpy<2 while CaliciTyper
runs numpy 2.5. See `docs/INSTALL_TOOLS.md`.

## Known limitations

- RdRp cannot flag a novel sapovirus P-type (above).
- 14 of 45 P-types have no reference carrying a complete RdRp, so queries of
  those types cannot be assigned; the tool says so rather than guessing.
- Norovirus genotyping is validated for GI and GII only.
- Without OpenRDP installed, 13 tests covering the recombination scanners are
  silently deselected — a green test run does not mean that path was exercised.
