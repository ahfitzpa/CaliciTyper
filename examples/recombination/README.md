# Worked recombination examples

Two synthetic datasets with a breakpoint planted at a known position, and the
output CaliciTyper produces for each. They exist to show what a positive call
looks like and, in particular, that the genome model changes how the same
breakpoint is described.

Both queries are named `QUERY` in their alignment.

## Running them

```bash
calicityper recomb --input norovirus_clear_recombinant.fasta \
    --query QUERY --virus norovirus -o norovirus_output

calicityper recomb --input sapovirus_clear_recombinant.fasta \
    --query QUERY --virus sapovirus -o sapovirus_output
```

OpenRDP must be on `PATH`, or `CALICITYPER_OPENRDP` must point at it. Without
it the run still completes but reports `not_assessed`, not `none` — see
"Reading the output" below. Each run takes roughly 20–50 minutes.

## What they show

| | norovirus | sapovirus |
|---|---|---|
| planted breakpoint | ORF1/ORF2 junction | position 5,139 |
| status | `strong` | `likely` |
| score | 0.898 | 0.787 |
| localised to | 5,060–5,350 | 4,949–5,415 |
| supported by | BOOTSCAN, GENECONV, 3SEQ | BOOTSCAN, GENECONV, 3SEQ, RDP |
| region reported | ORF1–ORF2 junction | RdRp–VP1 junction |

The last row is the point of the pair. Caliciviruses recombine at the junction
upstream of the capsid, but that junction is between ORFs in norovirus and
*within* ORF1 in sapovirus, whose non-structural proteins and VP1 are one
polyprotein. `--virus` selects the genome model, and the same event is
therefore named correctly for each genus. Passing the wrong model does not
suppress the call, but it mislabels where the breakpoint falls.

Both runs are deterministic: repeating either with the same seed returns an
identical status, score and event count.

## Reading the output

`recombination.tsv` is one row per event; `recombination.json` additionally
carries the full provenance block — the parameters, the scanners that ran, the
scanners that were skipped and why, and the reference-set identity.

The `status` column distinguishes three different things, and the difference
matters:

- `none` — scanners ran and found no support for recombination.
- `weak` / `possible` / `likely` / `strong` — increasing evidence.
- `not_assessed`, `failed`, `insufficient_sequence` — **no answer was
  obtained.** A scanner timed out, crashed, or was never found, or the
  sequence was too short. These are not negative results and must not be read
  as "no recombination". `recombination.json` names the reason under
  `provenance.tools_skipped`.

A `score` of 0.0 alongside `not_assessed` means nothing was measured, not that
nothing was found.

## Caveats

These are synthetic datasets with a single clean breakpoint and well-separated
parents, so they demonstrate the mechanics rather than establish sensitivity
on real data. Discordance between the VP1 and RdRp type of a real sequence is
a reason to run this analysis, not evidence of recombination on its own.
