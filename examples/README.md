# Example run

Four public GenBank genomes, typed with the default settings:

```bash
calicityper type examples/example.fasta -o example_typed.tsv --threads 4 --report example
```

| Accession | RIVM call | Typing | Novelty screen |
|---|---|---|---|
| [MW305639.1](https://www.ncbi.nlm.nih.gov/nuccore/MW305639.1) | GII.4 Sydney[P16] | ORF2 `GII.4_Sydney` (100%); ORF1 **unassigned**: `GII.P16` reached 65%, below the 70% threshold | both `known` |
| [MT451991.1](https://www.ncbi.nlm.nih.gov/nuccore/MT451991.1) | GI.3[P13] | `GI.P13` and `GI.3`, both 100% | both `known` |
| [MN509084.1](https://www.ncbi.nlm.nih.gov/nuccore/MN509084.1) | sapovirus GII.2 | genogroup only | `not_assessable`: no sapovirus reference sets bundled yet |
| [AB039780.1](https://www.ncbi.nlm.nih.gov/nuccore/AB039780.1) | GII.8[P8] | ORF2 `GII.8` (100%); ORF1 **unassigned**, no P-type drew any support | ORF1 **flagged** `putative_new_genotype` |

Each row is there to show something:

- **MW305639.1**: a near-threshold call. ORF1 is left `unassigned`, the
  runner-up column shows the best candidate, and `note` says why. The novelty
  screen still reads it as `known`: its average distance to GII.P16 (0.078) is
  inside that P-type's bound. Typing asks *which type does the tree support*;
  the screen asks *does it fall outside every known type*. They can disagree,
  and the disagreement is informative.
- **MN509084.1**: a non-norovirus. Detected to genogroup, nothing more.
- **AB039780.1**: a flag. Its RdRp is on average 0.190 (nt) from the nearest
  P-type, GII.P6, well outside that cluster's bound (0.127). RIVM calls it
  GII.P8, but no GII.P8 reference in the bundled set is close to it, so the
  flag here points at a gap in the reference set rather than a new virus. That
  is the right reading of any flag: *look at this one*, not *this is new*.

## Files

| File | Written by |
|---|---|
| `example.fasta` | input |
| `example_typed.tsv` | `-o` |
| `example.json` | `--report` |
| `example.md` | `--report` |
| `example.html` | `--report`, built here with `--html-cdn` |

`example.html` has the calls table with flagged rows highlighted, the novelty
screen, one alignment panel per typed region (the query against its closest
reference, with only the differences drawn), and a tree per region. It was
built with `--html-cdn`, which loads Bokeh from its CDN: the file is 190 KB
instead of about 5 MB, but needs an internet connection to draw the figures.
GitHub shows HTML as source, so download it and open it in a browser.

Bootstrapping is seeded (`--seed`, default 1), so the same version gives the
same numbers.
