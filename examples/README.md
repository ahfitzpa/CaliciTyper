# Example run

Four public GenBank genomes, typed with the default settings:

```bash
calicityper type examples/example.fasta -o example_typed.tsv --threads 4 --report example
```

| Accession | RIVM call | Typing | Putative novel ORF1/ORF2 type |
|---|---|---|---|
| [MW305639.1](https://www.ncbi.nlm.nih.gov/nuccore/MW305639.1) | GII.4 Sydney[P16] | `GII.P16` and `GII.4_Sydney`, both 100% | both `known` |
| [MT451991.1](https://www.ncbi.nlm.nih.gov/nuccore/MT451991.1) | GI.3[P13] | `GI.P13` and `GI.3`, both 100% | both `known` |
| [MN509084.1](https://www.ncbi.nlm.nih.gov/nuccore/MN509084.1) | sapovirus GII.2 | genogroup only | `not_assessable`: no sapovirus reference sets bundled yet |
| [AB242258.2](https://www.ncbi.nlm.nih.gov/nuccore/AB242258.2) | GII.3[P21] | `GII.P21` and `GII.3`, both 100% | ORF1 **flagged** `putative_novel_type` (putative novel RdRp P-type), borderline |

What each row shows:

- **MW305639.1** and **MT451991.1**: ordinary genomes, typed as RIVM types
  them, and inside their clusters on both ORFs.
- **MN509084.1**: a non-norovirus. Detected to genogroup, nothing more.
- **AB242258.2**: a flag, and a borderline one. Its RdRp is on average 0.047
  (nt) from the GII.P21 references, just outside that cluster's 2x SD bound of
  0.044, although it types GII.P21 at 100%. GII.P21's references are very
  close to each other, so its bound is tight. Typing asks *which type does the
  tree support*; the screen asks *is it outside every known cluster*. A flag
  this close to the bound, on a sequence that types cleanly, is the screen's
  main false-positive mode. Read any flag as *look at this one*, not *this is
  new*.

## Files

| File | Written by |
|---|---|
| `example.fasta` | input |
| `example_typed.tsv` | `-o` |
| `example.json` | `--report` |
| `example.md` | `--report` |
| `example.html` | `--report`, built here with `--html-cdn` |

`example.html` has the calls table with flagged rows highlighted, the list of
putative novel ORF1/ORF2 types, one alignment panel per typed region (the query
against its closest reference, with only the differences drawn), and a tree per
region. It was built with `--html-cdn`, which loads Bokeh from its CDN: the
file is about 190 KB instead of about 5 MB, but needs an internet connection to
draw the figures. GitHub shows HTML as source, so download it and open it in a
browser.

Bootstrapping is seeded (`--seed`, default 1), so the same version gives the
same numbers.
