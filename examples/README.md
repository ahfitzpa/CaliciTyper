# Example run

Three public GenBank genomes, typed with the default settings:

```bash
calicityper type examples/example.fasta -o example_typed.tsv --threads 4 --report example
```

| Accession | What it is (RIVM call) | What calicityper reports |
|---|---|---|
| [MW305639.1](https://www.ncbi.nlm.nih.gov/nuccore/MW305639.1) | norovirus GII.4 Sydney[P16] | ORF2 `GII.4_Sydney` at 100%; ORF1 **unassigned**: `GII.P16` reached 65%, below the 70% threshold |
| [MT451991.1](https://www.ncbi.nlm.nih.gov/nuccore/MT451991.1) | norovirus GI.3[P13] | `GI.P13` and `GI.3`, both at 100% |
| [MN509084.1](https://www.ncbi.nlm.nih.gov/nuccore/MN509084.1) | sapovirus GII.2 | genogroup only (`Caliciviridae_Sapovirus_GII`); sapovirus genotyping is not implemented |

The first row is deliberate. It shows a call that does not clear the
threshold: the type is left `unassigned`, the runner-up column still shows
what the best candidate was, and `note` says why.

## Files

| File | Written by |
|---|---|
| `example.fasta` | input |
| `example_typed.tsv` | `-o` |
| `example.json` | `--report` |
| `example.md` | `--report` |

`--report` also writes `example.html`, a self-contained page with alignment
panels and trees. It is about 5 MB, so it is not stored here; rerun the
command above to get it.

Bootstrapping is seeded (`--seed`, default 1), so the same version gives the
same numbers.
