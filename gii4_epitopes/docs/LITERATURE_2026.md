# Four 2017–2026 papers: evidence register for the antigenic annotation

These papers were read on 2026-11-19. They are the evidence behind the additive
`calicityper antigenic` subcommand. Nothing in this file changes a typing call.

| Paper | DOI / PMCID | Licence | Tool-relevant finding |
|---|---|---|---|
| Parra GI et al. 2026, *J Virol* 100(2):e01804-25 "Variable and conserved B cell epitopes of GII.4 human noroviruses" | 10.1128/jvi.01804-25 / PMC12911863 | CC BY 4.0 | GII.4 variable antigenic sites A,C,D,E,G in P2; conserved A1431 footprint 402/403/504/506; VP1 length 530–550 aa across human genotypes -> residue numbers must be applied via alignment columns, never raw indices on another genotype. |
| Kendra JA et al. 2021, *PNAS* 118(11):e2015874118 "Antigenic cartography reveals complexities of genetic determinants…" | 10.1073/pnas.2015874118 / PMC7980451 | CC BY-NC-ND | **≥18 VP1 aa changes** required for antigenic distance ≥4 (≥16-fold blockade difference); **~22% of ≥18-change pairs were antigenically similar**, 10 pairs reached ≥4 with <18 changes; shifts need simultaneous changes at A/D/E/G; site C (339–341,376,377) had minimal effect. Table 1 coevolving pairs: 297–372, 309–459, 357–372, 357–397, 395–504, 397–448, 465–483, 119–145. |
| Tohma K et al. 2017, *mSphere* 2(3):e00187-17 "Phylogenetic analyses suggest factors other than the capsid… GII.2" | 10.1128/mSphereDirect.00187-17 / PMC5437133 | CC BY 4.0 | GII.2 VP1 >99.5% identical in aa over ~40 y while nt evolved clock-like via **synonymous third-codon changes**; GII.P16-GII.2 re-emergence associated with 4 conservative **RdRp surface substitutions** (S293T near residue 291 known to alter polymerase kinetics). |
| Tohma K et al. 2022, *Cell Reports* 39(6):110738 "Dynamic immunodominance hierarchy of neutralizing antibody responses to evolving GII.4 noroviruses" | 10.1016/j.celrep.2022.00447-8 | open access | Single mAbs map to **multiple antigenic sites**; immunodominance hierarchy changed through GII.4 evolution; neutralisation and HBGA blockade titres strongly correlated. |

## What was actually implemented (additive only)

- `calicityper/antigenic.py`: static GII.4 site residues, conserved residues,
  coevolving pairs, verified immunotype lookup, and a `candidate_antigenic_shift`
  screen (≥18 aa differences **and** ≥2 of A/D/E/G changed). The ~22% false
  positive rate is printed in every output's `caveat` field.
- **Residue correction (2026-11-19):** site residues are the FULL sets from the
  J Virol 2026 paper, not the sparse examples in Kendra 2021 Table 1: A
  294–298/368/372/373, C 339–341/375–378, D 393–398, E 407–414 (variable 407
  and 412–414), G 352/355–357/359/364; conserved A1431 footprint 402/403/504/506
  plus site F 254–257/327/404. The earlier version ended site C at 377 — fixed.
- `calicityper antigenic` CLI subcommand: annotates an aligned GII.4 VP1 AA
  alignment; writes TSV + provenance sidecar. Does not touch `type`, `novel`,
  `recomb`, their thresholds or schemas.
- `tests/test_antigenic.py`: 5 tests, all passing locally.

## What was deliberately NOT implemented

- Antigenic cartography (needs HBGA blockade/neutralisation data).
- Immunodominance ranking (needs mAb binding data).
- Polymerase-kinetics prediction from RdRp substitutions (needs enzyme assays).
- GII.2/GII.P16 RdRp surface substitution flag: recorded as evidence only;
  would need an RdRp AA alignment and validation before becoming a flag.
- Immunotype lookup rows marked `TO_VERIFY_FROM_FIG5C` in
  `../inputs/immunotype_regime.tsv` are NOT shipped in the tool — only the 8
  supported rows are.

## Licence caution

PNAS is **CC BY-NC-ND**: findings are cited, tables are not embedded verbatim.
J Virol and mSphere are CC BY 4.0; residue sets are re-derived with attribution.
