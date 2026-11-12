# Attribution

## Upstream project

`calicityper` is a fork and partial rewrite of
[**zhuzhanji/NorovirusGenotyping**](https://github.com/zhuzhanji/NorovirusGenotyping),
a University College Dublin diploma project (GitHub:
[@zhuzhanji](https://github.com/zhuzhanji)).

What came from upstream:

- The curated norovirus reference sequence sets and BLAST database now in
  `calicityper/data/`, and the taxonomy tables in `validation/`.
- The RIVM comparison call set, `validation/results_rivm.csv`.

What is new in this fork (the `calicityper` package):

- A reimplementation as an installable, path-independent Python package with a
  command-line interface.
- MAFFT profile alignment onto fixed reference columns, replacing per-query
  ClustalW calls.
- Neighbour-joining with 100 bootstrap replicates and a >= 70% support
  threshold, in place of the original 5-replicate maximum-likelihood tree.
- Leave-one-out validation harness and the RIVM comparison script.

## Method

The typing method replicates the published approach of the **RIVM Norovirus
Typing Tool** (National Institute for Public Health and the Environment,
the Netherlands), available at <https://www.rivm.nl/mpf/typingtool/norovirus/>.

Please cite:

> Kroneman A, Vennema H, Deforche K, van der Avoort H, Peñaranda S, Oberste MS,
> Vinjé J, Koopmans M. An automated genotyping tool for enteroviruses and
> noroviruses. *Journal of Clinical Virology* 2011; 51(2): 121-125.
> doi:10.1016/j.jcv.2011.03.006

Reference sequence sets follow the 2019 norovirus classification:

> Chhabra P, de Graaf M, Parra GI, Chan MC, Green K, Martella V, Wang Q,
> White PA, Katayama K, Vennema H, Koopmans MPG, Vinjé J. Updated
> classification of norovirus genogroups and genotypes.
> *Journal of General Virology* 2019; 100(10): 1393-1406.
> doi:10.1099/jgv.0.001318

`calicityper` is an independent reimplementation. It is **not** produced,
endorsed by, or affiliated with RIVM.

## Third-party software

BLAST+, MAFFT, ClustalW, Biopython, NumPy, pandas and Numba are used under
their own licences.
