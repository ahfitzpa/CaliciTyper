import argparse, os, sys, tempfile
from collections import defaultdict
from multiprocessing import Pool
from Bio import SeqIO
import pandas as pd
from . import typing as ty, align, tree, __version__

_ARGS = None
def _work(chunk):
    """chunk = (region, group, [(row_index, seq), ...]): ONE MAFFT call for the chunk, then bootstrap each query."""
    region, g, lst = chunk
    a = _ARGS
    ref, aln = align.align_queries(region, g, [(str(i), s) for i, s in lst], threads=1)
    return [((i, region), ty.type_aligned(ref, qa, nrep=a.bootstraps, seed=a.seed, min_support=a.min_support,
                                          max_refs_per_type=a.max_refs_per_type))
            for (i, _), qa in zip(lst, aln)]

def _init(a):
    global _ARGS
    _ARGS = a

def main(argv=None):
    p = argparse.ArgumentParser(prog="norotyper", description="Norovirus genotyping (RIVM-style: BLAST, profile "
                                "alignment, NJ tree, bootstrap).")
    p.add_argument("fasta"); p.add_argument("-o", "--out", default="-")
    p.add_argument("--bootstraps", type=int, default=100, help="bootstrap replicates (RIVM: 100)")
    p.add_argument("--min-support", type=float, default=70.0, help="minimum bootstrap %% to assign a type (RIVM: 70)")
    p.add_argument("--max-refs-per-type", type=int, default=2, help="reference sequences per type used in the tree (RIVM uses a representative subset; 2)")
    p.add_argument("--threads", type=int, default=1); p.add_argument("--seed", type=int, default=1)
    p.add_argument("--version", action="version", version=__version__)
    a = p.parse_args(argv)
    recs = [(r.id, str(r.seq)) for r in SeqIO.parse(a.fasta, "fasta")]
    with tempfile.TemporaryDirectory(prefix="norotyper_") as td:
        q = os.path.join(td, "all.fa")
        SeqIO.write([r for r in SeqIO.parse(a.fasta, "fasta")], q, "fasta")
        hits = ty.blast_genogroup(q, a.threads)
    rows, todo = {}, defaultdict(list)             # todo[(region, group)] -> [(index, seq)]
    for i, (rid, s) in enumerate(recs):
        row, jobs = ty.plan_query(rid, s, hits.get(rid))
        rows[i] = row
        for region, g, seq in jobs:
            todo[(region, g)].append((i, seq))
    tree.warmup()
    chunks = []                                     # one aligner call per chunk of queries from one region/genogroup
    for (region, g), lst in todo.items():
        size = max(8, -(-len(lst) // (2 * a.threads)))
        chunks += [(region, g, lst[k:k + size]) for k in range(0, len(lst), size)]
    if a.threads > 1 and len(chunks) > 1:
        with Pool(a.threads, initializer=_init, initargs=(a,)) as pool:
            results = [r for part in pool.map(_work, chunks, chunksize=1) for r in part]
    else:
        _init(a); results = [r for c in chunks for r in _work(c)]
    for (i, region), res in sorted(results, key=lambda r: (r[0][0], r[0][1])):
        ty.apply_result(rows[i], region, res, a.min_support)
    df = pd.DataFrame([rows[i] for i in range(len(recs))])
    (df.to_csv(sys.stdout, sep="\t", index=False) if a.out == "-" else df.to_csv(a.out, sep="\t", index=False))

if __name__ == "__main__":
    main()
