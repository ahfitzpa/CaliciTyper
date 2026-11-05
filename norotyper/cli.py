import argparse, os, sys, tempfile
from multiprocessing import Pool
from Bio import SeqIO
import pandas as pd
from . import typing as ty, __version__

_ARGS = None
def _work(item):
    rid, seq, hit = item
    a = _ARGS
    return ty.type_query(rid, seq, hit, nrep=a.bootstraps, seed=a.seed, min_support=a.min_support,
                         max_refs_per_type=a.max_refs_per_type)

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
        hits = ty.blast_genogroup(q)
    items = [(i, s, hits.get(i)) for i, s in recs]
    with Pool(a.threads, initializer=_init, initargs=(a,)) as pool:
        rows = pool.map(_work, items, chunksize=1)
    df = pd.DataFrame(rows)
    (df.to_csv(sys.stdout, sep="\t", index=False) if a.out == "-" else df.to_csv(a.out, sep="\t", index=False))

if __name__ == "__main__":
    main()
