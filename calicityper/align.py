"""Add query sequences to a pre-computed reference alignment.

backend 'mafft'    : mafft --addfragments --keeplength, ALL queries of one (region, genogroup) in a single call
                     (default; the reference alignment is never recomputed and columns stay fixed)
(the original one-query-at-a-time ClustalW path is still available as typing.type_region(..., backend="clustalw"))
"""
import os, subprocess, tempfile
from functools import lru_cache
import numpy as np
from Bio import AlignIO, SeqIO
from . import config as cg
from . import tree


class RefSet:
    """Reference alignment of one region/genogroup, loaded once."""

    def __init__(self, region, group):
        d = cg.REF_DIR[region]
        aln = d / f"{group}.aln"
        orf = d / f"{group}.orf"
        if not aln.exists():
            aln, orf = d / "all.aln", d / "all.orf"
        self.aln_path = str(aln)
        self.types_by_acc = {r.id.split("|")[2]: r.id.split("|")[1] for r in SeqIO.parse(orf, "fasta")}
        a = AlignIO.read(aln, "clustal")
        self.names = [r.id.split("|")[-1] for r in a]          # bare accession
        self.seqs = [str(r.seq).upper() for r in a]
        self.length = len(self.seqs[0])
        keep = [i for i, n in enumerate(self.names) if n in self.types_by_acc]
        self.names = [self.names[i] for i in keep]
        self.seqs = [self.seqs[i] for i in keep]
        self.types = [self.types_by_acc[n] for n in self.names]
        self.codes = tree.encode(self.seqs)
        self.occupancy = (self.codes < 4).mean(0)             # fraction of references with a base, per column

    def fasta(self, exclude=()):
        ex = set(exclude)
        return "".join(f">{n}\n{s}\n" for n, s in zip(self.names, self.seqs) if n not in ex)


@lru_cache(maxsize=None)
def refset(region, group):
    return RefSet(region, group)


def _mafft(ref, queries, exclude, threads):
    """queries: list of (id, seq). Returns list of aligned strings of reference-alignment length."""
    with tempfile.TemporaryDirectory(prefix="calicityper_") as td:
        rf, qf = os.path.join(td, "ref.fa"), os.path.join(td, "q.fa")
        open(rf, "w").write(ref.fasta(exclude))
        with open(qf, "w") as fh:
            for i, (_, s) in enumerate(queries):
                fh.write(f">q{i}\n{str(s).upper().replace('U', 'T')}\n")
        cmd = [cg.exe("mafft"), "--quiet", "--thread", str(max(1, threads)), "--keeplength", "--6merpair",
               "--addfragments", qf, rf]
        out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    got = {}
    name, buf = None, []
    for line in out.splitlines():
        if line.startswith(">"):
            if name is not None:
                got[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    got[name] = "".join(buf).upper()
    return [got[f"q{i}"] for i in range(len(queries))]


def align_queries(region, group, queries, exclude=(), threads=1):
    """Align (id, seq) queries to the reference set; returns (RefSet, [aligned strings of reference length])."""
    ref = refset(region, group)
    return ref, (_mafft(ref, queries, exclude, threads) if queries else [])
