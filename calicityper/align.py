"""Add query sequences to a pre-computed reference alignment.

backend 'mafft'    : mafft --addfragments --keeplength, ALL queries of one (region, genogroup) in a single call
                     (default; the reference alignment is never recomputed and columns stay fixed)
(the original one-query-at-a-time ClustalW path is still available as typing.type_region(..., backend="clustalw"))
"""
import hashlib, os, subprocess, tempfile
from functools import lru_cache
from pathlib import Path
import numpy as np
from . import config as cg
from . import tree


def _refset_cache_dir():
    """Return a writable directory for pre-encoded RefSet .npz files, or None."""
    env = os.environ.get("CALICITYPER_CACHE")
    candidates = [Path(env) / "refset_cache"] if env else []
    base = Path(os.environ.get("XDG_CACHE_HOME") or (Path.home() / ".cache"))
    candidates += [base / "calicityper" / "refset_cache", cg.DATA / ".cache" / "refset_cache"]
    for d in candidates:
        try:
            d.mkdir(parents=True, exist_ok=True)
            probe = d / f".probe_{os.getpid()}"
            probe.touch()
            probe.unlink()
            return d
        except OSError:
            continue
    return None


def _pair_sha16(aln_path, orf_path):
    h = hashlib.sha256()
    for p in (aln_path, orf_path):
        h.update(Path(p).read_bytes())
        h.update(b"\0")
    return h.hexdigest()[:16]


class RefSet:
    """Reference alignment of one region/genogroup, loaded once (with .npz disk cache)."""

    def __init__(self, region, group, aln_path=None, orf_path=None):
        if aln_path is not None and orf_path is not None:
            aln, orf = Path(aln_path), Path(orf_path)
        else:
            paths = cg.ref_paths(region, group)
            if paths is None:
                raise FileNotFoundError(f"no {region} reference set for {group}")
            aln, orf = paths
        self.aln_path = str(aln)
        cache_file = None
        cdir = _refset_cache_dir()
        if cdir is not None:
            try:
                cache_file = cdir / f"refset_{_pair_sha16(aln, orf)}.npz"
            except OSError:
                cache_file = None
        if cache_file is not None and cache_file.exists():
            try:
                with np.load(cache_file, allow_pickle=False) as z:
                    self.names = [str(x) for x in z["names"].tolist()]
                    self.seqs = [str(x) for x in z["seqs"].tolist()]
                    self.types = [str(x) for x in z["types"].tolist()]
                    self.codes = np.array(z["codes"], dtype=np.int8, copy=True)
                    self.occupancy = np.array(z["occupancy"], dtype=np.float64, copy=True)
                self.types_by_acc = dict(zip(self.names, self.types))
                self.length = len(self.seqs[0]) if self.seqs else 0
                return
            except (OSError, ValueError, KeyError, EOFError) as exc:
                import logging
                logging.getLogger(__name__).debug("ignoring unreadable refset cache %s: %s", cache_file, exc)
        from Bio import AlignIO, SeqIO
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
        if cache_file is not None:
            tmp = cache_file.with_name(f".{cache_file.name}.{os.getpid()}.tmp.npz")
            try:
                np.savez(
                    tmp,
                    names=np.asarray(self.names, dtype=str),
                    seqs=np.asarray(self.seqs, dtype=str),
                    types=np.asarray(self.types, dtype=str),
                    codes=self.codes,
                    occupancy=self.occupancy,
                )
                os.replace(tmp, cache_file)
            except OSError:
                try:
                    tmp.unlink(missing_ok=True)
                except OSError:
                    pass

    def fasta(self, exclude=()):
        ex = set(exclude)
        return "".join(f">{n}\n{s}\n" for n, s in zip(self.names, self.seqs) if n not in ex)


def warm_cache(data_dir=None):
    """Pre-encode all available .aln/.orf reference pairs into the .npz cache."""
    base = Path(data_dir) if data_dir is not None else cg.DATA
    dirs = [
        ("orf1", base / "reference_sequences_RdRp"),
        ("orf2", base / "reference_sequences_VP1"),
        ("orf1", base / "reference_sequences_sapovirus_RdRp"),
        ("orf2", base / "reference_sequences_sapovirus_VP1"),
    ]
    built = 0
    for region, d in dirs:
        if not d.is_dir():
            continue
        for aln in sorted(d.glob("*.aln")):
            orf = aln.with_suffix(".orf")
            if orf.exists():
                RefSet(region, aln.stem, aln_path=aln, orf_path=orf)
                built += 1
    return built


@lru_cache(maxsize=None)
def refset(region, group):
    return RefSet(region, group)


def _mafft(ref, queries, exclude, threads):
    """queries: list of (id, seq). Returns list of aligned strings of reference-alignment length."""
    with tempfile.TemporaryDirectory(prefix="calicityper_") as td:
        rf, qf = os.path.join(td, "ref.fa"), os.path.join(td, "q.fa")
        Path(rf).write_text(ref.fasta(exclude), encoding="utf-8")
        with open(qf, "w", encoding="utf-8") as fh:
            for i, (_, s) in enumerate(queries):
                fh.write(f">q{i}\n{str(s).upper().replace('U', 'T')}\n")
        cmd = [cg.exe("mafft"), "--quiet", "--thread", str(max(1, threads)), "--keeplength", "--6merpair",
               "--addfragments", qf, rf]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            err = (p.stderr or "").strip()[-2000:] or "(no stderr)"
            raise RuntimeError(f"MAFFT failed (exit {p.returncode}) aligning {len(queries)} queries "
                               f"to {ref.aln_path}:\n{err}")
        out = p.stdout
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
