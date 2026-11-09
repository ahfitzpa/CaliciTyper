import argparse, os, subprocess, sys, tempfile
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
from Bio import SeqIO
import pandas as pd
from . import typing as ty, align, tree, __version__
from . import recomb_config as _rc

_ARGS = None
def _work(chunk):
    """chunk = (region, group, [(row_index, seq), ...]): ONE MAFFT call for the chunk, then bootstrap each query."""
    region, g, lst = chunk
    a = _ARGS
    ref, aln = align.align_queries(region, g, [(str(i), s) for i, s in lst], threads=1)
    out = []
    for (i, _), qa in zip(lst, aln):
        res = ty.type_aligned(ref, qa, nrep=a.bootstraps, seed=a.seed, min_support=a.min_support,
                              max_refs_per_type=a.max_refs_per_type)
        # keep the aligned query only when a report will use it: for thousands of
        # queries these strings are the bulk of what crosses back from a worker.
        out.append(((i, region), res, (g, qa) if a._want_aln else None))
    return out

def _init(a):
    global _ARGS
    _ARGS = a
    # Compile/load the numba kernel in THIS worker. Python 3.14+ defaults to the
    # "forkserver" start method on Linux, so workers are not forked from the
    # parent and a warmup() done there does not carry over -- every worker has to
    # load the kernel itself. njit(cache=True) makes this a disk-cache load
    # rather than a recompile, but it still has to happen here.
    tree.warmup()

def _add_type_args(p):
    p.add_argument("fasta"); p.add_argument("-o", "--out", default="-")
    p.add_argument("--bootstraps", type=int, default=100, help="bootstrap replicates (RIVM: 100)")
    p.add_argument("--min-support", type=float, default=70.0, help="minimum bootstrap %% to assign a type (RIVM: 70)")
    p.add_argument("--max-refs-per-type", type=int, default=2, help="reference sequences per type used in the tree (RIVM uses a representative subset; 2)")
    p.add_argument("--threads", type=int, default=1); p.add_argument("--seed", type=int, default=1)
    p.add_argument("--version", action="version", version=__version__)
    r = p.add_argument_group("reports (the TSV stays the primary output)")
    r.add_argument("--report", metavar="PREFIX",
                   help="write PREFIX.md, PREFIX.json and PREFIX.html as well as the TSV")
    r.add_argument("--json", metavar="PATH", help="write results as JSON (runner-ups structured)")
    r.add_argument("--md", metavar="PATH", help="write a Markdown report")
    r.add_argument("--html", metavar="PATH",
                   help="write a self-contained HTML report with Bokeh alignment panels and trees")
    r.add_argument("--report-max-panels", type=int, default=20, metavar="N",
                   help="most alignment panels to draw, highest-support calls first [default: 20]")
    r.add_argument("--no-tree", action="store_true", help="skip the per-region trees in the HTML report")
    return p


def _add_recomb_args(p):
    p.add_argument("--alignment", help="aligned FASTA containing the query and a reference panel")
    p.add_argument("--input", help="unaligned FASTA; aligned with MAFFT before scanning")
    p.add_argument("--query", action="append", default=[],
                   help="query sequence id; repeat for several, or omit for all non-reference ids")
    p.add_argument("-o", "--output", default="recomb_out", metavar="DIR",
                   help="output directory for the JSON and TSV [default: recomb_out]")
    p.add_argument("--reference-metadata", help="TSV/CSV/JSON mapping reference ids to lineages")
    p.add_argument("--config", help="YAML or TOML configuration file")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--min-length", type=int, help="override min_sequence_length")
    p.add_argument("--pvalue-threshold", type=float, help="override pvalue_threshold")
    p.add_argument("--max-breakpoint-distance", type=int, help="override max_breakpoint_distance")
    p.add_argument("--methods", nargs="+", metavar="M",
                   help=f"subset of {' '.join(_rc.ALL_METHODS)}")
    p.add_argument("--skip-gard", action="store_true", help="do not run HyPhy GARD")
    p.add_argument("--skip-rdp", action="store_true", help="do not run the OpenRDP scanners")
    p.add_argument("--no-segment-trees", action="store_true",
                   help="skip the segment phylogeny check (faster, loses a scoring component)")
    p.add_argument("--verbose", action="store_true")
    p.add_argument("--dry-run", action="store_true",
                   help="validate inputs, report which tools would run, then stop")
    return p


def main(argv=None):
    top = argparse.ArgumentParser(
        prog="norotyper",
        description="Norovirus genotyping (RIVM-style: BLAST, profile alignment, NJ tree, "
                    "bootstrap) and recombination analysis.")
    top.add_argument("--version", action="version", version=__version__)
    sub = top.add_subparsers(dest="command", required=True, metavar="{type,recomb}")
    _add_type_args(sub.add_parser("type", help="assign genogroup, ORF1 P-type and ORF2 capsid type"))
    _add_recomb_args(sub.add_parser("recomb", help="detect and score recombination (norovirus only)"))

    a = top.parse_args(argv)
    if a.command == "recomb":
        return _run_recomb(a)
    return _run_type(a)


def _run_type(a):
    if a.report:
        a.md = a.md or f"{a.report}.md"
        a.json = a.json or f"{a.report}.json"
        a.html = a.html or f"{a.report}.html"
    a._want_aln = bool(a.html)          # alignments/trees are only needed for the HTML
    ty.cg.check_data()          # clear error if the bundled reference data is missing
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
    tree.warmup()          # primes the on-disk kernel cache the workers load from
    chunks = []                                     # one aligner call per chunk of queries from one region/genogroup
    for (region, g), lst in todo.items():
        size = max(8, -(-len(lst) // (2 * a.threads)))
        chunks += [(region, g, lst[k:k + size]) for k in range(0, len(lst), size)]
    if a.threads > 1 and len(chunks) > 1:
        with Pool(a.threads, initializer=_init, initargs=(a,)) as pool:
            results = [r for part in pool.map(_work, chunks, chunksize=1) for r in part]
    else:
        _init(a); results = [r for c in chunks for r in _work(c)]
    alns = {}                                       # (row index, region) -> (genogroup, aligned query)
    for (i, region), res, aln in sorted(results, key=lambda r: (r[0][0], r[0][1])):
        ty.apply_result(rows[i], region, res, a.min_support)
        if aln is not None:
            alns[(i, region)] = aln
    ordered = [rows[i] for i in range(len(recs))]
    df = pd.DataFrame(ordered)
    (df.to_csv(sys.stdout, sep="\t", index=False) if a.out == "-" else df.to_csv(a.out, sep="\t", index=False))
    if a.md or a.json or a.html:
        _write_reports(a, ordered, alns)

def _support(row, region):
    try:
        return float(row.get(f"{region}_support") or 0.0)
    except (TypeError, ValueError):
        return 0.0


def _panels(a, rows, alns):
    """One alignment panel per call: the query against the references of the type
    it was assigned, plus the references of the runner-up type for contrast."""
    from . import report
    want = [(k, _support(rows[k[0]], k[1])) for k in alns
            if (rows[k[0]].get(f"{k[1]}_type") or "") not in ("", "unassigned")]
    want.sort(key=lambda kv: -kv[1])
    figs = []
    for (i, region), _ in want[:max(0, a.report_max_panels)]:
        g, qa = alns[(i, region)]
        ref = align.refset(region, g)
        row = rows[i]
        called = row.get(f"{region}_type")
        runner = next((t for t, _ in (_parse_runners(row.get(f"{region}_runners")) or [])
                       if t != called), None)
        names, seqs = [str(row["query"])[:38]], [qa]
        for t, tag in ((called, ""), (runner, " (runner-up)")):
            if not t:
                continue
            for nm, sq, ty_ in zip(ref.names, ref.seqs, ref.types):
                if ty_ == t and len(names) < 12:
                    names.append(f"{nm} {t}{tag}")
                    seqs.append(sq)
        if len(seqs) < 2:
            continue
        figs.append(report.alignment_figure(
            names, seqs,
            f"{row['query']} \u2014 {region.upper()} \u2192 {called} "
            f"({_support(row, region):.0f}% support)",
            highlight=names[0]))
    return [f for f in figs if f is not None]


def _parse_runners(s):
    out = []
    for part in (s or "").split("|"):
        if ":" in part:
            t, v = part.rsplit(":", 1)
            try:
                out.append((t, float(v)))
            except ValueError:
                pass
    return out


def _trees(a, rows, alns):
    """One neighbour-joining tree per region/genogroup: queries among references."""
    from . import report
    import numpy as np
    by_group = defaultdict(list)
    for (i, region), (g, qa) in alns.items():
        by_group[(region, g)].append((i, qa))
    figs, notes = [], []
    for (region, g), lst in sorted(by_group.items()):
        ref = align.refset(region, g)
        seen, keep = {}, []
        for k, t in enumerate(ref.types):             # same subset the tree uses when typing
            seen[t] = seen.get(t, 0) + 1
            if seen[t] <= (a.max_refs_per_type or 10 ** 9):
                keep.append(k)
        room = report.MAX_TREE_TAXA - len(keep)
        if room < 2:
            notes.append(f"{region.upper()}/{g}: tree skipped, the reference set alone "
                         f"exceeds {report.MAX_TREE_TAXA} taxa.")
            continue
        chosen = lst[:room]
        if len(chosen) < len(lst):
            notes.append(f"{region.upper()}/{g}: tree shows {len(chosen)} of {len(lst)} queries "
                         f"(capped at {report.MAX_TREE_TAXA} taxa for legibility).")
        qnames = [str(rows[i]["query"])[:38] for i, _ in chosen]
        names = qnames + [f"{ref.names[k]} {ref.types[k]}" for k in keep]
        codes = np.vstack([tree.encode([qa for _, qa in chosen]), ref.codes[keep]]) \
            if chosen else ref.codes[keep]
        D = tree.tn93(codes)
        figs.append(report.tree_figure(names, D, f"{region.upper()} / {g} \u2014 "
                                                 f"{len(chosen)} queries among references",
                                       highlight=set(qnames)))
    return [f for f in figs if f is not None], notes


def _write_reports(a, rows, alns):
    from . import report
    meta = dict(bootstraps=a.bootstraps, min_support=a.min_support,
                max_refs_per_type=a.max_refs_per_type, seed=a.seed, input=a.fasta)
    if a.json:
        report.to_json(rows, a.json, meta=meta)
    if a.md:
        report.to_markdown(rows, a.md, meta=meta)
    if a.html:
        figs, notes = ([], [])
        if not a.no_tree:
            figs, notes = _trees(a, rows, alns)
        figs = _panels(a, rows, alns) + figs
        if not report.HAVE_BOKEH:
            notes.append("bokeh is not installed, so no figures were drawn "
                         "(pip install bokeh).")
        report.to_html(rows, a.html, figures=figs, meta=meta, notes=notes)


def _run_recomb(a):
    """`norotyper recomb`: scan an alignment for recombination. Norovirus only."""
    from . import recomb as R

    cfg = _rc.load(a.config)
    if a.min_length is not None:
        cfg.min_sequence_length = a.min_length
    if a.pvalue_threshold is not None:
        cfg.pvalue_threshold = a.pvalue_threshold
    if a.max_breakpoint_distance is not None:
        cfg.max_breakpoint_distance = a.max_breakpoint_distance
    if a.methods:
        cfg.methods = tuple(m.lower() for m in a.methods)
    if a.skip_rdp:
        cfg.methods = ()
    if a.skip_gard:
        cfg.run_gard = False
    if a.no_segment_trees:
        cfg.build_segment_trees = False
    cfg.validate()

    if not a.alignment and not a.input:
        sys.exit("norotyper recomb: give --alignment (aligned FASTA) or --input (unaligned FASTA)")

    outdir = Path(a.output)
    outdir.mkdir(parents=True, exist_ok=True)

    aln_path = a.alignment
    if not aln_path:
        # --input is unaligned: align it first, because every scanner needs a
        # column-consistent alignment and silently mis-scanning unaligned input
        # would be worse than spending the MAFFT call.
        aln_path = str(outdir / "aligned.fasta")
        if a.verbose:
            print(f"[recomb] aligning {a.input} -> {aln_path}", file=sys.stderr)
        with open(aln_path, "w") as fh:
            subprocess.run([ty.cg.exe("mafft"), "--quiet", "--auto",
                            "--thread", str(max(1, a.threads)), a.input],
                           stdout=fh, check=True)

    aln = R.read_alignment(aln_path)
    queries = a.query or [n for n in aln.names]
    missing = [q for q in queries if q not in aln.names]
    if missing:
        sys.exit(f"norotyper recomb: query id(s) not in the alignment: {missing[:5]}")

    if a.dry_run:
        exe = cfg.openrdp()
        print(f"alignment        : {aln_path}")
        print(f"  sequences      : {len(aln.names)}")
        print(f"  columns        : {aln.length}")
        print(f"  masked to N    : {getattr(aln, 'masked_chars', 0)}")
        print(f"  informative    : {R.informative_sites(aln)}")
        print(f"queries          : {len(queries)}")
        print(f"methods requested: {' '.join(cfg.methods) or '(none)'}")
        print(f"OpenRDP          : {exe or 'NOT FOUND (set $NOROTYPER_OPENRDP)'}")
        print(f"segment trees    : {'yes' if cfg.build_segment_trees else 'no'}")
        print(f"GARD             : {'requested' if cfg.run_gard else 'off'}")
        for q in queries[:10]:
            ok, why, st = R.qc_query(aln, q, cfg)
            print(f"  {q}: {'OK' if ok else 'WOULD SKIP - ' + why}")
        return 0

    if a.verbose:
        print(f"[recomb] {len(queries)} quer{'y' if len(queries)==1 else 'ies'}, "
              f"{len(aln.names)} sequences, {aln.length} columns", file=sys.stderr)

    results, provenance = R.analyse(aln_path, queries, cfg,
                                    metadata_path=a.reference_metadata,
                                    workdir=str(outdir / "work"))
    js = R.write_json(results, provenance, outdir / "recombination.json")
    tsv = R.write_tsv(results, outdir / "recombination.tsv")

    for t in provenance["tools_skipped"]:
        print(f"[recomb] skipped {t['tool']}: {t['reason']}", file=sys.stderr)
    for r in results:
        rec = r["recombination"]
        print(f"{r['sequence_id']}\t{rec['status']}\t{rec['score']}\t"
              f"{len(rec['events'])} event(s)")
    print(f"[recomb] wrote {js} and {tsv}", file=sys.stderr)
    # exit 0 even when nothing is recombinant: "no recombination" is a result
    return 0


if __name__ == "__main__":
    main()
