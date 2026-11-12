import argparse, os, subprocess, sys, tempfile
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
from Bio import SeqIO
import pandas as pd
from . import typing as ty
from . import config as cg, align, tree, __version__
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
    r.add_argument("--report-max-panels", type=int, default=None, metavar="N",
                   help="cap the number of alignment panels, in input order "
                        "[default: one per typed region, every call]")
    r.add_argument("--no-tree", action="store_true", help="skip the per-region trees in the HTML report")
    r.add_argument("--html-cdn", action="store_true",
                   help="link BokehJS from a CDN instead of embedding it: a much smaller "
                        "HTML file that needs internet to render [default: embed, stands alone]")
    p.add_argument("--screen-all", action="store_true",
                   help="run the putative novel screen on every call, not only on calls the "
                        "bootstrap did not accept [default: only unaccepted calls]")
    return p


def _add_novel_args(p):
    """`novel`: the 2x SD screen as a command of its own, norovirus only."""
    p.add_argument("fasta"); p.add_argument("-o", "--out", default="-")
    p.add_argument("--region", choices=["orf1", "orf2", "both"], default="both",
                   help="which region to screen: orf2 needs a complete VP1 reading frame, "
                        "orf1 needs the 762 nt 3' RdRp typing region [both]")
    p.add_argument("--bootstraps", type=int, default=100)
    p.add_argument("--min-support", type=float, default=70.0)
    p.add_argument("--max-refs-per-type", type=int, default=2)
    p.add_argument("--threads", type=int, default=1); p.add_argument("--seed", type=int, default=1)
    p.add_argument("--version", action="version", version=__version__)
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
    p.add_argument("--virus", choices=sorted(_rc.REGIONS_BY_VIRUS), default="norovirus",
                   help="which genome model to score breakpoints against. Norovirus puts the "
                        "hotspot at the ORF1/ORF2 junction; sapovirus encodes VP1 inside ORF1, "
                        "so its hotspot is the RdRp/VP1 junction [default: norovirus]")
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


def _add_setup_args(p):
    p.add_argument("--source", help="URL or local tarball of the reference set "
                                    "[default: download_url from the manifest]")
    p.add_argument("--data-dir", help="where to install [default: a versioned "
                                      "directory under ~/.cache/calicityper]")
    p.add_argument("--in-package", action="store_true",
                   help="install into the package's own data directory instead")
    p.add_argument("--manifest", help="alternative reference_manifest.json")
    p.add_argument("--verify-only", action="store_true",
                   help="check an existing data directory against the manifest and stop")
    p.add_argument("--force", action="store_true", help="replace an existing data directory")
    return p


def main(argv=None):
    top = argparse.ArgumentParser(
        prog="calicityper",
        description="Norovirus genotyping (RIVM-style: BLAST, profile alignment, NJ tree, "
                    "bootstrap) and recombination analysis.")
    top.add_argument("--version", action="version", version=__version__)
    sub = top.add_subparsers(dest="command", required=True, metavar="{type,novel,recomb,setup}")
    _add_type_args(sub.add_parser("type", help="assign genogroup, ORF1 P-type and ORF2 capsid type"))
    _add_novel_args(sub.add_parser(
        "novel", help="screen for putative novel norovirus ORF1/ORF2 types (Chhabra et al. 2019, 2x SD)"))
    _add_recomb_args(sub.add_parser("recomb", help="detect and score recombination (norovirus or sapovirus)"))
    _add_setup_args(sub.add_parser("setup", help="fetch and verify the reference data set"))

    a = top.parse_args(argv)
    if a.command == "recomb":
        return _run_recomb(a)
    if a.command == "setup":
        from . import setup_data
        return setup_data.main(a)
    if a.command == "novel":
        return _run_novel(a)
    return _run_type(a)


def _align_and_type(a):
    """BLAST, plan, align and type every record. Returns (recs, rows, alns, todo).

    Shared by `type` and `novel`: the putative novel screen needs exactly the
    same alignment to the genogroup reference set that typing produces, so it
    runs off one pipeline rather than a parallel one that could drift from it.
    """
    a._want_aln = True                  # the putative novel type screen needs the alignment, not only the HTML
    ty.cg.check_data()          # clear error if the bundled reference data is missing
    recs = [(r.id, str(r.seq)) for r in SeqIO.parse(a.fasta, "fasta")]
    with tempfile.TemporaryDirectory(prefix="calicityper_") as td:
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
    _apply_window_rules(rows, alns)
    return recs, rows, alns, todo


def _run_type(a):
    if a.report:
        a.md = a.md or f"{a.report}.md"
        a.json = a.json or f"{a.report}.json"
        a.html = a.html or f"{a.report}.html"
    recs, rows, alns, todo = _align_and_type(a)
    _screen_putative_novel(rows, alns, todo, only_unaccepted=not a.screen_all)
    ordered = [rows[i] for i in range(len(recs))]
    df = pd.DataFrame(ordered)
    (df.to_csv(sys.stdout, sep="\t", index=False) if a.out == "-" else df.to_csv(a.out, sep="\t", index=False))
    if a.md or a.json or a.html:
        _write_reports(a, ordered, alns)

def _run_novel(a):
    """`calicityper novel`: the Chhabra et al. 2019 2x SD screen, on its own.

    Norovirus only, and deliberately so. The criterion is a norovirus
    criterion: sapovirus has its own published cut-offs (Zhao et al. 2024) on
    nucleotide distances, and applying norovirus bounds to a sapovirus would
    produce a number that means nothing.

    The paper defines the rule on what it calls complete sequence, and this
    command enforces that rather than screening whatever it is given:

      ORF2 capsid  complete VP1 as AMINO ACIDS (an ATG..stop reading frame of
                   reference length), distances on the translation
      ORF1 P-type  the 762 nt 3' RdRp typing region as NUCLEOTIDES, >= 95%
                   covered -- the paper's P-types are nucleotide-based, not
                   amino acid

    Anything that does not meet its region's requirement is reported
    `not_assessable` with the reason, never screened against a partial
    sequence. A flag is tentative by construction: one sequence cannot
    designate a genotype, which needs two non-identical sequences from
    geographically distinct locations (the paper's own NA rule).
    """
    from . import putative_novel
    a.report = a.md = a.json = a.html = None
    recs, rows, alns, todo = _align_and_type(a)
    subs = {(i, region): seq for (region, g), lst in todo.items() for i, seq in lst}
    regions = ("orf1", "orf2") if a.region == "both" else (a.region,)
    out = []
    for i, (rid, _) in enumerate(recs):
        row = rows[i]
        g = row.get("genogroup", "unassigned")
        for region in regions:
            rec = dict(query=rid, genogroup=g, region=region,
                       basis=putative_novel.REGION[region]["basis"],
                       type=row.get(f"{region}_type", ""),
                       support=row.get(f"{region}_support", ""),
                       accepted=row.get(f"{region}_accepted", ""),
                       nearest=row.get(f"{region}_nearest", ""),
                       nearest_identity=row.get(f"{region}_nearest_identity", ""))
            if g in ("", "unassigned") or g not in ty.cg.NOROVIRUS_GROUPS:
                rec.update(flag=putative_novel.NA, closest_type="", distance="", bound="",
                           detail=putative_novel.NOT_NOROVIRUS if g not in ("", "unassigned")
                           else "no norovirus genogroup assigned")
            elif (i, region) not in alns:
                rec.update(flag=putative_novel.NA, closest_type="", distance="", bound="",
                           detail=f"the {putative_novel.REGION[region]['label']} region is not "
                                  f"present in this query")
            else:
                _, qa = alns[(i, region)]
                res = putative_novel.assess(region, g, qa, subs.get((i, region), ""))
                rec.update(flag=res["flag"], closest_type=res.get("closest_type", ""),
                           distance=res.get("distance", ""), bound=res.get("bound", ""),
                           detail=res["detail"])
            out.append(rec)
    df = pd.DataFrame(out)
    (df.to_csv(sys.stdout, sep="\t", index=False) if a.out == "-"
     else df.to_csv(a.out, sep="\t", index=False))
    n = int((df.flag.isin(putative_novel.FLAGGED)).sum()) if len(df) else 0
    print(f"{n} of {len(df)} region calls flagged as putative novel "
          f"(tentative: one sequence cannot designate a genotype)", file=sys.stderr)
    return 0


def _apply_window_rules(rows, alns):
    """Record each query's coverage of the RIVM typing windows.

    This used to withhold a call below 80% coverage. The truncation
    leave-one-out harness showed that gate was wrong for this tool: RIVM needs
    window coverage because RIVM types ON the window, while this tool types on
    every aligned column, so a query displaced off the window still carries the
    signal -- withheld calls were 99-100% correct when allowed through. The
    coverages are therefore reported and acceptance is decided by bootstrap
    support and aligned sites instead (see typing.apply_result).

    CALICITYPER_RIVM_WINDOW_MODE=1 restores the gate, for replicating RIVM's
    own decisions rather than maximising correct calls.
    """
    for i, row in rows.items():
        for region in ("orf1", "orf2"):
            row[f"{region}_window_cov"] = ""
            if (i, region) not in alns:
                continue
            g, qa = alns[(i, region)]
            covered, width, frac = ty.rivm_window_coverage(
                region, g, qa, genotype=row.get(f"{region}_type"))
            row[f"{region}_window_cov"] = f"{covered}/{width}"
            if frac >= cg.RIVM_WINDOW_MIN_COVERAGE or not cg.RIVM_WINDOW_MODE:
                continue
            was = row.get(f"{region}_type") or ""
            if was and was != "unassigned":
                row[f"{region}_type"] = "unassigned"
            row[f"{region}_accepted"] = (f"rejected: {covered}/{width} nt of the RIVM window "
                                         f"({frac:.0%} < {cg.RIVM_WINDOW_MIN_COVERAGE:.0%})")
            extra = (f"{region}: {covered}/{width} nt of the RIVM window covered "
                     f"({frac:.0%} < {cg.RIVM_WINDOW_MIN_COVERAGE:.0%})"
                     + (f", withheld {was}" if was and was != "unassigned" else ""))
            row["note"] = "; ".join(x for x in (row.get("note"), extra) if x)
        _apply_gii4_variant_window(row, i, alns)


def _apply_gii4_variant_window(row, i, alns):
    """Record coverage of the GII.4 variant window (P2 to the VP1 terminus).

    Region C carries almost no variant signal (4 of 135 diagnostic columns), so
    a variant called without this window rests on a handful of sites, and the
    coverage belongs beside the call. It is reported, not enforced: at an 80%
    floor the harness found 86 of 103 correct variant assignments withheld and
    no error prevented. Under RIVM window mode the suffix is dropped below the
    floor and the genotype stem returned, as RIVM does.
    """
    call = row.get("orf2_type") or ""
    row["gii4_window_cov"] = ""
    if "_" not in call or (i, "orf2") not in alns:
        return
    g, qa = alns[(i, "orf2")]
    covered, width, frac = ty.gii4_variant_window_coverage(qa, g)
    row["gii4_window_cov"] = f"{covered}/{width}"
    if frac >= cg.RIVM_WINDOW_MIN_COVERAGE:
        return
    note = (f"GII.4 variant {call}: {covered}/{width} nt of the variant window covered "
            f"({frac:.0%})")
    if cg.RIVM_WINDOW_MODE:
        row["orf2_type"] = call.split("_")[0]
        note += f", withheld -> {call.split('_')[0]}"
    row["note"] = "; ".join(x for x in (row.get("note"), note) if x)


def _screen_putative_novel(rows, alns, todo, only_unaccepted=True):
    """Add putative novel ORF1/ORF2 type flags (2x SD screen); see putative_novel.py.

    By default the screen runs only where the bootstrap did NOT accept the
    call. A well-supported type is a known type by definition, so screening it
    spends the expensive amino-acid distance pass to confirm what the tree
    already said; an unsupported call is exactly the case where "this may be
    something not in the reference set" is the useful answer.
    """
    from . import putative_novel
    subs = {(i, region): seq for (region, g), lst in todo.items() for i, seq in lst}
    for i, row in rows.items():
        for region in ("orf1", "orf2"):
            row[f"{region}_putative_novel"], row[f"{region}_putative_novel_detail"] = "", ""
            accepted = str(row.get(f"{region}_accepted", "")).startswith("accepted")
            if only_unaccepted and accepted:
                row[f"{region}_putative_novel"] = putative_novel.KNOWN
                row[f"{region}_putative_novel_detail"] = (
                    "not screened: the bootstrap accepted this type")
                continue
            if (i, region) in alns:
                g, qa = alns[(i, region)]
                res = putative_novel.assess(region, g, qa, subs.get((i, region), ""))
                row[f"{region}_putative_novel"], row[f"{region}_putative_novel_detail"] = res["flag"], res["detail"]
            elif row.get("genogroup") not in ("", "unassigned") \
                    and row["genogroup"] not in ty.cg.NOROVIRUS_GROUPS:
                row[f"{region}_putative_novel"] = putative_novel.NA
                row[f"{region}_putative_novel_detail"] = putative_novel.NOT_NOROVIRUS
            elif row.get(f"{region}_type") in ("", None) and row.get("genogroup") != "unassigned":
                row[f"{region}_putative_novel"] = putative_novel.NA
                row[f"{region}_putative_novel_detail"] = "region not present in the query"


def _support(row, region):
    try:
        return float(row.get(f"{region}_support") or 0.0)
    except (TypeError, ValueError):
        return 0.0


def _panels(a, rows, alns):
    """One panel per typed region, assigned or not: the query against its single
    closest reference (uncorrected p-distance), drawing only the differences.

    Every region that reached the alignment step gets one, including
    `unassigned` calls -- those are the ones most worth inspecting. Order is
    input order, ORF1 before ORF2. --report-max-panels caps the count if set.
    """
    from . import report
    keys = sorted(alns, key=lambda k: (k[0], k[1]))
    if a.report_max_panels is not None:
        keys = keys[:max(0, a.report_max_panels)]
    figs = []
    for i, region in keys:
        g, qa = alns[(i, region)]
        ref = align.refset(region, g)
        best = report.closest_reference(qa, ref.seqs)
        if best is None:
            continue
        k, pdist, n, ndiff = best
        row = rows[i]
        called = row.get(f"{region}_type") or "untyped"
        sup = (f" ({_support(row, region):.0f}% support)"
               if row.get(f"{region}_support") not in (None, "") else "")
        nov = row.get(f"{region}_putative_novel") or ""
        tag = (f"\u2691 {report.novel_label(region, nov).upper()} \u00b7 "
               if nov in report.NOVEL_FLAGS else "")
        title = (f"{tag}{row['query']} \u2014 {region.upper()} \u2192 {called}{sup} \u00b7 "
                 f"closest {ref.names[k]} ({ref.types[k]}): {ndiff} differences "
                 f"in {n} nt, p = {pdist:.3f}")
        figs.append(report.difference_figure(
            str(row["query"])[:38], qa, f"{ref.names[k]} {ref.types[k]}"[:38],
            ref.seqs[k], title))
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
        report.to_html(rows, a.html, figures=figs, meta=meta, notes=notes,
                       embed_js=not a.html_cdn)


def _run_recomb(a):
    """`calicityper recomb`: scan an alignment for recombination.

    Norovirus and sapovirus; --virus selects the genome model used to judge
    whether a breakpoint falls in a biologically plausible place.
    """
    from . import recomb as R

    cfg = _rc.load(a.config)
    if a.virus:
        cfg.virus = a.virus
        cfg.regions = {}            # re-resolved from `virus` by validate()
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
        sys.exit("calicityper recomb: give --alignment (aligned FASTA) or --input (unaligned FASTA)")

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
        sys.exit(f"calicityper recomb: query id(s) not in the alignment: {missing[:5]}")

    if a.dry_run:
        exe = cfg.openrdp()
        print(f"alignment        : {aln_path}")
        print(f"  sequences      : {len(aln.names)}")
        print(f"  columns        : {aln.length}")
        print(f"  masked to N    : {getattr(aln, 'masked_chars', 0)}")
        print(f"  informative    : {R.informative_sites(aln)}")
        print(f"queries          : {len(queries)}")
        print(f"methods requested: {' '.join(cfg.methods) or '(none)'}")
        print(f"OpenRDP          : {exe or 'NOT FOUND (set $CALICITYPER_OPENRDP)'}")
        print(f"virus model      : {cfg.virus}")
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
