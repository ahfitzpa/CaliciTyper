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
    p.add_argument("fasta", help="input FASTA file of contigs/genomes (required; stdin not supported)")
    p.add_argument("-o", "--out", default="-",
                   help="output TSV path, or '-' for stdout [default: -]. "
                        "When a file path is used, a provenance sidecar (<out>.provenance.json) is written.")
    p.add_argument("--bootstraps", type=int, default=100, help="bootstrap replicates [default: 100]")
    p.add_argument("--min-support", type=float, default=70.0, help="minimum bootstrap %% to accept a type [default: 70]")
    p.add_argument("--max-refs-per-type", type=int, default=2, help="reference sequences per type used in the tree [default: 2]")
    p.add_argument("--threads", type=int, default=1,
                   help="worker processes for typing/alignment chunks [default: 1]. "
                        "BLAS/OpenMP threads are pinned per process by default to avoid oversubscription.")
    p.add_argument("--seed", type=int, default=1, help="random seed for reproducibility [default: 1]")
    p.add_argument("--version", action="version", version=__version__)
    r = p.add_argument_group("reports (TSV is the primary machine-readable table)")
    r.add_argument("--report", metavar="PREFIX",
                   help="write PREFIX.md, PREFIX.json, and PREFIX.html (predictable naming) in addition to the TSV")
    r.add_argument("--json", metavar="PATH", help="write JSON report/results (includes provenance metadata)")
    r.add_argument("--md", metavar="PATH", help="write Markdown report")
    r.add_argument("--html", metavar="PATH",
                   help="write self-contained HTML report with alignment panels and trees (requires bokeh for figures)")
    r.add_argument("--report-max-panels", type=int, default=None, metavar="N",
                   help="cap alignment panels in reports (input order). Recommended for large batches "
                        "[default: one per typed region]")
    r.add_argument("--no-tree", action="store_true",
                   help="skip per-region trees in HTML reports (faster/smaller for large batches)")
    r.add_argument("--html-cdn", action="store_true",
                   help="link BokehJS from CDN instead of embedding it (smaller HTML, needs internet)")
    p.add_argument("--screen-all", action="store_true",
                   help="run putative-novel screening on every call (default: only calls not accepted by bootstrap)")
    adv = p.add_argument_group("advanced (affect correctness/reproducibility; prefer CLI flags over env vars)")
    adv.add_argument("--min-sites", type=int, default=None, metavar="N",
                     help="minimum aligned informative sites to accept a type call "
                          f"[default: {cg.MIN_INFORMATIVE_SITES} or $CALICITYPER_MIN_SITES]")
    adv.add_argument("--novel-distance", choices=["ml", "nj", "p"], default=None,
                     help="distance measure for putative-novel screening: ml, nj, or p [default: ml]")
    adv.add_argument("--type-all-groups", action=argparse.BooleanOptionalAction, default=None,
                     help="type all norovirus genogroups (GIII-GX) instead of restricting to validated GI/GII "
                          "(not validated). Overrides $CALICITYPER_TYPE_ALL_GROUPS.")
    adv.add_argument("--rivm-window-mode", action=argparse.BooleanOptionalAction, default=None,
                     help="apply RIVM-style typing-window coverage gates (for comparison/replication). "
                          "Overrides $CALICITYPER_RIVM_WINDOW_MODE.")
    p.epilog = (
        "Examples:\n"
        "  calicityper type contigs.fasta -o typed.tsv --threads 8\n"
        "  calicityper type contigs.fasta -o typed.tsv --report run01 --threads 8\n"
        "  calicityper type contigs.fasta -o typed.tsv --report run01 --report-max-panels 200 --no-tree\n\n"
        "Notes: stdin is not supported; use '-' with -o/--out to write TSV to stdout (stderr carries status).\n"
        "Validated/supported typing scope: norovirus GI/GII. Requires a verified reference set (run 'calicityper setup')."
    )
    return p


def _add_novel_args(p):
    """`novel`: the 2x SD screen as a command of its own, norovirus only."""
    p.add_argument("fasta", help="input FASTA file of contigs/genomes (required; stdin not supported)")
    p.add_argument("-o", "--out", default="-",
                   help="output TSV path, or '-' for stdout [default: -]. "
                        "When a file path is used, a provenance sidecar (<out>.provenance.json) is written.")
    p.add_argument("--region", choices=["orf1", "orf2", "both"], default="both",
                   help="region(s) to screen [default: both]. ORF2 requires a complete VP1 reading frame; "
                        "ORF1 requires >=95%% of the 762 nt 3' RdRp typing region.")
    p.add_argument("--bootstraps", type=int, default=100, help="bootstrap replicates [default: 100]")
    p.add_argument("--min-support", type=float, default=70.0, help="minimum bootstrap %% for typing context [default: 70]")
    p.add_argument("--max-refs-per-type", type=int, default=2, help="reference sequences per type [default: 2]")
    p.add_argument("--threads", type=int, default=1,
                   help="worker processes for typing/alignment chunks [default: 1]")
    p.add_argument("--seed", type=int, default=1, help="random seed for reproducibility [default: 1]")
    p.add_argument("--version", action="version", version=__version__)
    adv = p.add_argument_group("advanced (affect correctness/reproducibility; prefer CLI flags over env vars)")
    adv.add_argument("--min-sites", type=int, default=None, metavar="N",
                     help="minimum aligned informative sites to accept a type call "
                          f"[default: {cg.MIN_INFORMATIVE_SITES} or $CALICITYPER_MIN_SITES]")
    adv.add_argument("--novel-distance", choices=["ml", "nj", "p"], default=None,
                     help="distance measure for putative-novel screening: ml, nj, or p [default: ml]")
    adv.add_argument("--type-all-groups", action=argparse.BooleanOptionalAction, default=None,
                     help="type all norovirus genogroups (not validated). Overrides $CALICITYPER_TYPE_ALL_GROUPS.")
    adv.add_argument("--rivm-window-mode", action=argparse.BooleanOptionalAction, default=None,
                     help="apply RIVM-style typing-window coverage gates (comparison/replication).")
    p.epilog = (
        "Examples:\n"
        "  calicityper novel sequences.fasta -o novel.tsv\n"
        "  calicityper novel sequences.fasta --region orf2 -o novel_vp1.tsv --threads 4\n\n"
        "Notes: norovirus-only, tentative screening (not genotype designation). "
        "stdin is not supported; use '-' with -o/--out for stdout TSV."
    )
    return p


def _add_recomb_args(p):
    p.add_argument("--alignment", help="aligned FASTA containing the query and a reference panel")
    p.add_argument("--input", help="unaligned FASTA; aligned with MAFFT before scanning (written to output dir)")
    p.add_argument("--query", action="append", default=[],
                   help="query sequence id; repeat for several, or omit to analyze all sequences in the alignment")
    p.add_argument("-o", "--output", default="recomb_out", metavar="DIR",
                   help="output directory for JSON/TSV/provenance and work artifacts [default: recomb_out]")
    p.add_argument("--reference-metadata", help="TSV/CSV/JSON mapping reference ids to lineages (improves parent labels)")
    p.add_argument("--config", help="YAML or TOML configuration file")
    p.add_argument("--virus", choices=sorted(_rc.REGIONS_BY_VIRUS), default="norovirus",
                   help="genome model for breakpoint plausibility (norovirus ORF1/ORF2 junction; "
                        "sapovirus RdRp/VP1 junction within ORF1) [default: norovirus]")
    p.add_argument("--threads", type=int, default=1,
                   help="threads for MAFFT when --input is used [default: 1]")
    p.add_argument("--min-length", type=int, help="override min_sequence_length")
    p.add_argument("--pvalue-threshold", type=float, help="override pvalue_threshold")
    p.add_argument("--max-breakpoint-distance", type=int, help="override max_breakpoint_distance")
    p.add_argument("--methods", nargs="+", metavar="M",
                   help=f"subset of methods to request: {' '.join(_rc.ALL_METHODS)}")
    p.add_argument("--skip-gard", action="store_true", help="do not run HyPhy GARD")
    p.add_argument("--skip-rdp", action="store_true", help="do not run OpenRDP scanners")
    p.add_argument("--no-segment-trees", action="store_true",
                   help="skip segment phylogeny checks (faster, less evidence)")
    p.add_argument("--verbose", action="store_true", help="print progress/info to stderr")
    p.add_argument("--dry-run", action="store_true",
                   help="validate inputs and report what would run (including tool availability), then stop")
    p.epilog = (
        "Examples:\n"
        "  calicityper recomb --alignment panel.fasta --query SAMPLE_01 -o recomb_out/\n"
        "  calicityper recomb --input samples.fasta -o recomb_out/ --threads 8\n"
        "  calicityper recomb --alignment panel.fasta --query SAMPLE_01 --dry-run\n\n"
        "Outputs: recombination.json, recombination.tsv, recombination.tsv.provenance.json, and work/ artifacts."
    )
    return p


def _add_setup_args(p):
    p.add_argument("--source", help="URL or local path to a manifest-matching reference-set tarball "
                                    "[default: download_url from the shipped manifest]")
    p.add_argument("--data-dir", help="install/verify location [default: versioned cache dir, or package dir with --in-package]")
    p.add_argument("--in-package", action="store_true",
                   help="install into the package's own data directory instead of the cache")
    p.add_argument("--manifest", help="alternative reference_manifest.json")
    p.add_argument("--verify-only", action="store_true",
                   help="verify an existing data directory against the manifest and exit")
    p.add_argument("--force", action="store_true", help="replace/reinstall even if data appears present")
    p.epilog = (
        "Examples:\n"
        "  calicityper setup\n"
        "  calicityper setup --verify-only\n"
        "  calicityper setup --source /shared/ref.tar.gz --data-dir /shared/calicityper-ref\n\n"
        "Offline/cluster: install once from a local manifest-matching archive, verify, then set CALICITYPER_DATA."
    )
    return p


def main(argv=None):
    top = argparse.ArgumentParser(
        prog="calicityper",
        description="Caliciviridae genogroup detection and norovirus typing/recombination analysis. "
                    "Validated/supported genotyping scope: norovirus GI/GII. Other norovirus "
                    "genogroups and non-norovirus Caliciviridae may be detected but genotyping "
                    "is not validated; sapovirus genotyping is not implemented.")
    top.add_argument("--version", action="version", version=__version__)
    sub = top.add_subparsers(dest="command", required=True, metavar="{type,novel,recomb,setup}")
    _add_type_args(sub.add_parser(
        "type",
        help="assign genogroup and norovirus GI/GII ORF1/ORF2 types (validated scope)",
        description="Norovirus GI/GII dual genotyping is the validated/supported typing scope. "
                    "Other norovirus genogroups and non-norovirus Caliciviridae may be detected "
                    "to genogroup, but genotyping is not validated for them. Sapovirus genotyping "
                    "is not implemented. Putative-novel screening is tentative screening, not a "
                    "genotype designation.",
        formatter_class=argparse.RawDescriptionHelpFormatter))
    _add_novel_args(sub.add_parser(
        "novel",
        help="screen for putative novel norovirus ORF1/ORF2 types (tentative, norovirus only)",
        description="Chhabra et al. 2019 2x SD screen applied to one query at a time "
                    "(tentative screening, not genotype designation). Norovirus only.",
        formatter_class=argparse.RawDescriptionHelpFormatter))
    _add_recomb_args(sub.add_parser(
        "recomb",
        help="detect and score recombination (norovirus or sapovirus)",
        description="Evidence-weighted recombination analysis for norovirus/sapovirus alignments. "
                    "Use --dry-run to validate inputs and tool availability.",
        formatter_class=argparse.RawDescriptionHelpFormatter))
    _add_setup_args(sub.add_parser(
        "setup",
        help="fetch and verify the reference data set",
        description="Install or verify the reference data set required for typing. "
                    "Refuses to install data that does not match the shipped manifest.",
        formatter_class=argparse.RawDescriptionHelpFormatter))

    a = top.parse_args(argv)
    try:
        if a.command == "recomb":
            return _run_recomb(a)
        if a.command == "setup":
            from . import setup_data
            return setup_data.main(a)
        if a.command == "novel":
            return _run_novel(a)
        return _run_type(a)
    except FileNotFoundError as e:
        msg = str(e)
        if "reference data not found" in msg or "Run `calicityper setup`" in msg:
            print(f"Error: {msg}", file=sys.stderr)
            print("Hint: run 'calicityper setup' (or provide a verified local archive) and/or set CALICITYPER_DATA.", file=sys.stderr)
            return 2
        print(f"Error: {msg}", file=sys.stderr)
        if any(x in msg for x in ("blastn", "mafft", "clustalw", "iqtree")) or "not found on PATH" in msg:
            print("Hint: install required external tools (e.g. conda install -c bioconda blast mafft) or set CALICITYPER_* paths.", file=sys.stderr)
        elif getattr(a, "fasta", None):
            print("Hint: check that the input FASTA path exists and is readable.", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        return 130


def _align_and_type(a):
    """BLAST, plan, align and type every record. Returns (recs, rows, alns, todo).

    Shared by `type` and `novel`: the putative novel screen needs exactly the
    same alignment to the genogroup reference set that typing produces, so it
    runs off one pipeline rather than a parallel one that could drift from it.
    """
    _apply_runtime_overrides(a)
    a._want_aln = True                  # the putative novel type screen needs the alignment, not only the HTML
    cmd = getattr(a, "command", "type")
    print(f"[{cmd}] checking reference data and input...", file=sys.stderr)
    ty.cg.check_data()          # clear error if the bundled reference data is missing
    recs = [(r.id, str(r.seq)) for r in SeqIO.parse(a.fasta, "fasta")]
    if not recs:
        print(f"[{cmd}] no sequences found in {a.fasta}", file=sys.stderr)
    else:
        print(f"[{cmd}] loaded {len(recs)} sequence(s) from {a.fasta}", file=sys.stderr)
    with tempfile.TemporaryDirectory(prefix="calicityper_") as td:
        q = os.path.join(td, "all.fa")
        SeqIO.write([r for r in SeqIO.parse(a.fasta, "fasta")], q, "fasta")
        print(f"[{cmd}] running genogroup detection (BLAST)...", file=sys.stderr)
        hits = ty.blast_genogroup(q, a.threads)
    rows, todo = {}, defaultdict(list)             # todo[(region, group)] -> [(index, seq)]
    for i, (rid, s) in enumerate(recs):
        row, jobs = ty.plan_query(rid, s, hits.get(rid))
        rows[i] = row
        for region, g, seq in jobs:
            todo[(region, g)].append((i, seq))
    n_jobs = sum(len(lst) for lst in todo.values())
    print(f"[{cmd}] planned {n_jobs} region typing job(s) across {len(todo)} region/genogroup group(s)", file=sys.stderr)
    tree.warmup()          # primes the on-disk kernel cache the workers load from
    chunks = []                                     # one aligner call per chunk of queries from one region/genogroup
    for (region, g), lst in todo.items():
        size = max(8, -(-len(lst) // (2 * a.threads)))
        chunks += [(region, g, lst[k:k + size]) for k in range(0, len(lst), size)]
    if chunks:
        mode = f"{a.threads} worker process(es)" if a.threads > 1 and len(chunks) > 1 else "1 process"
        print(f"[{cmd}] aligning/typing {len(chunks)} chunk(s) with {mode}...", file=sys.stderr)
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
    print(f"[{cmd}] completed typing for {len(recs)} sequence(s)", file=sys.stderr)
    return recs, rows, alns, todo


def _run_type(a):
    if a.report:
        a.md = a.md or f"{a.report}.md"
        a.json = a.json or f"{a.report}.json"
        a.html = a.html or f"{a.report}.html"
    recs, rows, alns, todo = _align_and_type(a)
    print("[type] running putative-novel screening...", file=sys.stderr)
    _screen_putative_novel(rows, alns, todo, only_unaccepted=not a.screen_all)
    ordered = [rows[i] for i in range(len(recs))]
    df = pd.DataFrame(ordered)
    print("Scope: norovirus GI/GII genotyping validated/supported; other norovirus genogroups and "
          "non-norovirus Caliciviridae may be detected but genotyping is not validated. "
          "Sapovirus genotyping not implemented. Putative-novel flags are tentative screening.",
          file=sys.stderr)
    if a.out == "-":
        df.to_csv(sys.stdout, sep="\t", index=False)
    else:
        df.to_csv(a.out, sep="\t", index=False)
        from . import report
        prov_path = _provenance_sidecar_path(a.out)
        report.write_type_provenance(prov_path, len(ordered), meta=_type_meta(a, command="type"))
        print(f"[type] wrote {a.out} and {prov_path}", file=sys.stderr)
    if a.md or a.json or a.html:
        _write_reports(a, ordered, alns)
        written = [p for p in (a.md, a.json, a.html) if p]
        if written:
            print(f"[type] wrote reports: {', '.join(written)}", file=sys.stderr)

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
    from . import report
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
                       distance=putative_novel.effective_measure(),
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
    print("[novel] assessing putative-novel status for screened region(s)...", file=sys.stderr)
    df = pd.DataFrame(out)
    if a.out == "-":
        df.to_csv(sys.stdout, sep="\t", index=False)
    else:
        df.to_csv(a.out, sep="\t", index=False)
        meta = _type_meta(a, command="novel")
        meta.update(region=a.region, screen_regions=list(regions))
        prov_path = _provenance_sidecar_path(a.out)
        report.write_novel_provenance(prov_path, len(out), meta=meta)
        print(f"[novel] wrote {a.out} and {prov_path}", file=sys.stderr)
    n = int((df.flag.isin(putative_novel.FLAGGED)).sum()) if len(df) else 0
    print(f"[novel] {n} of {len(df)} region call(s) flagged as putative novel "
          f"(tentative screening; not genotype designation)", file=sys.stderr)
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


def _type_meta(a, command="type"):
    from . import putative_novel
    from . import provenance
    return dict(bootstraps=a.bootstraps, min_support=a.min_support,
                max_refs_per_type=a.max_refs_per_type, seed=a.seed, input=a.fasta,
                min_sites=cg.MIN_INFORMATIVE_SITES,
                typeable_groups=list(cg.TYPEABLE_GROUPS),
                rivm_window_mode=bool(cg.RIVM_WINDOW_MODE),
                novel_distance=putative_novel.effective_measure(),
                threads=getattr(a, "threads", None),
                command=command,
                reference_set=provenance.reference_set_info())


def _provenance_sidecar_path(out_path):
    p = Path(out_path)
    return str(p.with_suffix(p.suffix + ".provenance.json"))


def _write_reports(a, rows, alns):
    from . import report
    meta = _type_meta(a, command="type")
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
        sys.exit("calicityper recomb: provide --alignment <aligned.fasta> or --input <unaligned.fasta>")

    outdir = Path(a.output)
    outdir.mkdir(parents=True, exist_ok=True)

    aln_path = a.alignment
    if not aln_path:
        # --input is unaligned: align it first, because every scanner needs a
        # column-consistent alignment and silently mis-scanning unaligned input
        # would be worse than spending the MAFFT call.
        aln_path = str(outdir / "aligned.fasta")
        print(f"[recomb] aligning input FASTA with MAFFT: {a.input} -> {aln_path}", file=sys.stderr)
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
    R.write_provenance_sidecar(provenance, outdir / "recombination.tsv.provenance.json",
                               n_results=len(results), n_events=sum(len(r["recombination"].get("events", [])) for r in results))

    for t in provenance.get("tools_skipped", []):
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
