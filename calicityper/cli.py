import argparse, csv, os, subprocess, sys, tempfile
from pathlib import Path
from collections import defaultdict
from . import typing as ty
from . import config as cg, align, tree, __version__
from . import recomb_config as _rc

_ARGS = None
def _work(chunk):
    """chunk = (region, group, [(row_index, seq), ...]): ONE MAFFT call for the chunk, then bootstrap each query."""
    region, g, lst = chunk
    a = _ARGS
    ref, aln = align.align_queries(region, g, [(str(i), s) for i, s in lst], threads=1)
    fast_boot = bool(getattr(a, "fast_boot", False))
    out = []
    for (i, _), qa in zip(lst, aln):
        res = ty.type_aligned(ref, qa, nrep=a.bootstraps, seed=a.seed, min_support=a.min_support,
                              max_refs_per_type=a.max_refs_per_type, fast_boot=fast_boot)
        # keep the aligned query only when a report will use it: for thousands of
        # queries these strings are the bulk of what crosses back from a worker.
        out.append(((i, region), res, (g, qa) if a._want_aln else None))
    return out

def _apply_runtime_overrides(a):
    """Apply the advanced command-line settings over their $CALICITYPER_* defaults.

    A flag left unset (None) keeps the environment/default value. Runs in the
    parent and again in every worker (_init): under the forkserver start method
    workers do not inherit module state changed at run time, so setting it in
    the parent alone would only take effect with --threads 1.
    """
    preset = getattr(a, "preset", None)
    if preset is not None:
        cfg = cg.PRESETS.get(preset, {})
        cg.MIN_INFORMATIVE_SITES_ORF1 = cfg.get("min_sites_orf1")
        cg.MIN_INFORMATIVE_SITES_ORF2 = cfg.get("min_sites_orf2")
    if getattr(a, "min_sites", None) is not None:
        cg.MIN_INFORMATIVE_SITES = a.min_sites
    if getattr(a, "min_sites_orf1", None) is not None:
        cg.MIN_INFORMATIVE_SITES_ORF1 = a.min_sites_orf1
    if getattr(a, "min_sites_orf2", None) is not None:
        cg.MIN_INFORMATIVE_SITES_ORF2 = a.min_sites_orf2
    if getattr(a, "rivm_window_mode", None) is not None:
        cg.RIVM_WINDOW_MODE = a.rivm_window_mode
    if getattr(a, "type_all_groups", None) is not None:
        cg.TYPEABLE_GROUPS = list(cg.NOROVIRUS_GROUPS if a.type_all_groups
                                  else cg.DEFAULT_TYPEABLE_GROUPS)
    if getattr(a, "type_sapovirus", None) is not None:
        cg.TYPE_SAPOVIRUS = a.type_sapovirus
    if getattr(a, "novel_distance", None) is not None:
        from . import putative_novel
        putative_novel.DISTANCE = a.novel_distance
    if getattr(a, "fast_novel", False):
        from . import putative_novel
        putative_novel.FAST_NOVEL = True

def _init(a):
    global _ARGS
    _ARGS = a
    _apply_runtime_overrides(a)     # workers don't inherit the parent's overrides
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
    p.add_argument("--fast-boot", action="store_true",
                   help="opt-in sequential binomial early stopping at 25/100 replicates when support "
                        "is 25/25 (100%%) or <=8/25 (<=32%%); borderline calls still run all replicates")
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
    adv.add_argument("--preset", choices=sorted(cg.PRESETS), default=None,
                     help="amplicon threshold preset: 'polcap' sets --min-sites-orf1 170 --min-sites-orf2 250 "
                          "for 570/579 bp dual PolCap amplicons")
    adv.add_argument("--min-sites", type=int, default=None, metavar="N",
                     help="minimum aligned informative sites to accept a type call "
                          f"[default: {cg.MIN_INFORMATIVE_SITES} or $CALICITYPER_MIN_SITES]")
    adv.add_argument("--min-sites-orf1", type=int, default=None, metavar="N",
                     help="minimum aligned informative sites for ORF1/RdRp only (overrides --min-sites / --preset)")
    adv.add_argument("--min-sites-orf2", type=int, default=None, metavar="N",
                     help="minimum aligned informative sites for ORF2/VP1 only (overrides --min-sites / --preset)")
    adv.add_argument("--novel-distance", choices=["ml", "nj", "p"], default=None,
                     help="distance measure for putative-novel screening: ml, nj, or p [default: ml]")
    adv.add_argument("--fast-novel", action="store_true",
                     help="use cached ML reference backbone constraint trees (IQ-TREE -g) during putative-novel screening")
    adv.add_argument("--type-all-groups", action=argparse.BooleanOptionalAction, default=None,
                     help="type all norovirus genogroups (GIII-GX) instead of restricting to validated GI/GII "
                          "(not validated). Overrides $CALICITYPER_TYPE_ALL_GROUPS.")
    adv.add_argument("--type-sapovirus", action=argparse.BooleanOptionalAction, default=None,
                     help="sapovirus VP1/RdRp dual typing by the Zhao et al. 2024 nucleotide "
                          "distance cut-offs [default: on]. --no-type-sapovirus reports genogroup "
                          "only. Overrides $CALICITYPER_TYPE_SAPOVIRUS.")
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
    p.add_argument("--fast-boot", action="store_true",
                   help="opt-in sequential binomial early stopping at 25/100 replicates for typing context")
    p.add_argument("--min-support", type=float, default=70.0, help="minimum bootstrap %% for typing context [default: 70]")
    p.add_argument("--max-refs-per-type", type=int, default=2, help="reference sequences per type [default: 2]")
    p.add_argument("--threads", type=int, default=1,
                   help="worker processes for typing/alignment chunks [default: 1]")
    p.add_argument("--seed", type=int, default=1, help="random seed for reproducibility [default: 1]")
    p.add_argument("--version", action="version", version=__version__)
    adv = p.add_argument_group("advanced (affect correctness/reproducibility; prefer CLI flags over env vars)")
    adv.add_argument("--preset", choices=sorted(cg.PRESETS), default=None,
                     help="amplicon threshold preset: 'polcap' sets --min-sites-orf1 170 --min-sites-orf2 250")
    adv.add_argument("--min-sites", type=int, default=None, metavar="N",
                     help="minimum aligned informative sites to accept a type call "
                          f"[default: {cg.MIN_INFORMATIVE_SITES} or $CALICITYPER_MIN_SITES]")
    adv.add_argument("--min-sites-orf1", type=int, default=None, metavar="N",
                     help="minimum aligned informative sites for ORF1/RdRp only")
    adv.add_argument("--min-sites-orf2", type=int, default=None, metavar="N",
                     help="minimum aligned informative sites for ORF2/VP1 only")
    adv.add_argument("--novel-distance", choices=["ml", "nj", "p"], default=None,
                     help="distance measure for putative-novel screening: ml, nj, or p [default: ml]")
    adv.add_argument("--fast-novel", action="store_true",
                     help="use cached ML reference backbone constraint trees (IQ-TREE -g) during screening")
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
    p.add_argument("--openrdp-timeout", type=int, default=7200, metavar="S",
                   help="seconds to allow the combined OpenRDP run; on failure the methods "
                        "are retried one at a time, each with a share of this budget "
                        "[default: 7200]")
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


def _antigenic_provenance(a, reference, n_seqs, entropy, path):
    """Record what was run, against what, and under which published rules."""
    import json as _json
    from . import antigenic as ag
    conserved = sorted(r for r, h in entropy.items()
                       if h is not None and h <= ag.CONSERVED_ENTROPY)
    rec = {
        "tool": "calicityper", "tool_version": __version__,
        "subcommand": "antigenic", "alignment": a.alignment,
        "reference": reference, "genotype": a.genotype,
        "n_sequences": n_seqs,
        "site_definitions": "Parra et al. 2026 J Virol doi:10.1128/jvi.01804-25 (CC BY 4.0)",
        "shift_rule": {
            "min_aa_differences": ag.MIN_AA_CHANGES,
            "min_sites_changed": 2, "sites": list(ag.SHIFT_SITES),
            "site_C_excluded": True,
            "source": "Kendra et al. 2021 PNAS doi:10.1073/pnas.2015874118"},
        "entropy": {"conserved_at_or_below": ag.CONSERVED_ENTROPY,
                    "source": "Kendra et al. 2021 (1,572 GII.4 VP1 sequences)",
                    "measured_over": "the supplied alignment, not that corpus",
                    "n_site_residues_conserved": len(conserved),
                    "conserved_residues": conserved},
        "caveat": ag.CAVEAT,
    }
    with open(path, "w") as fh:
        _json.dump(rec, fh, indent=1)


def _add_antigenic_args(p):
    p.add_argument("alignment", nargs="?",
                   help="aligned VP1 amino acids (FASTA), or unaligned ORF2 NT / VP1 AA when --unaligned is passed")
    p.add_argument("-o", "--out", default="antigenic.tsv", help="output TSV")
    p.add_argument("--genotype", default="GII.4",
                   help="genotype the site coordinates belong to (default GII.4)")
    p.add_argument("--unaligned", action="store_true",
                   help="translate (if ORF2 NT) and align input queries against bundled GII.4 variant references "
                        "(using packaged data/gii4_variant_references.tsv + GII.orf and MAFFT)")
    p.add_argument("--threads", type=int, default=1,
                   help="MAFFT threads when --unaligned is used [default: 1]")
    p.add_argument("--reference", metavar="NAME",
                   help="sequence id to compare against; default is the first record (or nearest GII.4 variant under --unaligned)")
    p.add_argument("--anchor", metavar="NAME",
                   help="optional 540-aa GII.4 coordinate-anchor id in the alignment "
                        "(e.g. JX459908.1|Sydney) when --reference is a 539-aa pre-2002 variant")
    p.add_argument("--skip-input-check", action="store_true",
                   help="annotate even when a query is not a norovirus GII VP1. The site "
                        "coordinates are GII VP1 numbering, so the residues reported for "
                        "another genogroup or genus are not the published sites; for "
                        "deliberate cross-genotype work and synthetic input only")
    p.add_argument("--show-sites", action="store_true",
                   help="print the residue definitions and exit")


def _run_antigenic(a):
    import csv as _csv
    from . import antigenic as ag

    if not a.show_sites and not a.alignment:
        print("Error: an alignment is required unless --show-sites is given.",
              file=sys.stderr)
        return 2
    if a.show_sites:
        for r in ag.sites_table():
            print(f"{r['site']:9s} {r['kind']:13s} {r['residues']}")
        return 0
    if a.genotype.split("_")[0] != "GII.4":
        print(f"Error: the site coordinates are GII.4 VP1 numbering. VP1 length varies "
              f"between 530 and 550 residues across human genotypes, so they cannot be "
              f"applied to {a.genotype} without structural mapping.", file=sys.stderr)
        return 2

    recs = ty.read_fasta(a.alignment)
    all_res = ag.all_tracked_residues()
    rows = []

    if getattr(a, "unaligned", False):
        if not recs:
            print("Error: no sequences found in input FASTA.", file=sys.stderr)
            return 2
        try:
            aln_refs, aln_qrys, anchor = ag.align_unaligned_queries(recs, threads=a.threads)
            if getattr(a, "anchor", None):
                ahit = [r for r in aln_refs if r[0] == a.anchor or r[0].split("|")[0] == a.anchor]
                if not ahit:
                    print(f"Error: coordinate anchor {a.anchor!r} is not in bundled GII.4 references.",
                          file=sys.stderr)
                    return 2
                anchor = ahit[0]
            cols = ag.residue_columns(anchor[1], all_res)
            pool_seqs = [s for _, s in aln_refs] + [s for _, s in aln_qrys]
            ent = ag.site_entropy(pool_seqs, cols)
            fixed_ref = None
            if a.reference:
                rhit = [r for r in aln_refs if r[0] == a.reference or r[0].split("|")[0] == a.reference]
                if not rhit:
                    print(f"Error: reference {a.reference!r} is not in bundled GII.4 references.",
                          file=sys.stderr)
                    return 2
                fixed_ref = rhit[0]

            for qname, qseq in aln_qrys:
                # norovirus GII VP1 only: see antigenic.check_gii_input
                if not getattr(a, "skip_input_check", False):
                    ag.check_gii_input(qname, qseq.replace("-", "").replace(".", ""))
                chosen_ref = fixed_ref or min(
                    aln_refs, key=lambda r: ag.total_aa_differences(qseq, r[1])
                )
                rows.append(ag.to_row(ag.annotate(
                    qname, qseq, chosen_ref[0], chosen_ref[1],
                    a.genotype, entropy=ent,
                    anchor_row=anchor[1], anchor_name=anchor[0],
                )))
            prov_ref_label = fixed_ref[0] if fixed_ref else f"nearest_of_{len(aln_refs)}_GII.4_variants(anchor={anchor[0]})"
            n_total_aln = len(pool_seqs)
        except ag.AntigenicError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 2
    else:
        if len(recs) < 2:
            print(f"Error: need at least two aligned sequences, found {len(recs)}.",
                  file=sys.stderr)
            return 2
        if a.reference:
            hit = [r for r in recs if r[0] == a.reference]
            if not hit:
                print(f"Error: reference {a.reference!r} is not in the alignment.",
                      file=sys.stderr)
                return 2
            ref = hit[0]
        else:
            ref = recs[0]

        if getattr(a, "anchor", None):
            ahit = [r for r in recs if r[0] == a.anchor]
            if not ahit:
                print(f"Error: coordinate anchor {a.anchor!r} is not in the alignment.",
                      file=sys.stderr)
                return 2
            anchor = ahit[0]
        else:
            anchor = ref

        try:
            cols = ag.residue_columns(anchor[1], all_res)
            ent = ag.site_entropy([s for _, s in recs], cols)
            for name, seq in recs:
                if name == ref[0]:
                    continue
                # norovirus GII VP1 only: see antigenic.check_gii_input
                if not getattr(a, "skip_input_check", False):
                    ag.check_gii_input(name, seq.replace("-", "").replace(".", ""))
                rows.append(ag.to_row(ag.annotate(name, seq, ref[0], ref[1],
                                                  a.genotype, entropy=ent,
                                                  anchor_row=anchor[1],
                                                  anchor_name=anchor[0])))
            prov_ref_label = ref[0]
            n_total_aln = len(recs)
        except ag.AntigenicError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 2

    with open(a.out, "w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=ag.TSV_COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    _antigenic_provenance(a, prov_ref_label, n_total_aln, ent, a.out + ".provenance.json")
    flagged = sum(r["candidate_antigenic_shift"] == "true" for r in rows)
    print(f"annotated {len(rows)} sequence(s) against {prov_ref_label} -> {a.out}")
    print(f"  {flagged} flagged candidate_antigenic_shift")
    print(f"  {ag.CAVEAT}")
    return 0


def main(argv=None):
    top = argparse.ArgumentParser(
        prog="calicityper",
        description="Caliciviridae genogroup detection and norovirus typing/recombination analysis. "
                    "Validated/supported genotyping scope: norovirus GI/GII. Other norovirus "
                    "genogroups and non-norovirus Caliciviridae may be detected but genotyping "
                    "is not validated; sapovirus genotyping (--type-sapovirus) is experimental and "
                    "not yet benchmarked.")
    top.add_argument("--version", action="version", version=__version__)
    sub = top.add_subparsers(dest="command", required=True, metavar="{type,novel,recomb,antigenic,setup}")
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
    _add_antigenic_args(sub.add_parser(
        "antigenic",
        help="annotate GII.4 antigenic sites on an aligned VP1 amino-acid alignment",
        description="Annotation only. Reports which published GII.4 antigenic-site "
                    "residues differ from a reference, and screens for the change "
                    "combination associated with antigenic shift. Does not type, does "
                    "not change any call, and cannot measure antigenic distance, which "
                    "needs blockade or neutralisation titres.",
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
        if a.command == "antigenic":
            return _run_antigenic(a)
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


def _write_tsv(rows, out_path):
    """Write a list of row dicts to TSV using stdlib csv (matches pandas.to_csv output)."""
    fields = []
    seen = set()
    for row in rows:
        for k in row.keys():
            if k not in seen:
                seen.add(k)
                fields.append(k)

    def _emit(fh):
        if not fields:
            fh.write("\n")
            return
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore", restval="", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: ("" if v is None else v) for k, v in row.items()})

    if out_path == "-":
        _emit(sys.stdout)
    else:
        with open(out_path, "w", newline="", encoding="utf-8") as fh:
            _emit(fh)


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
    recs = ty.read_fasta(a.fasta)
    if not recs:
        print(f"[{cmd}] no sequences found in {a.fasta}", file=sys.stderr)
    else:
        print(f"[{cmd}] loaded {len(recs)} sequence(s) from {a.fasta}", file=sys.stderr)
    blast_recs = [
        (rid, s) for rid, s in recs
        if not ty.invalid_characters(s)
        and len(str(s).replace("-", "").replace("N", "")) >= cg.MIN_INPUT_NT
    ]
    if blast_recs:
        with tempfile.TemporaryDirectory(prefix="calicityper_") as td:
            q = os.path.join(td, "all.fa")
            with open(q, "w", encoding="utf-8") as fh:
                for rid, s in blast_recs:
                    fh.write(f">{rid}\n{s}\n")
            print(f"[{cmd}] running genogroup detection (BLAST)...", file=sys.stderr)
            hits = ty.blast_genogroup(q, a.threads)
    else:
        hits = {}
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
    use_pool = a.threads > 1 and len(chunks) > 1 and n_jobs > 4
    if chunks:
        mode = f"{min(a.threads, len(chunks))} worker process(es)" if use_pool else "1 process"
        print(f"[{cmd}] aligning/typing {len(chunks)} chunk(s) with {mode}...", file=sys.stderr)
    if use_pool:
        from multiprocessing import Pool
        with Pool(min(a.threads, len(chunks)), initializer=_init, initargs=(a,)) as pool:
            results = [r for part in pool.map(_work, chunks, chunksize=1) for r in part]
    else:
        _init(a); results = [r for c in chunks for r in _work(c)]
    alns = {}                                       # (row index, region) -> (genogroup, aligned query)
    for (i, region), res, aln in sorted(results, key=lambda r: (r[0][0], r[0][1])):
        ty.apply_result(rows[i], region, res, a.min_support)
        if aln is not None:
            alns[(i, region)] = aln
    _apply_window_rules(rows, alns)
    _apply_dual_type(rows)
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
    print("Scope: norovirus GI/GII genotyping validated/supported; other norovirus genogroups and "
          "non-norovirus Caliciviridae may be detected but genotyping is not validated. "
          + ("Sapovirus genotyping ON: VP1 genotype and RdRp P-type, benchmarked; "
             "RdRp cannot flag a novel P-type, so P-types are provisional. "
             if cg.TYPE_SAPOVIRUS
             else "Sapovirus: genogroup only (CALICITYPER_TYPE_SAPOVIRUS=1 enables typing). ")
          + "Putative-novel flags are tentative screening.",
          file=sys.stderr)
    _write_tsv(ordered, a.out)
    if a.out != "-":
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
            if ty.cg.is_sapovirus(g) and ty.cg.TYPE_SAPOVIRUS:
                # sapovirus VP1 is screened by the published cut-offs, which the
                # typing step has already applied; RdRp is not screened at all
                if region == "orf2":
                    flag, detail = putative_novel.sapovirus_vp1_flag(row)
                    rec.update(flag=flag, closest_type=row.get("orf2_nearest_type", ""),
                               distance=row.get("orf2_distance", ""),
                               bound=ty.cg.SAPOVIRUS_CUTOFFS["orf2"][0], detail=detail)
                else:
                    rec.update(flag=putative_novel.NA, closest_type="", distance="", bound="",
                               detail=putative_novel.SAPOVIRUS_RDRP_NA)
            elif g in ("", "unassigned") or g not in ty.cg.NOROVIRUS_GROUPS:
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
    _write_tsv(out, a.out)
    if a.out != "-":
        meta = _type_meta(a, command="novel")
        meta.update(region=a.region, screen_regions=list(regions))
        prov_path = _provenance_sidecar_path(a.out)
        report.write_novel_provenance(prov_path, len(out), meta=meta)
        print(f"[novel] wrote {a.out} and {prov_path}", file=sys.stderr)
    n = sum(r.get("flag") in putative_novel.FLAGGED for r in out)
    print(f"[novel] {n} of {len(out)} region call(s) flagged as putative novel "
          f"(tentative screening; not genotype designation)", file=sys.stderr)
    return 0


def _apply_dual_type(rows):
    """Compose the dual-typing designation and flag genogroup incongruence.

    Discordance is reported as discordance, not as recombination: a capsid and
    a polymerase placed in different groups may reflect a recombination event,
    but also differential evolutionary rates or simple phylogenetic
    uncertainty. Zhao et al. make the same distinction for the sea lion strain,
    whose VP1/RdRp incongruence dedicated software found NOT to be
    recombination. `calicityper recomb` is what tests the recombination
    hypothesis; this column only says the two regions disagree.
    """
    for row in rows.values():
        gt, pt = row.get("orf2_type", ""), row.get("orf1_type", "")
        row["dual_type"] = ty.dual_type(gt, pt)
        row["dual_type_discordant"] = ""
        if not row["dual_type"]:
            continue
        gg = str(gt).split(".")[0]
        pg = str(pt).partition(".P")[0]
        if pg and gg and pg != gg:
            row["dual_type_discordant"] = (
                f"capsid {gg} with polymerase {pg}: regions placed in different groups; "
                f"potential recombination or evolutionary incongruence, test with `calicityper recomb`")


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
            if cg.is_sapovirus(g):
                continue                  # RIVM windows are norovirus typing windows
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
            if ty.cg.is_sapovirus(row.get("genogroup", "")):
                # VP1: report the cut-off outcome the type call already reached,
                # rather than leaving the column blank. RdRp: not screened.
                if region == "orf2" and ty.cg.TYPE_SAPOVIRUS:
                    row[f"{region}_putative_novel"], row[f"{region}_putative_novel_detail"] = \
                        putative_novel.sapovirus_vp1_flag(row)
                else:
                    row[f"{region}_putative_novel"] = putative_novel.NA
                    row[f"{region}_putative_novel_detail"] = (
                        putative_novel.SAPOVIRUS_RDRP_NA if region == "orf1" else
                        "sapovirus dual typing is off (--no-type-sapovirus): genogroup only")
                continue
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
    return dict(bootstraps=a.bootstraps, fast_boot=bool(getattr(a, "fast_boot", False)),
                min_support=a.min_support,
                max_refs_per_type=a.max_refs_per_type, seed=a.seed, input=a.fasta,
                preset=getattr(a, "preset", None),
                min_sites=cg.MIN_INFORMATIVE_SITES,
                min_sites_orf1=cg.min_sites_for("orf1"),
                min_sites_orf2=cg.min_sites_for("orf2"),
                typeable_groups=list(cg.TYPEABLE_GROUPS),
                rivm_window_mode=bool(cg.RIVM_WINDOW_MODE),
                type_sapovirus=bool(cg.TYPE_SAPOVIRUS),
                novel_distance=putative_novel.effective_measure(),
                fast_novel=bool(putative_novel.FAST_NOVEL),
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
                                    workdir=str(outdir / "work"),
                                    openrdp_timeout=a.openrdp_timeout)
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
    # propagate the exit code: main() returns 2 on a usage error and
    # discarding it made `python -m calicityper.cli` always exit 0,
    # which hid every failure from any script that checked.
    raise SystemExit(main())
