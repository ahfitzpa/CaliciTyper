#!/usr/bin/env python3
"""Figures for the CaliciTyper benchmarking manuscript.

Reads only benchmark outputs; writes PNG + PDF into figures/. Run through
calicityper_runs/scripts/manuscript_figures.sh (SLURM), not on the login node.

Figure design follows two models:
  * the tool-by-target dot heatmap and the threshold panel with a table
    beneath, after Figs 3 and 4 of the HLA genotyping benchmark (BMC Genomics
    2023, PMC10170851);
  * the pairwise-distance frequency distributions with cluster cut-offs
    marked, after the 2 x SD framework of Chhabra et al. 2019, which
    Zhao et al. 2024 reuse for sapovirus.

    python make_figures.py --out-dir ../figures
"""
from __future__ import annotations

import argparse
import itertools
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from Bio import AlignIO, SeqIO

HOME = Path(os.path.expanduser("~/scratch"))
REPO = HOME / "calicityper"
RUNS = HOME / "calicityper_runs"
BENCH070 = RUNS / "benchmarking_wip_0.7.0" / "results"
SAPO = REPO / "results" / "benchmark_sapovirus"
DEV = RUNS / "rivm_genogroup_dev"
sys.path.insert(0, str(REPO))

# ggsci journal palettes, read out of ggsci 5.2.0 (ggsci_db) rather than retyped.
PALETTES = {
    "bmj": ["#2A6EBB", "#F0AB00", "#C50084", "#7D5CC6", "#E37222",
            "#69BE28", "#00B2A9", "#CD202C", "#747678"],
    "jama": ["#374E55", "#DF8F44", "#00A1D5", "#B24745", "#79AF97", "#6A6599", "#80796B"],
    "nejm": ["#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91"],
}
PALETTE = os.environ.get("PALETTE", "bmj").lower()
if PALETTE not in PALETTES:
    raise SystemExit(f"PALETTE must be one of {sorted(PALETTES)}")
P = PALETTES[PALETTE]
# Categorical series take the palette IN ORDER, as scale_colour_bmj() does:
# picking favourites out of the middle of a journal palette produces a figure
# that is not recognisably that palette.
C_TOOL, C_RIVM, C_ZHAO = P[0], P[1], P[2]
# Outcome colours are also palette members, chosen for their meaning: the
# BMJ palette ends in a red and a grey, which suit "discordant" and
# "untypeable". Other palettes have no grey, so one is taken from their end.
if PALETTE == "bmj":
    C_OK, C_WRONG, C_UNTYPED = P[0], P[7], P[8]        # blue, red, grey
    C_MID, C_INK, C_NEUTRAL = P[1], "#2B2B2B", "#D9DBDC"
elif PALETTE == "jama":
    C_OK, C_WRONG, C_UNTYPED = P[0], P[3], P[6]
    C_MID, C_INK, C_NEUTRAL = P[1], P[0], "#D9DBDC"
else:
    C_OK, C_WRONG, C_UNTYPED = P[1], P[0], "#B4B8BB"
    C_MID, C_INK, C_NEUTRAL = P[2], "#2B2B2B", "#D9DBDC"

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.bbox": "tight", "axes.titlesize": 9,
                     "axes.titleweight": "bold", "legend.frameon": False,
                     "axes.labelcolor": C_INK, "text.color": C_INK,
                     "axes.edgecolor": "#8C9196", "xtick.color": C_INK, "ytick.color": C_INK})


def save(fig, out_dir, name):
    for ext in ("png", "pdf"):
        fig.savefig(Path(out_dir) / f"{name}.{ext}")
    plt.close(fig)
    print(f"  wrote {name}.png/.pdf")


def panel(ax, letter):
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight="bold",
            va="bottom", ha="right")


def pdist_matrix(codes):
    """Uncorrected pairwise p-distance with pairwise deletion of gaps."""
    n = len(codes)
    D = np.full((n, n), np.nan)
    for i, j in itertools.combinations(range(n), 2):
        a, b = codes[i], codes[j]
        both = (a < 4) & (b < 4)
        m = int(both.sum())
        if m >= 50:
            D[i, j] = D[j, i] = float((a[both] != b[both]).sum()) / m
    return D


# ====== Fig 1: accuracy against the simulated amplicon benchmark (ground truth)
def figure_amplicons(out_dir):
    f = REPO / "results" / "benchmark_zenodo" / "zenodo_comparison.csv"
    if not f.exists():
        print(f"  SKIP: {f} missing")
        return
    c = pd.read_csv(f).fillna("")
    s = c[c.expected_genotype != ""].copy()
    s["typed"] = ~s.called_genotype.isin(["", "unassigned"])
    s["ok"] = s.called_genotype == s.expected_genotype
    s["acc"] = s.accepted == "accepted"
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.4), gridspec_kw={"width_ratios": [1, 1.5, 1.2]})

    ax = axes[0]
    counts = [int((s.ok & s.acc).sum()), int((s.ok & ~s.acc).sum()),
              int((~s.ok & s.typed & ~s.acc).sum()), int((~s.ok & s.typed & s.acc).sum()),
              int((~s.typed).sum())]
    labels = ["Correct,\naccepted", "Correct,\nnot accepted", "Wrong,\ndeclined",
              "Wrong,\naccepted", "Not typed"]
    cols = [C_OK, C_MID, C_UNTYPED, C_WRONG, C_NEUTRAL]
    bars = ax.bar(labels, counts, color=cols, width=0.68)
    for b, v in zip(bars, counts):
        ax.text(b.get_x() + b.get_width() / 2, v + max(counts) * 0.015, str(v),
                ha="center", fontsize=8)
    ax.set_ylabel(f"Amplicons (n={len(s)})")
    ax.set_title("Outcome against known composition")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=6.5, rotation=30, ha="right")
    ax.set_ylim(0, max(counts) * 1.16)
    panel(ax, "A")

    ax = axes[1]
    per = (s.groupby("expected_genotype")
             .agg(n=("query", "size"), correct=("ok", "sum")).reset_index())
    per = per[per.n >= 5].copy()
    per["acc"] = per.correct / per.n * 100
    per = per.sort_values("acc")
    bars = ax.barh(per.expected_genotype, per.acc,
                   color=[C_WRONG if a < 90 else C_OK for a in per.acc], height=0.68)
    for b, r in zip(bars, per.itertuples()):
        ax.text(min(r.acc + 1.5, 101), b.get_y() + b.get_height() / 2, f"n={r.n}",
                va="center", fontsize=6.5)
    ax.set_xlim(0, 118)
    ax.set_xlabel("Correctly typed (%)")
    ax.set_title("Accuracy per expected genotype ($n\\geq$5)")
    ax.tick_params(axis="y", labelsize=6.5)
    panel(ax, "B")

    ax = axes[2]
    sup = pd.to_numeric(s.support, errors="coerce")
    bins = np.arange(0, 105, 5)
    ax.hist(sup[s.ok].dropna(), bins=bins, color=C_OK, alpha=0.85, label="Correct", edgecolor="white")
    ax.hist(sup[~s.ok & s.typed].dropna(), bins=bins, color=C_WRONG, alpha=0.95,
            label="Discordant", edgecolor="white")
    ax.axvline(70, color=C_INK, ls="--", lw=1.2, label="Acceptance threshold")
    ax.set_yscale("log")
    ax.set_xlabel("Bootstrap support (%)")
    ax.set_ylabel("Amplicons (log scale)")
    ax.set_title("Support separates right from wrong")
    ax.legend(fontsize=7.5, loc="upper left")
    panel(ax, "C")
    fig.tight_layout()
    save(fig, out_dir, "figure7_amplicon_accuracy")


# ============================= Fig 2: norovirus vs RIVM (secondary comparison)
def figure1(out_dir):
    f = BENCH070 / "rivm_comparison_all_0.7.0.csv"
    if not f.exists():
        print(f"  SKIP: {f} missing")
        return
    d = pd.read_csv(f)
    d = d[d.region != "Genogroup (norovirus only)"].reset_index(drop=True)
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.2), gridspec_kw={"width_ratios": [1.5, 1.2, 1]})

    ax, x, w = axes[0], np.arange(len(d)), 0.26
    for i, (m, c) in enumerate((("precision", C_TOOL), ("recall", C_RIVM), ("f1", C_ZHAO))):
        ax.bar(x + (i - 1) * w, d[m], w, label="F1" if m == "f1" else m.capitalize(), color=c)
    ax.set_xticks(x)
    ax.set_xticklabels([r.replace(" ", "\n") for r in d.region], fontsize=8)
    ax.set_ylim(85, 104)
    ax.set_ylabel("Concordance with RIVM (%)")
    ax.set_title("Agreement by typing level")
    ax.legend(ncol=3, loc="upper center", fontsize=7.5)
    ax.grid(axis="y", alpha=0.25)
    panel(ax, "A")

    ax = axes[1]
    ax.barh(x, d.both_assigned, color=C_NEUTRAL, label="Typed by both")
    ax.barh(x, d.both_assigned * d.agree_of_both / 100, color=C_TOOL, label="Concordant")
    ax.set_yticks(x)
    ax.set_yticklabels([r.replace(" ", "\n") for r in d.region], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Sequences")
    ax.set_title("Sequences typed and concordant")
    ax.legend(fontsize=8, loc="lower right")
    for i, r in d.iterrows():
        ax.text(r.both_assigned * 1.01, i, f"{r.agree_of_both:.1f}%", va="center", fontsize=7.5)
    panel(ax, "B")

    ax = axes[2]
    ax.bar(x - 0.18, d.rivm_only, 0.36, color=C_RIVM, label="RIVM only")
    ax.bar(x + 0.18, d.calicityper_only, 0.36, color=C_TOOL, label="CaliciTyper only")
    ax.set_xticks(x)
    ax.set_xticklabels([r.split()[0] for r in d.region], fontsize=8)
    ax.set_ylabel("Sequences typed by one tool only")
    ax.set_title("Calls unique to one tool")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.25)
    panel(ax, "C")
    fig.tight_layout()
    save(fig, out_dir, "figure1_norovirus_rivm_concordance")


# ============================================= Fig 2: norovirus leave-one-out
def figure2(out_dir):
    f = BENCH070 / "loo_0.7.0.csv"
    if not f.exists():
        print(f"  SKIP: {f} missing")
        return
    d = pd.read_csv(f)
    d["correct"] = d.call == d.truth
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.2))

    ax = axes[0]
    g = d.groupby(["region", "group"]).correct.agg(["mean", "size"]).reset_index()
    labs = [f"{r.region.upper()}\n{r.group}" for _, r in g.iterrows()]
    ax.bar(labs, g["mean"] * 100, color=[C_TOOL if r == "orf1" else C_RIVM for r in g.region])
    for i, r in g.iterrows():
        ax.text(i, r["mean"] * 100 + 0.6, f"n={r['size']}", ha="center", fontsize=7.5)
    ax.set_ylim(80, 103)
    ax.set_ylabel("Leave-one-out accuracy (%)")
    ax.set_title("Reference recovery by region")
    ax.grid(axis="y", alpha=0.25)
    panel(ax, "A")

    ax = axes[1]
    sup = pd.to_numeric(d.support, errors="coerce").dropna()
    ax.hist(sup, bins=np.arange(0, 105, 5), color=C_TOOL, edgecolor="white")
    ax.axvline(70, color=C_WRONG, ls="--", lw=1.2, label="Acceptance threshold (70%)")
    ax.set_xlabel("Bootstrap support (%)")
    ax.set_ylabel("Reference sequences")
    ax.set_title("Support for recovered types")
    ax.legend(fontsize=8, loc="upper left")
    panel(ax, "B")

    ax = axes[2]
    sites = pd.to_numeric(d.n_sites, errors="coerce")
    sup2 = pd.to_numeric(d.support, errors="coerce")
    ok = d.correct.values
    ax.scatter(sites[ok], sup2[ok], s=12, alpha=0.6, color=C_TOOL, label="Correct")
    if (~ok).any():
        ax.scatter(sites[~ok], sup2[~ok], s=28, color=C_WRONG, marker="x", label="Declined")
    ax.set_xlabel("Aligned informative sites (nt)")
    ax.set_ylabel("Bootstrap support (%)")
    ax.set_title("Support versus alignable sequence")
    ax.legend(fontsize=8, loc="lower right")
    panel(ax, "C")
    fig.tight_layout()
    save(fig, out_dir, "figure2_norovirus_leave_one_out")


# ================= Fig 3: per-genotype concordance dot heatmap (HLA benchmark Fig 3)
def _noro_pergenotype():
    """Per-type concordance, using the SAME label normalisation as the headline
    comparison (scripts/validate_vs_rivm.py). Comparing the raw strings instead
    scores 'GII.P31 (GII.Pe)' against our 'GII.P31', and RIVM's plain 'GII.4'
    against our variant-level 'GII.4_Sydney', as disagreements when they are
    nomenclature differences."""
    ours, rivm = DEV / "all_typed.tsv", REPO / "validation" / "results_rivm.csv"
    if not (ours.exists() and rivm.exists()):
        return None
    sys.path.insert(0, str(REPO / "scripts"))
    import validate_vs_rivm as V
    o = pd.read_csv(ours, sep="\t", dtype=str).fillna("")
    r = pd.read_csv(rivm, dtype=str).fillna("").set_index("name")
    rows = []
    for level, ourcol, rivmcol in (("ORF2\ncapsid", "orf2_type", "capsid type"),
                                   ("ORF1\nP-type", "orf1_type", "polymerase type")):
        for _, q in o.iterrows():
            if q["query"] not in r.index:
                continue
            truth = V.genotype(r.at[q["query"], rivmcol])
            call = V.genotype(q.get(ourcol, ""))
            if truth in ("", "unassigned") or call in ("", "unassigned"):
                continue
            rows.append(dict(level=level, ref=truth, agree=call == truth))
    if not rows:
        return None
    g = pd.DataFrame(rows).groupby(["level", "ref"]).agree.agg(["mean", "size"]).reset_index()
    return g[g["size"] >= 5]


def _sapo_pergenotype():
    f = SAPO / "rivm_comparison.csv"
    if not f.exists():
        return None
    d = pd.read_csv(f).fillna("")
    d = d[(d.rivm_vp1 != "") & ~d.our_vp1.isin(["", "unassigned"])]
    if not len(d):
        return None
    sys.path.insert(0, str(REPO / "scripts"))
    import validate_vs_rivm as V
    d = d.assign(agree=[V.genotype(a) == V.genotype(b) for a, b in zip(d.our_vp1, d.rivm_vp1)],
                 level="Sapovirus\nVP1")
    g = d.groupby(["level", "rivm_vp1"]).agree.agg(["mean", "size"]).reset_index()
    return g.rename(columns={"rivm_vp1": "ref"})


def figure3(out_dir):
    noro, sapo = _noro_pergenotype(), _sapo_pergenotype()
    blocks = [(t, p) for t, p in (("Norovirus", noro), ("Sapovirus", sapo)) if p is not None]
    if not blocks:
        print("  SKIP: no per-genotype comparison data")
        return
    heights = [max(4, 0.17 * p.ref.nunique()) for _, p in blocks]
    fig, axes = plt.subplots(1, len(blocks), figsize=(4.8 * len(blocks), max(heights)),
                             gridspec_kw={"width_ratios": [1.35, 1][:len(blocks)]})
    axes = np.atleast_1d(axes)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("conc", [C_WRONG, C_MID, C_OK])
    norm = matplotlib.colors.Normalize(vmin=0, vmax=100)

    def ssize(n):
        return float(np.interp(np.sqrt(n), [np.sqrt(5), np.sqrt(600)], [22, 400]))

    for ax, (title, sub), letter in zip(axes, blocks, "AB"):
        levels = list(sub.level.unique())
        order = sub.groupby("ref")["mean"].mean().sort_values().index.tolist()
        ypos = {t: i for i, t in enumerate(order)}
        xpos = {l: i for i, l in enumerate(levels)}
        for _, r in sub.iterrows():
            ax.scatter(xpos[r.level], ypos[r.ref], s=ssize(r["size"]),
                       c=[cmap(norm(r["mean"] * 100))], edgecolor="white", linewidth=0.5, zorder=3)
        ax.set_xticks(range(len(levels)))
        ax.set_xticklabels(levels, fontsize=8)
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels(order, fontsize=7)
        ax.set_xlim(-0.65, len(levels) - 0.35)
        ax.set_ylim(-1, len(order))
        ax.grid(axis="y", alpha=0.18, zorder=0)
        ax.set_title(f"{title}: concordance per reference type")
        panel(ax, letter)
    axes[0].set_ylabel("Type assigned by the Norovirus Typing Tool\n(rows ordered by mean concordance)")

    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=list(axes), fraction=0.022, pad=0.02)
    cb.set_label("Concordance with RIVM (%)")
    handles = [Line2D([], [], marker="o", ls="", markersize=np.sqrt(ssize(s)),
                      markerfacecolor="#9AA0A6", markeredgecolor="white", label=str(s))
               for s in (10, 50, 200, 600)]
    axes[-1].legend(handles=handles, title="Sequences", loc="lower right", fontsize=7,
                    title_fontsize=7, labelspacing=1.2, borderpad=0.9)
    save(fig, out_dir, "figure3_concordance_per_genotype")


# ====== Fig 4: pairwise distance distributions with cut-offs (Chhabra 2019 framework)
def _distance_classes(ref_dir, cap=200):
    """(same type, different type same genogroup, different genogroup) distances."""
    from calicityper import tree
    within, between, across = [], [], []
    for orf in sorted(Path(ref_dir).glob("*.orf")):
        if orf.stem == "all" or not orf.with_suffix(".aln").exists():
            continue
        types = {r.id.split("|")[2]: r.id.split("|")[1] for r in SeqIO.parse(orf, "fasta")}
        aln = AlignIO.read(orf.with_suffix(".aln"), "clustal")
        if len(aln) < 2:
            continue
        names = [r.id.split("|")[-1] for r in aln]
        D = pdist_matrix(tree.encode([str(r.seq).upper() for r in aln]))
        lab = [types.get(n, "") for n in names]
        for i, j in itertools.combinations(range(len(names)), 2):
            if np.isfinite(D[i, j]) and lab[i] and lab[j]:
                (within if lab[i] == lab[j] else between).append(D[i, j])
    allp, orfp = Path(ref_dir) / "all.aln", Path(ref_dir) / "all.orf"
    if allp.exists() and orfp.exists():
        grp = {r.id.split("|")[2]: r.id.split("|")[0] for r in SeqIO.parse(orfp, "fasta")}
        aln = AlignIO.read(allp, "clustal")
        names = [r.id.split("|")[-1] for r in aln]
        keep = [i for i, n in enumerate(names) if n in grp][:cap]
        if len(keep) > 1:
            codes = tree.encode([str(aln[i].seq).upper() for i in keep])
            D = pdist_matrix(codes)
            g = [grp[names[i]] for i in keep]
            for i, j in itertools.combinations(range(len(keep)), 2):
                if np.isfinite(D[i, j]) and g[i] != g[j]:
                    across.append(D[i, j])
    return np.array(within), np.array(between), np.array(across)


def figure4(out_dir):
    from calicityper import config as cg
    panels = [("Norovirus VP1 (ORF2)", cg.REF_DIR["orf2"], None),
              ("Norovirus RdRp (ORF1)", cg.REF_DIR["orf1"], None),
              ("Sapovirus VP1", cg.SAPOVIRUS_REF_DIR["orf2"], cg.SAPOVIRUS_CUTOFFS["orf2"]),
              ("Sapovirus RdRp", cg.SAPOVIRUS_REF_DIR["orf1"], cg.SAPOVIRUS_CUTOFFS["orf1"])]
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.4))
    bins = np.arange(0, 0.92, 0.02)
    for ax, (title, refdir, cuts), letter in zip(axes, panels, "ABCD"):
        if not Path(refdir).exists():
            ax.set_visible(False)
            continue
        w, b, a = _distance_classes(refdir)
        for vals, lab, c in ((w, "Same type", C_TOOL),
                             (b, "Different type,\nsame genogroup", C_RIVM),
                             (a, "Different genogroup", C_ZHAO)):
            if len(vals):
                ax.hist(vals, bins=bins, alpha=0.6, color=c, density=True, edgecolor="none",
                        label=f"{lab} (n={len(vals)})")
        if cuts:
            top = ax.get_ylim()[1]
            for v, lab in ((cuts[0], "genotype"), (cuts[1], "genogroup")):
                ax.axvline(v, color=C_INK, ls="--", lw=1.1)
                ax.text(v, top * 0.98, f" {lab} {v}", rotation=90, fontsize=6.5,
                        va="top", ha="left")
        ax.set_xlabel("Pairwise nucleotide distance")
        ax.set_title(title)
        ax.legend(fontsize=6.2, loc="upper right")
        ax.set_xlim(0, 0.9)
        panel(ax, letter)
    axes[0].set_ylabel("Density of pairwise comparisons")
    fig.tight_layout()
    save(fig, out_dir, "figure4_distance_distributions")


# ====== Fig 5: truncation and acceptance threshold, with a table (HLA benchmark Fig 4)
def figure5(out_dir):
    f = REPO / "results" / "truncation_loo_sapovirus.csv"
    ff = REPO / "results" / "truncation_loo_sapovirus.floors.csv"
    if not f.exists():
        print(f"  SKIP: {f} missing")
        return
    d = pd.read_csv(f)
    known = d[~d.singleton]
    fig = plt.figure(figsize=(11.5, 4.9))
    gs = fig.add_gridspec(2, 3, height_ratios=[3, 1.1], hspace=0.55, wspace=0.32)

    ax = fig.add_subplot(gs[0, 0])
    for region, c, m in (("VP1", C_TOOL, "o"), ("RdRp", C_RIVM, "s")):
        s = known[known.region == region]
        if not len(s):
            continue
        g = s.groupby("length").apply(lambda x: (x.outcome == "correct").mean() * 100,
                                      include_groups=False)
        lens = [l if l else int(s.region_sites.median()) for l in g.index]
        o = np.argsort(lens)
        ax.plot(np.array(lens)[o], g.values[o], marker=m, color=c, label=region, lw=1.6, ms=4)
    ax.set_xlabel("Fragment length submitted (nt)")
    ax.set_ylabel("Correctly typed (%)")
    ax.set_ylim(70, 102)
    ax.set_title("Accuracy against fragment length")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.25)
    panel(ax, "A")

    axb = fig.add_subplot(gs[0, 1])
    fl = pd.read_csv(ff) if ff.exists() else None
    if fl is not None:
        for region, c, m in (("VP1", C_TOOL, "o"), ("RdRp", C_RIVM, "s")):
            s = fl[fl.region == region]
            if not len(s):
                continue
            axb.plot(s.floor, s.precision * 100, marker=m, color=c, lw=1.6, ms=4,
                     label=f"{region} precision")
            axb.plot(s.floor, s.recall * 100, marker=m, color=c, lw=1.2, ms=4, ls="--",
                     alpha=0.6, label=f"{region} recall")
        axb.set_xlabel("Required coverage of the typing region")
        axb.set_ylabel("Precision / recall (%)")
        axb.set_title("Effect of an acceptance threshold")
        axb.legend(fontsize=7, ncol=2, loc="lower left")
        axb.grid(alpha=0.25)
    panel(axb, "B")

    axc = fig.add_subplot(gs[0, 2])
    rows, labels = [], []
    for region in ("VP1", "RdRp"):
        for grp, lab in ((known[known.region == region], f"{region}\nknown"),
                         (d[(d.region == region) & d.singleton], f"{region}\nsingleton")):
            g = grp[grp.length == 0]
            if not len(g):
                continue
            rows.append([(g.outcome == o).mean() * 100 for o in ("correct", "untypeable", "wrong")])
            labels.append(lab)
    rows = np.array(rows)
    bottom = np.zeros(len(rows))
    for j, (o, c) in enumerate((("Correct", C_OK), ("Untypeable", C_UNTYPED),
                                ("Discordant", C_WRONG))):
        axc.bar(labels, rows[:, j], bottom=bottom, color=c, label=o, width=0.6)
        bottom += rows[:, j]
    axc.set_ylabel("Complete-region calls (%)")
    axc.set_title("Outcome at full length")
    axc.legend(fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.13),
               columnspacing=1.0, handlelength=1.2)
    axc.tick_params(axis="x", labelsize=8)
    panel(axc, "C")

    axt = fig.add_subplot(gs[1, :])
    axt.axis("off")
    if fl is not None:
        show = fl[fl.region == "VP1"].sort_values("floor")
        cells = [[f"{int(v)}" for v in show.accepted],
                 [f"{v:.3f}" for v in show.precision],
                 [f"{v:.3f}" for v in show.recall]]
        t = axt.table(cellText=cells, rowLabels=["Calls accepted", "Precision", "Recall"],
                      colLabels=[f"{v:.2f}" for v in show.floor], loc="center", cellLoc="center")
        t.auto_set_font_size(False)
        t.set_fontsize(7)
        t.scale(1, 1.2)
        axt.set_title("VP1 calls retained at each coverage threshold (columns match panel B)",
                      fontsize=8, pad=10)
    save(fig, out_dir, "figure6_truncation_and_threshold")


# ====== Fig 6: sapovirus against the published classification, and RIVM differences
def figure6(out_dir):
    f, g = SAPO / "distance_ranges.csv", SAPO / "rivm_comparison.csv"
    if not f.exists():
        print(f"  SKIP: {f} missing")
        return
    d = pd.read_csv(f)
    d = d[d.published != "NA"].copy()
    d["pub_max"] = d.published.str.split("-").str[-1].astype(float)
    d["ours_max"] = pd.to_numeric(d.ours_max, errors="coerce")
    d = d.dropna(subset=["ours_max"])
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), gridspec_kw={"width_ratios": [1, 1, 1.3]})

    ax = axes[0]
    for region, c, m in (("VP1", C_TOOL, "o"), ("RdRp", C_RIVM, "s")):
        s = d[d.region == region]
        ax.scatter(s.pub_max, s.ours_max, s=26, alpha=0.85, color=c, marker=m,
                   label=f"{region} (n={len(s)})")
    lim = max(d.pub_max.max(), d.ours_max.max()) * 1.1
    ax.plot([0, lim], [0, lim], color="#888", ls="--", lw=1, label="Exact agreement")
    ax.set_xlabel("Published maximum (Zhao et al.)")
    ax.set_ylabel("Observed maximum (this study)")
    ax.set_title("Within-type distances")
    ax.legend(fontsize=7.5, loc="upper left")
    panel(ax, "A")

    ax = axes[1]
    agree = d.groupby("region").max_within_published.apply(
        lambda s: (s.astype(str) == "True").mean() * 100)
    cnt = d.groupby("region").size()
    bars = ax.bar(agree.index, agree.values, color=[C_RIVM, C_TOOL], width=0.55)
    for b, r in zip(bars, agree.index):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 1.5,
                f"{int(round(agree[r]/100*cnt[r]))}/{cnt[r]}", ha="center", fontsize=8)
    ax.set_ylim(0, 108)
    ax.set_ylabel("Types inside the published range (%)")
    ax.set_title("Agreement with the classification")
    ax.grid(axis="y", alpha=0.25)
    panel(ax, "B")

    ax = axes[2]
    if g.exists():
        c = pd.read_csv(g).fillna("")
        s = c[(c.rivm_vp1 != "") & ~c.our_vp1.isin(["", "unassigned", "untypeable"])]
        diff = s[s.our_vp1 != s.rivm_vp1]
        pairs = diff.groupby(["rivm_vp1", "our_vp1"]).size().sort_values()
        if len(pairs):
            labs = [f"{a}  $\\rightarrow$  {b}" for a, b in pairs.index]
            ax.barh(labs, pairs.values, color=C_RIVM, height=0.6)
            for i, v in enumerate(pairs.values):
                ax.text(v + 0.15, i, str(v), va="center", fontsize=8)
            ax.set_xlabel("Sequences")
            ax.set_xlim(0, pairs.values.max() * 1.25)
    ax.set_title("RIVM label $\\rightarrow$ Zhao et al. genotype")
    ax.tick_params(axis="y", labelsize=8)
    panel(ax, "C")
    fig.tight_layout()
    save(fig, out_dir, "figure5_sapovirus_classification")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(Path(__file__).resolve().parent.parent / "figures"))
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    print(f"writing figures to {out}  (palette: {PALETTE})")
    for fn in (figure_amplicons, figure1, figure2, figure3, figure4, figure5, figure6):
        print(fn.__name__)
        try:
            fn(out)
        except Exception as e:
            import traceback
            print(f"  FAILED {fn.__name__}: {type(e).__name__}: {e}")
            traceback.print_exc(limit=3)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
