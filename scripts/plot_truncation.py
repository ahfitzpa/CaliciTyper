#!/usr/bin/env python3
"""Figure: what sequence length does to genotyping, and where CaliciTyper and
RIVM/NoroNet actually diverge."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
TRUNC = "/home/people/fitzpatria/scratch/calicityper/results/truncation_loo.csv"
OUT = "/home/people/fitzpatria/scratch/calicityper/docs/figures/truncation_f1.png"

df = pd.read_csv(TRUNC).drop_duplicates(subset=["region", "group", "accession", "nt", "offset"])
live = df[df.rejected_at_input == 0].copy()
live["ts"] = live.truth.astype(str).str.split("_").str[0]
live["cs"] = live.call.astype(str).str.split("_").str[0]


def f1(d, gate=None, support=70.0):
    a = (~d.call.isin(["unassigned", "rejected"])) & (d.support >= support)
    if gate is not None:
        a = a & (d.window_cov >= gate)
    c = a & (d.cs == d.ts)
    n, na = len(d), int(a.sum())
    if n == 0:
        return np.nan
    p = 100 * c.sum() / na if na else 0.0
    r = 100 * c.sum() / n
    return 2 * p * r / (p + r) if (p + r) else 0.0


def style(ax, title, sub=None):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8.5, length=3)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.set_title(title, color=INK, fontsize=10.5, fontweight="semibold", loc="left", pad=14 if sub else 6)
    if sub:
        ax.text(0, 1.015, sub, transform=ax.transAxes, color=INK2, fontsize=8.5, va="bottom")


fig = plt.figure(figsize=(12.6, 8.6), facecolor=SURF)
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.08], hspace=0.52, wspace=0.2,
                      left=0.065, right=0.975, top=0.845, bottom=0.085)

# --- panels 1 and 2: F1 against submitted length, amplicon on target -------
LENGTHS = [172, 200, 252, 267, 300, 400, 600, 900, 1200]
for k, (region, label, window) in enumerate((
        ("orf1", "ORF1 P-type  ·  RdRp 3' window, 172 nt", 172),
        ("orf2", "ORF2 capsid  ·  region C window, 252-267 nt", 252))):
    ax = fig.add_subplot(gs[0, k])
    part = live[(live.region == region) & (live.offset == 0)]
    xs, no_gate, gated = [], [], []
    for L in LENGTHS:
        p = part[part.length == L]
        if len(p) < 20:
            continue
        xs.append(L)
        no_gate.append(f1(p))
        g = f1(p, gate=0.8)
        gated.append(0.0 if np.isnan(g) else g)
    ax.plot(xs, no_gate, color=BLUE, lw=2, marker="o", ms=5.5, label="no window gate", zorder=3)
    ax.plot(xs, gated, color=ORANGE, lw=2, marker="s", ms=5.5, label="80% window gate (current default)", zorder=3)
    ax.axvline(window, color=INK2, lw=1, ls=":", zorder=1)
    ax.text(window, 46, f" typing window\n {window} nt", color=INK2, fontsize=7.5, va="bottom")
    for x, y in ((xs[0], no_gate[0]), (xs[-1], no_gate[-1])):
        ax.annotate(f"{y:.0f}", (x, y), textcoords="offset points", xytext=(0, 9),
                    ha="center", fontsize=8, color=INK)
    if gated[0] == 0:
        ax.annotate("no call: shorter\nthan the window", (xs[0], 0),
                    textcoords="offset points", xytext=(-6, 16), fontsize=8, color=ORANGE)
    ax.set_xscale("log")
    ticks = [t for t in xs if t != 267]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.minorticks_off()
    ax.set_ylim(-4, 104)
    ax.set_xlabel("nucleotides submitted, amplicon centred on the typing window", color=INK2, fontsize=9)
    ax.set_ylabel("F1 (genotype, leave-one-out)", color=INK2, fontsize=9)
    style(ax, label, "on-target fragments: longer sequence types better, and the gate costs little")
    ax.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="center right")

# --- panel 3: displacement is what the gate actually punishes ---------------
ax = fig.add_subplot(gs[1, 0])
groups = [("orf1", 0), ("orf1", 100), ("orf1", 300), ("orf2", 0), ("orf2", 100), ("orf2", 300)]
x = np.arange(len(groups), dtype=float)
x[3:] += 0.6
ng = [f1(live[(live.region == r) & (live.offset == o)]) for r, o in groups]
gt = [f1(live[(live.region == r) & (live.offset == o)], gate=0.8) for r, o in groups]
gt = [0.0 if np.isnan(v) else v for v in gt]
w = 0.38
b1 = ax.bar(x - w / 2, ng, w, color=BLUE, label="no window gate", zorder=3)
b2 = ax.bar(x + w / 2, gt, w, color=ORANGE, label="80% window gate", zorder=3)
for rect, v in list(zip(b1, ng)) + list(zip(b2, gt)):
    ax.annotate("none" if v == 0 else f"{v:.0f}", (rect.get_x() + rect.get_width() / 2, v),
                textcoords="offset points", xytext=(0, 3), ha="center", fontsize=8,
                color=INK if v else ORANGE)
ax.set_xticks(x)
ax.set_xticklabels([f"{o} nt" for _, o in groups], fontsize=8.5)
ax.text(1, -0.17, "ORF1 P-type", ha="center", color=INK, fontsize=9, fontweight="semibold",
        transform=ax.get_xaxis_transform())
ax.text(4.6, -0.17, "ORF2 capsid", ha="center", color=INK, fontsize=9, fontweight="semibold",
        transform=ax.get_xaxis_transform())
ax.set_ylim(0, 125)
ax.set_xlabel("how far the fragment is displaced off the typing window", color=INK2, fontsize=9, labelpad=38)
ax.set_ylabel("F1 (genotype, leave-one-out)", color=INK2, fontsize=9)
style(ax, "What the gate costs is displaced sequence",
      "off-window fragments still type at 99-100% precision, but the gate withholds every call")
ax.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="upper center", ncol=2,
          bbox_to_anchor=(0.5, 1.02))

# --- panel 4: why RIVM calls where we abstain ------------------------------
ax = fig.add_subplot(gs[1, 1])
reasons = ["bootstrap support\n< 70%", "region not covered\n(<100 nt overlap)",
           "genotyping restricted\nto GI/GII", "RIVM window\ncoverage < 80%"]
orf2 = [82, 27, 23, 4]
orf1 = [14, 10, 21, 10]
y = np.arange(len(reasons))[::-1].astype(float)
h = 0.36
for off, vals, col, lab in ((+h / 2, orf2, BLUE, "ORF2 capsid (136 records)"),
                            (-h / 2, orf1, AQUA, "ORF1 P-type (55 records)")):
    bars = ax.barh(y + off, vals, h, color=col, label=lab, zorder=3)
    for rect, v in zip(bars, vals):
        ax.annotate(str(v), (v, rect.get_y() + rect.get_height() / 2),
                    textcoords="offset points", xytext=(4, 0), va="center", fontsize=8.5, color=INK)
ax.set_yticks(y)
ax.set_yticklabels(reasons, fontsize=8.5)
ax.set_xlim(0, 95)
ax.grid(axis="x", color=GRID, lw=0.8)
ax.grid(axis="y", visible=False)
ax.set_xlabel("records RIVM typed and CaliciTyper declined", color=INK2, fontsize=9)
style(ax, "Where the two tools actually diverge",
      "when both call, they agree on 100% of 3,956 capsid and 3,763 P-type records")
ax.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="lower right")

fig.suptitle("Genotyping accuracy against sequence length, and the real source of the RIVM gap",
             x=0.065, y=0.975, ha="left", color=INK, fontsize=14, fontweight="semibold")
fig.text(0.065, 0.905,
         "Leave-one-out over the bundled GI/GII reference sets (7,698 queries, 100 bootstraps each), scored against each "
         "reference's published type.  Fragments are cut relative to\neach record's own typing window, so a length at or "
         "above the window covers it fully and displacement, not length, is what varies coverage.",
         ha="left", color=INK2, fontsize=9)
fig.savefig(OUT, dpi=170, facecolor=SURF)
print("wrote", OUT)
