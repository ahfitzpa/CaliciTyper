"""Human-readable reports: Markdown, JSON, and a self-contained HTML page.

The TSV stays the machine-readable primary output. These are for looking at a
run: the same table, plus a Bokeh alignment panel per call showing the query
against the references of the type it was assigned to, plus a
neighbour-joining tree per region placing the queries among the references.

The point of the alignment panel is that a bootstrap percentage is not
inspectable. Seeing where a query actually differs from its assigned type, and
from the runner-up, is.

Bokeh is optional. Without it the Markdown and JSON reports still work and the
HTML is written without the figures.
"""
from __future__ import annotations

import html
import json
from datetime import datetime, timezone

import numpy as np

try:
    from bokeh.embed import components
    from bokeh.layouts import column as bk_column
    from bokeh.models import ColumnDataSource, HoverTool, Range1d
    from bokeh.plotting import figure
    from bokeh.resources import CDN, INLINE
    HAVE_BOKEH = True
except ImportError:                                     # pragma: no cover
    HAVE_BOKEH = False

# Colour-blind-safe base colours (Okabe-Ito derived); gaps/ambiguity stay grey.
BASE_COLOURS = {"A": "#0072B2", "C": "#E69F00", "G": "#009E73", "T": "#D55E00",
                "U": "#D55E00", "-": "#E8E8E8", "N": "#BBBBBB"}
OTHER_COLOUR = "#BBBBBB"
MAX_ALN_WIDTH = 1200        # px
MAX_TREE_TAXA = 120         # beyond this a tree is unreadable; it is skipped with a note


def _colour(ch):
    return BASE_COLOURS.get(ch.upper(), OTHER_COLOUR)


def alignment_figure(names, seqs, title, highlight=None):
    """Bokeh alignment panel: one row per sequence, one coloured cell per base.

    Follows the approach of Farrell's Bokeh sequence aligner
    (https://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner):
    a rect glyph per residue, with letters drawn only when zoomed in far enough
    to read them.
    """
    if not HAVE_BOKEH or not seqs:
        return None
    n, L = len(seqs), len(seqs[0])
    text = [ch for s in seqs for ch in s]
    colours = [_colour(ch) for ch in text]
    x = np.tile(np.arange(L) + 0.5, n)
    y = np.repeat(np.arange(n)[::-1] + 0.5, L)
    source = ColumnDataSource(dict(x=x, y=y, text=text, colours=colours))

    width = min(MAX_ALN_WIDTH, max(420, L * 7))
    height = max(120, 18 * n + 60)
    view_cols = min(L, 180)                 # initial window; the panel pans/zooms
    p = figure(title=title, width=width, height=height,
               x_range=Range1d(0, view_cols, bounds=(0, L)),
               y_range=[*reversed(names)], tools="xpan,xwheel_zoom,reset,save",
               active_scroll="xwheel_zoom", toolbar_location="right")
    p.rect("x", "y", 1.0, 1.0, source=source, fill_color="colours",
           line_color=None, fill_alpha=0.9)
    p.text("x", "y", text="text", source=source, text_font_size="8pt",
           text_align="center", text_baseline="middle", text_color="#222222")
    p.grid.visible = False
    p.xaxis.axis_label = "alignment column"
    p.yaxis.major_label_text_font_size = "8pt"
    if highlight and highlight in names:
        p.yaxis.major_label_overrides = {highlight: f"▶ {highlight}"}
    return p


_BASES = set("ACGTU")
MATCH_COLOUR = "#D9D9D9"     # covered and identical to the closest reference
GAP_COLOUR = "#444444"       # deletion in the query relative to the reference


def _runs(mask):
    """[start, end) runs of True in a boolean array."""
    if not mask.any():
        return []
    d = np.diff(np.concatenate([[0], mask.astype(np.int8), [0]]))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


def closest_reference(qa, ref_seqs):
    """Index of the reference nearest the aligned query, by uncorrected p-distance.

    Compared over columns where both query and reference carry an unambiguous
    base, so partial queries are not penalised for what they do not cover.
    Returns (index, p_distance, n_compared, n_different), or None when no
    reference shares a single comparable column with the query.
    """
    q = np.array(list(qa.upper()))
    qok = np.isin(q, list(_BASES))
    best = None
    for k, s in enumerate(ref_seqs):
        r = np.array(list(s.upper()))
        both = qok & np.isin(r, list(_BASES))
        n = int(both.sum())
        if n == 0:
            continue
        diff = int((q[both] != r[both]).sum())
        p = diff / n
        if best is None or (p, -n) < (best[1], -best[2]):
            best = (k, p, n, diff)
    return best


def difference_figure(qname, qseq, rname, rseq, title):
    """Query against its closest reference, drawing only what differs.

    Identical stretches are one grey bar per row; each differing base is a
    coloured, labelled cell in both rows, and a deletion in the query is a dark
    cell. Columns the query does not cover are left blank. Because only the
    differences become glyphs, a panel stays small enough to draw one for every
    call in a large run, and the whole region reads at a glance: zoom in to
    read the bases, or hover a cell for its column.
    """
    if not HAVE_BOKEH or not qseq:
        return None
    q = np.array(list(qseq.upper()))
    r = np.array(list(rseq.upper()))
    L = len(q)
    qb = np.isin(q, list(_BASES))
    rb = np.isin(r, list(_BASES))
    covered = np.zeros(L, bool)
    if qb.any():
        lo, hi = np.flatnonzero(qb)[[0, -1]]
        covered[lo:hi + 1] = True                   # query span, including internal gaps
    mism = qb & rb & (q != r)
    qdel = covered & (q == "-") & rb

    rows = [qname, rname]
    p = figure(title=title, width=MAX_ALN_WIDTH, height=150,
               x_range=Range1d(0, L, bounds=(0, L)), y_range=list(reversed(rows)),
               tools="xpan,xwheel_zoom,reset,save", active_scroll="xwheel_zoom",
               toolbar_location="right")
    bars = dict(y=[], left=[], right=[])
    for row, mask in ((qname, covered & qb & ~mism), (rname, covered & rb & ~mism & ~qdel)):
        for a, b in _runs(mask):
            bars["y"].append(row); bars["left"].append(int(a)); bars["right"].append(int(b))
    if bars["y"]:
        p.hbar(y="y", left="left", right="right", height=0.6, source=ColumnDataSource(bars),
               fill_color=MATCH_COLOUR, line_color=None)
    cols = np.flatnonzero(mism | qdel)
    if len(cols):
        qch, rch = list(q[cols]), list(r[cols])
        text = qch + rch
        src = ColumnDataSource(dict(
            x=np.concatenate([cols, cols]) + 0.5, y=[qname] * len(cols) + [rname] * len(cols),
            text=text, colours=[GAP_COLOUR if ch == "-" else _colour(ch) for ch in text],
            col=np.concatenate([cols, cols]) + 1, q=qch * 2, r=rch * 2))
        cells = p.rect("x", "y", 1.0, 0.8, source=src, fill_color="colours", line_color=None)
        p.text("x", "y", text="text", source=src, text_font_size="7pt",
               text_align="center", text_baseline="middle", text_color="#FFFFFF")
        p.add_tools(HoverTool(renderers=[cells],
                              tooltips=[("column", "@col"), ("query", "@q"), ("reference", "@r")]))
    p.grid.visible = False
    p.xaxis.axis_label = "alignment column (reference coordinates)"
    p.yaxis.major_label_text_font_size = "8pt"
    p.yaxis.major_label_overrides = {qname: f"▶ {qname}"}
    return p


def _nj_tree(names, D):
    """Neighbour-joining tree from a square distance matrix, via Biopython."""
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    lower = [[float(D[i][j]) for j in range(i + 1)] for i in range(len(names))]
    return DistanceTreeConstructor().nj(DistanceMatrix(list(names), lower))


def _layout(tree):
    """Rectangular cladogram coordinates: {clade: (x, y)} plus the leaf order."""
    leaves = tree.get_terminals()
    ypos = {lf: i for i, lf in enumerate(leaves)}
    xpos = {}

    def depth(clade, acc=0.0):
        xpos[clade] = acc
        for ch in clade.clades:
            depth(ch, acc + (ch.branch_length or 0.0))

    depth(tree.root)

    def yof(clade):
        if clade in ypos:
            return ypos[clade]
        v = sum(yof(c) for c in clade.clades) / len(clade.clades)
        ypos[clade] = v
        return v

    yof(tree.root)
    return xpos, ypos, leaves


def tree_figure(names, D, title, highlight=()):
    """Bokeh neighbour-joining tree. Queries are drawn in a contrasting colour."""
    if not HAVE_BOKEH or len(names) < 3:
        return None
    if len(names) > MAX_TREE_TAXA:
        return None
    tree = _nj_tree(names, D)
    tree.root_at_midpoint()
    xpos, ypos, leaves = _layout(tree)

    xs, ys = [], []
    for clade in tree.find_clades(order="level"):
        for ch in clade.clades:                        # elbow: vertical then horizontal
            xs += [[xpos[clade], xpos[clade]], [xpos[clade], xpos[ch]]]
            ys += [[ypos[clade], ypos[ch]], [ypos[ch], ypos[ch]]]

    hl = set(highlight)
    lx = [xpos[lf] for lf in leaves]
    ly = [ypos[lf] for lf in leaves]
    lab = [lf.name for lf in leaves]
    col = ["#D55E00" if nm in hl else "#444444" for nm in lab]

    height = max(220, 14 * len(leaves) + 70)
    p = figure(title=title, width=760, height=height,
               tools="pan,wheel_zoom,reset,save", toolbar_location="right")
    p.multi_line(xs=xs, ys=ys, line_color="#888888", line_width=1)
    p.scatter(lx, ly, size=5, fill_color=col, line_color=None)
    p.text(x=lx, y=ly, text=[f" {n}" for n in lab], text_font_size="8pt",
           text_baseline="middle", text_color=col)
    p.yaxis.visible = False
    p.grid.visible = False
    p.xaxis.axis_label = "substitutions per site"
    return p


# ---------------------------------------------------------------- serialisers

def _runners_to_list(s):
    out = []
    for part in (s or "").split("|"):
        if ":" in part:
            t, v = part.rsplit(":", 1)
            try:
                out.append({"type": t, "support": float(v)})
            except ValueError:
                pass
    return out


def to_json(rows, path, meta=None):
    """Full results as JSON, with runner-up types as structured lists."""
    recs = []
    for r in rows:
        rec = dict(r)
        for region in ("orf1", "orf2"):
            rec[f"{region}_runners"] = _runners_to_list(r.get(f"{region}_runners"))
        recs.append(rec)
    doc = {"tool": "calicityper",
           "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
           "n_queries": len(rows),
           "meta": meta or {},
           "results": recs}
    with open(path, "w") as fh:
        json.dump(doc, fh, indent=2)
    return path


CALL_COLUMNS = ["query", "genogroup", "orf1_type", "orf1_support", "orf1_runners", "orf1_putative_novel",
                "orf2_type", "orf2_support", "orf2_runners", "orf2_putative_novel", "note"]
NOVEL_FLAGS = ("putative_novel_type", "putative_novel_group")
NOVEL_LABEL = {("orf1", "putative_novel_type"): "putative novel RdRp P-type",
               ("orf1", "putative_novel_group"): "putative novel RdRp P-group",
               ("orf2", "putative_novel_type"): "putative novel VP1 genotype",
               ("orf2", "putative_novel_group"): "putative novel VP1 genogroup"}
NOVEL_HEADING = "Putative novel ORF1/ORF2 types"
NOVEL_NONE = "No query fell outside every known ORF1 (RdRp) or ORF2 (VP1) type."
NOVEL_NOTE = ("Putative novel ORF1/ORF2 types are screened with the 2x SD criterion of "
              "Chhabra et al. 2019 (complete VP1 amino acids for ORF2; partial RdRp "
              "nucleotides for ORF1) against the bundled references. A flag marks one "
              "sequence outside every known cluster: worth a closer look, tentative (NA) "
              "at most, never a new designation.")


def novel_label(region, flag):
    return NOVEL_LABEL.get((region, flag), "")


def flagged(rows):
    """[(query, region, flag, detail)] for every putative novel ORF1/ORF2 type."""
    out = []
    for r in rows:
        for region in ("orf1", "orf2"):
            f = r.get(f"{region}_putative_novel") or ""
            if f in NOVEL_FLAGS:
                out.append((r["query"], region, f, r.get(f"{region}_putative_novel_detail") or ""))
    return out


def _md_cell(x):
    """Escape a value for a Markdown table cell.

    Pipes are not decoration here: FASTA ids look like "GII|NC_029646.1" and the
    runner-up column is "GII.P12:93.0|GII.P31:1.0". Unescaped, each one ends the
    cell early and the whole table collapses into the wrong columns.
    """
    return str(x if x is not None else "").replace("|", "\\|").replace("\n", " ")


def _md_table(rows, columns):
    out = ["| " + " | ".join(columns) + " |",
           "| " + " | ".join("---" for _ in columns) + " |"]
    for r in rows:
        out.append("| " + " | ".join(_md_cell(r.get(c, "")) for c in columns) + " |")
    return "\n".join(out)


def to_markdown(rows, path, meta=None):
    cols = CALL_COLUMNS
    gg = {}
    for r in rows:
        gg[r.get("genogroup") or "unassigned"] = gg.get(r.get("genogroup") or "unassigned", 0) + 1
    lines = [
        "# calicityper results",
        "",
        f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')}  ",
        f"Queries: **{len(rows)}**",
        "",
        "## Genogroups",
        "",
        _md_table([{"genogroup": k, "n": v} for k, v in sorted(gg.items())], ["genogroup", "n"]),
        "",
        f"## {NOVEL_HEADING}",
        "",
        *([f"- **{q}** {reg.upper()}: {novel_label(reg, f)} -- {d}" for q, reg, f, d in flagged(rows)]
          or [NOVEL_NONE]),
        "",
        f"> {NOVEL_NOTE}",
        "",
        "## Calls",
        "",
        _md_table(rows, cols),
        "",
        "> Genotyping is validated for GI and GII only. Calls for other genogroups",
        "> are reported but unverified. Bootstrap support is not a guarantee of",
        "> correctness; check the runner-up column for near-threshold calls.",
        "",
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
    return path


_CSS = """
:root { --fg:#1a1a1a; --bg:#ffffff; --muted:#5c5c5c; --line:#e2e2e2; --accent:#0072B2;
  --flag:#B3261E; --flagbg:#FCEEEE; }
@media (prefers-color-scheme: dark) { :root:not([data-theme=light]) {
  --fg:#e8e8e8; --bg:#161616; --muted:#a8a8a8; --line:#333; --accent:#5bb0e8;
  --flag:#FF8A80; --flagbg:#3A1F1F; } }
* { box-sizing:border-box }
body { margin:0; padding:24px 16px; background:var(--bg); color:var(--fg);
  font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; }
.wrap { max-width:1240px; margin:0 auto }
h1 { font-size:1.6rem; margin:0 0 4px } h2 { font-size:1.15rem; margin:32px 0 10px }
.sub { color:var(--muted); font-size:.9rem; margin-bottom:8px }
table { border-collapse:collapse; width:100%; font-size:.85rem; margin:8px 0 }
th,td { border-bottom:1px solid var(--line); padding:6px 8px; text-align:left;
  white-space:nowrap; font-variant-numeric:tabular-nums }
th { font-weight:600; color:var(--muted); font-size:.78rem; text-transform:uppercase;
  letter-spacing:.04em }
tbody tr:hover { background:color-mix(in srgb, var(--accent) 7%, transparent) }
.tblwrap { overflow-x:auto; -webkit-overflow-scrolling:touch }
.note { color:var(--muted); font-size:.85rem }
tr.flagged td { background:var(--flagbg) }
td.flag { color:var(--flag); font-weight:600; white-space:nowrap }
.flagbox { border-left:3px solid var(--flag); padding:6px 14px; margin:12px 0; background:var(--flagbg) }
.warn { border-left:3px solid #E69F00; padding:10px 14px; margin:16px 0;
  background:color-mix(in srgb,#E69F00 8%,transparent); font-size:.88rem; border-radius:0 4px 4px 0 }
.panel { border:1px solid var(--line); border-radius:6px; padding:10px; margin:14px 0 }
.pill { display:inline-block; padding:1px 7px; border-radius:10px; font-size:.75rem;
  background:color-mix(in srgb,var(--accent) 14%,transparent); color:var(--accent) }
"""


def to_html(rows, path, figures=None, meta=None, notes=(), embed_js=True):
    """HTML report: table, alignment panels, trees.

    embed_js=True inlines BokehJS so the file genuinely stands alone -- it opens
    with no network, which is what an archived surveillance report needs. It
    costs roughly 3 MB. embed_js=False links the Bokeh CDN instead, giving a
    much smaller file that needs internet to render.
    """
    cols = CALL_COLUMNS
    figures = [f for f in (figures or []) if f is not None]

    script, divs = "", []
    if HAVE_BOKEH and figures:
        script, divs = components(bk_column(*figures, sizing_mode="scale_width")) \
            if len(figures) > 1 else components(figures[0])
        divs = [divs] if isinstance(divs, str) else list(divs)

    gg = {}
    for r in rows:
        k = r.get("genogroup") or "unassigned"
        gg[k] = gg.get(k, 0) + 1

    def esc(x):
        return html.escape(str(x if x is not None else ""))

    head = "".join(f"<th>{esc(c)}</th>" for c in cols)
    def cell(r, c):
        v = r.get(c, "")
        if c.endswith("_putative_novel") and v in NOVEL_FLAGS:
            return f'<td class="flag">&#9873; {esc(novel_label(c[:4], v))}</td>'
        return f"<td>{esc(v)}</td>"

    def is_flagged(r):
        return any((r.get(f"{g}_putative_novel") or "") in NOVEL_FLAGS for g in ("orf1", "orf2"))

    flag_items = "".join(f"<li><strong>{esc(q)}</strong> {reg.upper()}: "
                         f"{esc(novel_label(reg, f))} &mdash; {esc(d)}</li>"
                         for q, reg, f, d in flagged(rows))
    novel_html = (f'<div class="flagbox"><ul>{flag_items}</ul></div>' if flag_items else
                  f'<p class="note">{esc(NOVEL_NONE)}</p>')
    body = "".join(('<tr class="flagged">' if is_flagged(r) else "<tr>")
                   + "".join(cell(r, c) for c in cols) + "</tr>"
                   for r in rows)
    ggrows = "".join(f"<tr><td>{esc(k)}</td><td>{v}</td></tr>" for k, v in sorted(gg.items()))
    bokeh_css = ((INLINE if embed_js else CDN).render()
                 if (HAVE_BOKEH and figures) else "")
    notes_html = "".join(f'<div class="warn">{esc(n)}</div>' for n in notes)

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>calicityper report</title>
{bokeh_css}
<style>{_CSS}</style></head><body><div class="wrap">
<h1>calicityper report</h1>
<div class="sub">{esc(datetime.now(timezone.utc).isoformat(timespec='seconds'))}
 &middot; <span class="pill">{len(rows)} queries</span></div>
<div class="warn">Genotyping is validated for <strong>GI and GII only</strong>. Calls for
other genogroups are reported but unverified. Bootstrap support is not a guarantee of
correctness &mdash; check the runner-up column for near-threshold calls.</div>
{notes_html}
<h2>Genogroups</h2><div class="tblwrap"><table><thead><tr><th>genogroup</th><th>n</th></tr>
</thead><tbody>{ggrows}</tbody></table></div>
<h2>{NOVEL_HEADING}</h2>
{novel_html}
<p class="note">{esc(NOVEL_NOTE)}</p>
<h2>Calls</h2><div class="tblwrap"><table><thead><tr>{head}</tr></thead>
<tbody>{body}</tbody></table></div>
<h2>Alignments and trees</h2>
{"".join(f'<div class="panel">{d}</div>' for d in divs)
 if divs else '<p class="note">No figures: install bokeh for alignment panels and trees.</p>'}
{script}
</div></body></html>"""
    with open(path, "w") as fh:
        fh.write(doc)
    return path
