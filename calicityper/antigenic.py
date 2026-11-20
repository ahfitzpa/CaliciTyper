"""GII.4 antigenic-site annotation of aligned VP1 amino acids.

This is an annotation layer. It does not type a sequence, does not change a
call, and nothing here feeds the genotype, P-type, novelty or recombination
paths. It reports which published antigenic-site residues differ from the
nearest reference variant, and screens for a combination the literature
associates with antigenic change.

Residue numbers are GII.4 VP1 numbering and are resolved to *alignment
columns* through a named coordinate-anchor row (canonical 540-aa post-2002
GII.4 VP1, e.g. JX459908.1 GII.4_Sydney), never applied as raw indices. VP1
length varies between 530 and 550 residues across human genotypes (Parra et
al. 2026), so a residue number carried from GII.4 to another genotype points
at the wrong amino acid.

Sites and residues, verbatim from Parra et al. 2026 J Virol (CC BY 4.0,
doi:10.1128/jvi.01804-25):

    A  294-298, 368, 372, 373      variable, apex, near the HBGA pocket
    C  339-341, 375-378            variable, lateral
    D  393-398                     variable
    E  407-414, variable 407 and 412-414
    G  352, 355-357, 359, 364      variable
    F  254-257, 327, 404           conserved
    I  402, 403, 504-506           conserved; the mAb A1431 footprint

Informational secondary motifs (annotation only; never trigger shift):

    B          333, 382                       Lindesmith et al. 2012 dimer-interface site
    H          309, 310                       Lindesmith / Chhabra et al. 2025 Site H
    TOHMA      352, 357, 368, 378             Tohma et al. 2019 mBio emergence quartet
    HBGA_CORE  343, 344, 345, 374, 442-444    Conserved HBGA pocket (ONT frameshift/QC check)

The `candidate_antigenic_shift` screen follows Kendra et al. 2021 PNAS: at
least 18 VP1 amino-acid differences were needed for an antigenic distance of
4 or more, and shifts involved concurrent change at several of A, D, E and G.
Site C is excluded from the site criterion, and the ~22% false-positive caveat
travels with every output row.
"""
from __future__ import annotations

import math
from collections import Counter
from functools import lru_cache
from pathlib import Path
from dataclasses import dataclass, field

#: Primary variable antigenic sites, GII.4 540-aa VP1 numbering (Parra et al. 2026)
SITES = {
    "A": (294, 295, 296, 297, 298, 368, 372, 373),
    "C": (339, 340, 341, 375, 376, 377, 378),
    "D": (393, 394, 395, 396, 397, 398),
    "E": (407, 408, 409, 410, 411, 412, 413, 414),
    "G": (352, 355, 356, 357, 359, 364),
}
#: Variable subset of each primary site; Site E is only half variable
VARIABLE = dict(SITES, E=(407, 412, 413, 414))

#: Contiguous loop segments per primary site where internal insertions belong to the site
#: (e.g. GII.4 San Francisco 1-aa insertion inside Epitope A loop 294-298; Chhabra et al. 2025)
SITE_LOOPS = {
    "A": ((294, 298), (372, 373)),
    "C": ((339, 341), (375, 378)),
    "D": ((393, 398),),
    "E": ((407, 414),),
    "G": ((355, 357),),
}

#: Conserved antibody footprints, reported but never used to call a shift
CONSERVED = {
    "F": (254, 255, 256, 257, 327, 404),
    "I": (402, 403, 504, 505, 506),
}

#: Secondary / informational sites and QC motifs (annotation only; never trigger shift)
INFORMATIONAL = {
    "B": (333, 382),
    "H": (309, 310),
    "TOHMA": (352, 357, 368, 378),
    "HBGA_CORE": (343, 344, 345, 374, 442, 443, 444),
}

#: Sites whose concurrent change the literature associates with antigenic shift
SHIFT_SITES = ("A", "D", "E", "G")
#: Kendra et al. 2021 Table 1 epistatically coevolving residue pairs
COEVOLVING = (
    (297, 372),
    (309, 459),
    (357, 372),
    (357, 397),
    (395, 504),
    (397, 448),
    (465, 483),
    (119, 145),
)

CANONICAL_GII4_VP1_AA = 540
MIN_AA_CHANGES = 18

#: Input guard. The site coordinates are norovirus GII VP1 numbering, resolved
#: through a GII.4 anchor row, so this command accepts norovirus GII VP1 and
#: nothing else: a GI, animal-genogroup or sapovirus VP1 aligned against the
#: GII.4 anchor differs from it at far more than MIN_AA_CHANGES residues, so it
#: would be reported as an antigenic-change candidate on divergence alone.
#: A query is accepted when its nearest bundled VP1 reference, by shared
#: amino-acid k-mer containment, is a norovirus GII reference.
GUARD_K = 6
#: Floor for the nearest reference's containment. Below this the query does not
#: contain a bundled VP1 (truncated, frameshifted, wrong reading frame, or not
#: a calicivirus VP1) and the nearest-genogroup label carries no information.
#:
#: 0.3 from the panel check of 2026-11-20 (outputs/antigenic_guard_panel.tsv).
#: Across complete-VP1 records, true GII sat at a median containment of 0.955,
#: but the lowest true GII was KF931178.2 at 0.492, while the non-GII records
#: whose nearest reference happened to be GII sat at 0.093 or below. 0.05 let
#: those through and 0.5 refused KF931178.2, so the floor sits between them.
GUARD_MIN_CONTAINMENT = 0.3
CONSERVED_ENTROPY = 0.3
MAX_RESIDUE = max(
    max(r)
    for r in (
        list(SITES.values())
        + list(CONSERVED.values())
        + list(INFORMATIONAL.values())
        + [list(p) for p in COEVOLVING]
    )
)

CAVEAT = (
    "Screen only. Kendra et al. 2021: ~22% of pairs with >=18 amino-acid "
    "changes were not antigenically different, and 10 pairs were antigenically "
    "different with fewer than 18. Antigenic distance needs blockade or "
    "neutralisation titres and is not inferable from sequence."
)


class AntigenicError(ValueError):
    """Raised when the input cannot support a trustworthy annotation."""


def all_tracked_residues() -> list[int]:
    """Sorted union of all residues needed for primary, conserved, informational, and coevolving checks."""
    return sorted(
        {r for rs in SITES.values() for r in rs}
        | {r for rs in CONSERVED.values() for r in rs}
        | {r for rs in INFORMATIONAL.values() for r in rs}
        | {r for pair in COEVOLVING for r in pair}
    )


@dataclass
class Annotation:
    query: str
    reference: str
    genotype: str
    aa_differences: int
    compared_sites: int
    site_changes: dict = field(default_factory=dict)
    conserved_changes: dict = field(default_factory=dict)
    informational_changes: dict = field(default_factory=dict)
    coevolving_pairs_changed: list = field(default_factory=list)
    candidate_antigenic_shift: bool = False
    hbga_pocket_qc_alert: bool = False
    coordinate_warning: str = ""
    reason: str = ""
    site_entropy: dict = field(default_factory=dict)
    caveat: str = CAVEAT

    def haplotype(self) -> str:
        """A:TTH...|C:...|D:...|E:...|G:... from the per-site residues (including loop indels)."""
        return "|".join(
            f"{s}:{self.site_changes.get(s, {}).get('query_residues', '')}"
            for s in sorted(SITES)
        )


def shannon_entropy(column) -> float | None:
    """Shannon entropy in bits over the amino acids present in one column."""
    aa = [c for c in column if c not in "-.X"]
    if not aa:
        return None
    n = len(aa)
    h = -sum((k / n) * math.log2(k / n) for k in Counter(aa).values())
    return round(h, 4) + 0.0


def site_entropy(rows, columns) -> dict[int, float | None]:
    """{residue: entropy} over every sequence in the alignment."""
    out = {}
    for res, col in columns.items():
        out[res] = shannon_entropy([r[col] for r in rows if col < len(r)])
    return out


def _resolved_span(row: str) -> tuple[int, int] | None:
    idxs = [i for i, ch in enumerate(row.upper()) if ch not in "-.X?"]
    if not idxs:
        return None
    return idxs[0], idxs[-1]


def residue_columns(reference_row: str, residues, min_length: int | None = None) -> dict[int, int]:
    """{residue number: alignment column} for a gapped full-length VP1 coordinate anchor."""
    ungapped = sum(1 for c in reference_row if c not in "-.")
    need = min_length if min_length is not None else MAX_RESIDUE
    if ungapped < need:
        raise AntigenicError(
            f"the reference row has {ungapped} residues; annotation needs at "
            f"least {need}, because the published site coordinates are "
            f"full-length VP1 numbering. A trimmed or P2-only alignment "
            f"cannot be annotated: its residue 1 is not VP1 residue 1."
        )
    out, n = {}, 0
    wanted = set(residues)
    for col, ch in enumerate(reference_row):
        if ch in "-.":
            continue
        n += 1
        if n in wanted:
            out[n] = col
    missing = sorted(wanted - set(out))
    if missing:
        raise AntigenicError(f"reference does not reach residue(s) {missing}")
    return out


def compare(
    query_row: str,
    reference_row: str,
    columns: dict[int, int],
    include_internal_indels: bool = True,
) -> tuple[dict[int, tuple[str, str]], int]:
    """{residue: (ref_state, query_state)} for scorable differences (substitutions + internal indels)."""
    q_up, r_up = query_row.upper(), reference_row.upper()
    q_span = _resolved_span(q_up)
    r_span = _resolved_span(r_up)
    overlap = None
    if q_span is not None and r_span is not None:
        overlap = (max(q_span[0], r_span[0]), min(q_span[1], r_span[1]))

    diffs, scored = {}, 0
    for res, col in columns.items():
        if col >= len(q_up) or col >= len(r_up):
            continue
        q, r = q_up[col], r_up[col]
        if q == "X" or r == "X":
            continue
        q_gap = q in "-."
        r_gap = r in "-."
        if not q_gap and not r_gap:
            scored += 1
            if q != r:
                diffs[res] = (r, q)
        elif include_internal_indels and overlap is not None and overlap[0] <= col <= overlap[1]:
            if q_gap and not r_gap:
                scored += 1
                diffs[res] = (r, "del")
            elif not q_gap and r_gap:
                scored += 1
                diffs[res] = ("ins", q)
    return diffs, scored


def _loop_insertions(
    query_row: str,
    reference_row: str,
    anchor_row: str,
    columns: dict[int, int],
    loops: tuple[tuple[int, int], ...],
) -> list[tuple[int, int, str, str]]:
    """Detect internal insertions sitting in alignment columns where the coordinate anchor has a gap."""
    q_up, r_up, a_up = query_row.upper(), reference_row.upper(), anchor_row.upper()
    q_span = _resolved_span(q_up)
    r_span = _resolved_span(r_up)
    if q_span is None or r_span is None:
        return []
    ov_start, ov_end = max(q_span[0], r_span[0]), min(q_span[1], r_span[1])

    events = []
    for res_start, res_end in loops:
        if res_start not in columns or res_end not in columns:
            continue
        for res in range(res_start, res_end):
            if res not in columns or (res + 1) not in columns:
                continue
            c_left, c_right = columns[res], columns[res + 1]
            for col in range(c_left + 1, c_right):
                if not (ov_start <= col <= ov_end):
                    continue
                if a_up[col] not in "-.":
                    continue
                q, r = q_up[col], r_up[col]
                if q == "X" or r == "X" or q == r:
                    continue
                if q not in "-." and r in "-.":
                    events.append((res, col, "ins", q))
                elif q in "-." and r not in "-.":
                    events.append((res, col, r, "del"))
                elif q not in "-." and r not in "-." and q != r:
                    events.append((res, col, r, q))
    return events


def total_aa_differences(query_row: str, reference_row: str) -> int:
    """Amino-acid differences over internal overlapping columns (substitutions + internal indels)."""
    q_up, r_up = query_row.upper(), reference_row.upper()
    q_span = _resolved_span(q_up)
    r_span = _resolved_span(r_up)
    if q_span is None or r_span is None:
        return 0
    start = max(q_span[0], r_span[0])
    end = min(q_span[1], r_span[1])
    if start > end:
        return 0

    n = 0
    for col in range(start, end + 1):
        q, r = q_up[col], r_up[col]
        if q == "X" or r == "X":
            continue
        if (q in "-.") and (r in "-."):
            continue
        if q != r:
            n += 1
    return n


def _format_change(res: int, a: str, b: str) -> str:
    if a == "ins":
        return f"{res}ins{b}"
    if b == "del":
        return f"{res}{a}>del"
    return f"{res}{a}>{b}"


def annotate(
    query_name: str,
    query_row: str,
    reference_name: str,
    reference_row: str,
    genotype: str = "GII.4",
    entropy: dict[int, float | None] | None = None,
    anchor_row: str | None = None,
    anchor_name: str | None = None,
) -> Annotation:
    """Annotate one aligned query against one aligned reference."""
    if len(query_row) != len(reference_row):
        raise AntigenicError(
            f"query and reference are not from one alignment: "
            f"{len(query_row)} vs {len(reference_row)} columns"
        )
    coord_row = anchor_row if anchor_row is not None else reference_row
    if len(coord_row) != len(reference_row):
        raise AntigenicError(
            f"coordinate anchor and reference are not from one alignment: "
            f"{len(coord_row)} vs {len(reference_row)} columns"
        )

    all_res = all_tracked_residues()
    cols = residue_columns(coord_row, all_res)

    ungapped_anchor = sum(1 for c in coord_row if c not in "-.")
    coord_warn = ""
    if ungapped_anchor != CANONICAL_GII4_VP1_AA:
        who = anchor_name or reference_name
        coord_warn = (
            f"coordinate anchor {who} has {ungapped_anchor} ungapped aa "
            f"(canonical post-2002 GII.4 VP1 is {CANONICAL_GII4_VP1_AA} aa); "
            f"pass --anchor JX459908.1 to avoid indel-induced numbering shifts"
        )

    ann = Annotation(
        query=query_name,
        reference=reference_name,
        genotype=genotype,
        aa_differences=total_aa_differences(query_row, reference_row),
        compared_sites=0,
        coordinate_warning=coord_warn,
    )

    changed_shift_sites = []
    for site, residues in SITES.items():
        sub = {r: cols[r] for r in residues}
        diffs, scored = compare(query_row, reference_row, sub)
        loop_ins = _loop_insertions(
            query_row, reference_row, coord_row, cols, SITE_LOOPS.get(site, ())
        )
        scored += len(loop_ins)
        ann.compared_sites += scored

        change_tokens = [_format_change(r, a, b) for r, (a, b) in sorted(diffs.items())]
        for flank_res, _col, a, b in loop_ins:
            if a == "ins":
                change_tokens.append(f"{flank_res}_{flank_res + 1}ins{b}")
            elif b == "del":
                change_tokens.append(f"{flank_res}_{flank_res + 1}del{a}")
            else:
                change_tokens.append(f"{flank_res}_{flank_res + 1}{a}>{b}")

        ins_after: dict[int, list[str]] = {}
        for flank_res, col, _a, _b in loop_ins:
            q_ch = query_row[col].upper()
            if q_ch not in "-.":
                ins_after.setdefault(flank_res, []).append(q_ch)

        q_chars = []
        for r in residues:
            q_chars.append(query_row[cols[r]].upper() if cols[r] < len(query_row) else "?")
            if r in ins_after:
                q_chars.append(f"[+{''.join(ins_after[r])}]")

        variable = set(VARIABLE[site])
        has_var_change = any(r in variable for r in diffs) or any(
            (flank_res in variable or (flank_res + 1) in variable)
            for flank_res, *_ in loop_ins
        )
        ann.site_changes[site] = {
            "n_changed": len(diffs) + len(loop_ins),
            "n_scored": scored,
            "changes": ",".join(change_tokens),
            "query_residues": "".join(q_chars),
        }
        if site in SHIFT_SITES and has_var_change:
            changed_shift_sites.append(site)

    for site, residues in CONSERVED.items():
        diffs, _ = compare(query_row, reference_row, {r: cols[r] for r in residues})
        ann.conserved_changes[site] = ",".join(
            _format_change(r, a, b) for r, (a, b) in sorted(diffs.items())
        )

    for site, residues in INFORMATIONAL.items():
        diffs, _ = compare(query_row, reference_row, {r: cols[r] for r in residues})
        ann.informational_changes[site] = ",".join(
            _format_change(r, a, b) for r, (a, b) in sorted(diffs.items())
        )
    ann.hbga_pocket_qc_alert = bool(ann.informational_changes.get("HBGA_CORE"))

    coevolving_res = sorted({r for pair in COEVOLVING for r in pair})
    coev_diffs, _ = compare(query_row, reference_row, {r: cols[r] for r in coevolving_res})
    changed_coev = set(coev_diffs.keys())
    ann.coevolving_pairs_changed = [
        f"{a}-{b}" for a, b in COEVOLVING if a in changed_coev and b in changed_coev
    ]

    if entropy:
        ann.site_entropy = {
            s: round(
                sum(entropy[r] for r in SITES[s] if entropy.get(r) is not None)
                / max(1, sum(1 for r in SITES[s] if entropy.get(r) is not None)),
                4,
            )
            for s in SITES
        }

    enough_changes = ann.aa_differences >= MIN_AA_CHANGES
    enough_sites = len(changed_shift_sites) >= 2
    ann.candidate_antigenic_shift = enough_changes and enough_sites
    reason = (
        f"{ann.aa_differences} aa differences "
        f"({'>=' if enough_changes else '<'}{MIN_AA_CHANGES}); "
        f"variable changes at {len(changed_shift_sites)} of A/D/E/G"
        + (f" ({', '.join(changed_shift_sites)})" if changed_shift_sites else "")
        + "; site C excluded from the site criterion"
    )
    if ann.hbga_pocket_qc_alert:
        reason += f"; QC alert: conserved HBGA pocket change ({ann.informational_changes['HBGA_CORE']})"
    if ann.coordinate_warning:
        reason += f"; NOTE: {ann.coordinate_warning}"
    ann.reason = reason
    return ann


def sites_table() -> list[dict[str, str]]:
    """The residue definitions, for ``--show-sites`` and documentation."""
    rows = []
    for site, res in sorted(SITES.items()):
        rows.append({
            "site": site,
            "kind": "variable",
            "residues": ",".join(map(str, res)),
            "variable": ",".join(map(str, VARIABLE[site])),
        })
    for site, res in sorted(CONSERVED.items()):
        rows.append({
            "site": site,
            "kind": "conserved",
            "residues": ",".join(map(str, res)),
            "variable": "",
        })
    for site, res in sorted(INFORMATIONAL.items()):
        rows.append({
            "site": site,
            "kind": "informational",
            "residues": ",".join(map(str, res)),
            "variable": "",
        })
    return rows


TSV_COLUMNS = [
    "query",
    "reference",
    "genotype",
    "aa_differences",
    "residues_compared",
    "epitope_haplotype",
    "mean_site_entropy",
    "site_A_changes",
    "site_C_changes",
    "site_D_changes",
    "site_E_changes",
    "site_G_changes",
    "conserved_F_changes",
    "conserved_I_changes",
    "site_B_changes",
    "site_H_changes",
    "tohma_quartet_changes",
    "hbga_pocket_qc",
    "coevolving_pairs_changed",
    "candidate_antigenic_shift",
    "reason",
    "caveat",
]


def to_row(ann: Annotation) -> dict[str, object]:
    r = {
        "query": ann.query,
        "reference": ann.reference,
        "genotype": ann.genotype,
        "aa_differences": ann.aa_differences,
        "residues_compared": ann.compared_sites,
        "epitope_haplotype": ann.haplotype(),
        "mean_site_entropy": ";".join(
            f"{s}={ann.site_entropy[s]}" for s in sorted(ann.site_entropy)
        ),
        "site_B_changes": ann.informational_changes.get("B", ""),
        "site_H_changes": ann.informational_changes.get("H", ""),
        "tohma_quartet_changes": ann.informational_changes.get("TOHMA", ""),
        "hbga_pocket_qc": ann.informational_changes.get("HBGA_CORE", ""),
        "coevolving_pairs_changed": " ".join(ann.coevolving_pairs_changed),
        "candidate_antigenic_shift": str(ann.candidate_antigenic_shift).lower(),
        "reason": ann.reason,
        "caveat": ann.caveat,
    }
    for s in SITES:
        r[f"site_{s}_changes"] = ann.site_changes.get(s, {}).get("changes", "")
    for s in CONSERVED:
        r[f"conserved_{s}_changes"] = ann.conserved_changes.get(s, "")
    return r


# --------------------------------------------------------------- immunotypes
#: Fallback matching calicityper/data/immunotype_regime.tsv (Parra et al. 2017 rebuilt tree)
_DEFAULT_IMMUNOTYPE = {
    "GI.1": "A",
    "GI.2": "B", "GI.4": "B", "GI.5": "B", "GI.6": "B",
    "GI.3": "C", "GI.7": "C", "GI.8": "C", "GI.9": "C",
    "GII.1": "D", "GII.12": "D", "GII.16": "D",
    "GII.2": "E", "GII.5": "E", "GII.10": "E",
    "GII.3": "F",
    "GII.4": "G", "GII.20": "G",
    "GII.6": "H",
    "GII.7": "I", "GII.8": "I", "GII.9": "I", "GII.14": "I",
    "GII.13": "J", "GII.17": "J", "GII.21": "J",
    "GII.15": "K",
    "GII.22": "L",
}
EVOLVING = ("GII.4",)
CANONICAL_ANCHOR_ACC = "JX459908.1"


def _load_packaged_immunotypes() -> dict[str, str]:
    import csv
    from pathlib import Path

    tsv = Path(__file__).resolve().parent / "data" / "immunotype_regime.tsv"
    if not tsv.exists():
        return dict(_DEFAULT_IMMUNOTYPE)
    out: dict[str, str] = {}
    with tsv.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            imm = (r.get("immunotype") or "").strip()
            gt = (r.get("genotype") or "").strip()
            if gt and len(imm) == 1 and imm.isalpha():
                out[gt] = imm
    return out or dict(_DEFAULT_IMMUNOTYPE)


IMMUNOTYPE = _load_packaged_immunotypes()


def immunotype(genotype: str | None) -> tuple[str | None, str | None]:
    """(immunotype letter or None, 'evolving'|'static'|None) for a genotype."""
    g = (genotype or "").split("_")[0]
    letter = IMMUNOTYPE.get(g)
    if letter is None:
        return None, None
    return letter, ("evolving" if g in EVOLVING else "static")


# ----------------------------------------- turnkey reference loading & MAFFT
def _to_vp1_aa(seq: str) -> str:
    """Auto-convert nucleotide ORF2 CDS or amino-acid VP1 string into ungapped AA."""
    clean = seq.replace("-", "").replace(".", "").strip().upper()
    if not clean:
        return ""
    nt_frac = sum(1 for c in clean if c in "ACGTUN") / len(clean)
    if nt_frac >= 0.90 and len(clean) >= MAX_RESIDUE * 3:
        from Bio.Seq import Seq

        usable = len(clean) - (len(clean) % 3)
        aa = str(Seq(clean[:usable]).translate(to_stop=False))
        if aa.endswith("*"):
            aa = aa[:-1]
        return aa
    return clean


def _kmers(aa: str, k: int = GUARD_K) -> frozenset:
    return frozenset(aa[i:i + k] for i in range(len(aa) - k + 1))


@lru_cache(maxsize=None)
def vp1_genogroup_panel(data_dir=None):
    """[(virus, genogroup, accession, kmers)] from the bundled VP1 reference sets.

    Both genera are loaded because the guard has to tell a sapovirus VP1 from a
    norovirus one, and sapovirus genogroups reuse the names `GI`, `GII`, ... so
    they are distinguished by which reference set they came from, never by
    label. Headers are `genogroup|type|accession` and the sequences are
    nucleotide ORFs, translated here.
    """
    from Bio import SeqIO

    from . import config as cg

    panel = []
    sources = (
        ("norovirus", cg.REF_DIR["orf2"], "reference_sequences_VP1"),
        ("sapovirus", cg.SAPOVIRUS_REF_DIR["orf2"], "reference_sequences_sapovirus_VP1"),
    )
    for virus, ref_dir, leaf in sources:
        base = Path(data_dir) / leaf if data_dir else Path(ref_dir)
        orf = base / "all.orf"
        if not orf.exists():
            continue
        for rec in SeqIO.parse(str(orf), "fasta"):
            parts = rec.id.split("|")
            if len(parts) < 3:
                continue
            gg, _ty, acc = parts[:3]
            aa = _to_vp1_aa(str(rec.seq))
            if len(aa) >= GUARD_K:
                panel.append((virus, gg, acc, _kmers(aa)))
    return tuple(panel)


def nearest_vp1_reference(aa: str, data_dir=None):
    """(containment, virus, genogroup, accession) for the closest bundled VP1.

    Containment is the fraction of the REFERENCE's amino-acid k-mers present in
    the query, which separates genogroups without needing an alignment.

    Normalised by the reference and not by the query on purpose: dividing by the
    query's k-mer count penalises length, so a whole-genome translation scored
    about 0.001 against the VP1 it contains and the guard refused real GII
    sequences (2,817 of 3,747 in the panel check of 2026-11-20). Reference
    normalisation asks "is a known VP1 inside this sequence", which is the
    question the guard needs answered.
    """
    q = _kmers(aa)
    best = (0.0, None, None, None)
    if not q:
        return best
    for virus, gg, acc, ref_k in vp1_genogroup_panel(data_dir):
        if not ref_k:
            continue
        c = len(q & ref_k) / len(ref_k)
        if c > best[0]:
            best = (c, virus, gg, acc)
    return best


def check_gii_input(query_name: str, aa: str, data_dir=None) -> tuple[str, float]:
    """Accept norovirus GII VP1 only; raise AntigenicError otherwise.

    Returns (nearest accession, containment) so the caller can record what the
    query was matched against.
    """
    cont, virus, gg, acc = nearest_vp1_reference(aa, data_dir)
    if virus is None or cont < GUARD_MIN_CONTAINMENT:
        raise AntigenicError(
            f"query {query_name!r} matches no bundled VP1 reference "
            f"(best {GUARD_K}-mer containment {cont:.3f} < {GUARD_MIN_CONTAINMENT}); "
            f"`calicityper antigenic` takes a norovirus GII VP1 (complete ORF2)."
        )
    if virus != "norovirus" or gg.split(".")[0] != "GII":
        where = f"{virus} {gg}"
        raise AntigenicError(
            f"query {query_name!r} is nearest {where} reference {acc} "
            f"(containment {cont:.3f}), not norovirus GII. The site coordinates are "
            f"GII VP1 numbering resolved through a GII.4 anchor, so a {where} VP1 "
            f"would be scored as antigenic change on divergence alone. "
            f"Type the sequence first (`calicityper type`) and pass GII VP1 only."
        )
    return acc, cont


def load_gii4_variant_references(data_dir=None) -> list[tuple[str, str]]:
    """Load translated full-length GII.4 VP1 references using packaged gii4_variant_references.tsv.

    Places canonical 540-aa ``JX459908.1|Sydney`` first so residue-to-column mapping
    uses the canonical coordinate anchor automatically.
    """
    import csv
    from pathlib import Path
    from . import config as cg

    pkg_data = Path(__file__).resolve().parent / "data" / "gii4_variant_references.tsv"
    if not pkg_data.exists():
        raise AntigenicError(f"packaged variant table missing: {pkg_data}")

    root = Path(data_dir) if data_dir is not None else cg.DATA
    vp1_dir = root / "reference_sequences_VP1"
    nt_file = next(
        (p for p in (vp1_dir / "GII.orf", vp1_dir / "GII.fasta") if p.exists()),
        None,
    )
    if nt_file is None:
        raise FileNotFoundError(
            f"calicityper reference data not found in {vp1_dir}; run `calicityper setup` or set CALICITYPER_DATA"
        )

    nt_by_acc: dict[str, str] = {}
    hdr, buf = None, []
    with nt_file.open(encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    s = "".join(buf)
                    for tok in hdr.split("|"):
                        nt_by_acc.setdefault(tok, s)
                        nt_by_acc.setdefault(tok.split(".")[0], s)
                hdr, buf = line[1:].split()[0], []
            else:
                buf.append(line)
        if hdr is not None:
            s = "".join(buf)
            for tok in hdr.split("|"):
                nt_by_acc.setdefault(tok, s)
                nt_by_acc.setdefault(tok.split(".")[0], s)

    with pkg_data.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    rows.sort(key=lambda r: (0 if r["accession"].strip() == CANONICAL_ANCHOR_ACC else 1))

    refs: list[tuple[str, str]] = []
    for r in rows:
        acc, var = r["accession"].strip(), r["variant"].strip()
        nt = nt_by_acc.get(acc) or nt_by_acc.get(acc.split(".")[0])
        if not nt:
            continue
        aa = _to_vp1_aa(nt)
        if len(aa) >= MAX_RESIDUE:
            refs.append((f"{acc}|{var}", aa))
    if not refs:
        raise AntigenicError(f"no full-length GII.4 VP1 references recovered from {nt_file}")
    return refs


def align_unaligned_queries(
    queries: list[tuple[str, str]],
    threads: int = 1,
    data_dir=None,
) -> tuple[list[tuple[str, str]], list[tuple[str, str]], tuple[str, str]]:
    """Translate (if NT) and align queries together with bundled GII.4 variant references via MAFFT.

    Returns ``(aligned_refs, aligned_queries, anchor_record)``.
    """
    import subprocess
    import tempfile
    from pathlib import Path
    from . import config as cg

    refs = load_gii4_variant_references(data_dir=data_dir)
    conv_queries: list[tuple[str, str]] = []
    for qid, raw_seq in queries:
        aa = _to_vp1_aa(raw_seq)
        if len(aa) < MAX_RESIDUE:
            raise AntigenicError(
                f"query {qid!r} translates/spans {len(aa)} aa (< {MAX_RESIDUE} aa); "
                f"`calicityper antigenic` requires full-length VP1 (complete ORF2)."
            )
        conv_queries.append((qid, aa))

    with tempfile.TemporaryDirectory(prefix="calicityper_ag_") as td:
        in_fa = Path(td) / "in.faa"
        with in_fa.open("w", encoding="utf-8") as fh:
            for i, (_rid, rseq) in enumerate(refs):
                fh.write(f">__ref_{i}\n{rseq}\n")
            for j, (_qid, qseq) in enumerate(conv_queries):
                fh.write(f">__qry_{j}\n{qseq}\n")
        cmd = [cg.exe("mafft"), "--quiet", "--amino", "--thread", str(max(1, threads)), "--auto", str(in_fa)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"MAFFT failed aligning GII.4 VP1 amino acids:\n{proc.stderr}")

    parsed: dict[str, str] = {}
    cur, buf = None, []
    for line in proc.stdout.splitlines():
        if line.startswith(">"):
            if cur is not None:
                parsed[cur] = "".join(buf).upper()
            cur, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if cur is not None:
        parsed[cur] = "".join(buf).upper()

    aln_refs = [(refs[i][0], parsed[f"__ref_{i}"]) for i in range(len(refs))]
    aln_qrys = [(conv_queries[j][0], parsed[f"__qry_{j}"]) for j in range(len(conv_queries))]
    anchor = aln_refs[0]
    return aln_refs, aln_qrys, anchor

