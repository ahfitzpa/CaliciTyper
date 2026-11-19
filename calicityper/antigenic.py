"""GII.4 antigenic-site annotation of aligned VP1 amino acids.

This is an annotation layer. It does not type a sequence, does not change a
call, and nothing here feeds the genotype, P-type, novelty or recombination
paths. It reports which published antigenic-site residues differ from the
nearest reference variant, and screens for a combination the literature
associates with antigenic change.

Residue numbers are GII.4 VP1 numbering and are resolved to *alignment
columns* through a named reference row, never applied as raw indices. VP1
length varies between 530 and 550 residues across human genotypes (Parra et
al. 2026), so a residue number carried from GII.4 to another genotype points
at the wrong amino acid. The same reason the typing windows are column
ranges applies here.

Sites and residues, verbatim from Parra et al. 2026 J Virol (CC BY 4.0,
doi:10.1128/jvi.01804-25):

    A  294-298, 368, 372, 373      variable, apex, near the HBGA pocket
    C  339-341, 375-378            variable, lateral
    D  393-398                     variable
    E  407-414, variable 407 and 412-414
    G  352, 355-357, 359, 364      variable
    F  254-257, 327, 404           conserved
    I  402, 403, 504-506           conserved; the mAb A1431 footprint

The `candidate_antigenic_shift` screen follows Kendra et al. 2021 PNAS: at
least 18 VP1 amino-acid differences were needed for an antigenic distance of
4 or more, and shifts involved concurrent change at several of A, D, E and G.
Two things keep this honest. About 22% of pairs with 18 or more changes were
*not* antigenically different, and ten pairs reached antigenic distance 4
with fewer than 18 changes, so the screen is neither specific nor sensitive;
and site C was enriched among pairs that were not antigenically different, so
site C changes do not count toward the site criterion. Antigenic distance is
measured with blockade and neutralisation titres. It cannot be computed from
sequence, and this module does not claim to.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

#: residue sets, GII.4 VP1 numbering. See the module docstring for the source.
SITES = {
    "A": (294, 295, 296, 297, 298, 368, 372, 373),
    "C": (339, 340, 341, 375, 376, 377, 378),
    "D": (393, 394, 395, 396, 397, 398),
    "E": (407, 408, 409, 410, 411, 412, 413, 414),
    "G": (352, 355, 356, 357, 359, 364),
}
#: the variable subset of each site; site E is only half variable
VARIABLE = dict(SITES, E=(407, 412, 413, 414))
#: conserved sites, reported but never used to call a shift
CONSERVED = {"F": (254, 255, 256, 257, 327, 404),
             "I": (402, 403, 504, 505, 506)}
#: sites whose concurrent change the literature associates with antigenic shift
SHIFT_SITES = ("A", "D", "E", "G")
#: Kendra et al. 2021 Table 1
COEVOLVING = ((297, 372), (309, 459), (357, 372), (357, 397),
              (395, 504), (397, 448), (465, 483), (119, 145))

MIN_AA_CHANGES = 18
#: the lowest residue the annotation needs to be able to locate
MAX_RESIDUE = max(max(r) for r in list(SITES.values()) + list(CONSERVED.values()))

CAVEAT = (
    "Screen only. Kendra et al. 2021: ~22% of pairs with >=18 amino-acid "
    "changes were not antigenically different, and 10 pairs were antigenically "
    "different with fewer than 18. Antigenic distance needs blockade or "
    "neutralisation titres and is not inferable from sequence."
)


class AntigenicError(ValueError):
    """Raised when the input cannot support a trustworthy annotation."""


@dataclass
class Annotation:
    query: str
    reference: str
    genotype: str
    aa_differences: int
    compared_sites: int
    site_changes: dict = field(default_factory=dict)
    conserved_changes: dict = field(default_factory=dict)
    coevolving_pairs_changed: list = field(default_factory=list)
    candidate_antigenic_shift: bool = False
    reason: str = ""
    caveat: str = CAVEAT

    def haplotype(self):
        """A:TTH...|C:...|D:...|E:...|G:... from the per-site residues."""
        return "|".join(f"{s}:{self.site_changes.get(s, {}).get('query_residues', '')}"
                        for s in sorted(SITES))


def residue_columns(reference_row, residues, min_length=None):
    """{residue number: alignment column} for a gapped reference row.

    The reference must be full-length VP1. A P2-only alignment, such as one
    starting at residue 279, silently renumbers every residue: its first
    column would be read as residue 1 and the site coordinates would point at
    the wrong amino acids.
    """
    ungapped = sum(1 for c in reference_row if c not in "-.")
    need = min_length if min_length is not None else MAX_RESIDUE
    if ungapped < need:
        raise AntigenicError(
            f"the reference row has {ungapped} residues; annotation needs at "
            f"least {need}, because the published site coordinates are "
            f"full-length VP1 numbering. A trimmed or P2-only alignment "
            f"cannot be annotated: its residue 1 is not VP1 residue 1.")
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


def compare(query_row, reference_row, columns):
    """{residue: (ref aa, query aa)} for residues that differ and are scorable."""
    diffs, scored = {}, 0
    for res, col in columns.items():
        if col >= len(query_row):
            continue
        q, r = query_row[col].upper(), reference_row[col].upper()
        if q in "-.X" or r in "-.X":
            continue
        scored += 1
        if q != r:
            diffs[res] = (r, q)
    return diffs, scored


def total_aa_differences(query_row, reference_row):
    """Amino-acid differences over columns both sequences resolve."""
    n = 0
    for q, r in zip(query_row.upper(), reference_row.upper()):
        if q in "-.X" or r in "-.X":
            continue
        n += q != r
    return n


def annotate(query_name, query_row, reference_name, reference_row,
             genotype="GII.4"):
    """Annotate one aligned query against one aligned reference."""
    if len(query_row) != len(reference_row):
        raise AntigenicError(
            f"query and reference are not from one alignment: "
            f"{len(query_row)} vs {len(reference_row)} columns")

    all_res = sorted({r for rs in SITES.values() for r in rs}
                     | {r for rs in CONSERVED.values() for r in rs})
    cols = residue_columns(reference_row, all_res)

    ann = Annotation(query=query_name, reference=reference_name,
                     genotype=genotype,
                     aa_differences=total_aa_differences(query_row, reference_row),
                     compared_sites=0)

    changed_shift_sites = []
    for site, residues in SITES.items():
        sub = {r: cols[r] for r in residues}
        diffs, scored = compare(query_row, reference_row, sub)
        ann.compared_sites += scored
        variable = set(VARIABLE[site])
        ann.site_changes[site] = {
            "n_changed": len(diffs),
            "n_scored": scored,
            "changes": ",".join(f"{r}{a}>{b}" for r, (a, b) in sorted(diffs.items())),
            "query_residues": "".join(
                query_row[cols[r]].upper() if cols[r] < len(query_row) else "?"
                for r in residues),
        }
        if diffs and site in SHIFT_SITES and any(r in variable for r in diffs):
            changed_shift_sites.append(site)

    for site, residues in CONSERVED.items():
        diffs, _ = compare(query_row, reference_row,
                           {r: cols[r] for r in residues})
        ann.conserved_changes[site] = ",".join(
            f"{r}{a}>{b}" for r, (a, b) in sorted(diffs.items()))

    changed = {r for s in SITES for r in
               (int(x[:-3]) for x in ann.site_changes[s]["changes"].split(",") if x)}
    ann.coevolving_pairs_changed = [f"{a}-{b}" for a, b in COEVOLVING
                                    if a in changed and b in changed]

    enough_changes = ann.aa_differences >= MIN_AA_CHANGES
    enough_sites = len(changed_shift_sites) >= 2
    ann.candidate_antigenic_shift = enough_changes and enough_sites
    ann.reason = (
        f"{ann.aa_differences} aa differences "
        f"({'>=' if enough_changes else '<'}{MIN_AA_CHANGES}); "
        f"variable changes at {len(changed_shift_sites)} of A/D/E/G"
        + (f" ({', '.join(changed_shift_sites)})" if changed_shift_sites else "")
        + "; site C excluded from the site criterion")
    return ann


def sites_table():
    """The residue definitions, for `--show-sites` and for documentation."""
    rows = []
    for site, res in sorted(SITES.items()):
        rows.append({"site": site, "kind": "variable",
                     "residues": ",".join(map(str, res)),
                     "variable": ",".join(map(str, VARIABLE[site]))})
    for site, res in sorted(CONSERVED.items()):
        rows.append({"site": site, "kind": "conserved",
                     "residues": ",".join(map(str, res)), "variable": ""})
    return rows


TSV_COLUMNS = ["query", "reference", "genotype", "aa_differences",
               "residues_compared", "epitope_haplotype",
               "site_A_changes", "site_C_changes", "site_D_changes",
               "site_E_changes", "site_G_changes",
               "conserved_F_changes", "conserved_I_changes",
               "coevolving_pairs_changed", "candidate_antigenic_shift",
               "reason", "caveat"]


def to_row(ann):
    r = {"query": ann.query, "reference": ann.reference, "genotype": ann.genotype,
         "aa_differences": ann.aa_differences,
         "residues_compared": ann.compared_sites,
         "epitope_haplotype": ann.haplotype(),
         "coevolving_pairs_changed": " ".join(ann.coevolving_pairs_changed),
         "candidate_antigenic_shift": str(ann.candidate_antigenic_shift).lower(),
         "reason": ann.reason, "caveat": ann.caveat}
    for s in SITES:
        r[f"site_{s}_changes"] = ann.site_changes.get(s, {}).get("changes", "")
    for s in CONSERVED:
        r[f"conserved_{s}_changes"] = ann.conserved_changes.get(s, "")
    return r
