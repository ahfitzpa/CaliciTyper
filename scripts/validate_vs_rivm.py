#!/usr/bin/env python3
"""Compare calicityper output with RIVM Norovirus Typing Tool calls.

Needs validation/results_rivm.csv (not redistributed here -- see
validation/README.md).

    calicityper rivm_test_set.fasta -o calicityper_out.tsv --threads 16
    python scripts/validate_vs_rivm.py calicityper_out.tsv

Label conventions differ and must be reconciled before anything is comparable:

* RIVM splits GII.4 across two columns: `capsid type` is plain "GII.4" for all
  1,868 GII.4 records, and the variant lives in `capsid subtype`
  ("Sydney_2012", "Den_Haag_2006b", ...).
* calicityper carries the variant inside the type itself ("GII.4_Sydney"),
  because its ORF2 reference set is labelled at variant level.

Comparing those two fields literally scores every GII.4 sequence as a
disagreement -- 42% of the test set. So genotype agreement is measured on the
genotype stem (GII.4_Sydney -> GII.4), and the GII.4 variant is reported
separately against RIVM's subtype column.
"""
import argparse
import re

import pandas as pd

UNASSIGNED = {"", "could not assign", "unassigned", "nan", "none"}


def norm(x):
    """RIVM/calicityper label -> comparable genotype, or 'unassigned'."""
    if pd.isna(x) or str(x).strip().lower() in UNASSIGNED:
        return "unassigned"
    s = re.sub(r"\s*\(.*\)$", "", str(x).strip())     # 'GI.P10 (GI.Pa)' -> 'GI.P10'
    return s


def genotype(x):
    """Drop a variant suffix to leave the genotype: 'GII.4_Sydney' -> 'GII.4'."""
    s = norm(x)
    return s if s == "unassigned" else re.sub(r"_.*$", "", s)


def variant(x):
    """Just the variant part, lowercased: 'GII.4_Sydney' -> 'sydney'."""
    s = norm(x)
    if s == "unassigned" or "_" not in s:
        return ""
    return s.split("_", 1)[1].strip().lower()


def rivm_variant(x):
    """RIVM subtype -> comparable variant: 'Sydney_2012' -> 'sydney'.

    RIVM appends the emergence year; calicityper's reference labels do not, so the
    year is stripped before comparing.
    """
    if pd.isna(x) or str(x).strip().lower() in UNASSIGNED:
        return ""
    s = str(x).strip().lower()
    return re.sub(r"_\d{4}[a-z]?$", "", s)


NOROVIRUS_RE = re.compile(r"^caliciviridae[ _]norovirus[ _](g[ivx]+)$", re.I)


def genogroup_label(x):
    """Reduce either tool's genogroup to one comparable label.

    RIVM reports it in the `BLAST` column as a lineage string
    ("Caliciviridae Norovirus GII", "Caliciviridae Sapovirus GI",
    "Caliciviridae Lagovirus unclassified"). calicityper reports the subject
    prefix of its BLAST database, which is a bare genogroup for norovirus
    ("GII") and an underscored lineage for everything else
    ("Caliciviridae_Sapovirus_GI"). Both collapse to the same form here:
    norovirus -> "GII"; anything else -> "Sapovirus GI", "Lagovirus unclassified".
    """
    if pd.isna(x) or str(x).strip().lower() in UNASSIGNED:
        return "unassigned"
    s = re.sub(r"[_]+", " ", str(x).strip())
    s = re.sub(r"\s+", " ", s)
    m = NOROVIRUS_RE.match(s)
    if m:
        return m.group(1).upper()
    if re.fullmatch(r"g[ivx]+", s, re.I):             # calicityper's bare norovirus genogroup
        return s.upper()
    return re.sub(r"^caliciviridae ", "", s, flags=re.I)


def summarise(d, label):
    """Agreement plus a micro-averaged F1 over the types RIVM assigns.

    Treating "unassigned" as a label inflates agreement when both tools decline
    the same sequence, and `agree_of_both` hides every sequence one tool typed
    and the other did not. F1 counts each typed call once:

      TP  both assign and agree        FP  we assign, RIVM differs or declines
      FN  RIVM assigns and we decline or differ

    A disagreement where both assign is therefore both an FP and an FN, which
    is what makes F1 the honest single number here.
    """
    both = d[(d.nt != "unassigned") & (d.rivm != "unassigned")]
    tp = int((both.nt == both.rivm).sum())
    fp = int(((d.nt != "unassigned") & (d.nt != d.rivm)).sum())
    fn = int(((d.rivm != "unassigned") & (d.nt != d.rivm)).sum())
    prec = tp / (tp + fp) if tp + fp else None
    rec = tp / (tp + fn) if tp + fn else None
    f1 = (2 * prec * rec / (prec + rec)) if prec and rec else None
    return dict(
        region=label,
        n=len(d),
        both_assigned=len(both),
        agree_of_both=round(100 * (both.nt == both.rivm).mean(), 2) if len(both) else None,
        rivm_only=int(((d.nt == "unassigned") & (d.rivm != "unassigned")).sum()),
        calicityper_only=int(((d.nt != "unassigned") & (d.rivm == "unassigned")).sum()),
        agree_overall=round(100 * (d.nt == d.rivm).mean(), 2) if len(d) else None,
        precision=round(100 * prec, 2) if prec is not None else None,
        recall=round(100 * rec, 2) if rec is not None else None,
        f1=round(100 * f1, 2) if f1 is not None else None,
    )


def _accepted_only(d, type_col, accepted_col):
    """Blank the type wherever CaliciTyper did not accept it.

    Not a filter: the row stays, so a rejected call is scored as a decline
    (rivm_only) exactly as a pre-0.7.0 withheld call was, rather than vanishing
    from the denominator.
    """
    d = d.copy()
    if accepted_col not in d.columns:
        return d
    ok = d[accepted_col].fillna("").str.startswith("accepted")
    d.loc[~ok, type_col] = None
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("calicityper_tsv")
    ap.add_argument("--rivm", default="validation/results_rivm.csv")
    ap.add_argument("--out", default="results/rivm_comparison.csv")
    ap.add_argument("--accepted-only", action="store_true",
                    help="Score only calls CaliciTyper accepted. Since 0.7.0 the type is "
                         "always reported, with *_accepted saying whether it was accepted, so "
                         "the default scoring credits reported-but-rejected calls as assigned. "
                         "That is the right default (it scores what the tool found), but it is "
                         "not comparable with the pre-0.7.0 tables, which could only score "
                         "accepted calls because nothing else was written out.")
    a = ap.parse_args()

    mine = pd.read_csv(a.calicityper_tsv, sep="\t", dtype=str)
    mine["acc"] = mine["query"].str.split(".").str[0].str.split("|").str[-1]
    rv = pd.read_csv(a.rivm, dtype=str)
    rv["acc"] = rv["name"].str.split(".").str[0]
    m = mine.merge(rv, on="acc", suffixes=("_nt", "_rivm"))
    print(f"matched {len(m)} of {len(mine)} calicityper rows to RIVM records")

    # Genotyping is a NOROVIRUS comparison and must exclude everything else.
    # RIVM's norovirus tool also emits sapovirus genotypes, in nomenclature
    # that collides with norovirus ("GI.1", "GII.3"), so scoring them against a
    # norovirus typer counts 294 correctly-identified sapoviruses as capsid
    # disagreements and drags the capsid figure from ~99% to 87.6%.
    noro = m[m["BLAST"].map(genogroup_label).str.fullmatch(r"G[IVX]+", na=False)]
    print(f"genotyping scored on {len(noro)} records RIVM calls norovirus "
          f"({len(m) - len(noro)} non-norovirus excluded)")

    rows = []
    # Genotype-level agreement: variant suffixes stripped from BOTH sides.
    for label, mc, rc, ms, rs in (
        ("ORF2 capsid", "orf2_type", "capsid type", "orf2_support", "capsid type support"),
        ("ORF1 P-type", "orf1_type", "polymerase type", "orf1_support", "polymerase type support"),
    ):
        d = noro.copy()
        if a.accepted_only:
            d = _accepted_only(d, mc, mc.split("_")[0] + "_accepted")
        d["nt"] = d[mc].map(genotype)
        d["rivm"] = d[rc].map(genotype)
        d = d[(d.nt != "unassigned") | (d.rivm != "unassigned")]
        rows.append(summarise(d, label))
        d[d.nt != d.rivm][["query", "nt", "rivm", ms, rs]].to_csv(
            a.out.replace(".csv", f"_{label.split()[0]}_discordant.csv"), index=False)

    # Genogroup detection, scored against RIVM's BLAST lineage column. This
    # covers every matched record, including the non-norovirus caliciviruses,
    # so it is a wider test than genotyping (which needs a curated alignment).
    d = m.copy()
    d["nt"] = d["genogroup"].map(genogroup_label)
    d["rivm"] = d["BLAST"].map(genogroup_label)
    d = d[(d.nt != "unassigned") | (d.rivm != "unassigned")]
    rows.insert(0, summarise(d, "Genogroup"))
    d[d.nt != d.rivm][["query", "nt", "rivm"]].to_csv(
        a.out.replace(".csv", "_Genogroup_discordant.csv"), index=False)

    # Genogroup restricted to the sequences RIVM calls norovirus at all.
    nv = d[d.rivm.str.fullmatch(r"G[IVX]+", na=False)]
    if len(nv):
        rows.insert(1, summarise(nv, "Genogroup (norovirus only)"))

    # GII.4 variant, compared against RIVM's separate subtype column.
    g = noro[noro["capsid type"].map(genotype) == "GII.4"].copy()
    if len(g):
        if a.accepted_only:
            g = _accepted_only(g, "orf2_type", "orf2_accepted")
        g["nt"] = g["orf2_type"].map(variant).replace("", "unassigned")
        g["rivm"] = g["capsid subtype"].map(rivm_variant).replace("", "unassigned")
        g = g[(g.nt != "unassigned") | (g.rivm != "unassigned")]
        rows.append(summarise(g, "GII.4 variant"))
        g[g.nt != g.rivm][["query", "nt", "rivm", "orf2_support", "capsid subtype support"]].to_csv(
            a.out.replace(".csv", "_GII4variant_discordant.csv"), index=False)

    s = pd.DataFrame(rows)
    s.to_csv(a.out, index=False)
    print(s.to_string(index=False))


if __name__ == "__main__":
    main()
