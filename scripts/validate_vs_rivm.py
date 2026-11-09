#!/usr/bin/env python3
"""Compare norotyper output with RIVM Norovirus Typing Tool calls.

Needs validation/results_rivm.csv (not redistributed here -- see
validation/README.md).

    norotyper rivm_test_set.fasta -o norotyper_out.tsv --threads 16
    python scripts/validate_vs_rivm.py norotyper_out.tsv

Label conventions differ and must be reconciled before anything is comparable:

* RIVM splits GII.4 across two columns: `capsid type` is plain "GII.4" for all
  1,868 GII.4 records, and the variant lives in `capsid subtype`
  ("Sydney_2012", "Den_Haag_2006b", ...).
* norotyper carries the variant inside the type itself ("GII.4_Sydney"),
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
    """RIVM/norotyper label -> comparable genotype, or 'unassigned'."""
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

    RIVM appends the emergence year; norotyper's reference labels do not, so the
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
    "Caliciviridae Lagovirus unclassified"). norotyper reports the subject
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
    if re.fullmatch(r"g[ivx]+", s, re.I):             # norotyper's bare norovirus genogroup
        return s.upper()
    return re.sub(r"^caliciviridae ", "", s, flags=re.I)


def summarise(d, label):
    both = d[(d.nt != "unassigned") & (d.rivm != "unassigned")]
    return dict(
        region=label,
        n=len(d),
        both_assigned=len(both),
        agree_of_both=round(100 * (both.nt == both.rivm).mean(), 2) if len(both) else None,
        rivm_only=int(((d.nt == "unassigned") & (d.rivm != "unassigned")).sum()),
        norotyper_only=int(((d.nt != "unassigned") & (d.rivm == "unassigned")).sum()),
        agree_overall=round(100 * (d.nt == d.rivm).mean(), 2) if len(d) else None,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("norotyper_tsv")
    ap.add_argument("--rivm", default="validation/results_rivm.csv")
    ap.add_argument("--out", default="results/rivm_comparison.csv")
    a = ap.parse_args()

    mine = pd.read_csv(a.norotyper_tsv, sep="\t", dtype=str)
    mine["acc"] = mine["query"].str.split(".").str[0].str.split("|").str[-1]
    rv = pd.read_csv(a.rivm, dtype=str)
    rv["acc"] = rv["name"].str.split(".").str[0]
    m = mine.merge(rv, on="acc", suffixes=("_nt", "_rivm"))
    print(f"matched {len(m)} of {len(mine)} norotyper rows to RIVM records")

    rows = []
    # Genotype-level agreement: variant suffixes stripped from BOTH sides.
    for label, mc, rc, ms, rs in (
        ("ORF2 capsid", "orf2_type", "capsid type", "orf2_support", "capsid type support"),
        ("ORF1 P-type", "orf1_type", "polymerase type", "orf1_support", "polymerase type support"),
    ):
        d = m.copy()
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
    g = m[m["capsid type"].map(genotype) == "GII.4"].copy()
    if len(g):
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
