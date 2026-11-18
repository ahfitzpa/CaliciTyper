#!/usr/bin/env python3
"""End-to-end leave-one-out validation of the sapovirus dual type.

The components have been validated separately: VP1 genotyping and RdRp P-type
assignment each have a leave-one-out figure. The composed designation,
GI.4[P1], has not. Its only test was a five-case unit test of the string
format, and the cross-check against Zhao et al. Table 2 compares the genotype
and the P-type in separate columns, never the pair.

This types each reference in BOTH regions with that reference removed from
both alignments and both trees, composes the dual type from the two calls, and
scores it. A dual type is only as good as its weaker half, so the quantity of
interest is the joint outcome, not either region alone.

Only strains present in both reference sets can be dual-typed at all: 51 of
144 VP1 references. Of those, a strain that is the sole member of its type in
either region cannot have that half recovered once it is removed, so the
correct outcome there is a declined call, scored separately as for the
single-region harness.

Scored against two truths:
  * the reference tables, available for every dual-typeable strain;
  * the dual type published in Table 2, for the representatives it names.

    python scripts/dual_type_loo.py --out results/dual_type_loo.csv --threads 16
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from calicityper import align, config as cg, typing as ty   # noqa: E402
from calicityper.typing import dual_type                    # noqa: E402

REGIONS = {"orf1": "RdRp", "orf2": "VP1"}
A = None


def references(region):
    """{accession: (group label, type, ungapped region sequence)}."""
    out = {}
    d = cg.SAPOVIRUS_REF_DIR[region]
    for orf in sorted(d.glob("*.orf")):
        if orf.stem == "all":
            continue
        for r in SeqIO.parse(orf, "fasta"):
            _, typ, acc = r.id.split("|")
            out[acc] = (cg.SAPOVIRUS_PREFIX + orf.stem, typ,
                        str(r.seq).upper().replace("-", ""))
    return out


def published_table2():
    """{accession: published dual type} from Zhao et al. Table 2."""
    p = ROOT / "validation" / "zhao_table2_dual_types.csv"
    out = {}
    if p.exists():
        for r in csv.DictReader(l for l in open(p) if not l.startswith("#")):
            d = (r.get("dual_type") or "").strip()
            if d:
                out[r["accession"].split(".")[0]] = d
    return out


def init(a):
    global A
    A = a


def call_one(region, group, acc, seq):
    """Type one reference against its genogroup with itself removed."""
    ref, alns = align.align_queries(region, group, [(acc, seq)],
                                    exclude=[acc], threads=1)
    res = ty.type_aligned(ref, alns[0], exclude=[acc], nrep=A.bootstraps, seed=1,
                          min_support=70.0, max_refs_per_type=A.max_refs_per_type)
    row = {"query": acc, "genogroup": group, "note": ""}
    ty.apply_result(row, region, res, min_support=70.0)
    return row.get(f"{region}_type", "") or "unassigned"


def work(job):
    acc, vp1, rdrp, singles = job
    g2, t2, s2 = vp1
    g1, t1, s1 = rdrp
    try:
        call2 = call_one("orf2", g2, acc, s2)
        call1 = call_one("orf1", g1, acc, s1)
    except Exception as e:                       # one bad reference must not sink the run
        return dict(accession=acc, error=f"{type(e).__name__}: {e}")

    declined = {"untypeable", "unassigned", ""}
    truth = dual_type(t2, t1)
    ours = dual_type(call2, call1)
    sing2, sing1 = singles

    # A region whose type has no other member cannot be recovered; declining is
    # the correct outcome there, so such strains are scored apart.
    if sing2 or sing1:
        outcome = "singleton (declining is correct)" if ours == "" else "singleton given a label"
    elif ours == "":
        outcome = "no dual type (a region declined)"
    elif ours == truth:
        outcome = "correct"
    else:
        outcome = "wrong"

    return dict(accession=acc, vp1_truth=t2, rdrp_truth=t1, dual_truth=truth,
                vp1_call=call2, rdrp_call=call1, dual_call=ours,
                vp1_ok=call2 == t2, rdrp_ok=call1 == t1,
                vp1_declined=call2 in declined, rdrp_declined=call1 in declined,
                singleton_vp1=sing2, singleton_rdrp=sing1, outcome=outcome)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default="results/dual_type_loo.csv")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--bootstraps", type=int, default=100)
    ap.add_argument("--max-refs-per-type", type=int, default=None)
    a = ap.parse_args()

    v, p = references("orf2"), references("orf1")
    both = sorted(set(v) & set(p))
    cv = Counter(t for _, t, _ in v.values())
    cp = Counter(t for _, t, _ in p.values())
    print(f"VP1 references {len(v)}, RdRp {len(p)}, in both {len(both)}")
    print(f"  dual-typeable fraction of the VP1 set: {100*len(both)/len(v):.1f}%")

    jobs = [(acc, v[acc], p[acc], (cv[v[acc][1]] == 1, cp[p[acc][1]] == 1))
            for acc in both]
    rows = []
    with ProcessPoolExecutor(max_workers=a.threads, initializer=init,
                             initargs=(a,)) as ex:
        for i, r in enumerate(ex.map(work, jobs), 1):
            rows.append(r)
            if i % 10 == 0:
                print(f"  {i}/{len(jobs)}", flush=True)

    errs = [r for r in rows if r.get("error")]
    rows = [r for r in rows if not r.get("error")]
    for r in errs:
        print(f"  ERROR {r['accession']}: {r['error']}")

    scored = [r for r in rows if not (r["singleton_vp1"] or r["singleton_rdrp"])]
    sing = [r for r in rows if r["singleton_vp1"] or r["singleton_rdrp"]]
    c = Counter(r["outcome"] for r in scored)
    produced = [r for r in scored if r["dual_call"]]

    print(f"\n{'=' * 72}\nDUAL-TYPE LEAVE-ONE-OUT\n{'=' * 72}")
    print(f"strains testable (non-singleton in both regions): {len(scored)}")
    print(f"  dual type produced : {len(produced)} "
          f"({100*len(produced)/len(scored):.1f}% call rate)" if scored else "")
    print(f"  correct            : {c['correct']}")
    print(f"  wrong              : {c['wrong']}")
    print(f"  declined by a region: {c['no dual type (a region declined)']}")
    if produced:
        print(f"  accuracy among produced dual types: "
              f"{100*c['correct']/len(produced):.1f}%")
    # which half is the limiting factor
    v_bad = sum(1 for r in scored if not r["vp1_ok"])
    p_bad = sum(1 for r in scored if not r["rdrp_ok"])
    print(f"\n  of {len(scored)}: VP1 half not recovered {v_bad}, "
          f"RdRp half not recovered {p_bad}")
    if sing:
        cs = Counter(r["outcome"] for r in sing)
        print(f"\nsingleton in one or both regions: {len(sing)}")
        for k, n in cs.items():
            print(f"  {k}: {n}")

    t2 = published_table2()
    named = [r for r in rows if r["accession"] in t2]
    if named:
        ok = sum(1 for r in named if r["dual_call"] == t2[r["accession"]])
        print(f"\nagainst Zhao Table 2's published dual type ({len(named)} "
              f"representatives held in both sets): {ok} match")
        for r in named:
            if r["dual_call"] != t2[r["accession"]]:
                print(f"  {r['accession']:10s} published {t2[r['accession']]:18s} "
                      f"ours {r['dual_call'] or '(declined)'}")

    outp = Path(a.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    cols = ["accession", "vp1_truth", "rdrp_truth", "dual_truth", "vp1_call",
            "rdrp_call", "dual_call", "vp1_ok", "rdrp_ok", "vp1_declined",
            "rdrp_declined", "singleton_vp1", "singleton_rdrp", "outcome"]
    with open(outp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    print(f"\nwrote {outp} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
