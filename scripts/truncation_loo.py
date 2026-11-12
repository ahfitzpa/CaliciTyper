#!/usr/bin/env python3
"""Truncation leave-one-out harness: how typing behaves on partial sequence.

Every record in the RIVM benchmark is at least 2,000 nt, so nothing in that
benchmark says what happens to a 300 nt amplicon -- which is most of what a
surveillance laboratory actually submits. This harness answers that on the
reference set, where the truth label is known independently of any call.

Each reference ORF is truncated to a series of lengths and then typed against
all OTHER references of its genogroup -- its own accession is dropped from both
the alignment and the tree, exactly as tests/loo_validation.py does.

Fragments are anchored on the TYPING WINDOW, in each record's own coordinates,
not on the ends of the record. A diagnostic amplicon is positioned by primers
that target the typing region, so that is the geometry to simulate; anchoring
on the record's ends instead measures how far each reference happens to extend
beyond the window, which is an artefact of what was submitted to GenBank. The
difference is not cosmetic: the shortest ORF1 references stop 5' of the window,
so a 3'-anchored 200 nt fragment of one covers 21% of the window and would be
scored as a typing failure at every length.

Offsets are included as well, so that long-but-displaced sequences appear in
the curve: a 1 kb fragment that starts 300 nt into region C must fail region C,
and only an offset sweep tests that.

Each job records the call, its bootstrap support, and the coverage of all three
windows WITHOUT applying any coverage threshold. The thresholds are then swept
over the recorded table in analyse(), so one expensive pass answers "what
should CALICITYPER_MIN_WINDOW_COVERAGE be?" rather than one pass per candidate
value.

Scoring is against the record's published type, not against the full-length
call, so a window that is too short to type shows up as a loss of recall and a
window that types wrongly shows up as a loss of precision. The full-length
column is kept alongside, because a record its own reference set cannot recover
at full length should not be counted against truncation.

    python scripts/truncation_loo.py --out results/truncation_loo.csv --threads 16
    python scripts/truncation_loo.py --analyse-only results/truncation_loo.csv
"""
from __future__ import annotations

import argparse
import collections
import sys
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from calicityper import align, config as cg, typing as ty   # noqa: E402

# Lengths to test, in nucleotides. The low end brackets MIN_INPUT_NT (172) and
# the operational windows (172-267); 300-600 is the common amplicon range;
# the high end is whole-ORF. 0 means "no truncation".
LENGTHS = (150, 172, 200, 252, 267, 300, 400, 600, 900, 1200, 0)
# Displacement of the fragment away from the window, in nucleotides: how far
# the amplicon's window-proximal end sits inside the window. A displaced
# fragment can be long and still miss most of what it needs.
OFFSETS = (0, 100, 300)
A = None


def window_anchors(region, group):
    """{accession: (start, end)} of the typing window in each reference's OWN
    ungapped coordinates, read off the reference alignment it is already in."""
    import numpy as np
    ref = align.refset(region, group)
    cols, _ = ty.rivm_window_columns(region, group)
    if len(cols) == 0:
        return {}
    out = {}
    lo, hi = int(cols[0]), int(cols[-1])
    for name, aligned in zip(ref.names, ref.seqs):
        a = np.frombuffer(aligned.encode(), dtype="S1")
        ungapped = np.cumsum(a != b"-") - 1          # column -> index in the record
        present = np.flatnonzero((a != b"-"))
        inside = present[(present >= lo) & (present <= hi)]
        if len(inside) == 0:
            continue                                 # this record does not reach its window
        out[name] = (int(ungapped[inside[0]]), int(ungapped[inside[-1]]) + 1)
    return out


def jobs(regions, groups, lengths, offsets, min_per_type=2):
    """One job per (region, group, record, length, offset)."""
    out = []
    for region in regions:
        for g in groups:
            orf = cg.REF_DIR[region] / f"{g}.orf"
            if not orf.exists():
                print(f"no reference set for {region} {g}, skipping", file=sys.stderr)
                continue
            recs = list(SeqIO.parse(orf, "fasta"))
            n_of_type = collections.Counter(r.id.split("|")[1] for r in recs)
            anchors = window_anchors(region, g)
            for r in recs:
                _, truth, acc = r.id.split("|")[:3]
                if n_of_type[truth] < min_per_type:
                    continue           # a type with one reference cannot be recovered when left out
                if acc not in anchors:
                    print(f"{acc} ({region} {g}) does not reach its typing window, skipping",
                          file=sys.stderr)
                    continue
                seq = str(r.seq).upper()
                for length in lengths:
                    for off in offsets:
                        if length == 0 and off:
                            continue                  # untruncated: offset is meaningless
                        out.append((region, g, acc, truth, seq, anchors[acc], length, off))
    return out


def truncate(seq, region, anchor, length, offset):
    """Cut `seq` to `length` nt, anchored on its typing window.

    `anchor` is (start, end) of the window in this record's coordinates. ORF2
    windows are at the 5' end of VP1, so the fragment starts at the window and
    runs 3-ward; ORF1 windows are at the 3' end of RdRp, so it ends at the
    window and runs 5-ward. `offset` displaces the fragment INTO the window, so
    a positive offset deliberately leaves part of the window uncovered.
    """
    if length == 0:
        return seq
    ws, we = anchor
    if region == "orf2":
        start = min(ws + offset, len(seq))
        return seq[start: start + length]
    end = max(0, we - offset)
    return seq[max(0, end - length): end]


def work(job):
    region, g, acc, truth, seq, anchor, length, off = job
    sub = truncate(seq, region, anchor, length, off)
    row = dict(region=region, group=g, accession=acc, truth=truth, length=length or len(seq),
               offset=off, nt=len(sub), truncated=int(length != 0 or off != 0))
    if len(sub) < cg.MIN_INPUT_NT:
        # the production path rejects these before BLAST; record that, do not type
        row.update(call="rejected", support=0.0, n_sites=0, rejected_at_input=1,
                   window_cov=0.0, window_nt=0, window_width=0, gii4_cov="", note="shorter than MIN_INPUT_NT")
        return row
    row["rejected_at_input"] = 0
    ref, (qaln,) = align.align_queries(region, g, [(acc, sub)], exclude=[acc], threads=1)
    res = ty.type_aligned(ref, qaln, exclude=[acc], nrep=A.bootstraps,
                          min_support=A.min_support, max_refs_per_type=A.max_refs_per_type)
    call = res["type"] or "unassigned"
    covered, width, frac = ty.rivm_window_coverage(region, g, qaln, genotype=res["type"])
    row.update(call=call, support=res["support"], n_sites=res["n_sites"],
               window_nt=covered, window_width=width, window_cov=frac,
               best=res.get("best_type", ""), note=res.get("note", ""))
    # the GII.4 variant window is a third window, scored only where a variant
    # label exists to recover (the suffix is part of the reference label)
    row["gii4_cov"] = ""
    if region == "orf2" and "_" in truth:
        c2, w2, f2 = ty.gii4_variant_window_coverage(qaln, g)
        row["gii4_cov"] = f2
        row["gii4_nt"] = c2
        row["gii4_width"] = w2
    return row


def init(a):
    global A
    A = a


# ------------------------------------------------------------------ analysis

def _apply(df, cov, gii4_cov, min_support):
    """Re-derive what the production path WOULD report at these thresholds."""
    out = df.copy()
    call = out["call"].astype(str)
    withheld = (out["window_cov"] < cov) | (out["support"] < min_support)
    call = call.where(~withheld, "unassigned")
    # variant suffix dropped below the variant-window floor, genotype stem kept
    stem = call.str.split("_").str[0]
    drop = (out["region"] == "orf2") & call.str.contains("_") & \
           (pd.to_numeric(out["gii4_cov"], errors="coerce").fillna(0) < gii4_cov)
    out["reported"] = call.where(~drop, stem)
    out["truth_stem"] = out["truth"].astype(str).str.split("_").str[0]
    return out


def _score(out, genotype_level=True):
    """Micro-averaged precision/recall/F1 over records with a call to make."""
    truth = out["truth"] if genotype_level else out["truth_stem"]
    reported = out["reported"]
    assigned = ~reported.isin(["unassigned", "rejected"])
    correct = assigned & (reported == truth)
    n = len(out)
    prec = 100 * correct.sum() / assigned.sum() if assigned.sum() else float("nan")
    rec = 100 * correct.sum() / n if n else float("nan")
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else float("nan")
    return dict(n=n, assigned=int(assigned.sum()), correct=int(correct.sum()),
                precision=round(prec, 2), recall=round(rec, 2), f1=round(f1, 2),
                wrong=int((assigned & ~correct).sum()))


def analyse(df, sweep=(0.5, 0.6, 0.7, 0.8, 0.9, 1.0), min_support=70.0):
    print("\n=== coverage of the typing window, by submitted length ===")
    t = df[df.rejected_at_input == 0].groupby(["region", "length"]).agg(
        n=("accession", "size"), window_cov=("window_cov", "mean"),
        full_window=("window_cov", lambda s: 100 * (s >= 1.0).mean()),
        ge80=("window_cov", lambda s: 100 * (s >= 0.8).mean())).round(3)
    print(t.to_string())

    print(f"\n=== threshold sweep, genotype level (support >= {min_support}) ===")
    rows = []
    for cov in sweep:
        for region, part in df.groupby("region"):
            s = _score(_apply(part, cov, cov, min_support))
            rows.append(dict(coverage=cov, region=region, **s))
    sw = pd.DataFrame(rows)
    print(sw.to_string(index=False))
    for region, part in sw.groupby("region"):
        best = part.sort_values("f1", ascending=False).iloc[0]
        print(f"  best F1 for {region}: {best.f1} at coverage {best.coverage} "
              f"(precision {best.precision}, recall {best.recall}, {best.wrong} wrong)")

    print("\n=== by length, at the current threshold "
          f"({cg.RIVM_WINDOW_MIN_COVERAGE:.0%}) ===")
    cur = _apply(df, cg.RIVM_WINDOW_MIN_COVERAGE, cg.RIVM_WINDOW_MIN_COVERAGE, min_support)
    per = []
    for (region, length), part in cur.groupby(["region", "length"]):
        per.append(dict(region=region, length=length, **_score(part)))
    print(pd.DataFrame(per).to_string(index=False))

    v = df[(df.region == "orf2") & (df.truth.astype(str).str.contains("_")) & (df.rejected_at_input == 0)]
    if len(v):
        print("\n=== GII.4 variant window: its own sweep ===")
        rows = []
        for cov in sweep:
            s = _score(_apply(v, 0.0, cov, min_support))      # window floor off, variant floor swept
            rows.append(dict(variant_coverage=cov, **s))
        print(pd.DataFrame(rows).to_string(index=False))
        stem = _apply(v, 0.0, 1.01, min_support)
        print(f"  stem-only (variant always withheld): recall {_score(stem)['recall']}")

    print("\n=== calls that are wrong, not merely withheld, at the current threshold ===")
    bad = cur[(~cur.reported.isin(["unassigned", "rejected"])) & (cur.reported != cur.truth)]
    if bad.empty:
        print("  none")
    else:
        print(bad[["region", "group", "accession", "truth", "reported", "support",
                   "length", "offset", "window_cov"]].sort_values(
                       ["region", "accession", "length"]).to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/truncation_loo.csv")
    ap.add_argument("--analyse-only", help="skip typing, analyse an existing CSV")
    ap.add_argument("--bootstraps", type=int, default=100)
    ap.add_argument("--min-support", type=float, default=70.0)
    ap.add_argument("--max-refs-per-type", type=int, default=2)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--regions", default="orf1,orf2")
    ap.add_argument("--groups", default="GI,GII")
    ap.add_argument("--lengths", default=",".join(str(x) for x in LENGTHS))
    ap.add_argument("--offsets", default=",".join(str(x) for x in OFFSETS))
    ap.add_argument("--limit", type=int, default=0, help="subsample jobs, for a smoke test")
    a = ap.parse_args()

    if a.analyse_only:
        analyse(pd.read_csv(a.analyse_only), min_support=a.min_support)
        return 0

    cg.check_data()
    js = jobs(a.regions.split(","), a.groups.split(","),
              tuple(int(x) for x in a.lengths.split(",")),
              tuple(int(x) for x in a.offsets.split(",")))
    if a.limit:
        js = js[:: max(1, len(js) // a.limit)][:a.limit]
    print(f"{len(js)} typing jobs ({a.bootstraps} bootstraps each) on {a.threads} threads", flush=True)
    if a.threads > 1:
        with Pool(a.threads, initializer=init, initargs=(a,)) as p:
            rows = p.map(work, js, chunksize=4)
    else:
        init(a)
        rows = [work(j) for j in js]
    df = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, index=False)
    print(f"wrote {a.out} ({len(df)} rows)")
    analyse(df, min_support=a.min_support)
    return 0


if __name__ == "__main__":
    sys.exit(main())
