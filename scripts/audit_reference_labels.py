#!/usr/bin/env python3
"""Audit every reference's label against the sequence it is attached to.

Two mislabelled references have already cost real accuracy: LC342058.1, a GI
genome filed into the GII RdRp set as GII.PNA4 because one cell of the
published table reads GII where the strain is GI, and before that two GII RdRp
references whose typing region came from the wrong CDS. Both were found by
accident. This looks for the rest on purpose.

Three checks per reference, all distance-based, so a wrong label shows up as a
sequence that does not belong where its label puts it:

  genogroup   the genogroup in the label matches the set the record sits in
  own type    its distance to the other members of its own type
  nearest     its distance to the nearest member of any OTHER type

A reference whose nearest neighbour of a different type is closer than its own
type's members is reported: either the label is wrong, the region was extracted
from the wrong place, or the type is genuinely paraphyletic. The report says
which references to look at; it does not change anything.

    python scripts/audit_reference_labels.py [--distance ml|nj|p]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from calicityper import align, config as cg, putative_novel as pn   # noqa: E402


def sets_to_audit():
    for region, label in (("orf1", "RdRp"), ("orf2", "VP1")):
        d = cg.REF_DIR[region]
        for orf in sorted(d.glob("*.orf")):
            if orf.stem == "all":
                continue
            yield region, label, orf.stem


def audit(region, group, measure):
    """Rows of findings for one region/genogroup reference set."""
    ref = align.refset(region, group)
    names = [n for n in ref.names]
    idx = [i for i, n in enumerate(names)]
    # only the records the .orf file assigns to THIS group
    own = {}
    for line in open(cg.REF_DIR[region] / f"{group}.orf"):
        if line.startswith(">"):
            g, t, acc = line[1:].strip().split("|")[:3]
            own[acc] = (g, t)
    keep = [i for i in idx if names[i] in own]
    if len(keep) < 2:
        return []
    codes = ref.codes[keep]
    labels = [own[names[i]] for i in keep]
    accs = [names[i] for i in keep]

    D = np.zeros((len(keep), len(keep)))
    for a in range(len(keep)):
        D[a], _ = pn.msa_distances(codes[a], codes)
    D = np.nan_to_num(D, nan=1.0)
    np.fill_diagonal(D, 0.0)
    if measure != "p":
        recs = [(accs[a], "".join(np.asarray(list(ref.seqs[keep[a]])))) for a in range(len(keep))]
        D = pn.as_patristic(D, measure, records=recs, alphabet="nt")

    rows = []
    for a, acc in enumerate(accs):
        g_label, t_label = labels[a]
        problems = []
        if g_label != group:
            problems.append(f"label says {g_label}, filed under {group}")
        if t_label.split(".")[0] != group:
            problems.append(f"type {t_label} is not a {group} type")
        same = [b for b in range(len(accs)) if b != a and labels[b][1] == t_label]
        other = [b for b in range(len(accs)) if labels[b][1] != t_label]
        d_own = float(np.min(D[a, same])) if same else None
        d_other = float(np.min(D[a, other])) if other else None
        nearest = accs[int(other[int(np.argmin(D[a, other]))])] if other else ""
        nearest_t = labels[accs.index(nearest)][1] if nearest else ""
        if d_own is not None and d_other is not None and d_other < d_own:
            problems.append(f"closer to {nearest_t} ({nearest}, {d_other:.3f}) than to its own "
                            f"{t_label} ({d_own:.3f})")
        if same and d_own is not None and d_own > 0.5:
            problems.append(f"{d_own:.3f} from the nearest member of its own type")
        if problems:
            rows.append(dict(region=region, set=group, accession=acc, label=t_label,
                             d_own=round(d_own, 4) if d_own is not None else "",
                             d_other=round(d_other, 4) if d_other is not None else "",
                             nearest=nearest, nearest_type=nearest_t,
                             problem="; ".join(problems)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--distance", default="p", choices=["ml", "nj", "p"],
                    help="distance measure; p is enough to spot a misfiled sequence [p]")
    a = ap.parse_args()
    findings = []
    for region, label, group in sets_to_audit():
        try:
            findings += audit(region, group, a.distance)
        except Exception as e:                                   # noqa: BLE001
            print(f"  {region} {group}: could not audit ({type(e).__name__}: {e})", file=sys.stderr)
    if not findings:
        print("no mislabelled references found")
        return 0
    print(f"{len(findings)} reference(s) to look at ({a.distance} distances):\n")
    for f in findings:
        print(f"  {f['region']:5} {f['set']:5} {f['accession']:13} {f['label']:14} {f['problem']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
