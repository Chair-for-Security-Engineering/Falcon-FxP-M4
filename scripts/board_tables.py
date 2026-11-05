#!/usr/bin/env python3
"""Board cycle tables from the logs of scripts/board.sh.

usage: scripts/board_tables.py results/board/<date>

Per scheme, one row per configuration: the sign-cycle median of every
pass, their mean, the spread over all samples, the verify median, and the
change against the ieee baseline (ieee-base for the assembly rows, ieee-c
for the portable C rows), as in the paper's table.
"""
import re
import sys
from pathlib import Path
from statistics import median

out = Path(sys.argv[1])
configs = [l.split("|") for l in (out / "configs.txt").read_text().split("\n") if l]
passes = sorted(p for p in out.glob("p*") if p.is_dir())
schemes = [s for s in ("512", "1024") if any(out.glob(f"p*/*-{s}-*.log"))]
# The ieee baseline row of each implementation.
BASE = {"ieee-m4f": "ieee-base", "fxp3264-cm4": "ieee-base",
        "tw-m4f": "ieee-base", "ieee-c": "ieee-c", "fxp3264-c": "ieee-c",
        "tw-c": "ieee-c"}


def samples(log):
    if not log.exists():
        return None, None, "missing"
    txt = log.read_text(errors="replace")
    sig = [int(x) for x in re.findall(r"^sign cycles:\n(\d+)", txt, re.M)]
    ver = [int(x) for x in re.findall(r"^verify cycles:\n(\d+)", txt, re.M)]
    note = ""
    if re.search(r"#error|Error \d|error:", txt) and not sig:
        note = "build error"
    elif not sig:
        note = "no samples"
    return sig, ver, note


for scheme in schemes:
    rows = []
    for impl, name, knobs in configs:
        meds, allsig, allver, notes = [], [], [], []
        for p in passes:
            sig, ver, note = samples(p / f"{impl}-{scheme}-{name}.log")
            if sig:
                meds.append(median(sig))
                allsig += sig
                allver += ver
            else:
                meds.append(None)
                notes.append(f"{p.name}: {note}")
        rows.append((impl, name, knobs, meds, allsig, allver, notes))

    def mean_med(r):
        m = [x for x in r[3] if x is not None]
        return sum(m) / len(m) if m else None

    ref = {r[1]: mean_med(r) for r in rows if r[1] in BASE.values()}

    print(f"\n## fndsa{scheme} — sign cycles (signbench medians)\n")
    hdr = ["impl", "config", "knobs"] + [p.name for p in passes] + \
        ["mean", "spread", "verify", "vs baseline", "note"]
    print("| " + " | ".join(hdr) + " |")
    print("|" + "---|" * len(hdr))
    for impl, name, knobs, meds, allsig, allver, notes in rows:
        mm = mean_med((impl, name, knobs, meds))
        cells = [impl, name, f"`{knobs}`" if knobs else "(defaults)"]
        cells += [f"{m:,.0f}" if m else "—" for m in meds]
        if mm:
            spread = (max(allsig) - min(allsig)) / median(allsig) * 100
            cells += [f"{mm:,.0f}", f"{spread:.1f} %",
                      f"{median(allver):,.0f}" if allver else "—"]
            b = ref.get(BASE.get(impl))
            cells.append(f"{(mm / b - 1) * 100:+.1f} %" if b else "—")
        else:
            cells += ["—"] * 4
        cells.append("; ".join(notes))
        print("| " + " | ".join(cells) + " |")
