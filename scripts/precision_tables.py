#!/usr/bin/env python3
"""The paper's precision table from the runs of scripts/precision.sh.

usage: scripts/precision_tables.py [--dir results/precision-paper]

One row per table row, Falcon-512 / Falcon-1024 in each cell. Exponents
are e = -log2(error) at the worst key, larger is better: delta_t is the
per-key RMSE of the center error, delta_sigma the per-key RMSE of the
relative sigma error. log2 Qs is Theorem 2 of ePrint 2026/1531,
  log2 Qs = 2 e_t - log2(4 n alpha (1 + c)),  alpha = 2 lambda - 1,
with lambda 128 / 256 and c the measured restart rate. Writes
DIR/precision-table.md and prints it.
"""

import argparse
import math
import os
import re

NAME_RE = re.compile(r"^(?P<target>ieee|fxp32\.64).*-logn(?P<logn>\d+)-r\d+"
                     r".*-row(?P<row>\d+)\.summary\.txt$")
PERKEY_RE = re.compile(r"^perkey (mu_abs|sigma_rel)\s+keys=\d+\s+rmse-e:"
                       r" min\(worst\)=([\d.]+)")
VERIFY_RE = re.compile(r"verify: ok=\d+ fail=(\d+).*attempts=(\d+)"
                       r" \(retries=(\d+)\)")
LAMBDA = {9: 128, 10: 256}
LABELS = [
    "baseline (non-split)", "Symplecroot", "split", "split, int. Gram",
    "split, Symplecroot", "+ int. Gram", "+ FFT-new",
    "non-split, fused reciprocal", "split, fused reciprocal",
    "+ fused sqrt", "+ wide 1/q", "+ Symplecroot",
    "+ int. Gram (recommended)",
    "non-split, fused recip./sqrt, wide 1/q, Symplecroot",
    "recommended, FFT-new",
]


def parse(path):
    r = {}
    with open(path) as f:
        for line in f:
            m = PERKEY_RE.match(line)
            if m:
                r["t" if m.group(1) == "mu_abs" else "s"] = float(m.group(2))
            m = VERIFY_RE.search(line)
            if m:
                r["fail"] = int(m.group(1))
                r["c"] = int(m.group(3)) / max(1, int(m.group(2)))
    return r


def qs(e_t, logn, c):
    n = 1 << logn
    return 2 * e_t - math.log2(4 * n * (2 * LAMBDA[logn] - 1) * (1 + c))


def cell(res, row, target, key):
    vals = []
    for logn in (9, 10):
        r = res.get((row, target, logn))
        if r is None or key not in r and key != "qs":
            vals.append("—")
        elif key == "qs":
            vals.append(f"{qs(r['t'], logn, r['c']):.1f}" if "t" in r else "—")
        else:
            vals.append(f"{r[key]:.1f}")
    return " / ".join(vals)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/precision-paper")
    args = ap.parse_args()
    res = {}
    for name in sorted(os.listdir(args.dir)):
        m = NAME_RE.match(name)
        if m:
            res[(int(m.group("row")), m.group("target"),
                 int(m.group("logn")))] = parse(os.path.join(args.dir, name))
    hdr = ["#", "configuration", "ieee -log2 δt", "ieee -log2 δσ",
           "fxp32.64 -log2 δt", "fxp32.64 -log2 δσ", "ieee log2 Qs",
           "fxp32.64 log2 Qs"]
    lines = ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    for i, label in enumerate(LABELS, 1):
        cells = [str(i), label]
        cells += [cell(res, i, t, k) for t in ("ieee", "fxp32.64")
                  for k in ("t", "s")]
        cells += [cell(res, i, t, "qs") for t in ("ieee", "fxp32.64")]
        lines.append("| " + " | ".join(cells) + " |")
    failed = sum(r.get("fail", 0) for r in res.values())
    lines.append("")
    lines.append(f"Falcon-512 / Falcon-1024; {len(res)} runs, "
                 f"{failed} failed verifications.")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(args.dir, "precision-table.md"), "w") as f:
        f.write(text)
    print(text, end="")


if __name__ == "__main__":
    main()
