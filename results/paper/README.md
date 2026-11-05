# Data behind paper tables

## gamma/

Key generation threshold (gamma) tuples for Falcon-512 and Falcon-1024,
from 10^5 key generations per degree with all four checks disabled, each
tuple replayed on that set. Computed with the experiments of
[fxp-falcon-py](https://github.com/fixed-point-fndsa/fxp-falcon-py), from
our pull request at commit [`0d6aa8b`](https://github.com/Stef-Halmans/fxp-falcon-py/tree/0d6aa8be8ca6c87dd26ed99328394ae568a43444) of our fork
(`experiments/simulate_ntrugen_dist.py` runs the key generations,
`experiments/ntrugen_gamma_tables.py` renders the tables).

- `gamma_tables.md`: the tables in Markdown.
- `gamma_tables.tex`: the same tables in LaTeX.

Per degree there are two tables:

- m budgets and acceptance per tuple: backs the paper's threshold-tuple
  table (thresholds, accepted keys, m budgets) and the m columns of the
  variable-budget table. "twk" and "std" are the split and non-split
  target budgets. The paper's tuples are "paper" at n = 512 and
  "512-literal", A, B and C at n = 1024; "512-parity" is not in the paper.
- rejections per 1000 accepted keys: the tuple B row at n = 1024 is the
  paper's rejection table, and the "attempts" column divided by 1000 gives
  the attempts per key of the threshold-tuple table ("stock" is unfiltered
  key generation).

Tuple B is the one used in `kgen_ntru.c`.

## fprbench/

Device output of `make pqm4-bench-run` (Falcon-512) with CONFIG=ieee-base
and with CONFIG=fxp-split L10_64=0 (round to nearest even, as in the
paper's per-operation table). Each section says how it measures. The
table takes mul from "one direct call from assembly"; fxp32.64 add and
add-sub from "inlined fpr operations", since signing inlines them; the
other columns from the FPC_MUL, FPR_RECIP_X, FPR_MUL_RECIP and
FPR_SQRT_SIGMA rows; and FFT from fpoly_FFT at logn=9.
