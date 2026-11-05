# Falcon_FxP_M4: fixed-point FN-DSA signing on the Cortex-M4

This repository is the code artifact of the paper "Flying, Not Floating:
Precise Fixed-Point Falcon on Embedded Processors"
([ePrint 2026/2431](https://eprint.iacr.org/2026/2431)). It signs FN-DSA
(Falcon-512 and Falcon-1024) in 32.64 fixed-point arithmetic, with
constant-time Cortex-M4 assembly, and also builds natively for tests and
precision measurements.

## Summary

Falcon signing needs floating-point arithmetic for the FFT, the LDL tree,
and the Gaussian sampler. Microcontrollers like the Cortex-M4 have no
double-precision FPU, so the floating-point operations are emulated in
software. We replace them with 32.64 fixed point: 96-bit values in three
32-bit words, with 32 integer and 64 fractional bits. Fixed-point
operations are cheaper on the M4, their errors are easier to bound, and
they suit masking better. Key generation rejects keys whose norms
exceed fixed bounds (the gamma checks of ePrint 2026/1531), which bounds
every value in signing, so 32 integer bits are enough.

The repository contains:

- The fixed-point signer, in C and Cortex-M4 assembly
  (`fpr/fxp32.64/`, pqm4 implementation `fxp3264-cm4`). It signs in
  about 45 % fewer cycles than the double-precision emulation of
  c-fn-dsa, at both degrees.
- The precise fixed-point operations (`fpr/fxp32.64/fxp_fused.h`). Most
  of the error in signing comes from small values that keep few
  significant bits in fixed point: the reciprocals in the LDL tree and
  the constants 1/q and 1/sigma. We keep these wider and round once after
  the multiplication, which makes fixed-point signing much more precise
  than double precision.
- The precision framework (`fpr/prec/`, `eval_precision.c`). It runs
  signing in the target arithmetic next to a 200-bit MPFR reference and
  records the error of every sampler input. We used it to find the main
  error sources and to compare all configurations of the paper.
- The signing techniques of the fixed-point Falcon paper
  ([ePrint 2026/1531](https://eprint.iacr.org/2026/1531)) and of
  c-fn-dsa-alt ([ePrint 2026/1915](https://eprint.iacr.org/2026/1915)) as
  compile-time knobs: the split target, Symplecroot and Symplectree, the
  integer Gram matrix, and FFT-new. They work on every backend, and they
  also bring double precision above the 2^64 target.
- Gamma bounds for Falcon-1024 (`kgen_ntru.c`, tables in
  `results/paper/gamma/`). They accept 72 % of keys, where reusing the
  Falcon-512 bounds accepts 8.5 %.
- The constant-time checks (`scripts/ct-*.sh`), and the KATs over the
  ieee, tw, and fxp32.64 backends and their knobs (`scripts/kat-matrix.sh`).

The cost is memory: the 12-byte values raise the signing RAM from 59n
bytes in c-fn-dsa to 83n bytes.

**Concurrent work.** Thomas Pornin's
[c-fn-dsa-fx32.64](https://github.com/pornin/c-fn-dsa-fx32.64), published
around the same time as the paper, puts the same 32.64 arithmetic on the
structure of c-fn-dsa-alt, at about a third of our signing RAM and with
somewhat fewer cycles on the Cortex-M4. The two were developed
independently, and the paper does not yet discuss it. Our contribution is
on the arithmetic side: the precise operations, identified error sources,
the measured precision of different implementation techniques, and
integer-bit budgets at both degrees.

The code is a fork of [c-fn-dsa](https://github.com/pornin/c-fn-dsa)
(Thomas Pornin) @ 3aa5048. We made its floating-point layer (`fpr`) a
compile-time backend, so all arithmetics run the same signing code. That
makes them directly comparable: in speed on a NUCLEO-L4R5ZI board through
[pqm4](https://github.com/mupq/pqm4), in precision against the MPFR
reference, and in correctness with the c-fn-dsa KATs.

## Backends

| backend | arithmetic | native | pqm4 dirs |
|---|---|---|---|
| `fxp32.64` | 32.64 fixed point, 12-byte words (ours) | `BACKEND=fxp32.64` (default) | `fxp3264-c`, `fxp3264-cm4` |
| `ieee` | IEEE-754 binary64 emulated in integers (c-fn-dsa) | `BACKEND=ieee` | `ieee-c`, `ieee-m4f` |
| `tw` | triple-word arithmetic, a port of [TWFalcon](https://github.com/NXP-Research/TWFalcon) | `BACKEND=tw` | `tw-c`, `tw-m4f` |
| `prec` | MPFR 200-bit reference next to one embedded target (`PREC_TARGET=fxp32.64\|ieee\|tw\|none`), for the precision measurements | `BACKEND=prec` | none |

The `-c` pqm4 directories are portable C. All Cortex-M4 directories share
c-fn-dsa's assembly for the NTT mod q, encoding, and the base sampler.
They also share the SHAKE code (on the Keccak of Huang et al. by default)
and the mod-q' signature code of c-fn-dsa-alt for the split target.

- `fxp3264-cm4` has our assembly:
  - `sign_fpr_cm4.S`, the fixed-point operations: multiplications,
    complex multiplications, the fused reciprocal and square root, the
    multiplication by 1/q, and the twiddle multiplications,
  - `sign_fpoly.S`, the FFT and inverse FFT butterflies with split and
    merge,
  - inline assembly in `fpr_fxp.h` for addition.

  The signing structure (ffLDL, ffSampling, SamplerZ) and key generation
  stay C, built on these operations.
- `ieee-m4f` has c-fn-dsa's assembly for its floating-point emulation
  and ffSampling.
- `tw-m4f` has TWFalcon's triple-word assembly.

## Requirements

The easiest way to get all tools is [Nix](https://nixos.org/download/),
which is all you need to install. The repository has two shells, both on
the nixpkgs 24.11 commit pinned in `nixpkgs.nix`:

- `shell.nix`, native: clang, gcc, make, MPFR, GMP, Boehm GC, Valgrind, and
  rsync,
- `shell-cm4.nix`, Cortex-M4: the Arm GNU toolchain 13.3 (binutils 2.42),
  OpenOCD, stlink, rsync, and Python 3 with pyserial and tqdm, plus the
  host gcc and gcov.

The Makefile enters the shells by itself: `shell.nix` when MPFR, GMP, or
the Boehm GC is missing (fxp32.64 and prec), and `shell-cm4.nix` for the
pqm4 targets when the Arm toolchain or the Python modules are missing. The
ieee and tw builds use clang from your PATH. Without it, run them inside
the native shell. The first run downloads nixpkgs and the toolchain,
which takes a while. To enter a shell by hand:

    nix-shell -I nixpkgs=./nixpkgs.nix
    nix-shell -I nixpkgs=./nixpkgs.nix shell-cm4.nix

The `-I` points nix-shell at the pinned nixpkgs. Without it, nix-shell
warns when no nixpkgs channel is set up.

Without Nix, install the same packages. The TWFalcon assembly needs
binutils 2.42, as 2.44 rejects it. The ieee and tw backends need only
clang (or another compiler, with `CC=gcc`), fxp32.64 also needs MPFR and
GMP for its operation test, and prec needs MPFR, GMP, and the Boehm GC.

The Python scripts in `scripts/` use only the standard library. The
board runs go through pqm4's host scripts, which need pyserial and tqdm
(both in `shell-cm4.nix`). The constant-time scripts also call
`objdump` and `addr2line` of the Arm toolchain, `gcov`, and `rsync`.
These scripts enter the right shell by themselves.

The board runs need a NUCLEO-L4R5ZI connected over USB, and access to
its ST-Link probe and serial port (`/dev/ttyACM*` on Linux, via the udev
rules of stlink or OpenOCD and usually the `dialout` group). Other pqm4
platforms build with `PLATFORM=<name>`, but we only ran the NUCLEO.

## Quick start

Run the native tests first. The default backend is fxp32.64 in the
recommended configuration:

    make test

This builds and runs two tests, in about 10 seconds:

- `test_fpr` checks every fpr operation against MPFR,
- `test_fndsa` runs the c-fn-dsa self-tests and KATs.

A failure stops `make` with an error, and a pass ends with
`Test KAT: ... done.`

The other backends:

    make BACKEND=ieee test
    make BACKEND=tw test
    make BACKEND=prec PREC_TARGET=fxp32.64 test

The build tracks the backend and the knobs, so switching between them
rebuilds what is needed. A knob goes on the command line
(`make SPLIT=0 test`), and a whole paper configuration is selected with
`CONFIG=<preset>` (`make CONFIG=fxp-nonsplit test`, see below).

On the board, build the pqm4 images and run the tests. `CONFIG` selects
one implementation and its knobs. Without it, every implementation is
built for both schemes:

    make pqm4 CONFIG=fxp-rec                              # about 1 minute
    make pqm4-test CONFIG=fxp-rec PQM4_UART=/dev/ttyACM0  # about 3 minutes

`pqm4-test` prints `PASS` for Falcon-512 and Falcon-1024.

Then measure:

    make pqm4-signbench-run CONFIG=fxp-rec KEYS=2 ITER=25   # sign cycles
    make pqm4-run CONFIG=fxp-rec BIN=speed ITER=10          # keygen/sign/verify
    make pqm4-run CONFIG=fxp-rec BIN=stack                  # stack use
    make pqm4-bench-run BENCH_IMPL=fxp3264-cm4              # per-operation cycles

Each of these takes a few minutes. The signbench and speed runs print
the cycles of every key generation, signature, and verification, not a
summary. For
`fxp-rec`, signing takes about 9.2 M cycles at Falcon-512 and 19.5 M at
Falcon-1024, and verification about 0.33 M. `BIN=stack` prints the stack
use of key generation, signing, and verification.

Notes:

- `PQM4_UART` defaults to `/dev/ttyACM0`.
- With several probes attached, also set `STLINK_SERIAL`.
- `SCHEME=512` or `SCHEME=1024` restricts the build to one scheme.
  `pqm4-run` runs Falcon-512 unless `SCHEME=1024` is set.
- The full target list is at the top of the `Makefile`.

## Default setup

The knob defaults are in `knobs.h`. Each knob is also a Makefile
variable that becomes a `-D` define.

- fxp32.64 builds the recommended configuration of the paper:
  - the split target with an integer Gram matrix,
  - Symplecroot,
  - the fused reciprocal and square root,
  - the wide 1/q,
  - L10 in 64 bits,
  - ties-up rounding,
  - the single-rounding complex multiplication.

  On `fxp3264-cm4` the FFT also runs in assembly (`FFT_ASM`).
- ieee and tw keep the structure of c-fn-dsa: non-split, without
  Symplecroot.
- Every backend runs the key generation gamma checks, and its Cortex-M4
  assembly directory uses the Keccak of Huang et al.

| knob | paper term | default fxp32.64 | default ieee, tw |
|---|---|---|---|
| `SPLIT` | split target, signature computed mod q' | 1 | 0 |
| `SYMP_ROOT` | Symplecroot, D11 = q^2/G00 at the root | 1 | 0 |
| `SYMP_TREE` | Symplectree (excludes `SYMP_ROOT`) | 0 | 0 |
| `INT_GRAM` | integer Gram (needs `SPLIT`) | = `SPLIT` | = `SPLIT` |
| `FFT_NEW` | FFT-new of c-fn-dsa-alt | 0 | 0 |
| `RECIP_FUSED` | fused reciprocal | 1 | n/a |
| `SQRT_FUSED` | fused square root with wide 1/sigma | 1 | n/a |
| `INVQ_WIDE` | wide 1/q | 1 | n/a |
| `L10_64` | L10 in 64 bits | 1 with fused reciprocal and Symplecroot/-tree | n/a |
| `TIES_UP` | ties-up rounding | 1 | n/a |
| `CMUL_FUSED` | single-rounding complex multiplication | 1 | n/a |
| `TWIDDLE_U64` | unsigned 64-bit twiddle table | 1 (0 under `prec`) | n/a |
| `FFT_ASM` | fused FFT butterflies in assembly | 1 on `fxp3264-cm4` | n/a |
| `SHAKE_HUANG` | Keccak of Huang et al.; 0 = c-fn-dsa's own | 1 on `fxp3264-cm4` | 1 on `ieee-m4f`, `tw-m4f` |
| `GAMMA` | key generation gamma checks | 1 | 1 |

Without the split, fxp32.64 multiplies the target by 1/q before the
basis (scale-first), so that it fits in 32 integer bits. ieee and tw keep
the multiply-first order of c-fn-dsa.

`CONFIG=<preset>` (from `configs.mk`) sets the backend, the pqm4
directory, and the knobs of one row of the paper's sign-cycle table. It
works for native builds, `make precision`, and the pqm4 targets.

| preset | paper row |
|---|---|
| `ieee-base-origshake` | ieee baseline, original SHAKE |
| `ieee-base` | ieee baseline (non-split) |
| `ieee-base-fftnew` | ieee baseline, FFT-new |
| `ieee-split` | ieee split, int. Gram |
| `ieee-rec` | ieee split, int. Gram, Symplecroot (recommended ieee) |
| `ieee-rec-fftnew` | `ieee-rec` with FFT-new |
| `fxp-split` | fxp32.64 split, fused recip./sqrt, wide 1/q, Symplecroot, int. Gram |
| `fxp-unfused` | `fxp-split` with the unfused reciprocal |
| `fxp-rec` | `fxp-split` with ties-up (recommended, the default build) |
| `fxp-rec-fftnew` | `fxp-rec` with FFT-new |
| `fxp-nonsplit` | fxp32.64 non-split, fused recip./sqrt, wide 1/q, Symplecroot |
| `fxp-nonsplit-tiesup` | `fxp-nonsplit` with ties-up |
| `fxp-nonsplit-tiesup-fftnew` | `fxp-nonsplit-tiesup` with FFT-new |
| `ieee-c`, `fxp-c` | portable C comparison in the text |
| `tw-m4f`, `tw-c` | triple-word port mentioned in the text |

## Reproducing the paper

**Per-operation cycles** (the per-operation table): the table comes from
the fpr benchmark (`fprbench.c`) of the ieee baseline and of the fxp32.64
split row with round to nearest even:

    make pqm4-bench-run CONFIG=ieee-base
    make pqm4-bench-run CONFIG=fxp-split L10_64=0

Each takes under a minute. Our output is in `results/paper/fprbench/`, and
`results/paper/README.md` says which rows the table uses. The values
should match exactly, as they are cycle counts of constant-time code.

**Sign cycles** (the sign-cycle table, 100 keys x 100 signatures): one
row is

    make pqm4-signbench-run CONFIG=fxp-rec KEYS=100 ITER=100

which takes about 2 hours at Falcon-512 and 4 at Falcon-1024 for an fxp
row, and about twice that for an ieee row. The whole table is

    scripts/board.sh
    scripts/board_tables.py results/board/<date>

`board.sh` runs every row of `scripts/board-configs.txt` at both degrees:
the rows of the table, plus the portable C and triple-word rows from the
text. That takes about a week on one board, and a re-run continues where
it stopped. Set `KEYS=2 ITER=25` for a quick pass, or `CONFIG_FILE` for
another list of rows. `board_tables.py` prints the medians and the change
against the ieee baseline. The medians vary by about 0.1 to 0.3 % between
runs, so they should match the paper to within that.

Key generation, verification, and stack come from pqm4's own binaries,
e.g. `make pqm4-run CONFIG=fxp-rec BIN=speed ITER=200` or `BIN=stack`.

**Precision** (the sampler-error table and Qs): a short run of one
configuration, which takes seconds, is

    make precision CONFIG=fxp-rec LOGN=9 RUNS=100 SIGS=10 JOBS=8

It writes `results/precision/<name>.summary.txt`, whose `perkey`
lines give the worst-key exponents: `mu_abs` is -log2(delta_t), and
`sigma_rel` is -log2(delta_sigma). This run gives about 59.5 and 64.3,
against 59.3 and 64.2 in the paper.

`scripts/precision.sh` runs every row of the paper's precision table, on
ieee and fxp32.64, for both degrees, and `scripts/precision_tables.py`
renders the table into `results/precision-paper/precision-table.md`. This
is a stripped-down version of our evaluation, with far fewer runs: by
default 10 keys x 10 signatures per row, which takes minutes. The paper
used 10^5 keys x 100 signatures, about 10^10 sampler calls per row and
degree, each also in 200-bit MPFR. On our server that is about 50
core-hours per row at Falcon-512 and 90 at Falcon-1024, so the paper's
size (`KEYS=100000 SIGS=100`) needs a large many-core machine. With fewer
keys the worst key is less extreme, so the exponents come out higher than
in the paper. The runs are seeded, so the same seed and sizes give the
same exponents.

**Constant time**: there are three checks.

- `scripts/ct-matrix.sh` runs TIMECOP: `test_ct.c` under Valgrind, with
  the signing key marked secret, on every sign-cycle configuration, at
  100 keys x 100 signatures per degree. The default build is gcc `-O3`,
  and `CT_OPT=-O0` gives the paper's second run.
- `scripts/ct-scan-matrix.sh` scans the disassembled Cortex-M4 signing
  images with `scripts/ct_scan.py`, our own heuristic scanner. It fails on
  any division, FPU arithmetic, or indirect call, and on any branch that
  is neither a C condition nor reviewed in `scripts/ct_scan.allow`.
  `scripts/ct_cover.py` checks that the TIMECOP run executed every C
  condition the scan relies on.
- `scripts/ct-board-matrix.sh` times every fpr operation and FFT-domain
  routine on the board on many inputs (`make pqm4-ctbench-run` per
  configuration).

**KATs**: `make test` checks the byte-exact c-fn-dsa KATs on ieee, tw, and
fxp32.64, while the prec backend runs a keygen, sign, and verify round
trip instead. For a few vectors the gamma checks reject the seeded key, so
their output differs from c-fn-dsa. These are listed as expected to differ
and must still verify. `scripts/kat-matrix.sh` runs `make test` over 54
backend and knob combinations. On the board, `make pqm4-fprtest-run`
compares every assembly fpr operation with the portable C bit for bit.

**Gamma tuples** (the Falcon-1024 threshold tables): computed with the
experiments of [fxp-falcon-py](https://github.com/fixed-point-fndsa/fxp-falcon-py).
Our gamma scripts are in a pull request to that repository, at commit
[`0d6aa8b`](https://github.com/Stef-Halmans/fxp-falcon-py/tree/0d6aa8be8ca6c87dd26ed99328394ae568a43444) of our fork:
`experiments/simulate_ntrugen_dist.py` runs the key generations, and
`experiments/ntrugen_gamma_tables.py` renders the tables. The resulting
tables are in `results/paper/gamma/`, and the chosen tuple B is in
`kgen_ntru.c`.

## Repository layout

    *.c, *.h, *_cm4.s   c-fn-dsa with fpr as a backend (fpr/fpr_api.h), knobs.h, configs.mk
    fpr/fxp32.64/       32.64 fixed point: C and Cortex-M4 assembly, fused
                        operations in fxp_fused.h, FFT and twiddles in fxp_fft.h
    fpr/ieee/           binary64 emulation, its Cortex-M4 assembly and ffSampling
    fpr/tw/             triple-word arithmetic (TWFalcon)
    fpr/prec/           MPFR reference backend
    pqm4/               vendored pqm4, our implementations in pqm4/crypto_sign/
    scripts/            board, precision, constant-time, and KAT drivers
    results/paper/      gamma tables and per-operation benchmark logs

## Use of AI

We used Anthropic's Claude models (Claude Opus 4.6 to 5.5 and Claude Fable
5 and 5.1) while developing this code. They helped us set up the build
system and the Nix environments, integrate the precision framework (the
`prec` backend and its drivers), and write the heuristic constant-time
evaluation of the Cortex-M4 binaries (`scripts/ct_scan.py`,
`scripts/ct_cover.py`). We also used them to find bugs, for assisted
coding, and for the documentation.

All ideas are our own. We reviewed and tested all AI-assisted code and
take full responsibility for it.

## Authors and licenses

Stef Halmans and Amit Kumar, Chair for Security Engineering,
Ruhr-Universität Bochum.

Our code is under the license of the Chair (non-commercial use), in
`LICENSES/LicenseRef-RUB-SecEng.txt`. The code we build on keeps its own
license, and `LICENSE` lists which files fall under which license.

- [c-fn-dsa](https://github.com/pornin/c-fn-dsa) (Thomas Pornin), commit
  3aa50485d5f9 (2026-06-02). This repository shares its git history, so
  upstream changes merge after
  `git remote add upstream https://github.com/pornin/c-fn-dsa` with
  `git fetch upstream && git merge upstream/main`. Unlicense.
- `fpr/tw/` comes from
  [TWFalcon](https://github.com/NXP-Research/TWFalcon) (NXP), commit
  6c3145b: `triple_float.{c,h}`, `triple_float_round.c`, and the assembly
  `triple_float_cm4_ct*.s` are copied, while `fpr_tw.h`, `ops_c.c`, and
  `fpr_const_tw.c` port its fpr code to our backend interface. MIT.
- `sign_mq2.c` and `sign_mq2_cm4.s` (signature mod q'), and the FFT-new
  code in `sign_fpoly.c`, come from
  [c-fn-dsa-alt](https://github.com/pornin/c-fn-dsa-alt) @ 1fcdede
  (Thomas Pornin). Unlicense.
- `sha3_huang_cm4.s` is SHAKE glue around the Keccak permutation of
  Huang et al., ["Revisiting Keccak and Dilithium Implementations on
  ARMv7-M"](https://doi.org/10.46586/tches.v2024.i2.1-24), TCHES 2024(2)
  ([code](https://github.com/UIC-ESLAS/Dilithium-Multi-Moduli)), which
  pqm4 ships in `pqm4/common/keccakf1600.S`.
- `pqm4/` is a pruned copy of [pqm4](https://github.com/mupq/pqm4) with
  its own licenses (see `pqm4/VENDORED`).
