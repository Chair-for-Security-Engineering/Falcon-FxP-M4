# FN-DSA scheme directories

`fndsa512/` and `fndsa1024/` hold the FN-DSA implementations in pqm4's
`crypto_sign/<scheme>/<impl>/` layout. Everything is committed; there is
no generation step.

- Sources are relative symlinks into the repository root. pqm4 only
  compiles files at the top level of an implementation directory, so
  every `.c`, `.s` and `.S` file is linked on its own, and a new source
  file needs a link in each directory that builds it. The `fpr` link is
  a directory link that keeps `#include "fpr/..."` working.
- Per scheme: `api.h` and `pqm4_glue.c`, the NIST API over `fndsa.h`
  (raw messages, attached signatures).
- Shared: `sysrng_fndsa.c` maps c-fn-dsa's `sysrng()` onto pqm4's
  `randombytes()`.
- Per implementation: `config.mk` sets the backend defines. The knobs
  come from the root Makefile (`make pqm4 CONFIG=...` or `SPLIT=...`),
  which `make pqm4` passes through `MUPQ_EXTRA_CPPFLAGS`.

| dir | backend |
|---|---|
| `ieee-c` | binary64 emulation, portable C |
| `ieee-m4f` | binary64 emulation, Cortex-M4 assembly; assembly ffSampling without Symplecroot/-tree |
| `tw-c` | triple-word (TWFalcon), portable C |
| `tw-m4f` | triple-word, Cortex-M4 assembly; ffSampling in C |
| `fxp3264-c` | 32.64 fixed point, portable C |
| `fxp3264-cm4` | 32.64 fixed point, Cortex-M4 assembly, incl. fused FFT butterflies |

The Cortex-M4 directories also link the representation-independent
assembly (SHAKE, mod q and mod q' arithmetic, codec, Gaussian table) and
the on-board fpr test (`fprtest_wrap.c`, `fprtest_cref_<rep>.c`). All six
link the fpr benchmark (`fpr_bench_wrap.c`, `fpr_bench_cm4.S`). The root
`make pqm4` finds every directory that has a `config.mk`.
