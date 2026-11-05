#ifndef FNDSA_KNOBS_H__
#define FNDSA_KNOBS_H__

/*
 * Backend selection and build knobs. Each one can be set with -D; the
 * defaults below are the upstream structure for ieee and tw, and the
 * recommended configuration for fxp32.64. Preprocessor only, so the
 * assembly includes it too.
 */

/* Representations of fpr. */
#define FNDSA_FPR_REP_IEEE      1   /* IEEE-754 binary64 in a uint64_t */
#define FNDSA_FPR_REP_TW        2   /* triple-word (TWFalcon) */
#define FNDSA_FPR_REP_FXP3264   3   /* fixed point, 32 int + 64 frac bits */
#define FNDSA_FPR_REP_PREC      4   /* MPFR reference + embedded target,
                                       native analysis only */

/* Implementations of the out-of-line operations. */
#define FNDSA_FPR_IMPL_C     1   /* portable C */
#define FNDSA_FPR_IMPL_CM4   2   /* ARM Cortex-M4 assembly */

#ifndef FNDSA_ASM_CORTEXM4
#define FNDSA_ASM_CORTEXM4   0
#endif

#ifndef FNDSA_FPR_REP
#define FNDSA_FPR_REP   FNDSA_FPR_REP_IEEE
#endif

#ifndef FNDSA_FPR_IMPL
#if FNDSA_ASM_CORTEXM4
#define FNDSA_FPR_IMPL   FNDSA_FPR_IMPL_CM4
#else
#define FNDSA_FPR_IMPL   FNDSA_FPR_IMPL_C
#endif
#endif

/* 1 when the arithmetic is fxp32.64: the backend itself, or the target
   that prec evaluates (FNDSA_PREC_TARGET, set by the Makefile for every
   unit of a prec build). */
#ifdef FNDSA_PREC_TARGET
#define FNDSA_FXP3264   (FNDSA_PREC_TARGET == FNDSA_FPR_REP_FXP3264)
#else
#define FNDSA_FXP3264   (FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264)
#endif

/* Split target: only t_frac enters ffSampling, and the signature is
   computed mod q' from z and the integer residues. */
#ifndef FNDSA_SIGN_SPLIT
#define FNDSA_SIGN_SPLIT   FNDSA_FXP3264
#endif

/* Symplectree: the right ffLDL subtree is derived from the left one. */
#ifndef FNDSA_SYMP_TREE
#define FNDSA_SYMP_TREE   0
#endif

/* Symplecroot: d11 = q^2/g00 at the tree root. */
#ifndef FNDSA_SYMP_ROOT
#define FNDSA_SYMP_ROOT   (FNDSA_FXP3264 && !FNDSA_SYMP_TREE)
#endif

/* Exact integer Gram matrix, via NTTs mod q and q'. */
#ifndef FNDSA_SIGN_INT_GRAM
#define FNDSA_SIGN_INT_GRAM   FNDSA_SIGN_SPLIT
#endif

/* Compact FFT on self-adjoint halves (c-fn-dsa-alt). */
#ifndef FNDSA_FFT_NEW
#define FNDSA_FFT_NEW   0
#endif

/* fxp32.64 only: fused reciprocal-multiply at the LDL divisions. */
#ifndef FNDSA_SIGN_RECIP_FUSED
#define FNDSA_SIGN_RECIP_FUSED   FNDSA_FXP3264
#endif

/* fxp32.64 only: sqrt(x)*INV_SIGMA at the sampler leaves, rounded once. */
#ifndef FNDSA_SIGN_SQRT_FUSED
#define FNDSA_SIGN_SQRT_FUSED   FNDSA_FXP3264
#endif

/* fxp32.64 only: 1/q as a normalized 64-bit word. */
#ifndef FNDSA_SIGN_INVQ_WIDE
#define FNDSA_SIGN_INVQ_WIDE   FNDSA_FXP3264
#endif

/* fxp32.64 only: non-root l10 stored in Q0.63 format (8 bytes). */
#ifndef FNDSA_SIGN_L10_64
#define FNDSA_SIGN_L10_64   (FNDSA_SIGN_RECIP_FUSED \
	&& (FNDSA_SYMP_ROOT || FNDSA_SYMP_TREE))
#endif

/* fxp32.64 only: products rounded half-up instead of ties to even. */
#ifndef FNDSA_FXP_TIES_UP
#define FNDSA_FXP_TIES_UP   FNDSA_FXP3264
#endif

/* fxp32.64 only: complex products rounded once. */
#ifndef FNDSA_FXP_CMUL_FUSED
#define FNDSA_FXP_CMUL_FUSED   1
#endif

/* fxp32.64 only: twiddles stored as unsigned 64-bit magnitudes. Not
   under prec, which reads the target's table as fpr values. */
#ifndef FNDSA_TWIDDLE_U64
#ifdef FNDSA_PREC_TARGET
#define FNDSA_TWIDDLE_U64   0
#else
#define FNDSA_TWIDDLE_U64   FNDSA_FXP3264
#endif
#endif

/* fxp32.64 on Cortex-M4: FFT butterflies in assembly (no effect in C). */
#ifndef FNDSA_FFT_ASM
#define FNDSA_FFT_ASM   (FNDSA_FXP3264 && FNDSA_TWIDDLE_U64 \
	&& FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4)
#endif

/* Cortex-M4: SHAKE on the Keccak permutation of Huang et al.
   (pqm4/common/keccakf1600.S); 0 uses c-fn-dsa's own sha3_cm4.s. */
#ifndef FNDSA_SHAKE_HUANG_BUF
#define FNDSA_SHAKE_HUANG_BUF   FNDSA_ASM_CORTEXM4
#endif

/* Keygen gamma checks (Algorithm 2 of ePrint 2026/1531). */
#ifndef FNDSA_KGEN_GAMMA
#define FNDSA_KGEN_GAMMA   1
#endif

/* The SIMD and RISC-V paths treat fpr as a native double: ieee only,
   not in prec builds and not under the Symplectree. Split and
   Symplecroot keep AVX2 and SSE2, but not NEON and RV64D. */
#if FNDSA_FPR_REP != FNDSA_FPR_REP_IEEE || defined FNDSA_PREC_TARGET \
	|| FNDSA_SYMP_TREE
#ifndef FNDSA_AVX2
#define FNDSA_AVX2   0
#endif
#ifndef FNDSA_SSE2
#define FNDSA_SSE2   0
#endif
#endif
#if FNDSA_FPR_REP != FNDSA_FPR_REP_IEEE || defined FNDSA_PREC_TARGET \
	|| FNDSA_SIGN_SPLIT || FNDSA_SYMP_ROOT || FNDSA_SYMP_TREE
#ifndef FNDSA_NEON
#define FNDSA_NEON   0
#endif
#ifndef FNDSA_RV64D
#define FNDSA_RV64D   0
#endif
#endif

#if FNDSA_SIGN_INT_GRAM && !FNDSA_SIGN_SPLIT
#error "FNDSA_SIGN_INT_GRAM requires FNDSA_SIGN_SPLIT"
#endif
#if FNDSA_SYMP_TREE && FNDSA_SYMP_ROOT
#error "FNDSA_SYMP_TREE and FNDSA_SYMP_ROOT exclude each other"
#endif
#if !FNDSA_FXP3264 && (FNDSA_SIGN_RECIP_FUSED || FNDSA_SIGN_SQRT_FUSED \
	|| FNDSA_SIGN_INVQ_WIDE || FNDSA_SIGN_L10_64 || FNDSA_FXP_TIES_UP \
	|| FNDSA_TWIDDLE_U64)
#error "this knob needs the fxp32.64 representation"
#endif
#if FNDSA_TWIDDLE_U64 && defined FNDSA_PREC_TARGET
#error "FNDSA_TWIDDLE_U64 is not supported by the prec backend"
#endif
#if FNDSA_FFT_ASM && !(FNDSA_FXP3264 && FNDSA_TWIDDLE_U64)
#error "FNDSA_FFT_ASM needs fxp32.64 with FNDSA_TWIDDLE_U64"
#endif
#if FNDSA_SHAKE_HUANG_BUF && !FNDSA_ASM_CORTEXM4
#error "FNDSA_SHAKE_HUANG_BUF requires Cortex-M4 assembly"
#endif

#endif
