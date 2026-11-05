/*
 * Benchmark targets for the fpr contract, timed by BENCH in
 * fpr_bench_cm4.S. Operands are passed in globals because an fpr does not
 * always fit in the argument registers. Each op has a baseline wrapper of
 * the same shape with a nop in place of the operation, so
 *     op cycles = BENCH(op wrapper) - BENCH(baseline)
 * leaves only the operation body, as in c-fn-dsa's bench_cm4.
 */

#include "sign_inner.h"

#define NOINLINE   __attribute__((noinline, used))

/* One 'bx lr' under several prototypes (fpr_bench_cm4.S). */
extern void fndsa_bench_nop(void);
extern fpr fndsa_bench_nop_f1(fpr x);
extern fpr fndsa_bench_nop_f2(fpr x, fpr y);
extern void fndsa_bench_nop_sqrt_sigma(uint32_t *d, const uint32_t *x,
	unsigned logn);
#ifndef FNDSA_FPR_RECIP_FALLBACK
extern void fndsa_bench_nop_p3(uint32_t *m, uint32_t *k, const uint32_t *x);
extern void fndsa_bench_nop_p4(uint32_t *z, const uint32_t *a,
	const uint32_t *m, uint32_t k);
#endif

/* volatile so the operation cannot be folded or hoisted. */
volatile fpr fndsa_bw_x, fndsa_bw_y;
volatile fpr fndsa_bw_sink_a, fndsa_bw_sink_b;
volatile int64_t fndsa_bw_sink_i;
volatile int32_t fndsa_bw_i;
volatile unsigned fndsa_bw_e;

/* Polynomial arguments for the routines that take more than four words. */
fpr *fndsa_bw_p00, *fndsa_bw_p01, *fndsa_bw_p10, *fndsa_bw_p11;

/* Written by the recip_x row, read by the mul_recip row (run in that
   order). Not volatile, as at the real call sites. */
fpr_rcp fndsa_bw_rcp;

/* --- calibration baselines ---------------------------------------------- */

NOINLINE void
fndsa_bw_none(void)
{
	fndsa_bench_nop();
}

NOINLINE void
fndsa_bw_base1(void)
{
	fndsa_bw_sink_a = fndsa_bench_nop_f1(fndsa_bw_x);
}

/* No-call variant for the inline ops, whose wrappers contain no call. */
NOINLINE void
fndsa_bw_base1nc(void)
{
	fndsa_bw_sink_a = fndsa_bw_x;
}

NOINLINE void
fndsa_bw_base2(void)
{
	fndsa_bw_sink_a = fndsa_bench_nop_f2(fndsa_bw_x, fndsa_bw_y);
}

/* No-call two-operand variant for inline binary operations. */
NOINLINE void
fndsa_bw_base2nc(void)
{
	fpr a = fndsa_bw_x;
	fpr b = fndsa_bw_y;
	(void)b;
	fndsa_bw_sink_a = a;
}

/* Two-store shape for the two-output macros. */
NOINLINE void
fndsa_bw_base2s2(void)
{
	fpr a = fndsa_bench_nop_f2(fndsa_bw_x, fndsa_bw_y);
	fndsa_bw_sink_a = a;
	fndsa_bw_sink_b = a;
}

/* No-call two-operand, two-store shape for inline two-output macros. */
NOINLINE void
fndsa_bw_base2s2nc(void)
{
	fpr a = fndsa_bw_x;
	fpr b = fndsa_bw_y;
	fndsa_bw_sink_a = a;
	fndsa_bw_sink_b = b;
}

/* With the ieee cm4 shim (addsub_cm4.h), copy its register pinning and
   clobbers so the spills around the call cancel out. */
NOINLINE void
fndsa_bw_base_addsub(void)
{
#ifdef FNDSA_FPR_IEEE_ADDSUB_CM4_H__
	fpr t_x = fndsa_bw_x;
	fpr t_y = fndsa_bw_y;
	register uint32_t t_x0 __asm__("r0") = (uint32_t)t_x;
	register uint32_t t_x1 __asm__("r1") = (uint32_t)(t_x >> 32);
	register uint32_t t_y0 __asm__("r2") = (uint32_t)t_y;
	register uint32_t t_y1 __asm__("r3") = (uint32_t)(t_y >> 32);
	__asm__(
		"bl	fndsa_bench_nop"
		: "+r" (t_x0), "+r" (t_x1), "+r" (t_y0), "+r" (t_y1)
		:
		: "r4", "r5", "r6", "r7", "r8",
		  "r10", "r11", "r12", "r14", "cc");
	fndsa_bw_sink_a = (uint64_t)t_x0 | ((uint64_t)t_x1 << 32);
	fndsa_bw_sink_b = (uint64_t)t_y0 | ((uint64_t)t_y1 << 32);
#else
	fpr a = fndsa_bench_nop_f2(fndsa_bw_x, fndsa_bw_y);
	fndsa_bw_sink_a = a;
	fndsa_bw_sink_b = a;
#endif
}

/* These mirror whichever expansion FPR_RECIP_X / FPR_MUL_RECIP have:
   fused (pointer arguments) or the composed fallback (by value). */
NOINLINE void
fndsa_bw_base_recip_x(void)
{
#ifdef FNDSA_FPR_RECIP_FALLBACK
	fpr x = fndsa_bw_x;
	fndsa_bw_rcp = fndsa_bench_nop_f1(x);
#else
	fpr x = fndsa_bw_x;
	fndsa_bench_nop_p3(fndsa_bw_rcp.m, &fndsa_bw_rcp.k, x.w);
#endif
}

NOINLINE void
fndsa_bw_base_mul_recip(void)
{
#ifdef FNDSA_FPR_RECIP_FALLBACK
	fpr a = fndsa_bw_x;
	fndsa_bw_sink_a = fndsa_bench_nop_f2(a, fndsa_bw_rcp);
#else
	fpr a = fndsa_bw_x;
	fpr d;
	fndsa_bench_nop_p4(d.w, a.w, fndsa_bw_rcp.m, fndsa_bw_rcp.k);
	fndsa_bw_sink_a = d;
#endif
}

/* Same for FPR_SQRT_MUL_INV_SIGMA; the INV_SIGMA load is counted as
   part of the operation. */
NOINLINE void
fndsa_bw_base_sqrt_sigma(void)
{
#ifdef FNDSA_FPR_SQRT_MUL_INV_SIGMA_FALLBACK
	fpr x = fndsa_bw_x;
	fpr t = fndsa_bench_nop_f1(x);
	fndsa_bw_sink_a = fndsa_bench_nop_f2(t, x);
#else
	fpr x = fndsa_bw_x;
	fpr d;
	fndsa_bench_nop_sqrt_sigma(d.w, x.w, fndsa_bw_e);
	fndsa_bw_sink_a = d;
#endif
}

/* --- out-of-line operations --------------------------------------------- */

NOINLINE void
fndsa_bw_add(void)
{
	fndsa_bw_sink_a = fpr_add(fndsa_bw_x, fndsa_bw_y);
}

NOINLINE void
fndsa_bw_mul(void)
{
	fndsa_bw_sink_a = fpr_mul(fndsa_bw_x, fndsa_bw_y);
}

NOINLINE void
fndsa_bw_div(void)
{
	fndsa_bw_sink_a = fpr_div(fndsa_bw_x, fndsa_bw_y);
}

NOINLINE void
fndsa_bw_sqrt(void)
{
	fndsa_bw_sink_a = fpr_sqrt(fndsa_bw_x);
}

NOINLINE void
fndsa_bw_sqrt_sigma(void)
{
	fpr x = fndsa_bw_x;
	fpr d;
	FPR_SQRT_MUL_INV_SIGMA(d, x, fndsa_bw_e);
	fndsa_bw_sink_a = d;
}

/* --- inline operations -------------------------------------------------- */

NOINLINE void
fndsa_bw_sub(void)
{
	fndsa_bw_sink_a = fpr_sub(fndsa_bw_x, fndsa_bw_y);
}

NOINLINE void
fndsa_bw_sqr(void)
{
	fndsa_bw_sink_a = fpr_sqr(fndsa_bw_x);
}

NOINLINE void
fndsa_bw_inv(void)
{
	fndsa_bw_sink_a = fpr_inv(fndsa_bw_x);
}

/* Reciprocal and multiply-by-reciprocal, as used in the LDL. */
NOINLINE void
fndsa_bw_recip_x(void)
{
	fpr x = fndsa_bw_x;
	FPR_RECIP_X(fndsa_bw_rcp, x);
}

NOINLINE void
fndsa_bw_mul_recip(void)
{
	fpr a = fndsa_bw_x;
	fpr d;
	FPR_MUL_RECIP(d, a, fndsa_bw_rcp);
	fndsa_bw_sink_a = d;
}

NOINLINE void
fndsa_bw_neg(void)
{
	fndsa_bw_sink_a = fpr_neg(fndsa_bw_x);
}

NOINLINE void
fndsa_bw_half(void)
{
	fndsa_bw_sink_a = fpr_half(fndsa_bw_x);
}

NOINLINE void
fndsa_bw_double(void)
{
	fndsa_bw_sink_a = fpr_double(fndsa_bw_x);
}

NOINLINE void
fndsa_bw_div2e(void)
{
	fndsa_bw_sink_a = fpr_div2e(fndsa_bw_x, fndsa_bw_e);
}

NOINLINE void
fndsa_bw_mul2e(void)
{
	fndsa_bw_sink_a = fpr_mul2e(fndsa_bw_x, fndsa_bw_e);
}

NOINLINE void
fndsa_bw_of(void)
{
	fndsa_bw_sink_a = fpr_of(fndsa_bw_i);
}

NOINLINE void
fndsa_bw_of32(void)
{
	fndsa_bw_sink_a = fpr_of32(fndsa_bw_i);
}

NOINLINE void
fndsa_bw_rint(void)
{
	fndsa_bw_sink_i = fpr_rint(fndsa_bw_x);
}

NOINLINE void
fndsa_bw_floor(void)
{
	fndsa_bw_sink_i = fpr_floor(fndsa_bw_x);
}

NOINLINE void
fndsa_bw_trunc(void)
{
	fndsa_bw_sink_i = fpr_trunc(fndsa_bw_x);
}

NOINLINE void
fndsa_bw_add_sub(void)
{
	fpr a, b;
	FPR_ADD_SUB(a, b, fndsa_bw_x, fndsa_bw_y);
	fndsa_bw_sink_a = a;
	fndsa_bw_sink_b = b;
}

NOINLINE void
fndsa_bw_fpc_mul(void)
{
	fpr x = fndsa_bw_x;
	fpr y = fndsa_bw_y;
	fpr d_re, d_im;
	FPC_MUL_NOALIAS(d_re, d_im, x, y, y, x);
	fndsa_bw_sink_a = d_re;
	fndsa_bw_sink_b = d_im;
}

/* Complex multiply by a u64 twiddle pair (fxp32.64 only). The b side is
   a pointer to a constant, so the base2s2 baseline slightly
   overcorrects. */
#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 && FNDSA_TWIDDLE_U64
/* GM pair 2 (cos, sin at 45 degrees). */
static const uint64_t fndsa_bw_gm_pair[2] = {
	0xB504F333F9DE6484ull, 0xB504F333F9DE6484ull
};

NOINLINE void
fndsa_bw_fpc_mul_gm(void)
{
	fpr d_re, d_im;
	fxp_fpc fpct_a, fpct_d;
	fpct_a.z.re = fndsa_bw_x;
	fpct_a.z.im = fndsa_bw_y;
	FNDSA_FXP96_CMUL_GM(fpct_d.w, fpct_a.w, fndsa_bw_gm_pair);
	d_re = fpct_d.z.re;
	d_im = fpct_d.z.im;
	fndsa_bw_sink_a = d_re;
	fndsa_bw_sink_b = d_im;
}

NOINLINE void
fndsa_bw_fpc_mul_gmj(void)
{
	fpr d_re, d_im;
	fxp_fpc fpct_a, fpct_d;
	fpct_a.z.re = fndsa_bw_x;
	fpct_a.z.im = fndsa_bw_y;
	FNDSA_FXP96_CMUL_GMJ(fpct_d.w, fpct_a.w, fndsa_bw_gm_pair);
	d_re = fpct_d.z.re;
	d_im = fpct_d.z.im;
	fndsa_bw_sink_a = d_re;
	fndsa_bw_sink_b = d_im;
}

/* Scalar twiddle multiply; one operand, so baseline base1. */
NOINLINE void
fndsa_bw_mul_gm(void)
{
	fndsa_bw_sink_a = FNDSA_FPR_MUL_GM(fndsa_bw_x,
		fndsa_bw_gm_pair[0]);
}
#endif

/* --- polynomial routines with more than four argument words ----------- */

NOINLINE void
fndsa_bw_gram_fft(unsigned logn)
{
	fpoly_gram_fft(logn, fndsa_bw_p00, fndsa_bw_p01,
		fndsa_bw_p10, fndsa_bw_p11);
}

NOINLINE void
fndsa_bw_mulconst(unsigned logn, fpr *a)
{
	fpoly_mulconst(logn, a, fndsa_bw_x);
}

NOINLINE void
fndsa_bw_LDL_fft(unsigned logn, const fpr *g00, fpr *g01, fpr *g11)
{
	/* l10 aliases g01, as in the sampler. */
	fpoly_LDL_fft(logn, g00, g01, g11, (l10_t *)g01);
}
