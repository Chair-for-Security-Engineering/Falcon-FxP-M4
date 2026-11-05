#ifndef FNDSA_FPR_TW_H__
#define FNDSA_FPR_TW_H__

/*
 * Copyright 2025 NXP
 * SPDX-License-Identifier: MIT
 *
 * The tw-specific content of this file is derived from TWFalcon
 * (NXP-Research/TWFalcon, MIT-licensed); the c-fn-dsa code it binds to
 * is public domain. See LICENSES/MIT-TWFalcon.txt.
 */

/* ==================================================================== */
/*
 * 'tw' representation: triple-word arithmetic (TWFalcon).
 *
 * A value is the unevaluated sum of three floats, x[0] + x[1] + x[2],
 * which gives about 66 bits of mantissa with only single-precision
 * hardware, as on the Cortex-M4F. The arithmetic is TWFalcon's
 * triple_float code, unchanged; this header maps it onto the fpr
 * contract, with the constant-time (_ct) variants of the arithmetic.
 */

#include "triple_float.h"

typedef tw_fpr fpr;

/* Three 4-byte floats. */
#define FPR_STRIDE   12

/* ==================================================================== */
/* Constants.                                                           */

/* Constant from three float limbs. */
#define FPR(x0, x1, x2)   ((tw_fpr){{(x0), (x1), (x2)}})

#define FPR_ZERO    ((tw_fpr){{0, 0, 0}})
#define FPR_NZERO   ((tw_fpr){{-0, -0, -0}})
#define FPR_ONE     ((tw_fpr){{1065353216, 0, 0}})

/* ==================================================================== */
/* Constant-time 64-bit shift helpers (same as in the ieee header).     */
static inline uint64_t
fpr_ursh(uint64_t x, int n)
{
#if FNDSA_64 || FNDSA_ASM_CORTEXM4
	return x >> n;
#else
	x ^= (x ^ (x >> 32)) & -(uint64_t)(n >> 5);
	return x >> (n & 31);
#endif
}

static inline int64_t
fpr_irsh(int64_t x, int n)
{
#if FNDSA_64 || FNDSA_ASM_CORTEXM4
	return x >> n;
#else
	x ^= (x ^ (x >> 32)) & -(int64_t)(n >> 5);
	return x >> (n & 31);
#endif
}

static inline uint64_t
fpr_ulsh(uint64_t x, int n)
{
#if FNDSA_64 || FNDSA_ASM_CORTEXM4
	return x << n;
#else
	x ^= (x ^ (x << 32)) & -(uint64_t)(n >> 5);
	return x << (n & 31);
#endif
}

/* ==================================================================== */
/* Conversions and rounding.                                            */

/* Given integer i and scale sc, return i*2^sc (used by tests only). */
#define fpr_scaled   fndsa_fpr_scaled
fpr fpr_scaled(int64_t i, int sc);

/* Integer to fpr. */
#define fpr_of(i)     tw_from_int(i)
/* 32-bit integer to fpr. */
#define fpr_of32(i)   tw_from_int(i)

/* Rounding. tw_*_fast are exact for the small values seen in signing.
   The cm4 rint/floor take the two high limbs in VFP registers. */
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
extern int64_t tw_rint_fast_ct_asm(float x0, float x1);
extern int32_t tw_floor_fast_asm(float x0, float x1);
extern int64_t tw_trunc_full_ct_asm(fpr x);
extern int32_t tw_trunc_fast_ct_asm(fpr x);
static inline int64_t fpr_rint(fpr x)  { return tw_rint_fast_ct_asm(x.x[0], x.x[1]); }
static inline int64_t fpr_floor(fpr x) { return (int64_t)tw_floor_fast_asm(x.x[0], x.x[1]); }
static inline int64_t fpr_trunc(fpr x) { return tw_trunc_full_ct_asm(x); }
static inline int64_t fpr_trunc_fast(fpr x) { return (int64_t)tw_trunc_fast_ct_asm(x); }
#else
static inline int64_t fpr_rint(fpr x)  { return tw_rint_fast(x); }
static inline int64_t fpr_floor(fpr x) { return tw_floor_fast(x); }
static inline int64_t fpr_trunc(fpr x) { return tw_trunc_full(x); }
static inline int64_t fpr_trunc_fast(fpr x) { return tw_trunc_fast(x); }
#endif

/* ==================================================================== */
/* Arithmetic.                                                          */

/* Addition: ops_c.c, or the assembly under cm4. */
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
#define fpr_add   tw_sum_ct_asm
#else
#define fpr_add   fndsa_fpr_add
#endif
fpr fpr_add(fpr x, fpr y);

/* Subtraction. Under cm4, negate the limbs and use the assembly add. */
static inline fpr
fpr_sub(fpr x, fpr y)
{
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
	y.x[0] = -y.x[0];
	y.x[1] = -y.x[1];
	y.x[2] = -y.x[2];
	return fpr_add(x, y);
#else
	return tw_sub_ct(x, y);
#endif
}

/* Combined addition/subtraction: a <- x + y, b <- x - y. */
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
/* tw_add_sub_ct_asm returns a in s0-s2 and b in s3-s5, which needs
   inline asm. */
#define FPR_ADD_SUB(a, b, xx, yy)   do { \
		fpr t_add_sub_x = (xx); \
		fpr t_add_sub_y = (yy); \
		register float t_add_sub_s0 __asm__("s0") = t_add_sub_x.x[0]; \
		register float t_add_sub_s1 __asm__("s1") = t_add_sub_x.x[1]; \
		register float t_add_sub_s2 __asm__("s2") = t_add_sub_x.x[2]; \
		register float t_add_sub_s3 __asm__("s3") = t_add_sub_y.x[0]; \
		register float t_add_sub_s4 __asm__("s4") = t_add_sub_y.x[1]; \
		register float t_add_sub_s5 __asm__("s5") = t_add_sub_y.x[2]; \
		__asm__ volatile( \
			"bl	tw_add_sub_ct_asm\n" \
			: "+t"(t_add_sub_s0), "+t"(t_add_sub_s1), \
			  "+t"(t_add_sub_s2), "+t"(t_add_sub_s3), \
			  "+t"(t_add_sub_s4), "+t"(t_add_sub_s5) \
			: \
			: "lr", "r0", "r1", "r2", "r3", "r4", "r5", "r6", \
			  "r7", "r8", "r9", "r10", "r11", "s6", "s7", "s8", \
			  "s9", "s10", "s11", "s12", "s13", "s14", "s15", "cc"); \
		(a) = (fpr){{t_add_sub_s0, t_add_sub_s1, t_add_sub_s2}}; \
		(b) = (fpr){{t_add_sub_s3, t_add_sub_s4, t_add_sub_s5}}; \
	} while (0)
#else
#define FPR_ADD_SUB(a, b, x, y)   do { \
		fpr t_add_sub_x = (x); \
		fpr t_add_sub_y = (y); \
		(a) = fpr_add(t_add_sub_x, t_add_sub_y); \
		(b) = fpr_sub(t_add_sub_x, t_add_sub_y); \
	} while (0)
#endif

/* Negation. */
static inline fpr
fpr_neg(fpr x)
{
	return (tw_fpr){{-x.x[0], -x.x[1], -x.x[2]}};
}

/* Halving. */
static inline fpr
fpr_half(fpr x)
{
	return (tw_fpr){{0.5f * x.x[0], 0.5f * x.x[1], 0.5f * x.x[2]}};
}

/* Doubling. */
static inline fpr
fpr_double(fpr x)
{
	return (tw_fpr){{2.0f * x.x[0], 2.0f * x.x[1], 2.0f * x.x[2]}};
}

/* Multiplication. */
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
#define fpr_mul   tw_prod_fast_ct_asm
#else
#define fpr_mul   fndsa_fpr_mul
#endif
fpr fpr_mul(fpr x, fpr y);

/* Squaring. */
static inline fpr
fpr_sqr(fpr x)
{
	return fpr_mul(x, x);
}

/* Division. */
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
#define fpr_div   tw_div_fast_ct_asm
#else
#define fpr_div   fndsa_fpr_div
#endif
fpr fpr_div(fpr x, fpr y);

/* Inversion. */
static inline fpr
fpr_inv(fpr x)
{
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
	extern tw_fpr tw_reci_fast_ct_asm(tw_fpr x);
	return tw_reci_fast_ct_asm(x);
#else
	return tw_reci_fast_ct(x);
#endif
}

/* Square root. */
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
#define fpr_sqrt   tw_sqrt_fast_ct_asm
#else
#define fpr_sqrt   fndsa_fpr_sqrt
#endif
fpr fpr_sqrt(fpr x);

/* Division by 2^e (scale each limb by 2^-e). */
static inline fpr
fpr_div2e(fpr x, unsigned e)
{
	funion_t xu;
	xu.i = (-(signed int)e + 127) << 23;
	return (tw_fpr){{xu.x * x.x[0], xu.x * x.x[1], xu.x * x.x[2]}};
}

/* Multiplication by 2^e (scale each limb by 2^e). */
static inline fpr
fpr_mul2e(fpr x, unsigned e)
{
	funion_t xu;
	xu.i = (e + 127) << 23;
	return (tw_fpr){{xu.x * x.x[0], xu.x * x.x[1], xu.x * x.x[2]}};
}

/* Complex multiplication (four real multiplications). */
#define FPC_MUL(d_re, d_im, a_re, a_im, b_re, b_im)   do { \
		fpr fpct_a_re = (a_re), fpct_a_im = (a_im); \
		fpr fpct_b_re = (b_re), fpct_b_im = (b_im); \
		fpr fpct_d_re = fpr_sub( \
			fpr_mul(fpct_a_re, fpct_b_re), \
			fpr_mul(fpct_a_im, fpct_b_im)); \
		fpr fpct_d_im = fpr_add( \
			fpr_mul(fpct_a_re, fpct_b_im), \
			fpr_mul(fpct_a_im, fpct_b_re)); \
		(d_re) = fpct_d_re; \
		(d_im) = fpct_d_im; \
	} while (0)

/* ==================================================================== */
/* Constants, from TWFalcon; the tables are in fpr_const_tw.c.         */

typedef union {
	fpr f;
} fpr_u;

/* 1/(2*(1.8205^2)) */
#define INV_2SQRSIGMA0   FPR(0x1.34f8bcp-3, 0x1.83bbc2p-31, 0x1.5d1caap-60)
/* log(2) */
#define LOG2       FPR(0x1.62e43p-1, -0x1.05c61p-29, -0x1.950d88p-54)
/* 1/log(2) */
#define INV_LOG2   FPR(0x1.715476p+0, 0x1.4ae0cp-26, -0x1.e8883p-52)
/* 1/q and -1/q (q = 12289), used to compute the target. */
#define INV_Q         FPR(0x1.554e3ap-14, -0x1.ed0b1p-39, 0x1.d7f62ap-66)
#define MINUS_INV_Q   FPR(-0x1.554e3ap-14, 0x1.ed0b1p-39, -0x1.d7f62ap-66)
/* q^2 = 151019521 (exact), the determinant of the Gram matrix. */
#define Q2            FPR(0x1.200cp+27, 0x1p+0, 0x0p+0)

#define GM          fndsa_fpr_GM
#define INV_SIGMA   fndsa_fpr_INV_SIGMA
#define SIGMA_MIN   fndsa_fpr_SIGMA_MIN
extern const fpr GM[];
extern const fpr_u INV_SIGMA[];
extern const fpr_u SIGMA_MIN[];

/* ==================================================================== */

#endif
