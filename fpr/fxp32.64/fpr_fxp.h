#ifndef FNDSA_FPR_FXP3264_H__
#define FNDSA_FPR_FXP3264_H__

/*
 * fxp32.64: signed 96-bit fixed point, value = raw * 2^-64, three 32-bit
 * limbs with w[0] lowest. Values outside [-2^31, 2^31) wrap; the key
 * generation gamma checks keep every signing value in range, and fpr_mul
 * relies on that by forming only the low 160 bits of the product.
 *
 * add/sub/neg are exact; mul rounds to nearest (ties to even, or ties up
 * with FNDSA_FXP_TIES_UP). Reciprocal and square root use Newton-Raphson
 * kernels whose wide result is rounded once by the consumer. There is no
 * signed zero, infinity or NaN; divisors must be in [1, 2^31).
 */

#include <stdint.h>

typedef struct {
	uint32_t w[3];
} fxp_fpr;

typedef fxp_fpr fpr;

typedef union {
	struct {
		fpr re;
		fpr im;
	} z;
	uint32_t w[6];
} fxp_fpc;

/* The signing buffer math uses FPR_STRIDE, so it must equal sizeof(fpr). */
#define FPR_STRIDE   12
typedef char fxp_stride_check[(sizeof(fxp_fpr) == FPR_STRIDE) ? 1 : -1];
typedef char fxp_fpc_stride_check[(sizeof(fxp_fpc) == 2 * FPR_STRIDE) ? 1 : -1];

/* Read by the MPFR test oracle and the prec backend. */
#define FXP_INT_BITS    32
#define FXP_FRAC_BITS   64

/* For a cm4 build with all fxp ops in C, pass
   -DFNDSA_FPR_IMPL=FNDSA_FPR_IMPL_C. */
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4 && !FNDSA_ASM_CORTEXM4
#error fxp32.64 cm4 fpr backend requires FNDSA_ASM_CORTEXM4=1
#endif

/* ==================================================================== */
/* Constants.                                                           */

/* Limbs in storage order: low fraction, high fraction, integer. */
#define FPR(w0, w1, w2)   ((fxp_fpr){{(w0), (w1), (w2)}})

#define FPR_ZERO    ((fxp_fpr){{0, 0, 0}})
/* No negative zero in two's complement. */
#define FPR_NZERO   ((fxp_fpr){{0, 0, 0}})
#define FPR_ONE     ((fxp_fpr){{0, 0, 1}})

/* ==================================================================== */
/* Constant-time 64-bit shift helpers.                                  */
/* Integer (not fpr) shifts for the sampler, as in the other backends. */
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
/* 96-bit limb helpers, in portable C with nothing wider than uint64_t.  */
/* The shifts are constant time and accept counts above 95.             */

#define FXP_N   3   /* limbs per fpr value */

/* Sign bit of a (1 if a < 0). */
static inline uint32_t
fxp_is_neg96(const uint32_t *a)
{
	return a[FXP_N - 1] >> 31;
}

/* r <- a << n (mod 2^96); n >= 96 yields zero. */
static inline void
fxp_shl96(uint32_t *r, const uint32_t *a, unsigned n)
{
	uint32_t bsh = n & 31;
	uint32_t w0, w1, w2, m;

	/* (v >> 1) >> (31 - bsh) is the funnel term, defined (and zero)
	   at bsh == 0 */
	w0 = a[0] << bsh;
	w1 = (a[1] << bsh) | ((a[0] >> 1) >> (31 - bsh));
	w2 = (a[2] << bsh) | ((a[1] >> 1) >> (31 - bsh));
	/* move up one limb if bit 5 of n */
	m = -(uint32_t)((n >> 5) & 1);
	w2 ^= (w2 ^ w1) & m;
	w1 ^= (w1 ^ w0) & m;
	w0 &= ~m;
	/* move up two limbs if bit 6 (n in [96,127] zeroes out here) */
	m = -(uint32_t)((n >> 6) & 1);
	w2 ^= (w2 ^ w0) & m;
	w1 &= ~m;
	w0 &= ~m;
	/* n >= 128: zero */
	m = (uint32_t)((0 - (uint64_t)(n >> 7)) >> 32);
	r[0] = w0 & ~m;
	r[1] = w1 & ~m;
	r[2] = w2 & ~m;
}

/* r <- a / 2^n, rounded to nearest, ties to even; n >= 96 yields zero. */
static inline void
fxp_sar_rn96(uint32_t *r, const uint32_t *a, unsigned n)
{
	uint32_t q[FXP_N], t[FXP_N];
	uint32_t sgn, bsh, rb, st, inc, m, s;

	/* clamp n at 96: every larger shift gives the same result (0) */
	m = 0 - (uint32_t)(((uint64_t)96 - n) >> 63);
	n ^= (n ^ 96u) & m;

	/* q = floor(a / 2^n) */
	sgn = -(uint32_t)fxp_is_neg96(a);
	bsh = n & 31;
	q[0] = (a[0] >> bsh) | ((a[1] << 1) << (31 - bsh));
	q[1] = (a[1] >> bsh) | ((a[2] << 1) << (31 - bsh));
	q[2] = (uint32_t)((int32_t)a[2] >> bsh);
	/* move down one limb if bit 5 of n */
	m = -(uint32_t)((n >> 5) & 1);
	q[0] ^= (q[0] ^ q[1]) & m;
	q[1] ^= (q[1] ^ q[2]) & m;
	q[2] ^= (q[2] ^ sgn) & m;
	/* move down two limbs if bit 6 */
	m = -(uint32_t)((n >> 6) & 1);
	q[0] ^= (q[0] ^ q[2]) & m;
	q[1] ^= (q[1] ^ sgn) & m;
	q[2] ^= (q[2] ^ sgn) & m;

	/* the dropped bits, shifted to the top, give round bit and sticky */
	fxp_shl96(t, a, 96 - n);
	rb = t[2] >> 31;
	s = (t[2] << 1) | t[1] | t[0];
	st = (uint32_t)((0 - (uint64_t)s) >> 63);
	inc = rb & (st | (q[0] & 1));
	uint64_t cc = inc;
	for (int i = 0; i < FXP_N; i ++) {
		cc += (uint64_t)q[i];
		r[i] = (uint32_t)cc;
		cc >>= 32;
	}
}

/* Integer part, floor(a / 2^64). */
static inline int64_t
fxp_int96(const uint32_t *a)
{
	return (int64_t)(int32_t)a[2];
}

/* ==================================================================== */
/* Conversions and rounding.                                            */

/* i*2^sc, rounded to nearest (ties to even). */
#define fpr_scaled   fndsa_fpr_scaled
fpr fpr_scaled(int64_t i, int sc);

/* Exact for |i| < 2^31; larger values wrap. */
static inline fpr
fpr_of(int64_t i)
{
	fpr d;

	d.w[0] = 0;
	d.w[1] = 0;
	d.w[2] = (uint32_t)(uint64_t)i;
	return d;
}

static inline fpr
fpr_of32(int32_t i)
{
	fpr d;

	d.w[0] = 0;
	d.w[1] = 0;
	d.w[2] = (uint32_t)i;
	return d;
}

/* Round to nearest integer (ties to even). */
static inline int64_t
fpr_rint(fpr x)
{
	fpr d;

	fxp_sar_rn96(d.w, x.w, 64);
	return (int64_t)(((uint64_t)d.w[1] << 32) | (uint64_t)d.w[0]);
}

/* floor is the integer limb. */
static inline int64_t
fpr_floor(fpr x)
{
	return fxp_int96(x.w);
}

/* Round toward zero: floor, plus one if negative with a fraction. */
static inline int64_t
fpr_trunc(fpr x)
{
	int64_t f = fxp_int96(x.w);
	uint32_t frac_nz = (uint32_t)((x.w[0] | x.w[1]) != 0);

	return f + (int64_t)(uint64_t)(fxp_is_neg96(x.w) & frac_nz);
}

/* floor(x * 2^63) for x in [0, 1): the fraction limbs shifted right
   once. The generic x * 2^63 would not fit 32 integer bits. */
static inline uint64_t
fpr_p63(fpr x)
{
	return (((uint64_t)x.w[1] << 32) | (uint64_t)x.w[0]) >> 1;
}
#define fpr_p63   fpr_p63

/* ==================================================================== */
/* Arithmetic.                                                          */

/* Exact (mod 2^96). */
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
#define FNDSA_FPR_ADD_INLINE   1
static inline fpr
fpr_add(fpr x, fpr y)
{
	uint32_t z0 = x.w[0];
	uint32_t z1 = x.w[1];
	uint32_t z2 = x.w[2];

	__asm__ (
		"adds   %[z0], %[z0], %[y0]\n\t"
		"adcs   %[z1], %[z1], %[y1]\n\t"
		"adcs   %[z2], %[z2], %[y2]"
		: [z0] "+&r" (z0), [z1] "+&r" (z1), [z2] "+&r" (z2)
		: [y0] "r" (y.w[0]), [y1] "r" (y.w[1]),
		  [y2] "r" (y.w[2])
		: "cc");

	return (fpr){{ z0, z1, z2 }};
}
#else
static inline fpr
fpr_add(fpr x, fpr y)
{
	fpr d;
	uint64_t cc = 0;

	for (int i = 0; i < FXP_N; i ++) {
		cc += (uint64_t)x.w[i] + (uint64_t)y.w[i];
		d.w[i] = (uint32_t)cc;
		cc >>= 32;
	}
	return d;
}
#endif

/* Exact (mod 2^96). */
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
#define FNDSA_FPR_SUB_INLINE   1
static inline fpr
fpr_sub(fpr x, fpr y)
{
	uint32_t z0 = x.w[0];
	uint32_t z1 = x.w[1];
	uint32_t z2 = x.w[2];

	__asm__ (
		"subs   %[z0], %[z0], %[y0]\n\t"
		"sbcs   %[z1], %[z1], %[y1]\n\t"
		"sbcs   %[z2], %[z2], %[y2]"
		: [z0] "+&r" (z0), [z1] "+&r" (z1), [z2] "+&r" (z2)
		: [y0] "r" (y.w[0]), [y1] "r" (y.w[1]),
		  [y2] "r" (y.w[2])
		: "cc");

	return (fpr){{ z0, z1, z2 }};
}
#else
static inline fpr
fpr_sub(fpr x, fpr y)
{
	fpr d;
	uint64_t cc = 1;

	for (int i = 0; i < FXP_N; i ++) {
		cc += (uint64_t)x.w[i] + (uint64_t)(~y.w[i]);
		d.w[i] = (uint32_t)cc;
		cc >>= 32;
	}
	return d;
}
#endif

/* a = x + y, b = x - y. */
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
#define FNDSA_FPR_ADD_SUB_INLINE   1
#define FPR_ADD_SUB(a, b, x, y)   do { \
		fpr t_add_sub_x = (x); \
		fpr t_add_sub_y = (y); \
		uint32_t t_add_sub_x0 = t_add_sub_x.w[0]; \
		uint32_t t_add_sub_x1 = t_add_sub_x.w[1]; \
		uint32_t t_add_sub_x2 = t_add_sub_x.w[2]; \
		uint32_t t_add_sub_b0, t_add_sub_b1, t_add_sub_b2; \
		__asm__ ( \
			"subs   %[b0], %[x0], %[y0]\n\t" \
			"sbcs   %[b1], %[x1], %[y1]\n\t" \
			"sbcs   %[b2], %[x2], %[y2]\n\t" \
			"adds   %[x0], %[x0], %[y0]\n\t" \
			"adcs   %[x1], %[x1], %[y1]\n\t" \
			"adcs   %[x2], %[x2], %[y2]" \
			: [x0] "+&r" (t_add_sub_x0), \
			  [x1] "+&r" (t_add_sub_x1), \
			  [x2] "+&r" (t_add_sub_x2), \
			  [b0] "=&r" (t_add_sub_b0), \
			  [b1] "=&r" (t_add_sub_b1), \
			  [b2] "=&r" (t_add_sub_b2) \
			: [y0] "r" (t_add_sub_y.w[0]), \
			  [y1] "r" (t_add_sub_y.w[1]), \
			  [y2] "r" (t_add_sub_y.w[2]) \
			: "cc"); \
		(a) = (fpr){{ t_add_sub_x0, t_add_sub_x1, t_add_sub_x2 }}; \
		(b) = (fpr){{ t_add_sub_b0, t_add_sub_b1, t_add_sub_b2 }}; \
	} while (0)
#else
#define FPR_ADD_SUB(a, b, x, y)   do { \
		fpr t_add_sub_x = (x); \
		fpr t_add_sub_y = (y); \
		(a) = fpr_add(t_add_sub_x, t_add_sub_y); \
		(b) = fpr_sub(t_add_sub_x, t_add_sub_y); \
	} while (0)
#endif

/* Exact (mod 2^96). */
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
static inline fpr
fpr_neg(fpr x)
{
	uint32_t z0 = x.w[0];
	uint32_t z1 = x.w[1];
	uint32_t z2 = x.w[2];
	uint32_t zero;

	__asm__ (
		"movs   %[zero], #0\n\t"
		"rsbs   %[z0], %[z0], #0\n\t"
		"sbcs   %[z1], %[zero], %[z1]\n\t"
		"sbcs   %[z2], %[zero], %[z2]"
		: [z0] "+&r" (z0), [z1] "+&r" (z1), [z2] "+&r" (z2),
		  [zero] "=&r" (zero)
		:
		: "cc");

	return (fpr){{ z0, z1, z2 }};
}
#else
static inline fpr
fpr_neg(fpr x)
{
	fpr d;

	uint64_t cc = 1;
	for (int i = 0; i < FXP_N; i ++) {
		cc += (uint64_t)(~x.w[i]);
		d.w[i] = (uint32_t)cc;
		cc >>= 32;
	}
	return d;
}
#endif

/* x/2, rounded to nearest (ties to even): add 1 when the low two raw
   bits are 11. */
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
static inline fpr
fpr_half(fpr x)
{
	uint32_t z0, z1, z2, t;

	__asm__ (
		"and    %[t], %[x0], %[x0], lsr #1\n\t"
		"and    %[t], %[t], #1\n\t"
		"asrs   %[z2], %[x2], #1\n\t"
		"rrxs   %[z1], %[x1]\n\t"
		"rrx    %[z0], %[x0]\n\t"
		"adds   %[z0], %[z0], %[t]\n\t"
		"adcs   %[z1], %[z1], #0\n\t"
		"adc    %[z2], %[z2], #0"
		: [z0] "=&r" (z0), [z1] "=&r" (z1), [z2] "=&r" (z2),
		  [t] "=&r" (t)
		: [x0] "r" (x.w[0]), [x1] "r" (x.w[1]), [x2] "r" (x.w[2])
		: "cc");

	return (fpr){{ z0, z1, z2 }};
}
#else
static inline fpr
fpr_half(fpr x)
{
	fpr d;
	uint32_t inc = x.w[0] & (x.w[0] >> 1) & 1;
	uint64_t cc;

	cc = (uint64_t)((x.w[0] >> 1) | (x.w[1] << 31)) + inc;
	d.w[0] = (uint32_t)cc;
	cc = (cc >> 32) + (uint64_t)((x.w[1] >> 1) | (x.w[2] << 31));
	d.w[1] = (uint32_t)cc;
	d.w[2] = (uint32_t)(cc >> 32) + (uint32_t)((int32_t)x.w[2] >> 1);
	return d;
}
#endif

/* Exact (mod 2^96). */
static inline fpr
fpr_double(fpr x)
{
	fpr d;

	d.w[0] = x.w[0] << 1;
	d.w[1] = (x.w[1] << 1) | (x.w[0] >> 31);
	d.w[2] = (x.w[2] << 1) | (x.w[1] >> 31);
	return d;
}

/* Rounded to nearest, ties to even. FNDSA_FXP_TIES_UP=1 rounds ties up
   in every multiply (mul, complex and twiddle products, reciprocal and
   sqrt tails); half, div2e, rint and scaled keep ties to even. */
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
void fndsa_fpr_mul(void);
static inline fpr
fpr_mul(fpr x, fpr y)
{
	register uint32_t x0 __asm__("r0") = x.w[0];
	register uint32_t x1 __asm__("r1") = x.w[1];
	register uint32_t x2 __asm__("r2") = x.w[2];
	register uint32_t y0 __asm__("r3") = y.w[0];
	register uint32_t y1 __asm__("r4") = y.w[1];
	register uint32_t y2 __asm__("r5") = y.w[2];

	__asm__ volatile (
		"@ fndsa_fpr_mul %0 %1 %2 %3 %4 %5\n\t"
		"bl	fndsa_fpr_mul"
		: "+r" (x0), "+r" (x1), "+r" (x2)
		: "r" (y0), "r" (y1), "r" (y2)
		: "r6", "r7", "r8", "r9", "r12", "r14", "cc");

	return (fpr){{ x0, x1, x2 }};
}
#else
#define fpr_mul   fndsa_fpr_mul
fpr fpr_mul(fpr x, fpr y);
#endif

static inline fpr
fpr_sqr(fpr x)
{
	return fpr_mul(x, x);
}

/* 1/x, correctly rounded, for 1 <= x < 2^31. */
#define fpr_inv   fndsa_fpr_inv
fpr fpr_inv(fpr x);

/* x * fpr_inv(y), rounded twice. y must be in the fpr_inv domain. */
#define fpr_div   fndsa_fpr_div
fpr fpr_div(fpr x, fpr y);

/* x * rsqrt(x) with a single rounding; zero for negative x. */
#define fpr_sqrt   fndsa_fpr_sqrt
fpr fpr_sqrt(fpr x);

/* Rounded to nearest, ties to even. */
static inline fpr
fpr_div2e(fpr x, unsigned e)
{
	fpr d;

	fxp_sar_rn96(d.w, x.w, e);
	return d;
}

/* Exact (mod 2^96). */
static inline fpr
fpr_mul2e(fpr x, unsigned e)
{
	fpr d;

	fxp_shl96(d.w, x.w, e);
	return d;
}

/* ==================================================================== */
/* Constants, generated by fpr/gen/gen_const.c.                         */

typedef union {
	fpr f;
} fpr_u;

/* INV_2SQRSIGMA0, LOG2, INV_LOG2, INV_Q, MINUS_INV_Q, Q2, INV_Q_Q77_HI/_LO. */
#include "fxp_const_gen.h"

#define GM          fndsa_fpr_GM
#define RNJ         fndsa_fpr_RNJ
#define INV_SIGMA   fndsa_fpr_INV_SIGMA
#define INV_SIGMA_Q71   fndsa_fpr_INV_SIGMA_Q71
#define SIGMA_MIN   fndsa_fpr_SIGMA_MIN


extern const fpr RNJ[];
extern const fpr_u INV_SIGMA[];
extern const uint64_t INV_SIGMA_Q71[];
extern const fpr_u SIGMA_MIN[];

/* ==================================================================== */

#endif
