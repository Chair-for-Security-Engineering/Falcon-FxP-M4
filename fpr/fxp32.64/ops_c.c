/*
 * Out-of-line fpr operations for fxp32.64, in portable 32-bit C (no
 * __int128). Most of them are left out of the cm4 build, which has them
 * in sign_fpr_cm4.S with the same results.
 *
 * fpr_mul keeps only the low 160 bits of the product, which is exact as
 * long as the result is in [-2^31, 2^31). Reciprocal and inverse square
 * root use a few narrow Newton-Raphson steps followed by one full-width
 * step, and their unrounded result is rounded once by the consumer.
 */

#include "sign_inner.h"
#include <inner.h>

/* Rounding increment for multiplies: rb = round bit, st = sticky (0/1),
   odd = quotient parity. Ties to even, or ties up with FNDSA_FXP_TIES_UP. */
#define FXP_RN_INC(rb, st, odd) \
	((rb) & (FNDSA_FXP_TIES_UP ? 1u : ((st) | (odd))))

/* ==================================================================== */
/* 192-bit unsigned helpers.                                            */

#define U_N   6   /* limbs in a 192-bit value */

/* r <- a * b, unsigned 96x96 to 192 bits. */
static void
u192_mul96(uint32_t *r, const uint32_t *a, const uint32_t *b)
{
	for (int i = 0; i < U_N; i ++) {
		r[i] = 0;
	}
	for (int i = 0; i < FXP_N; i ++) {
		uint64_t cc = 0;
		for (int j = 0; j < FXP_N; j ++) {
			cc += (uint64_t)a[i] * (uint64_t)b[j]
				+ (uint64_t)r[i + j];
			r[i + j] = (uint32_t)cc;
			cc >>= 32;
		}
		r[i + FXP_N] = (uint32_t)cc;
	}
}

static inline uint32_t
ct_mask_nonzero32(uint32_t x)
{
	return (uint32_t)(0 - ((x | (uint32_t)(0 - x)) >> 31));
}

static inline uint32_t
ct_mask_eq32(uint32_t x, uint32_t y)
{
	return ~ct_mask_nonzero32(x ^ y);
}

/* r <- low 96 bits of round(a / 2^n), for a 192-bit a and
   63 <= n <= 110, constant time. */
static inline void
u192_shr_rn96_ct(uint32_t *r, const uint32_t *a, unsigned n)
{
	uint32_t q0, q1, q2, rb, st, inc;
	uint32_t q01, q11, q21, q02, q12, q22, q03, q13, q23;
	uint32_t rb1, rb2, rb3, st1, st2, st3, m1, m2, m3, pm;
	uint64_t cc;
	unsigned wsh, bsh, rbw, rbb;

	wsh = n >> 5;
	bsh = n & 31;
	m1 = ct_mask_eq32((uint32_t)wsh, 1);
	m2 = ct_mask_eq32((uint32_t)wsh, 2);
	m3 = ~(m1 | m2);
	q01 = (a[1] >> bsh) | ((a[2] << 1) << (31 - bsh));
	q11 = (a[2] >> bsh) | ((a[3] << 1) << (31 - bsh));
	q21 = (a[3] >> bsh) | ((a[4] << 1) << (31 - bsh));
	q02 = q11;
	q12 = q21;
	q22 = (a[4] >> bsh) | ((a[5] << 1) << (31 - bsh));
	q03 = q21;
	q13 = q22;
	q23 = a[5] >> bsh;
	q0 = (q01 & m1) | (q02 & m2) | (q03 & m3);
	q1 = (q11 & m1) | (q12 & m2) | (q13 & m3);
	q2 = (q21 & m1) | (q22 & m2) | (q23 & m3);

	rbw = (n - 1) >> 5;
	rbb = (n - 1) & 31;
	pm = ((uint32_t)1 << rbb) - 1u;
	m1 = ct_mask_eq32((uint32_t)rbw, 1);
	m2 = ct_mask_eq32((uint32_t)rbw, 2);
	m3 = ~(m1 | m2);
	rb1 = (a[1] >> rbb) & 1;
	rb2 = (a[2] >> rbb) & 1;
	rb3 = (a[3] >> rbb) & 1;
	st1 = a[0] | (a[1] & pm);
	st2 = a[0] | a[1] | (a[2] & pm);
	st3 = a[0] | a[1] | a[2] | (a[3] & pm);
	rb = (rb1 & (m1 & 1u)) | (rb2 & (m2 & 1u)) | (rb3 & (m3 & 1u));
	st = (st1 & m1) | (st2 & m2) | (st3 & m3);
	st = (st | (0 - st)) >> 31;
	inc = FXP_RN_INC(rb, st, q0 & 1);
	cc = (uint64_t)q0 + inc;
	r[0] = (uint32_t)cc;
	cc = (uint64_t)q1 + (cc >> 32);
	r[1] = (uint32_t)cc;
	r[2] = q2 + (uint32_t)(cc >> 32);
}

/* ==================================================================== */
/* fpr_mul.                                                             */

#if FNDSA_FPR_IMPL != FNDSA_FPR_IMPL_CM4
/* p <- x * y (mod 2^160), signed: the unsigned limb product, then
   subtract sx*y + sy*x at the 2^96 column. A result in [-2^31, 2^31)
   needs no more than 160 bits, so x[2]*y[2] contributes only its low
   word. */
static void
smul96_160(uint32_t *p, const uint32_t *x, const uint32_t *y)
{
	uint64_t cc;
	uint32_t mx, my;

	for (int i = 0; i < 5; i ++) {
		p[i] = 0;
	}
	/* rows 0 and 1 */
	for (int i = 0; i < 2; i ++) {
		cc = 0;
		for (int j = 0; j < FXP_N; j ++) {
			cc += (uint64_t)x[i] * (uint64_t)y[j]
				+ (uint64_t)p[i + j];
			p[i + j] = (uint32_t)cc;
			cc >>= 32;
		}
		p[i + FXP_N] = (uint32_t)cc;
	}
	/* row 2: only the low word of x[2]*y[2] is below 2^160 */
	cc = (uint64_t)x[2] * (uint64_t)y[0] + (uint64_t)p[2];
	p[2] = (uint32_t)cc;
	cc >>= 32;
	cc += (uint64_t)x[2] * (uint64_t)y[1] + (uint64_t)p[3];
	p[3] = (uint32_t)cc;
	cc >>= 32;
	p[4] += (uint32_t)cc + x[2] * y[2];

	/* signed correction */
	mx = (uint32_t)((int32_t)x[2] >> 31);
	my = (uint32_t)((int32_t)y[2] >> 31);
	cc = (uint64_t)p[3] - (y[0] & mx) - (x[0] & my);
	p[3] = (uint32_t)cc;
	p[4] = p[4] - (y[1] & mx) - (x[1] & my) + (uint32_t)(cc >> 32);
}
#endif

/* r <- p / 2^64 (mod 2^96), rounded to nearest, for a 160-bit
   two's-complement p. */
#if FNDSA_FPR_IMPL != FNDSA_FPR_IMPL_CM4
static void
shr64_rn160(uint32_t *r, const uint32_t *p)
{
	uint32_t rb, st, inc;
	uint64_t cc;

	rb = p[1] >> 31;
	st = p[0] | (p[1] & 0x7FFFFFFF);
	st = (st | (0 - st)) >> 31;
	inc = FXP_RN_INC(rb, st, p[2] & 1);
	cc = (uint64_t)p[2] + inc;
	r[0] = (uint32_t)cc;
	cc = (uint64_t)p[3] + (cc >> 32);
	r[1] = (uint32_t)cc;
	r[2] = p[4] + (uint32_t)(cc >> 32);
}
#endif

/* ==================================================================== */
/* Newton-Raphson kernels.                                              */

/* 3 in Q94. */
static const uint32_t THREE_Q94[FXP_N] = { 0, 0, 0xC0000000 };

/* hi:lo <- a*b. */
static void
umul64_wide(uint64_t *hi, uint64_t *lo, uint64_t a, uint64_t b)
{
	uint64_t p0, p1, p2, p3, s;

	p0 = (uint64_t)(uint32_t)a * (uint32_t)b;
	p1 = (uint64_t)(uint32_t)(a >> 32) * (uint32_t)b;
	p2 = (uint64_t)(uint32_t)a * (uint32_t)(b >> 32);
	p3 = (uint64_t)(uint32_t)(a >> 32) * (uint32_t)(b >> 32);
	s = (p0 >> 32) + (uint32_t)p1 + (uint32_t)p2;
	*lo = (p0 & 0xFFFFFFFFu) | (s << 32);
	*hi = p3 + (p1 >> 32) + (p2 >> 32) + (s >> 32);
}

/* return (a*b) >> n, for n in {62, 63, 64}. */
static uint64_t
umul64_shr(uint64_t a, uint64_t b, unsigned n)
{
	uint64_t hi, lo;

	umul64_wide(&hi, &lo, a, b);
	if (n == 64) {
		return hi;
	}
	return (hi << (64 - n)) | (lo >> n);
}

/* r <- a*b, unsigned 96x64 to 160 bits. */
static void
umul96_64(uint32_t *r, const uint32_t *a, uint64_t b)
{
	uint32_t bw[2];
	uint64_t cc;

	bw[0] = (uint32_t)b;
	bw[1] = (uint32_t)(b >> 32);
	for (int i = 0; i < 5; i ++) {
		r[i] = 0;
	}
	for (int i = 0; i < FXP_N; i ++) {
		cc = 0;
		for (int j = 0; j < 2; j ++) {
			cc += (uint64_t)a[i] * (uint64_t)bw[j]
				+ (uint64_t)r[i + j];
			r[i + j] = (uint32_t)cc;
			cc >>= 32;
		}
		r[i + 2] = (uint32_t)cc;
	}
}

/* r <- a >> n, for a 160-bit source and n in {62, 63, 64}; low 96 bits. */
static void
u160_shr96(uint32_t *r, const uint32_t *a, unsigned n)
{
	unsigned wsh = n >> 5, bsh = n & 31;

	for (int i = 0; i < FXP_N; i ++) {
		uint32_t lo = a[(unsigned)i + wsh];
		uint32_t hi = ((unsigned)i + wsh + 1 < 5) ?
			a[(unsigned)i + wsh + 1] : 0;
		r[i] = lo >> bsh;
		if (bsh != 0) {
			r[i] |= hi << (32 - bsh);
		}
	}
}

/* r <- a << n (mod 2^96), for n < 96. */
static inline void
fxp_shl96_lt96(uint32_t *r, const uint32_t *a, unsigned n)
{
	uint32_t bsh, w0, w1, w2, m;

	bsh = n & 31;
	w0 = a[0] << bsh;
	w1 = (a[1] << bsh) | ((a[0] >> 1) >> (31 - bsh));
	w2 = (a[2] << bsh) | ((a[1] >> 1) >> (31 - bsh));
	m = -(uint32_t)((n >> 5) & 1);
	w2 ^= (w2 ^ w1) & m;
	w1 ^= (w1 ^ w0) & m;
	w0 &= ~m;
	m = -(uint32_t)((n >> 6) & 1);
	w2 ^= (w2 ^ w0) & m;
	w1 &= ~m;
	w0 &= ~m;
	r[0] = w0;
	r[1] = w1;
	r[2] = w2;
}

/* Bit length of a non-negative 96-bit value, with fixed control flow.
   Returns 0 for zero. */
static inline unsigned
fxp_bitlen96_ct(const uint32_t *a)
{
	uint32_t m2, m1, m0, v, base;

	m2 = ct_mask_nonzero32(a[2]);
	m1 = ct_mask_nonzero32(a[1]) & ~m2;
	m0 = ~(m2 | m1);
	v = (a[2] & m2) | (a[1] & m1) | (a[0] & m0);
	base = (m2 & 64u) | (m1 & 32u);
	return (unsigned)(base + 32u - (uint32_t)lzcnt(v));
}

#if FNDSA_FPR_IMPL != FNDSA_FPR_IMPL_CM4
/* Seed 48/17 - 32/17*d, as Q31 words mod 2^32 (M4 immediates). */
#define RECIP31_C1   ((uint32_t)0x69696969u)
#define RECIP31_C2   ((uint32_t)0xF0F0F0F0u)

/* r[0..3] <- 1/d in Q126, unrounded, for d normalized to [1/2, 1)
   (bit 95 set). Two Newton steps on 32-bit words, one on 64 bits, then
   one against the full 96-bit d, which also corrects the earlier
   truncations of d. */
static void
fxp_recip_fast(uint32_t *r, const uint32_t *d)
{
	uint64_t d64, y, t, cc;
	uint32_t yw[2], p[5], en[4], q[6];

	d64 = ((uint64_t)d[2] << 32) | (uint64_t)d[1];

	/* Seed and two steps on 32-bit words against d[2]. y stays below 2
	   (Newton approaches 1/d from below) and is kept in Q31. */
	{
		uint32_t d32 = d[2], z, w, t32;

		/* y0 (Q31); y1 = high32(y0*(2 - d*y0)) lands in Q30; y2 in
		   Q28 (with 2 in Q30 = 2^31); one shift brings it to Q31. */
		z = RECIP31_C1 - (uint32_t)(((uint64_t)RECIP31_C2 * d32) >> 32);
		t32 = (uint32_t)(((uint64_t)d32 * z) >> 32);
		w = (uint32_t)(((uint64_t)z * (0u - t32)) >> 32);
		t32 = (uint32_t)(((uint64_t)d32 * w) >> 32);
		z = (uint32_t)(((uint64_t)w * (0x80000000u - t32)) >> 32) << 3;
		y = (uint64_t)z << 32;
	}
	/* t = d*y in Q63, e = 2 - t (2^64 wraps to 0), then
	   y3 = (y*e) >> 64 in Q62. */
	t = umul64_shr(d64, y, 64);
	y = umul64_shr(y, 0 - t, 64);

	/* Final step against the full divisor: P = d * y3 in Q158, exact. */
	yw[0] = (uint32_t)y;
	yw[1] = (uint32_t)(y >> 32);
	for (int i = 0; i < 5; i ++) {
		p[i] = 0;
	}
	for (int i = 0; i < FXP_N; i ++) {
		cc = 0;
		for (int j = 0; j < 2; j ++) {
			cc += (uint64_t)d[i] * (uint64_t)yw[j]
				+ (uint64_t)p[i + j];
			p[i + j] = (uint32_t)cc;
			cc >>= 32;
		}
		p[i + 2] = (uint32_t)cc;
	}

	/* en = (2^159 - P) >> 32, i.e. E = 2 - d*y3 in Q126. */
	cc = 1 + (uint64_t)(uint32_t)~p[0];
	cc >>= 32;
	for (int i = 1; i < 4; i ++) {
		cc += (uint64_t)(uint32_t)~p[i];
		en[i - 1] = (uint32_t)cc;
		cc >>= 32;
	}
	cc += (uint64_t)0x80000000u + (uint64_t)(uint32_t)~p[4];
	en[3] = (uint32_t)cc;

	/* r = (y3 * E) >> 62, back to Q126. */
	for (int i = 0; i < 6; i ++) {
		q[i] = 0;
	}
	for (int i = 0; i < 2; i ++) {
		cc = 0;
		for (int j = 0; j < 4; j ++) {
			cc += (uint64_t)yw[i] * (uint64_t)en[j]
				+ (uint64_t)q[i + j];
			q[i + j] = (uint32_t)cc;
			cc >>= 32;
		}
		q[i + 4] = (uint32_t)cc;
	}
	for (int i = 0; i < 4; i ++) {
		r[i] = (q[i + 1] >> 30) | (q[i + 2] << 2);
	}
}
#endif

#if FNDSA_FPR_IMPL != FNDSA_FPR_IMPL_CM4
/* p[0..7] <- (a << k) * rc, exact and signed: the quotient in Q158. */
static inline void
mul_recip_prod(uint32_t *p, const uint32_t *a, const uint32_t *rc,
	unsigned k)
{
	uint32_t am[4], ms;
	uint64_t cc;

	/* a << k, sign-extended to 128 bits */
	ms = (uint32_t)((int32_t)a[2] >> 31);
	am[0] = a[0] << k;
	am[1] = (a[1] << k) | (a[0] >> (32 - k));
	am[2] = (a[2] << k) | (a[1] >> (32 - k));
	am[3] = (ms << k) | (a[2] >> (32 - k));

	for (int i = 0; i < 8; i ++) {
		p[i] = 0;
	}
	for (int i = 0; i < 4; i ++) {
		cc = 0;
		for (int j = 0; j < 4; j ++) {
			cc += (uint64_t)am[i] * (uint64_t)rc[j]
				+ (uint64_t)p[i + j];
			p[i + j] = (uint32_t)cc;
			cc >>= 32;
		}
		p[i + 4] = (uint32_t)cc;
	}

	/* signed correction: subtract rc at 2^128 if a < 0 */
	cc = 1;
	for (int j = 0; j < 4; j ++) {
		cc += (uint64_t)p[4 + j] + (uint32_t)~(rc[j] & ms);
		p[4 + j] = (uint32_t)cc;
		cc >>= 32;
	}
}

/* z <- round(a * rc * 2^(k-158)) = a/x, with (rc, k) from FPR_RECIP_X. */
static void
fxp_mul_recip(uint32_t *z, const uint32_t *a, const uint32_t *rc,
	unsigned k)
{
	uint32_t p[8];
	uint32_t rb, st, inc;
	uint64_t cc;

	mul_recip_prod(p, a, rc, k);

	/* round at bit 158 */
	rb = (p[4] >> 29) & 1;
	st = p[0] | p[1] | p[2] | p[3] | (p[4] & 0x1FFFFFFF);
	st = (st | (0 - st)) >> 31;
	for (int i = 0; i < FXP_N; i ++) {
		z[i] = (p[i + 4] >> 30) | (p[i + 5] << 2);
	}
	inc = FXP_RN_INC(rb, st, z[0] & 1);
	cc = (uint64_t)z[0] + inc;
	z[0] = (uint32_t)cc;
	cc = (uint64_t)z[1] + (cc >> 32);
	z[1] = (uint32_t)cc;
	z[2] += (uint32_t)(cc >> 32);
}

/* a/x rounded to a signed Q0.63 value, for the non-root l10. Requires
   |a/x| < 1. */
static uint64_t
fxp_mul_recip_p63(const uint32_t *a, const uint32_t *rc, unsigned k)
{
	uint32_t p[8], rb, st, inc;
	uint64_t q;

	mul_recip_prod(p, a, rc, k);

	/* round at bit 159 */
	rb = (p[4] >> 30) & 1;
	st = p[0] | p[1] | p[2] | p[3] | (p[4] & 0x3FFFFFFF);
	st = (st | (0 - st)) >> 31;
	q = ((uint64_t)((p[5] >> 31) | (p[6] << 1)) << 32)
		| ((p[4] >> 31) | (p[5] << 1));
	inc = FXP_RN_INC(rb, st, (uint32_t)q & 1);
	return q + inc;
}
#endif

#if FNDSA_FPR_IMPL != FNDSA_FPR_IMPL_CM4
void
fndsa_fxp96_recip_x(uint32_t *m, uint32_t *k, const uint32_t *x)
{
	uint32_t dn[FXP_N];
	unsigned kk;

	kk = lzcnt_nonzero(x[2]);
	/* 1 <= x < 2^31, so kk is in 1..31 and a funnel shift suffices */
	dn[0] = x[0] << kk;
	dn[1] = (x[1] << kk) | ((x[0] >> 1) >> (31 - kk));
	dn[2] = (x[2] << kk) | ((x[1] >> 1) >> (31 - kk));
	fxp_recip_fast(m, dn);
	*k = kk;
}
#endif

#if FNDSA_FPR_IMPL != FNDSA_FPR_IMPL_CM4
void
fndsa_fxp96_mul_recip(uint32_t *z, const uint32_t *a, const uint32_t *m,
	uint32_t k)
{
	fxp_mul_recip(z, a, m, k);
}

uint64_t
fndsa_fxp96_mul_recip_p63(const uint32_t *a, const uint32_t *m, uint32_t k)
{
	return fxp_mul_recip_p63(a, m, k);
}
#endif

/* r <- y*(3 - m*y^2)/2 in Q94, from y in Q62, against the full 96-bit m:
   a = (m*y) >> 64, b = (a*y) >> 62, r = (y*(3 - b)) >> 63. */
static void
fxp_rsqrt_step_full96_q62(uint32_t *r, const uint32_t *m, uint64_t y)
{
	uint32_t p[5], a[FXP_N], b[FXP_N], e[FXP_N];

	umul96_64(p, m, y);
	for (int j = 0; j < FXP_N; j ++) {
		a[j] = p[j + 2];
	}
	umul96_64(p, a, y);
	u160_shr96(b, p, 62);
	uint64_t cc = 1;
	for (int j = 0; j < FXP_N; j ++) {
		cc += (uint64_t)THREE_Q94[j] + (uint64_t)(~b[j]);
		e[j] = (uint32_t)cc;
		cc >>= 32;
	}
	umul96_64(p, e, y);
	u160_shr96(r, p, 63);
}

/* Seed 2.2118 - (4/3)*a, as Q31 words mod 2^32 (M4 immediates). */
#define RSQRT31_C3   ((uint32_t)0x1B1B1B1Bu)
#define RSQRT31_C4   ((uint32_t)0xAAAAAAAAu)
#define THREE_Q62    ((uint64_t)0xC000000000000000u)

/* r <- 1/sqrt(m/2^96) in Q94, unrounded, for m in [2^94, 2^96). Three
   Newton steps on 32-bit words, one on 64 bits, then one against the
   full 96-bit m. */
static void
fxp_rsqrt_fast(uint32_t *r, const uint32_t *m)
{
	uint64_t m64, y, a, b, e;

	m64 = ((uint64_t)m[2] << 32) | (uint64_t)m[1];
	/* y stays below 2 and is kept in Q31. */
	{
		uint32_t m32 = m[2], z, a32, b32;

		/* y0 Q31, y1 Q30, y2 Q27 << 4, y3 Q30 << 1 */
		z = RSQRT31_C3 - (uint32_t)(((uint64_t)RSQRT31_C4 * m32) >> 32);
		a32 = (uint32_t)(((uint64_t)m32 * z) >> 32);
		b32 = (uint32_t)(((uint64_t)a32 * z) >> 32);
		z = (uint32_t)(((uint64_t)z * (0xC0000000u - b32)) >> 32);
		a32 = (uint32_t)(((uint64_t)m32 * z) >> 32);
		b32 = (uint32_t)(((uint64_t)a32 * z) >> 32);
		z = (uint32_t)(((uint64_t)z * (0x30000000u - b32)) >> 32) << 4;
		a32 = (uint32_t)(((uint64_t)m32 * z) >> 32);
		b32 = (uint32_t)(((uint64_t)a32 * z) >> 32);
		z = (uint32_t)(((uint64_t)z * (0xC0000000u - b32)) >> 32) << 1;
		y = (uint64_t)z << 32;
	}
	/* a = m*y (Q63), b = a*y (Q62), e = 3 - b (Q62), y = y*e/2 (Q62). */
	a = umul64_shr(m64, y, 64);
	b = umul64_shr(a, y, 64);
	e = THREE_Q62 - b;
	y = umul64_shr(y, e, 64);

	fxp_rsqrt_step_full96_q62(r, m, y);
}

/* ==================================================================== */
/* The fpr operations.                                                  */

fpr
fpr_scaled(int64_t i, int sc)
{
	fpr d;
	uint32_t v[FXP_N];
	uint64_t u = (uint64_t)i;
	uint32_t sgn = (uint32_t)-(uint32_t)(i < 0);
	/* 64-bit so that an extreme sc cannot overflow */
	int64_t e = (int64_t)64 + (int64_t)sc;

	/* raw = i*2^(64+sc) */
	v[0] = (uint32_t)u;
	v[1] = (uint32_t)(u >> 32);
	v[2] = sgn;
	if (e >= 96 || e <= -96) {
		return FPR_ZERO;
	}
	if (e >= 0) {
		fxp_shl96(d.w, v, (unsigned)e);
	} else {
		fxp_sar_rn96(d.w, v, (unsigned)-e);
	}
	return d;
}

#if FNDSA_FPR_IMPL != FNDSA_FPR_IMPL_CM4
fpr
fpr_mul(fpr x, fpr y)
{
	fpr d;
	uint32_t p[5];

	smul96_160(p, x.w, y.w);
	shr64_rn160(d.w, p);
	return d;
}

/* Each component is summed exactly (mod 2^160) and rounded once. */
#if FNDSA_FXP_CMUL_FUSED
void
fndsa_fxp96_cmul_rne64(uint32_t *out, const uint32_t *a, const uint32_t *b)
{
	fpr a_re, a_im, b_re, b_im, d_re, d_im;

	for (int i = 0; i < FXP_N; i ++) {
		a_re.w[i] = a[i];
		a_im.w[i] = a[FXP_N + i];
		b_re.w[i] = b[i];
		b_im.w[i] = b[FXP_N + i];
	}
    uint32_t y[5], y2[5];
    uint64_t cc;


    smul96_160(y, a_re.w, b_re.w);
    smul96_160(y2, a_im.w, b_im.w);

    /* y <- y - y2 (mod 2^160), two's complement. */
    cc = 1;
    for (int i = 0; i < 5; i ++) {
        cc += (uint64_t)y[i] + (uint64_t)(uint32_t)~y2[i];
        y[i] = (uint32_t)cc;
        cc >>= 32;
    }

    shr64_rn160(d_re.w, y);

    smul96_160(y, a_re.w, b_im.w);
    smul96_160(y2, a_im.w, b_re.w);

    /* y <- y + y2 (mod 2^160), two's complement. */
    cc = 0;
    for (int i = 0; i < 5; i ++) {
        cc += (uint64_t)y[i] + (uint64_t)y2[i];
        y[i] = (uint32_t)cc;
        cc >>= 32;
    }

    shr64_rn160(d_im.w, y);

	for (int i = 0; i < FXP_N; i ++) {
		out[i] = d_re.w[i];
		out[FXP_N + i] = d_im.w[i];
	}
}
#endif /* FNDSA_FXP_CMUL_FUSED */
#endif

#if FNDSA_FPR_IMPL != FNDSA_FPR_IMPL_CM4
/* p <- a * b, exact, for signed 96-bit a and unsigned 64-bit b. */
static void
smul96x64(uint32_t *p, const uint32_t *a, uint64_t b)
{
	uint32_t bw[2], ma;
	uint64_t cc;

	bw[0] = (uint32_t)b;
	bw[1] = (uint32_t)(b >> 32);
	for (int i = 0; i < 5; i ++) {
		p[i] = 0;
	}
	for (int i = 0; i < FXP_N; i ++) {
		cc = 0;
		for (int j = 0; j < 2; j ++) {
			cc += (uint64_t)a[i] * (uint64_t)bw[j]
				+ (uint64_t)p[i + j];
			p[i + j] = (uint32_t)cc;
			cc >>= 32;
		}
		p[i + 2] = (uint32_t)cc;
	}

	/* signed correction: subtract b at 2^96 if a < 0 */
	ma = (uint32_t)((int32_t)a[2] >> 31);
	cc = (uint64_t)p[3] - (bw[0] & ma);
	p[3] = (uint32_t)cc;
	p[4] = p[4] - (bw[1] & ma) + (uint32_t)(cc >> 32);
}

/* x/q as x * INV_Q_Q77 / 2^77, rounded once. */
fpr
fndsa_fpr_mul_inv_q(fpr x)
{
	fpr d;
	uint32_t p[5], rb, st, inc;
	uint64_t cc;

	smul96x64(p, x.w, INV_Q_Q77);
	rb = (p[2] >> 12) & 1;
	st = p[0] | p[1] | (p[2] & 0xFFF);
	st = (st | (0 - st)) >> 31;
	d.w[0] = (p[2] >> 13) | (p[3] << 19);
	d.w[1] = (p[3] >> 13) | (p[4] << 19);
	d.w[2] = (uint32_t)((int32_t)p[4] >> 13);
	inc = FXP_RN_INC(rb, st, d.w[0] & 1);
	cc = (uint64_t)d.w[0] + inc;
	d.w[0] = (uint32_t)cc;
	cc = (uint64_t)d.w[1] + (cc >> 32);
	d.w[1] = (uint32_t)cc;
	d.w[2] += (uint32_t)(cc >> 32);
	return d;
}
#endif

#if FNDSA_TWIDDLE_U64 && FNDSA_FPR_IMPL != FNDSA_FPR_IMPL_CM4

/* out = a*(A+iB) (cmul_gm) or a*(A-iB) (cmul_gmj) for an unsigned Q0.64
   GM pair (A, B), each component rounded once. out may alias a. */
void
fndsa_fxp96_cmul_gm(uint32_t *out, const uint32_t *a, const uint64_t *gm)
{
	uint32_t p1[5], p2[5], p3[5];
	uint64_t cc;

	/* re = a_re*A - a_im*B */
	smul96x64(p1, a, gm[0]);
	smul96x64(p2, a + FXP_N, gm[1]);
	cc = 1;
	for (int i = 0; i < 5; i ++) {
		cc += (uint64_t)p1[i] + (uint64_t)(uint32_t)~p2[i];
		p1[i] = (uint32_t)cc;
		cc >>= 32;
	}

	/* im = a_re*B + a_im*A */
	smul96x64(p2, a, gm[1]);
	smul96x64(p3, a + FXP_N, gm[0]);
	cc = 0;
	for (int i = 0; i < 5; i ++) {
		cc += (uint64_t)p2[i] + (uint64_t)p3[i];
		p2[i] = (uint32_t)cc;
		cc >>= 32;
	}

	shr64_rn160(out, p1);
	shr64_rn160(out + FXP_N, p2);
}

void
fndsa_fxp96_cmul_gmj(uint32_t *out, const uint32_t *a, const uint64_t *gm)
{
	uint32_t p1[5], p2[5], p3[5];
	uint64_t cc;

	/* re = a_re*A + a_im*B */
	smul96x64(p1, a, gm[0]);
	smul96x64(p2, a + FXP_N, gm[1]);
	cc = 0;
	for (int i = 0; i < 5; i ++) {
		cc += (uint64_t)p1[i] + (uint64_t)p2[i];
		p1[i] = (uint32_t)cc;
		cc >>= 32;
	}

	/* im = a_im*A - a_re*B */
	smul96x64(p2, a + FXP_N, gm[0]);
	smul96x64(p3, a, gm[1]);
	cc = 1;
	for (int i = 0; i < 5; i ++) {
		cc += (uint64_t)p2[i] + (uint64_t)(uint32_t)~p3[i];
		p2[i] = (uint32_t)cc;
		cc >>= 32;
	}

	shr64_rn160(out, p1);
	shr64_rn160(out + FXP_N, p2);
}

/* x times one GM word. */
fpr
fndsa_fpr_mul_gm(fpr x, uint64_t c)
{
	fpr d;
	uint32_t p[5];

	smul96x64(p, x.w, c);
	shr64_rn160(d.w, p);
	return d;
}

/* x * (integer + frac/2^64), with the integer part added before the
   single rounding. */
fpr
fndsa_fpr_mul_rnj(fpr x, uint64_t frac, uint32_t integer)
{
	fpr d;
	uint32_t p[5], a0, a1, a2;
	uint64_t cc;

	smul96x64(p, x.w, frac);
	cc = (uint64_t)x.w[0] * integer;
	a0 = (uint32_t)cc;
	cc = (cc >> 32) + (uint64_t)x.w[1] * integer;
	a1 = (uint32_t)cc;
	cc = (cc >> 32) + (uint64_t)x.w[2] * integer;
	a2 = (uint32_t)cc;
	cc = (uint64_t)p[2] + a0;
	p[2] = (uint32_t)cc;
	cc = (cc >> 32) + (uint64_t)p[3] + a1;
	p[3] = (uint32_t)cc;
	p[4] += a2 + (uint32_t)(cc >> 32);
	shr64_rn160(d.w, p);
	return d;
}
#endif

/* 1/x for 1 <= x < 2^31: the FNDSA_FXP96_RECIP_RAW kernel rounded to the
   grid, then corrected to the correctly rounded value. Constant time. */
fpr
fpr_inv(fpr x)
{
	fpr d;
	uint32_t m[7], q[U_N], t[U_N];
	uint32_t k, up, dn, b;
	uint64_t ca, cb, cc;

	FNDSA_FXP96_RECIP_RAW(m, k, x.w);
	m[4] = 0;
	m[5] = 0;
	m[6] = 0;
	u192_shr_rn96_ct(d.w, m, 94 - k);

	/* d can be one ulp off near x = 1. With Q = x*d in raw units,
	   2^128/x is above d + 1/2 iff 2Q + x < 2^129, and below d - 1/2
	   iff 2Q - x > 2^129. */
	u192_mul96(q, x.w, d.w);
	for (int i = U_N - 1; i > 0; i --) {
		t[i] = (q[i] << 1) | (q[i - 1] >> 31);
	}
	t[0] = q[0] << 1;
	ca = 0;
	cb = 0;
	up = 0;
	dn = 0;
	for (int i = 0; i < U_N; i ++) {
		uint32_t xi = (i < FXP_N) ? x.w[i] : 0;
		uint32_t ci = (i == 4) ? 2 : 0;
		uint32_t ai, bi;

		/* a = 2Q + x and b = 2Q - x, limb by limb */
		ca += (uint64_t)t[i] + xi;
		ai = (uint32_t)ca;
		ca >>= 32;
		cb = (uint64_t)t[i] - xi - cb;
		bi = (uint32_t)cb;
		cb = (cb >> 63) & 1;
		/* up <- borrow of a - 2^129, dn <- borrow of 2^129 - b */
		cc = (uint64_t)ai - ci - up;
		up = (uint32_t)(cc >> 63) & 1;
		cc = (uint64_t)ci - bi - dn;
		dn = (uint32_t)(cc >> 63) & 1;
	}
	b = up - dn;
	cc = (uint64_t)d.w[0] + b;
	d.w[0] = (uint32_t)cc;
	cc = (uint64_t)d.w[1] + (uint32_t)-(b >> 31) + (cc >> 32);
	d.w[1] = (uint32_t)cc;
	d.w[2] += (uint32_t)-(b >> 31) + (uint32_t)(cc >> 32);
	return d;
}

fpr
fpr_div(fpr x, fpr y)
{
	return fpr_mul(x, fpr_inv(y));
}

fpr
fpr_sqrt(fpr x)
{
	fpr d;
	uint32_t sx[FXP_N], m[FXP_N], y[FXP_N], prod[U_N + 1];
	uint32_t neg;
	unsigned p, t, th, sh;

	/* m = raw_x << t, t even, lands in [2^94, 2^96); then
	   raw_x * rsqrt(m) >> (110 - t/2) = sqrt(raw_x) * 2^32, the raw
	   result. Negative x is masked to zero, and p is forced to 1 for
	   zero so the shift stays defined. */
	neg = (uint32_t)0 - fxp_is_neg96(x.w);
	for (int i = 0; i < FXP_N; i ++) {
		sx[i] = x.w[i] & ~neg;
	}
	p = fxp_bitlen96_ct(sx);
	p |= ct_mask_eq32((uint32_t)p, 0) & 1u;
	t = 95u - p;
	t += t & 1u;   /* even, p+t in {95,96} */
	fxp_shl96_lt96(m, sx, t);
	fxp_rsqrt_fast(y, m);
	th = t >> 1;
	sh = 110 - th;
	u192_mul96(prod, sx, y);
	prod[U_N] = 0;
	u192_shr_rn96_ct(d.w, prod, sh);
	return d;
}

/* ==================================================================== */
/* Q0.63 l10 helpers.                                                   */

#if FNDSA_SIGN_L10_64
/* a <- a * b with b stored in Q0.63. */
void
fpoly_mul_l10_fft(unsigned logn, fpr *a, const l10_t *b)
{
	size_t hn = (size_t)1 << (logn - 1);
	size_t im = L10_IM_OFF(logn);

	for (size_t i = 0; i < hn; i ++) {
		fpr b_re = l10_to_fpr(b[i]);
		fpr b_im = l10_to_fpr(b[i + im]);
		FPC_MUL_SELF(a[i], a[i + hn], b_re, b_im);
	}
}

void
fpoly_l10_neg(unsigned logn, l10_t *a)
{
	size_t hn = (size_t)1 << (logn - 1);
	size_t im = L10_IM_OFF(logn);

	for (size_t i = 0; i < hn; i ++) {
		a[i] = 0 - a[i];
		a[i + im] = 0 - a[i + im];
	}
}
#endif

#if FNDSA_SIGN_SQRT_FUSED
/* r <- low 96 bits of round(xy * inv_sigma_q71 / 2^n), for the 192-bit
   xy = raw_x * rsqrt and 134 <= n <= 181. Constant time. */
static inline void
sqrt_sigma_q71_shr_rn96_ct(uint32_t *r, const uint32_t *xy,
	uint64_t inv_sigma_q71,
	unsigned n)
{
	uint32_t b0, b1;
	uint32_t p0, p1, p2, p3, p4, p5, p6, p7;
	uint32_t q0, q1, q2, rb, st, inc, stlo;
	uint32_t q04, q14, q24, q05, q15, q25;
	uint32_t rb4, rb5, st4, st5, m4, m5, pm;
	uint64_t cc;
	unsigned wsh, bsh, rbw, rbb;

	b0 = (uint32_t)inv_sigma_q71;
	b1 = (uint32_t)(inv_sigma_q71 >> 32);

	p0 = p1 = p2 = 0;

#define UMUL192_64_STEP(ai, pi, pj, pk)   do { \
		cc = (uint64_t)(ai) * (uint64_t)b0 + (uint64_t)(pi); \
		(pi) = (uint32_t)cc; \
		cc >>= 32; \
		cc += (uint64_t)(ai) * (uint64_t)b1 + (uint64_t)(pj); \
		(pj) = (uint32_t)cc; \
		(pk) = (uint32_t)(cc >> 32); \
	} while (0)

	UMUL192_64_STEP(xy[0], p0, p1, p2);
	UMUL192_64_STEP(xy[1], p1, p2, p3);
	UMUL192_64_STEP(xy[2], p2, p3, p4);
	UMUL192_64_STEP(xy[3], p3, p4, p5);
	UMUL192_64_STEP(xy[4], p4, p5, p6);
	UMUL192_64_STEP(xy[5], p5, p6, p7);

#undef UMUL192_64_STEP

	wsh = n >> 5;
	bsh = n & 31;
	m4 = ct_mask_eq32((uint32_t)wsh, 4);
	m5 = ct_mask_eq32((uint32_t)wsh, 5);
	q04 = (p4 >> bsh) | ((p5 << 1) << (31 - bsh));
	q14 = (p5 >> bsh) | ((p6 << 1) << (31 - bsh));
	q24 = (p6 >> bsh) | ((p7 << 1) << (31 - bsh));
	q05 = (p5 >> bsh) | ((p6 << 1) << (31 - bsh));
	q15 = (p6 >> bsh) | ((p7 << 1) << (31 - bsh));
	q25 = p7 >> bsh;
	q0 = (q04 & m4) | (q05 & m5);
	q1 = (q14 & m4) | (q15 & m5);
	q2 = (q24 & m4) | (q25 & m5);

	rbw = (n - 1) >> 5;
	rbb = (n - 1) & 31;
	pm = ((uint32_t)1 << rbb) - 1u;
	m4 = ct_mask_eq32((uint32_t)rbw, 4);
	m5 = ct_mask_eq32((uint32_t)rbw, 5);
	rb4 = (p4 >> rbb) & 1;
	rb5 = (p5 >> rbb) & 1;
	stlo = p0 | p1 | p2 | p3;
	st4 = stlo | (p4 & pm);
	st5 = stlo | p4 | (p5 & pm);
	rb = (rb4 & (m4 & 1u)) | (rb5 & (m5 & 1u));
	st = (st4 & m4) | (st5 & m5);
	st = (st | (0 - st)) >> 31;
	inc = FXP_RN_INC(rb, st, q0 & 1);
	cc = (uint64_t)q0 + inc;
	r[0] = (uint32_t)cc;
	cc = (uint64_t)q1 + (cc >> 32);
	r[1] = (uint32_t)cc;
	r[2] = q2 + (uint32_t)(cc >> 32);
}

/* FPR_SQRT_MUL_INV_SIGMA without the gamma bounds: same domain as
   fpr_sqrt. */
void
fndsa_fxp96_sqrt_mul_inv_sigma(uint32_t *d, const uint32_t *x, unsigned logn)
{
	uint32_t sx[FXP_N], m[FXP_N], y[FXP_N], xy[U_N];
	uint64_t inv_sigma;
	uint32_t neg;
	unsigned p, t, th, sh;

	/* fpr_sqrt with the multiply by INV_SIGMA[logn] folded in:
	   raw_x * rsqrt(m) * INV_SIGMA_Q71 >> (181 - t/2)
	   = sqrt(raw_x) * INV_SIGMA[logn] * 2^32. */
	neg = (uint32_t)0 - fxp_is_neg96(x);
	for (int i = 0; i < FXP_N; i ++) {
		sx[i] = x[i] & ~neg;
	}
	p = fxp_bitlen96_ct(sx);
	p |= ct_mask_eq32((uint32_t)p, 0) & 1u;
	t = 95u - p;
	t += t & 1u;   /* even, p+t in {95,96} */
	fxp_shl96_lt96(m, sx, t);
	fxp_rsqrt_fast(y, m);
	u192_mul96(xy, sx, y);
	inv_sigma = INV_SIGMA_Q71[logn];
	th = t >> 1;
	sh = 181 - th;
	sqrt_sigma_q71_shr_rn96_ct(d, xy, inv_sigma, sh);
}
#endif
