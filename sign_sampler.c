/*
 * Gaussian sampling.
 */

#include "sign_inner.h"

/* fpr_u, the constants INV_2SQRSIGMA0, LOG2, INV_LOG2 and the tables
   INV_SIGMA and SIGMA_MIN come from the fpr representation. */

#if !FNDSA_ASM_CORTEXM4
/* Distribution for gaussian0() (this is the RCDT table from the
   specification, expressed in base 2^24). */
static const uint32_t GAUSS0[][3] = {
	{ 10745844,  3068844,  3741698 },
	{  5559083,  1580863,  8248194 },
	{  2260429, 13669192,  2736639 },
	{   708981,  4421575, 10046180 },
	{   169348,  7122675,  4136815 },
	{    30538, 13063405,  7650655 },
	{     4132, 14505003,  7826148 },
	{      417, 16768101, 11363290 },
	{       31,  8444042,  8086568 },
	{        1, 12844466,   265321 },
	{        0,  1232676, 13644283 },
	{        0,    38047,  9111839 },
	{        0,      870,  6138264 },
	{        0,       14, 12545723 },
	{        0,        0,  3104126 },
	{        0,        0,    28824 },
	{        0,        0,      198 },
	{        0,        0,        1 }
};
#endif

/* We access the PRNG through macros so that they can be overridden by some
   compatiblity tests with the original Falcon implementation. */
#ifndef prng_init
#if FNDSA_SHAKE256X4
#define prng_init       shake256x4_init
#define prng_next_u8    shake256x4_next_u8
#define prng_next_u64   shake256x4_next_u64
#else
#define prng_init(pc, seed, seed_len)   do { \
		shake_init(pc, 256); \
		shake_inject(pc, seed, seed_len); \
		shake_flip(pc); \
	} while (0)
#define prng_next_u8    shake_next_u8
#define prng_next_u64   shake_next_u64
#endif
#endif

/* see sign_inner.h */
void
sampler_init(sampler_state *ss, unsigned logn,
	const void *seed, size_t seed_len)
{
	prng_init(&ss->pc, seed, seed_len);
	ss->logn = logn;
}

#if FNDSA_ASM_CORTEXM4
int32_t fndsa_gaussian0_helper(uint64_t lo, uint32_t hi);
#endif

static inline int32_t
gaussian0(sampler_state *ss)
{
	/* Get a random 72-bit value, into three 24-bit limbs (v0..v2). */
	uint64_t lo = prng_next_u64(&ss->pc);
	uint32_t hi = prng_next_u8(&ss->pc);
#if FNDSA_ASM_CORTEXM4
	return fndsa_gaussian0_helper(lo, hi);
#else
	uint32_t v0 = (uint32_t)lo & 0xFFFFFF;
	uint32_t v1 = (uint32_t)(lo >> 24) & 0xFFFFFF;
	uint32_t v2 = (uint32_t)(lo >> 48) | (hi << 16);

	/* Sampled value is z such that v0..v2 is lower than the first
	   z elements of the table. */
	int32_t z = 0;
	for (size_t i = 0; i < (sizeof GAUSS0) / sizeof(GAUSS0[0]); i ++) {
		uint32_t cc;
		cc = (v0 - GAUSS0[i][2]) >> 31;
		cc = (v1 - GAUSS0[i][1] - cc) >> 31;
		cc = (v2 - GAUSS0[i][0] - cc) >> 31;
		z += (int32_t)cc;
	}
	return z;
#endif
}

#if FNDSA_SSE2
/* ========================= SSE2 IMPLEMENTATION ========================= */

/* Input: 0 <= x < log(2)
   Output: trunc(x*2^63) */
TARGET_SSE2
static inline int64_t
mtwop63(__m128d x)
{
#if FNDSA_64
	static const union {
		fpr f[2];
		__m128d x;
	} twop63 = { {
		FPR(4503599627370496, 11),
		FPR(4503599627370496, 11),
	} };

	return _mm_cvttsd_si64(_mm_mul_sd(x, twop63.x));
#else
	/* 32-bit x86 does not have an SSE2 opcode to convert floating-point
	   values to 64-bit integers, only 32-bit signed integers. We must
	   do the conversion in three steps with factor 2^21. */
	static const union {
		fpr f[2];
		__m128d x;
	} twop21 = { {
		FPR(4503599627370496, -31),
		FPR(4503599627370496, -31),
	} };
	x = _mm_mul_sd(x, twop21.x);
	int32_t z2 = _mm_cvttsd_si32(x);
	x = _mm_sub_sd(x, _mm_cvtsi32_sd(_mm_setzero_pd(), z2));
	x = _mm_mul_sd(x, twop21.x);
	int32_t z1 = _mm_cvttsd_si32(x);
	x = _mm_sub_sd(x, _mm_cvtsi32_sd(_mm_setzero_pd(), z1));
	x = _mm_mul_sd(x, twop21.x);
	int32_t z0 = _mm_cvttsd_si32(x);
	return ((int64_t)z2 << 42) + ((int64_t)z1 << 21) + (int64_t)z0;
#endif
}

/* Compute ccs*exp(-x)*2^63, rounded to an integer. This function assumes
   that 0 <= x < log(2), and 0 <= ccs <= 1. It returns a value in [0,2^63]. */
TARGET_SSE2
static inline uint64_t
expm_p63(__m128d x, __m128d ccs)
{
	/* The polynomial approximation of exp(-x) is from FACCT:
	      https://eprint.iacr.org/2018/1234
	   Specifically, the values are extracted from the implementation
	   referenced by the FACCT paper, available at:
	      https://github.com/raykzhao/gaussian  */
	static const uint64_t EXPM_COEFFS[] = {
		0x00000004741183A3,
		0x00000036548CFC06,
		0x0000024FDCBF140A,
		0x0000171D939DE045,
		0x0000D00CF58F6F84,
		0x000680681CF796E3,
		0x002D82D8305B0FEA,
		0x011111110E066FD0,
		0x0555555555070F00,
		0x155555555581FF00,
		0x400000000002B400,
		0x7FFFFFFFFFFF4800,
		0x8000000000000000
	};

	uint64_t y = EXPM_COEFFS[0];
	uint64_t z = (uint64_t)mtwop63(x) << 1;
	uint64_t w = (uint64_t)mtwop63(ccs) << 1;
#if FNDSA_64
	/* On 64-bit x86, we have 64x64->128 multiplication, then we can use
	   it, it's normally constant-time.
	   MSVC uses a different syntax for this operation. */
	for (size_t i = 1; i < (sizeof EXPM_COEFFS) / sizeof(uint64_t); i ++) {
#if defined _MSC_VER
		y = EXPM_COEFFS[i] - __umulh(z, y);
#else
		unsigned __int128 c =
			(unsigned __int128)z * (unsigned __int128)y;
		y = EXPM_COEFFS[i] - (uint64_t)(c >> 64);
#endif
	}
#if defined _MSC_VER
	y = __umulh(w, y);
#else
	y = (uint64_t)(((unsigned __int128)w * (unsigned __int128)y) >> 64);
#endif
#else
	/* On 32-bit x86, no 64x64->128 multiplication, we must use
	   four 32x32->64 multiplications. */
	uint32_t z0 = (uint32_t)z, z1 = (uint32_t)(z >> 32);
	uint32_t w0 = (uint32_t)w, w1 = (uint32_t)(w >> 32);

	for (size_t i = 1; i < (sizeof EXPM_COEFFS) / sizeof(uint64_t); i ++) {
		uint32_t y0 = (uint32_t)y, y1 = (uint32_t)(y >> 32);
		uint64_t f = (uint64_t)z0 * (uint64_t)y0;
		uint64_t a = (uint64_t)z0 * (uint64_t)y1 + (f >> 32);
		uint64_t b = (uint64_t)z1 * (uint64_t)y0;
		uint64_t c = (a >> 32) + (b >> 32)
			+ (((uint64_t)(uint32_t)a
			  + (uint64_t)(uint32_t)b) >> 32)
			+ (uint64_t)z1 * (uint64_t)y1;
		y = EXPM_COEFFS[i] - c;
	}
	uint32_t y0 = (uint32_t)y, y1 = (uint32_t)(y >> 32);
	uint64_t f = (uint64_t)w0 * (uint64_t)y0;
	uint64_t a = (uint64_t)w0 * (uint64_t)y1 + (f >> 32);
	uint64_t b = (uint64_t)w1 * (uint64_t)y0;
	y = (a >> 32) + (b >> 32)
		+ (((uint64_t)(uint32_t)a + (uint64_t)(uint32_t)b) >> 32)
		+ (uint64_t)w1 * (uint64_t)y1;
#endif
	return y;
}

/* Sample a bit with probability ccs*exp(-x) (for x >= 0). */
TARGET_SSE2
static inline int
ber_exp(sampler_state *ss, __m128d x, __m128d ccs)
{
	static union { fpr f[2]; __m128d x; }
		LOG2_u = { { LOG2, LOG2 } },
		INV_LOG2_u = { { INV_LOG2, INV_LOG2 } };

	/* Reduce x modulo log(2): x = s*log(2) + r, with s an integer,
	   and 0 <= r < log(2). We can use a truncating conversion because
	   x >= 0. Moreover, x is small, so we can stick to 32-bit values. */
	int32_t si = _mm_cvttsd_si32(_mm_mul_sd(x, INV_LOG2_u.x));
	__m128d r = _mm_sub_sd(x,
		_mm_mul_sd(_mm_cvtsi32_sd(_mm_setzero_pd(), si), LOG2_u.x));

	/* If s >= 64, sigma = 1.2, r = 0 and b = 1, then we get s >= 64
	   if the half-Gaussian produced z >= 13, which happens with
	   probability about 2^(-32). When s >= 64, ber_exp() will return
	   true with probability less than 2^(-64), so we can simply
	   saturate s at 63 (the bias introduced here is lower than 2^(-96),
	   and would require about 2^192 samplings to be detectable, which
	   is way beyond the formal bound of 2^64 signatures with the
	   same key. */
	uint32_t s = (uint32_t)si;
	s |= (uint32_t)(63 - s) >> 26;
	s &= 63;

	/* Compute ccs*exp(-x). Since x = s*log(2) + r, we compute
	   ccs*exp(-r)/2^s. We know that 0 <= r < log(2), so we can
	   use expm_p63(), which yields a result scaled by 63 bits. We
	   scale it up 1 bit further, then right-shift by s bits.

	   We subtract 1 to make sure that the value fits on 64 bits
	   (i.e. if r = 0 then we may get 2^64 and we prefer 2^64-1
	   in that case, to avoid the overflow). The bias is negligible
	   since expm_p63() has precision only 51 bits or so. */
	uint64_t z = fpr_ursh((expm_p63(r, ccs) << 1) - 1, s);

	/* Sample a bit. We lazily compare the value z with a uniform 64-bit
	   integer, consuming only as many bytes as necessary. Since the PRNG
	   is cryptographically strong, we leak no information from the
	   conditional jumps below. */
	for (int i = 56; i >= 0; i -= 8) {
		unsigned w = prng_next_u8(&ss->pc);
		unsigned bz = (unsigned)(z >> i) & 0xFF;
		if (w != bz) {
			return w < bz;
		}
	}
	return 0;
}

TARGET_SSE2
static int32_t
sampler_next_sse2(sampler_state *ss, __m128d mu, __m128d isigma)
{
	static union { fpr f[2]; __m128d x; }
		HALF_u = { {
			FPR(4503599627370496, -53),
			FPR(4503599627370496, -53)
		} },
		INV_2SQRSIGMA0_u = { {
			INV_2SQRSIGMA0, INV_2SQRSIGMA0
		} };

	/* Split center mu into s + r, for an integer s, and 0 <= r < 1. */
	int32_t s = _mm_cvttsd_si32(mu);
	s -= _mm_comilt_sd(mu, _mm_cvtsi32_sd(_mm_setzero_pd(), s));
	__m128d r = _mm_sub_sd(mu, _mm_cvtsi32_sd(_mm_setzero_pd(), s));

	/* dss = 1/(2*sigma^2) = 0.5*(isigma^2)  */
	__m128d dss = _mm_mul_sd(_mm_mul_sd(isigma, isigma), HALF_u.x);

	/* css = sigma_min / sigma = sigma_min * isigma  */
	__m128d ccs = _mm_mul_sd(isigma,
		_mm_load_sd((const double *)SIGMA_MIN + ss->logn));

	/* We sample on centre r. */
	for (;;) {
		/* Sample z for a Gaussian distribution (non-negative only),
		   then get a random bit b to turn the sampling into a
		   bimodal distribution (we use z+1 if b = 1, or -z
		   otherwise). */
		int32_t z0 = gaussian0(ss);
		int32_t b = prng_next_u8(&ss->pc) & 1;
		int32_t z = b + ((b << 1) - 1) * z0;

		/* Rejection sampling. We want a Gaussian centred on r,
		   but we sampled against a bimodal distribution (with
		   "centres" at 0 and 1). However, we know that z is
		   always in the range where our sampling distribution is
		   greater than the Gaussian distribution, so rejection works.

		   We got z from distribution:
		      G(z) = exp(-((z-b)^2)/(2*sigma0^2))
		   We target distribution:
		      S(z) = exp(-((z-r)^2)/(2*signa^2))
		   Rejection sampling works by keeping the value z with
		   probability S(z)/G(z), and starting again otherwise.
		   This requires S(z) <= G(z), which is the case here.
		   Thus, we simply need to keep our z with probability:
		      P = exp(-x)
		   where:
		      x = ((z-r)^2)/(2*sigma^2) - ((z-b)^2)/(2*sigma0^2)
		   Here, we scale up the Bernouilli distribution, which
		   makes rejection more probable, but also makes the
		   rejection rate sufficiently decorrelated from the Gaussian
		   centre and standard deviation, so that measurement of the
		   rejection rate does not leak enough usable information
		   to attackers (which is how the implementation can claim
		   to be "constant-time").  */
		__m128d x = _mm_sub_sd(_mm_cvtsi32_sd(_mm_setzero_pd(), z), r);
		x = _mm_mul_sd(_mm_mul_sd(x, x), dss);
		x = _mm_sub_sd(x, _mm_mul_sd(
			_mm_cvtsi32_sd(_mm_setzero_pd(), z0 * z0),
			INV_2SQRSIGMA0_u.x));
		if (ber_exp(ss, x, ccs)) {
			return s + z;
		}
	}
}

/* see sign_inner.h */
TARGET_SSE2
int32_t
sampler_next(sampler_state *ss, fpr mu, fpr isigma)
{
	return sampler_next_sse2(ss,
		_mm_load_sd((double *)&mu),
		_mm_load_sd((double *)&isigma));
}

#elif FNDSA_NEON
/* ========================= NEON IMPLEMENTATION ========================= */

/* Input: 0 <= x < log(2)
   Output: trunc(x*2^63) */
TARGET_NEON
static inline int64_t
mtwop63(float64x1_t x)
{
	static const fpr_u twop63 = { FPR(4503599627370496, 11) };
	return vget_lane_s64(vcvt_s64_f64(vmul_f64(x, twop63.v)), 0);
}

/* Compute ccs*exp(-x)*2^63, rounded to an integer. This function assumes
   that 0 <= x < log(2), and 0 <= ccs <= 1. It returns a value in [0,2^63]. */
TARGET_NEON
static inline uint64_t
expm_p63(float64x1_t x, float64x1_t ccs)
{
	/* The polynomial approximation of exp(-x) is from FACCT:
	      https://eprint.iacr.org/2018/1234
	   Specifically, the values are extracted from the implementation
	   referenced by the FACCT paper, available at:
	      https://github.com/raykzhao/gaussian  */
	static const uint64_t EXPM_COEFFS[] = {
		0x00000004741183A3,
		0x00000036548CFC06,
		0x0000024FDCBF140A,
		0x0000171D939DE045,
		0x0000D00CF58F6F84,
		0x000680681CF796E3,
		0x002D82D8305B0FEA,
		0x011111110E066FD0,
		0x0555555555070F00,
		0x155555555581FF00,
		0x400000000002B400,
		0x7FFFFFFFFFFF4800,
		0x8000000000000000
	};

	uint64_t y = EXPM_COEFFS[0];
	uint64_t z = (uint64_t)mtwop63(x) << 1;
	uint64_t w = (uint64_t)mtwop63(ccs) << 1;

	/* We assume here that 64x64->128 multiplications are constant-time,
	   which is not exactly true on some aarch64 systems (e.g. ARM
	   Cortex A53 and A55 return the result one cycle earlier when
	   the operands fit on 32 bits).
	   ARM compiler calls the 128-bit type '__uint128_t' while GCC
	   and Clang use 'unsigned __int128'. */
	for (size_t i = 1; i < (sizeof EXPM_COEFFS) / sizeof(uint64_t); i ++) {
#if defined __GNUC__ || defined __clang__
		unsigned __int128 c =
			(unsigned __int128)z * (unsigned __int128)y;
		y = EXPM_COEFFS[i] - (uint64_t)(c >> 64);
#else
		__uint128_t c = (__uint128_t)z * (__uint128_t)y;
		y = EXPM_COEFFS[i] - (uint64_t)(c >> 64);
#endif
	}
#if defined __GNUC__ || defined __clang__
	y = (uint64_t)(((unsigned __int128)w * (unsigned __int128)y) >> 64);
#else
	y = (uint64_t)(((__uint128_t)w * (__uint128_t)y) >> 64);
#endif
	return y;
}

/* Sample a bit with probability ccs*exp(-x) (for x >= 0). */
TARGET_NEON
static inline int
ber_exp(sampler_state *ss, float64x1_t x, float64x1_t ccs)
{
	static const fpr_u LOG2_u = { LOG2 };
	static const fpr_u INV_LOG2_u = { INV_LOG2 };

	/* Reduce x modulo log(2): x = s*log(2) + r, with s an integer,
	   and 0 <= r < log(2). We can use a truncating conversion because
	   x >= 0. Moreover, x is small, so we can stick to 32-bit values. */
	int32_t si = (int32_t)vget_lane_s64(
		vcvt_s64_f64(vmul_f64(x, INV_LOG2_u.v)), 0);
	float64x1_t r = vsub_f64(x,
		vmul_f64(vcvt_f64_s64(vcreate_s64(si)), LOG2_u.v));

	/* If s >= 64, sigma = 1.2, r = 0 and b = 1, then we get s >= 64
	   if the half-Gaussian produced z >= 13, which happens with
	   probability about 2^(-32). When s >= 64, ber_exp() will return
	   true with probability less than 2^(-64), so we can simply
	   saturate s at 63 (the bias introduced here is lower than 2^(-96),
	   and would require about 2^192 samplings to be detectable, which
	   is way beyond the formal bound of 2^64 signatures with the
	   same key. */
	uint32_t s = (uint32_t)si;
	s |= (uint32_t)(63 - s) >> 26;
	s &= 63;

	/* Compute ccs*exp(-x). Since x = s*log(2) + r, we compute
	   ccs*exp(-r)/2^s. We know that 0 <= r < log(2), so we can
	   use expm_p63(), which yields a result scaled by 63 bits. We
	   scale it up 1 bit further, then right-shift by s bits.

	   We subtract 1 to make sure that the value fits on 64 bits
	   (i.e. if r = 0 then we may get 2^64 and we prefer 2^64-1
	   in that case, to avoid the overflow). The bias is negligible
	   since expm_p63() has precision only 51 bits or so. */
	uint64_t z = fpr_ursh((expm_p63(r, ccs) << 1) - 1, s);

	/* Sample a bit. We lazily compare the value z with a uniform 64-bit
	   integer, consuming only as many bytes as necessary. Since the PRNG
	   is cryptographically strong, we leak no information from the
	   conditional jumps below. */
	for (int i = 56; i >= 0; i -= 8) {
		unsigned w = prng_next_u8(&ss->pc);
		unsigned bz = (unsigned)(z >> i) & 0xFF;
		if (w != bz) {
			return w < bz;
		}
	}
	return 0;
}

TARGET_NEON
static int32_t
sampler_next_neon(sampler_state *ss, float64x1_t mu, float64x1_t isigma)
{
	static const fpr_u HALF_u = { FPR(4503599627370496, -53) };
	static const fpr_u INV_2SQRSIGMA0_u = { INV_2SQRSIGMA0 };

	/* Split center mu into s + r, for an integer s, and 0 <= r < 1. */
	int32_t s = (int32_t)vcvtmd_s64_f64(mu);
	float64x1_t r = vsub_f64(mu, vcvt_f64_s64(vcreate_s64(s)));

	/* dss = 1/(2*sigma^2) = 0.5*(isigma^2)  */
	float64x1_t dss = vmul_f64(vmul_f64(isigma, isigma), HALF_u.v);

	/* css = sigma_min / sigma = sigma_min * isigma  */
	float64x1_t ccs = vmul_f64(isigma, SIGMA_MIN[ss->logn].v);

	/* We sample on centre r. */
	for (;;) {
		/* Sample z for a Gaussian distribution (non-negative only),
		   then get a random bit b to turn the sampling into a
		   bimodal distribution (we use z+1 if b = 1, or -z
		   otherwise). */
		int32_t z0 = gaussian0(ss);
		int32_t b = prng_next_u8(&ss->pc) & 1;
		int32_t z = b + ((b << 1) - 1) * z0;

		/* Rejection sampling. We want a Gaussian centred on r,
		   but we sampled against a bimodal distribution (with
		   "centres" at 0 and 1). However, we know that z is
		   always in the range where our sampling distribution is
		   greater than the Gaussian distribution, so rejection works.

		   We got z from distribution:
		      G(z) = exp(-((z-b)^2)/(2*sigma0^2))
		   We target distribution:
		      S(z) = exp(-((z-r)^2)/(2*signa^2))
		   Rejection sampling works by keeping the value z with
		   probability S(z)/G(z), and starting again otherwise.
		   This requires S(z) <= G(z), which is the case here.
		   Thus, we simply need to keep our z with probability:
		      P = exp(-x)
		   where:
		      x = ((z-r)^2)/(2*sigma^2) - ((z-b)^2)/(2*sigma0^2)
		   Here, we scale up the Bernouilli distribution, which
		   makes rejection more probable, but also makes the
		   rejection rate sufficiently decorrelated from the Gaussian
		   centre and standard deviation, so that measurement of the
		   rejection rate does not leak enough usable information
		   to attackers (which is how the implementation can claim
		   to be "constant-time").  */
		float64x1_t x = vsub_f64(vcvt_f64_s64(vcreate_s64(z)), r);
		x = vmul_f64(vmul_f64(x, x), dss);
		x = vsub_f64(x, vmul_f64(
			vcvt_f64_s64(vcreate_s64(z0 * z0)),
			INV_2SQRSIGMA0_u.v));
		if (ber_exp(ss, x, ccs)) {
			return s + z;
		}
	}
}

/* see sign_inner.h */
TARGET_NEON
int32_t
sampler_next(sampler_state *ss, fpr mu, fpr isigma)
{
	return sampler_next_neon(ss,
		vld1_f64((const float64_t *)&mu),
		vld1_f64((const float64_t *)&isigma));
}

#elif FNDSA_RV64D
/* ========================= RISC-V IMPLEMENTATION ======================= */

/* Input: 0 <= x < log(2)
   Output: trunc(x*2^63) */
TARGET_NEON
static inline int64_t
mtwop63(f64 x)
{
	static const fpr_u twop63 = { FPR(4503599627370496, 11) };
	return f64_trunc(f64_mul(x, twop63.v));
}

static inline uint64_t
expm_p63(f64 x, f64 ccs)
{
	/* The polynomial approximation of exp(-x) is from FACCT:
	      https://eprint.iacr.org/2018/1234
	   Specifically, the values are extracted from the implementation
	   referenced by the FACCT paper, available at:
	      https://github.com/raykzhao/gaussian  */
	static const uint64_t EXPM_COEFFS[] = {
		0x00000004741183A3,
		0x00000036548CFC06,
		0x0000024FDCBF140A,
		0x0000171D939DE045,
		0x0000D00CF58F6F84,
		0x000680681CF796E3,
		0x002D82D8305B0FEA,
		0x011111110E066FD0,
		0x0555555555070F00,
		0x155555555581FF00,
		0x400000000002B400,
		0x7FFFFFFFFFFF4800,
		0x8000000000000000
	};

	uint64_t y = EXPM_COEFFS[0];
	uint64_t z = (uint64_t)mtwop63(x) << 1;
	uint64_t w = (uint64_t)mtwop63(ccs) << 1;

	/* We assume here that 64x64->128 multiplications are constant-time
	   and that the compiler is GCC/Clang compatible (i.e. supports
	   the 'unsigned __int128' type). */
	for (size_t i = 1; i < (sizeof EXPM_COEFFS) / sizeof(uint64_t); i ++) {
		unsigned __int128 c =
			(unsigned __int128)z * (unsigned __int128)y;
		y = EXPM_COEFFS[i] - (uint64_t)(c >> 64);
	}
	y = (uint64_t)(((unsigned __int128)w * (unsigned __int128)y) >> 64);
	return y;
}

/* Sample a bit with probability ccs*exp(-x) (for x >= 0). */
static inline int
ber_exp(sampler_state *ss, f64 x, f64 ccs)
{
	static const fpr_u LOG2_u = { LOG2 };
	static const fpr_u INV_LOG2_u = { INV_LOG2 };

	/* Reduce x modulo log(2): x = s*log(2) + r, with s an integer,
	   and 0 <= r < log(2). We can use f64_trunc() because x >= 0. */
	int32_t si = (int32_t)f64_trunc(f64_mul(x, INV_LOG2_u.v));
	f64 r = f64_sub(x, f64_mul(f64_of(si), LOG2_u.v));

	/* If s >= 64, sigma = 1.2, r = 0 and b = 1, then we get s >= 64
	   if the half-Gaussian produced z >= 13, which happens with
	   probability about 2^(-32). When s >= 64, ber_exp() will return
	   true with probability less than 2^(-64), so we can simply
	   saturate s at 63 (the bias introduced here is lower than 2^(-96),
	   and would require about 2^192 samplings to be detectable, which
	   is way beyond the formal bound of 2^64 signatures with the
	   same key. */
	uint32_t s = (uint32_t)si;
	s |= (uint32_t)(63 - s) >> 26;
	s &= 63;

	/* Compute ccs*exp(-x). Since x = s*log(2) + r, we compute
	   ccs*exp(-r)/2^s. We know that 0 <= r < log(2), so we can
	   use expm_p63(), which yields a result scaled by 63 bits. We
	   scale it up 1 bit further, then right-shift by s bits.

	   We subtract 1 to make sure that the value fits on 64 bits
	   (i.e. if r = 0 then we may get 2^64 and we prefer 2^64-1
	   in that case, to avoid the overflow). The bias is negligible
	   since expm_p63() has precision only 51 bits or so. */
	uint64_t z = fpr_ursh((expm_p63(r, ccs) << 1) - 1, s);

	/* Sample a bit. We lazily compare the value z with a uniform 64-bit
	   integer, consuming only as many bytes as necessary. Since the PRNG
	   is cryptographically strong, we leak no information from the
	   conditional jumps below. */
	for (int i = 56; i >= 0; i -= 8) {
		unsigned w = prng_next_u8(&ss->pc);
		unsigned bz = (unsigned)(z >> i) & 0xFF;
		if (w != bz) {
			return w < bz;
		}
	}
	return 0;
}

static int32_t
sampler_next_rv64d(sampler_state *ss, f64 mu, f64 isigma)
{
	static const fpr_u INV_2SQRSIGMA0_u = { INV_2SQRSIGMA0 };

	/* Split center mu into s + r, for an integer s, and 0 <= r < 1. */
	int64_t s = f64_floor(mu);
	f64 r = f64_sub(mu, f64_of(s));

	/* dss = 1/(2*sigma^2) = 0.5*(isigma^2)  */
	f64 dss = f64_half(f64_sqr(isigma));

	/* css = sigma_min / sigma = sigma_min * isigma  */
	f64 ccs = f64_mul(isigma, SIGMA_MIN[ss->logn].v);

	/* We sample on centre r. */
	for (;;) {
		/* Sample z for a Gaussian distribution (non-negative only),
		   then get a random bit b to turn the sampling into a
		   bimodal distribution (we use z+1 if b = 1, or -z
		   otherwise). */
		int32_t z0 = gaussian0(ss);
		int32_t b = prng_next_u8(&ss->pc) & 1;
		int32_t z = b + ((b << 1) - 1) * z0;

		/* Rejection sampling. We want a Gaussian centred on r,
		   but we sampled against a bimodal distribution (with
		   "centres" at 0 and 1). However, we know that z is
		   always in the range where our sampling distribution is
		   greater than the Gaussian distribution, so rejection works.

		   We got z from distribution:
		      G(z) = exp(-((z-b)^2)/(2*sigma0^2))
		   We target distribution:
		      S(z) = exp(-((z-r)^2)/(2*signa^2))
		   Rejection sampling works by keeping the value z with
		   probability S(z)/G(z), and starting again otherwise.
		   This requires S(z) <= G(z), which is the case here.
		   Thus, we simply need to keep our z with probability:
		      P = exp(-x)
		   where:
		      x = ((z-r)^2)/(2*sigma^2) - ((z-b)^2)/(2*sigma0^2)
		   Here, we scale up the Bernouilli distribution, which
		   makes rejection more probable, but also makes the
		   rejection rate sufficiently decorrelated from the Gaussian
		   centre and standard deviation, so that measurement of the
		   rejection rate does not leak enough usable information
		   to attackers (which is how the implementation can claim
		   to be "constant-time").  */
		f64 x = f64_mul(f64_sqr(f64_sub(f64_of(z), r)), dss);
		x = f64_sub(x, f64_mul(f64_of(z0 * z0), INV_2SQRSIGMA0_u.v));
		if (ber_exp(ss, x, ccs)) {
			return (int32_t)s + z;
		}
	}
}

/* see sign_inner.h */
int32_t
sampler_next(sampler_state *ss, fpr mu, fpr isigma)
{
	return sampler_next_rv64d(ss, f64_from_raw(mu), f64_from_raw(isigma));
}

#else
/* ========================= PLAIN IMPLEMENTATION ======================== */

static inline uint64_t
expm_p63(fpr x, fpr ccs)
{
	/* The polynomial approximation of exp(-x) is from FACCT:
	      https://eprint.iacr.org/2018/1234
	   Specifically, the values are extracted from the implementation
	   referenced by the FACCT paper, available at:
	      https://github.com/raykzhao/gaussian  */
	static const uint64_t EXPM_COEFFS[] = {
		0x00000004741183A3,
		0x00000036548CFC06,
		0x0000024FDCBF140A,
		0x0000171D939DE045,
		0x0000D00CF58F6F84,
		0x000680681CF796E3,
		0x002D82D8305B0FEA,
		0x011111110E066FD0,
		0x0555555555070F00,
		0x155555555581FF00,
		0x400000000002B400,
		0x7FFFFFFFFFFF4800,
		0x8000000000000000
	};

	/* TODO: maybe use 64x64->128 multiplications if available? It
	   is a bit tricky to decide, because the plain code is used for
	   unknown architectures, and we do not know if the larger
	   multiplication is constant-time (it often happens that it
	   is not). */

	uint64_t y = EXPM_COEFFS[0];
	uint64_t z = fpr_p63(x) << 1;
	uint32_t z0 = (uint32_t)z, z1 = (uint32_t)(z >> 32);
#if FNDSA_ASM_CORTEXM4
#pragma GCC unroll 12
#endif
	for (size_t i = 1; i < (sizeof EXPM_COEFFS) / sizeof(uint64_t); i ++) {
		uint32_t y0 = (uint32_t)y, y1 = (uint32_t)(y >> 32);
#if FNDSA_ASM_CORTEXM4
		uint32_t tt, r0, r1;
		__asm__(
			"umull	%0, %2, %3, %5\n\t"
			"umull	%0, %1, %3, %6\n\t"
			"umaal	%2, %0, %4, %5\n\t"
			"umaal	%0, %1, %4, %6\n\t"
			: "=&r" (r0), "=&r" (r1), "=&r" (tt)
			: "r" (y0), "r" (y1), "r" (z0), "r" (z1));
		y = EXPM_COEFFS[i] - ((uint64_t)r0 | ((uint64_t)r1 << 32));
#else
		uint64_t f = (uint64_t)z0 * (uint64_t)y0;
		uint64_t a = (uint64_t)z0 * (uint64_t)y1 + (f >> 32);
		uint64_t b = (uint64_t)z1 * (uint64_t)y0;
		uint64_t c = (a >> 32) + (b >> 32)
			+ (((uint64_t)(uint32_t)a
			  + (uint64_t)(uint32_t)b) >> 32)
			+ (uint64_t)z1 * (uint64_t)y1;
		y = EXPM_COEFFS[i] - c;
#endif
	}

	/* The scaling factor must be applied at the end. Since y is now
	   in fixed-point notation, we have to convert the factor to the
	   same format, and we do an extra integer multiplication. */
	uint64_t w = fpr_p63(ccs) << 1;
	uint32_t w0 = (uint32_t)w, w1 = (uint32_t)(w >> 32);
	uint32_t y0 = (uint32_t)y, y1 = (uint32_t)(y >> 32);
#if FNDSA_ASM_CORTEXM4
	uint32_t tt, r0, r1;
	__asm__(
		"umull	%0, %2, %3, %5\n\t"
		"umull	%0, %1, %3, %6\n\t"
		"umaal	%2, %0, %4, %5\n\t"
		"umaal	%0, %1, %4, %6\n\t"
		: "=&r" (r0), "=&r" (r1), "=&r" (tt)
		: "r" (y0), "r" (y1), "r" (w0), "r" (w1));
	y = (uint64_t)r0 | ((uint64_t)r1 << 32);
#else
	uint64_t f = (uint64_t)w0 * (uint64_t)y0;
	uint64_t a = (uint64_t)w0 * (uint64_t)y1 + (f >> 32);
	uint64_t b = (uint64_t)w1 * (uint64_t)y0;
	y = (a >> 32) + (b >> 32)
		+ (((uint64_t)(uint32_t)a + (uint64_t)(uint32_t)b) >> 32)
		+ (uint64_t)w1 * (uint64_t)y1;
#endif
	return y;
}

/* Sample a bit with probability ccs*exp(-x) (for x >= 0). */
static inline int
ber_exp(sampler_state *ss, fpr x, fpr ccs)
{
	/* Reduce x modulo log(2): x = s*log(2) + r, with s an integer,
	   and 0 <= r < log(2). We can use fpr_trunc() because x >= 0
	   (fpr_trunc() is presumably a bit faster than fpr_floor()). */
	int32_t si = (int32_t)fpr_trunc(fpr_mul(x, INV_LOG2));
	fpr r = fpr_sub(x, fpr_mul(fpr_of32(si), LOG2));

	/* If s >= 64, sigma = 1.2, r = 0 and b = 1, then we get s >= 64
	   if the half-Gaussian produced z >= 13, which happens with
	   probability about 2^(-32). When s >= 64, ber_exp() will return
	   true with probability less than 2^(-64), so we can simply
	   saturate s at 63 (the bias introduced here is lower than 2^(-96),
	   and would require about 2^192 samplings to be detectable, which
	   is way beyond the formal bound of 2^64 signatures with the
	   same key. */
	uint32_t s = (uint32_t)si;
	s |= (uint32_t)(63 - s) >> 26;
	s &= 63;

	/* Compute ccs*exp(-x). Since x = s*log(2) + r, we compute
	   ccs*exp(-r)/2^s. We know that 0 <= r < log(2), so we can
	   use expm_p63(), which yields a result scaled by 63 bits. We
	   scale it up 1 bit further, then right-shift by s bits.

	   We subtract 1 to make sure that the value fits on 64 bits
	   (i.e. if r = 0 then we may get 2^64 and we prefer 2^64-1
	   in that case, to avoid the overflow). The bias is negligible
	   since expm_p63() has precision only 51 bits or so. */
	uint64_t z = fpr_ursh((expm_p63(r, ccs) << 1) - 1, s);

	/* Sample a bit. We lazily compare the value z with a uniform 64-bit
	   integer, consuming only as many bytes as necessary. Since the PRNG
	   is cryptographically strong, we leak no information from the
	   conditional jumps below: w == bz has probability 1/256 whatever
	   bz is, hence the declassification. */
	for (int i = 56; i >= 0; i -= 8) {
		unsigned w = prng_next_u8(&ss->pc);
		unsigned bz = (unsigned)(z >> i) & 0xFF;
		if (ct_public_int(w != bz)) {
			return w < bz;
		}
	}
	return 0;
}

#ifdef FNDSA_PREC_CHECKPOINTS
/* Hook of the prec backend: records the errors of mu and 1/sigma. */
void fndsa_prec_checkpoint_mu_isigma(fpr mu, fpr isigma);
#endif

/* see sign_inner.h */
int32_t
sampler_next(sampler_state *ss, fpr mu, fpr isigma)
{
#ifdef FNDSA_PREC_CHECKPOINTS
	fndsa_prec_checkpoint_mu_isigma(mu, isigma);
#endif
	/* Split center mu into s + r, for an integer s, and 0 <= r < 1. */
	int64_t s = fpr_floor(mu);
	fpr r = fpr_sub(mu, fpr_of32((int32_t)s));

	/* dss = 1/(2*sigma^2) = 0.5*(isigma^2)  */
	fpr dss = fpr_half(fpr_sqr(isigma));

	/* css = sigma_min / sigma = sigma_min * isigma  */
	fpr ccs = fpr_mul(isigma, SIGMA_MIN[ss->logn].f);

	/* We sample on centre r. */
	for (;;) {
		/* Sample z for a Gaussian distribution (non-negative only),
		   then get a random bit b to turn the sampling into a
		   bimodal distribution (we use z+1 if b = 1, or -z
		   otherwise). */
		int32_t z0 = gaussian0(ss);
		int32_t b = prng_next_u8(&ss->pc) & 1;
		int32_t z = b + ((b << 1) - 1) * z0;

		/* Rejection sampling. We want a Gaussian centred on r,
		   but we sampled against a bimodal distribution (with
		   "centres" at 0 and 1). However, we know that z is
		   always in the range where our sampling distribution is
		   greater than the Gaussian distribution, so rejection works.

		   We got z from distribution:
		      G(z) = exp(-((z-b)^2)/(2*sigma0^2))
		   We target distribution:
		      S(z) = exp(-((z-r)^2)/(2*signa^2))
		   Rejection sampling works by keeping the value z with
		   probability S(z)/G(z), and starting again otherwise.
		   This requires S(z) <= G(z), which is the case here.
		   Thus, we simply need to keep our z with probability:
		      P = exp(-x)
		   where:
		      x = ((z-r)^2)/(2*sigma^2) - ((z-b)^2)/(2*sigma0^2)
		   Here, we scale up the Bernouilli distribution, which
		   makes rejection more probable, but also makes the
		   rejection rate sufficiently decorrelated from the Gaussian
		   centre and standard deviation, so that measurement of the
		   rejection rate does not leak enough usable information
		   to attackers (which is how the implementation can claim
		   to be "constant-time").  */
		fpr x = fpr_mul(fpr_sqr(fpr_sub(fpr_of32(z), r)), dss);
		x = fpr_sub(x, fpr_mul(fpr_of32(z0 * z0), INV_2SQRSIGMA0));
		/* The acceptance rate does not depend on mu or sigma, so
		   each decision is public. */
		if (ct_public_int(ber_exp(ss, x, ccs))) {
			return (int32_t)s + z;
		}
	}
}
#endif

#if !FNDSA_FFSAMP_ASM
TARGET_SSE2 TARGET_NEON static
#endif
void
ffsamp_fft_deepest(sampler_state *ss, fpr *tmp)
{
	fpr *t0 = tmp;
	fpr *t1 = tmp + 2;
#if FNDSA_SIGN_L10_64
	/* Same band layout as ffsamp_fft_inner, at logn = 1. */
	fpr *g00 = (fpr *)((uint8_t *)tmp + 2 * (sizeof(fpr) << 1)
		+ (sizeof(l10_t) << 1));
	fpr *g11 = (fpr *)((uint8_t *)g00 + sizeof(fpr));
	fpr *g01 = (fpr *)((uint8_t *)g11 + sizeof(fpr));
#else
	fpr *g01 = tmp + 4;
	fpr *g00 = tmp + 6;
	fpr *g11 = tmp + 7;
#endif
#if FNDSA_SSE2
	static const union {
		fpr f[2];
		__m128d x;
	} one_u = { { FPR_ONE, FPR_ONE } };
	__m128d cz = _mm_castsi128_pd(
		_mm_setr_epi32(0, 0, 0, -0x80000000));

	/* Decompose G into LDL. g00 and g11 are self-adjoint,
	   thus only one (real) coefficient each. */
	__m128d g00_re = _mm_load_sd((double *)g00);
	__m128d g01_cc = _mm_loadu_pd((double *)g01);
	__m128d g11_re = _mm_load_sd((double *)g11);
	__m128d inv_g00_re = _mm_div_sd(one_u.x, g00_re);
	__m128d inv_g00 = _mm_shuffle_pd(inv_g00_re, inv_g00_re, 0);
	__m128d mu = _mm_mul_pd(g01_cc, inv_g00);
	__m128d zo = _mm_mul_pd(mu, g01_cc);
	__m128d zo_re = _mm_add_sd(zo, _mm_shuffle_pd(zo, zo, 1));
	__m128d d00_re = g00_re;
	__m128d l01 = _mm_xor_pd(cz, mu);
	__m128d d11_re = _mm_sub_sd(g11_re, zo_re);

	/* No split on d00 and d11, since they are one-coeff each. */

	/* The half-size Gram matrices for the recursive LDL tree
	   exploration are now:
	     - left sub-tree:   d00_re, zero, d00_re
	     - right sub-tree:  d11_re, zero, d11_re
	   t1 split is trivial. */
	__m128d w = _mm_loadu_pd((double *)t1);
	__m128d w0 = w;
	__m128d w1 = _mm_shuffle_pd(w, w, 3);
	__m128d leaf = _mm_mul_sd(
		_mm_sqrt_sd(_mm_setzero_pd(), d11_re),
		_mm_load_sd((const double *)INV_SIGMA + ss->logn));
	__m128d y0 = _mm_cvtsi32_sd(_mm_setzero_pd(),
		sampler_next_sse2(ss, w0, leaf));
	__m128d y1 = _mm_cvtsi32_sd(_mm_setzero_pd(),
		sampler_next_sse2(ss, w1, leaf));

	/* Merge is trivial, since logn = 1. */

	/* At this point:
	     t0 and t1 are unmodified; t1 is also [w0, w1]
	     l10 is in [l10_re, l10_im]
	     z1 is [y0, y1]
	   Compute tb0 = t0 + (t1 - z1)*l10  (into [x0, x1]).
	   z1 is moved into t1. */
	__m128d y = _mm_shuffle_pd(y0, y1, 0);
	__m128d a = _mm_sub_pd(w, y);
	__m128d b1 = _mm_mul_pd(a, _mm_xor_pd(cz, l01));
	__m128d b2 = _mm_mul_pd(a, _mm_shuffle_pd(l01, l01, 1));
	__m128d b = _mm_add_pd(
		_mm_shuffle_pd(b1, b2, 2),
		_mm_shuffle_pd(b1, b2, 1));
	__m128d x = _mm_add_pd(b, _mm_loadu_pd((double *)t0));
	_mm_storeu_pd((double *)t1, y);

	/* Second recursive invocation, on the split tb0, using
	   the left sub-tree. tb0 is [x0, x1], and the split is
	   trivial since logn = 1. */
	__m128d x0 = x;
	__m128d x1 = _mm_shuffle_pd(x, x, 3);
	leaf = _mm_mul_sd(
		_mm_sqrt_sd(_mm_setzero_pd(), d00_re),
		_mm_load_sd((const double *)INV_SIGMA + ss->logn));
	x0 = _mm_cvtsi32_sd(_mm_setzero_pd(),
		sampler_next_sse2(ss, x0, leaf));
	x1 = _mm_cvtsi32_sd(_mm_setzero_pd(),
		sampler_next_sse2(ss, x1, leaf));
	_mm_store_sd((double *)t0, x0);
	_mm_store_sd((double *)t0 + 1, x1);
#elif FNDSA_NEON
	static const fpr_u one_u = { FPR_ONE };
	static const union { fpr f[2]; float64x2_t x; }
		cz = { { FPR_ZERO, FPR_NZERO } };

	/* Decompose G into LDL. g00 and g11 are self-adjoint,
	   thus only one (real) coefficient each. */
	float64x1_t g00_re = vld1_f64((float64_t *)g00);
	float64x2_t g01_cc = vld1q_f64((float64_t *)g01);
	float64x1_t g11_re = vld1_f64((float64_t *)g11);
	float64x1_t inv_g00_re = vdiv_f64(one_u.v, g00_re);
	float64x2_t inv_g00 = vdupq_lane_f64(inv_g00_re, 0);
	float64x2_t mu = vmulq_f64(g01_cc, inv_g00);
	float64x2_t zo = vmulq_f64(mu, g01_cc);
	float64x1_t zo_re = vget_low_f64(vpaddq_f64(zo, zo));
	float64x1_t d00_re = g00_re;
	float64x2_t l01 = vreinterpretq_f64_u64(
		veorq_u64(cz.x, vreinterpretq_u64_f64(mu)));
	float64x1_t d11_re = vsub_f64(g11_re, zo_re);

	/* No split on d00 and d11, since they are one-coeff each. */

	/* The half-size Gram matrices for the recursive LDL tree
	   exploration are now:
	     - left sub-tree:   d00_re, zero, d00_re
	     - right sub-tree:  d11_re, zero, d11_re
	   t1 split is trivial. */
	float64x2_t w = vld1q_f64((const float64_t *)t1);
	float64x1_t w0 = vget_low_f64(w);
	float64x1_t w1 = vget_high_f64(w);
	float64x1_t leaf = vmul_f64(
		vsqrt_f64(d11_re),
		INV_SIGMA[ss->logn].v);
	float64x1_t y0 = vcvt_f64_s64(vcreate_s64(
		sampler_next_neon(ss, w0, leaf)));
	float64x1_t y1 = vcvt_f64_s64(vcreate_s64(
		sampler_next_neon(ss, w1, leaf)));

	/* Merge is trivial, since logn = 1. */

	/* At this point:
	     t0 and t1 are unmodified; t1 is also [w0, w1]
	     l10 is in [l10_re, l10_im]
	     z1 is [y0, y1]
	   Compute tb0 = t0 + (t1 - z1)*l10  (into [x0, x1]).
	   z1 is moved into t1. */
	float64x2_t y = vcombine_f64(y0, y1);
	float64x2_t a = vsubq_f64(w, y);
	float64x2_t b1 = vmulq_f64(a,
		vreinterpretq_f64_u64(veorq_u64(cz.x,
			vreinterpretq_u64_f64(l01))));
	float64x2_t b2 = vmulq_f64(a, vextq_f64(l01, l01, 1));
	float64x2_t b = vpaddq_f64(b1, b2);
	float64x2_t x = vaddq_f64(b, vld1q_f64((const float64_t *)t0));
	vst1q_f64((float64_t *)t1, y);

	/* Second recursive invocation, on the split tb0, using
	   the left sub-tree. tb0 is [x0, x1], and the split is
	   trivial since logn = 1. */
	float64x1_t x0 = vget_low_f64(x);
	float64x1_t x1 = vget_high_f64(x);
	leaf = vmul_f64(
		vsqrt_f64(d00_re),
		INV_SIGMA[ss->logn].v);
	x0 = vcvt_f64_s64(vcreate_s64(
		sampler_next_neon(ss, x0, leaf)));
	x1 = vcvt_f64_s64(vcreate_s64(
		sampler_next_neon(ss, x1, leaf)));
	vst1_f64((float64_t *)t0, x0);
	vst1_f64((float64_t *)t0 + 1, x1);
#elif FNDSA_RV64D
	/* Decompose G into LDL. g00 and g11 are self-adjoint,
	   thus only one (real) coefficient each. */
	f64 g00_re = ((const f64 *)g00)[0];
	f64 g01_re = ((const f64 *)g01)[0];
	f64 g01_im = ((const f64 *)g01)[1];
	f64 g11_re = ((const f64 *)g11)[0];
	f64 inv_g00_re = f64_inv(g00_re);
	f64 mu_re = f64_mul(g01_re, inv_g00_re);
	f64 mu_im = f64_mul(g01_im, inv_g00_re);
	f64 zo_re = f64_add(
		f64_mul(mu_re, g01_re),
		f64_mul(mu_im, g01_im));
	f64 d00_re = g00_re;
	f64 l01_re = mu_re;
	f64 l01_im = f64_neg(mu_im);
	f64 d11_re = f64_sub(g11_re, zo_re);

	/* No split on d00 and d11, since they are one-coeff each. */

	/* The half-size Gram matrices for the recursive LDL tree
	   exploration are now:
	     - left sub-tree:   d00_re, zero, d00_re
	     - right sub-tree:  d11_re, zero, d11_re
	   t1 split is trivial. */
	f64 w0 = ((const f64 *)t1)[0];
	f64 w1 = ((const f64 *)t1)[1];
	f64 leaf = f64_mul(f64_sqrt(d11_re), INV_SIGMA[ss->logn].v);
	f64 y0 = f64_of(sampler_next_rv64d(ss, w0, leaf));
	f64 y1 = f64_of(sampler_next_rv64d(ss, w1, leaf));

	/* Merge is trivial, since logn = 1. */

	/* At this point:
	     t0 and t1 are unmodified; t1 is also [w0, w1]
	     l10 is in [l10_re, l10_im]
	     z1 is [y0, y1]
	   Compute tb0 = t0 + (t1 - z1)*l10  (into [x0, x1]).
	   z1 is moved into t1. */
	f64 a_re = f64_sub(w0, y0);
	f64 a_im = f64_sub(w1, y1);
	f64 b_re = f64_sub(
		f64_mul(a_re, l01_re),
		f64_mul(a_im, l01_im));
	f64 b_im = f64_add(
		f64_mul(a_im, l01_re),
		f64_mul(a_re, l01_im));
	f64 x0 = f64_add(((const f64 *)t0)[0], b_re);
	f64 x1 = f64_add(((const f64 *)t0)[1], b_im);
	((f64 *)t1)[0] = y0;
	((f64 *)t1)[1] = y1;

	/* Second recursive invocation, on the split tb0, using
	   the left sub-tree. tb0 is [x0, x1], and the split is
	   trivial since logn = 1. */
	leaf = f64_mul(f64_sqrt(d00_re), INV_SIGMA[ss->logn].v);
	((f64 *)t0)[0] = f64_of(sampler_next_rv64d(ss, x0, leaf));
	((f64 *)t0)[1] = f64_of(sampler_next_rv64d(ss, x1, leaf));
#else
	/* Decompose G into LDL. g00 and g11 are self-adjoint,
	   thus only one (real) coefficient each. */
	fpr g00_re = g00[0];
	fpr g01_re = g01[0], g01_im = g01[1];
	fpr g11_re = g11[0];
	fpr_rcp inv_g00;
	fpr mu_re, mu_im;
	FPR_RECIP_X(inv_g00, g00_re);
	FPR_MUL_RECIP(mu_re, g01_re, inv_g00);
	FPR_MUL_RECIP(mu_im, g01_im, inv_g00);
	fpr zo_re = fpr_add(
		fpr_mul(mu_re, g01_re),
		fpr_mul(mu_im, g01_im));
	fpr d00_re = g00_re;
	fpr l01_re = mu_re;
	fpr l01_im = fpr_neg(mu_im);
	fpr d11_re = fpr_sub(g11_re, zo_re);

	/* No split on d00 and d11, since they are one-coeff each. */

	/* The half-size Gram matrices for the recursive LDL tree
	   exploration are now:
	     - left sub-tree:   d00_re, zero, d00_re
	     - right sub-tree:  d11_re, zero, d11_re
	   t1 split is trivial. */
	fpr w0 = t1[0];
	fpr w1 = t1[1];
	fpr leaf;
	FPR_SQRT_MUL_INV_SIGMA(leaf, d11_re, ss->logn);
	fpr y0 = fpr_of32(sampler_next(ss, w0, leaf));
	fpr y1 = fpr_of32(sampler_next(ss, w1, leaf));

	/* Merge is trivial, since logn = 1. */

	/* At this point:
	     t0 and t1 are unmodified; t1 is also [w0, w1]
	     l10 is in [l10_re, l10_im]
	     z1 is [y0, y1]
	   Compute tb0 = t0 + (t1 - z1)*l10  (into [x0, x1]).
	   z1 is moved into t1. */
	fpr a_re = fpr_sub(w0, y0);
	fpr a_im = fpr_sub(w1, y1);
	fpr b_re, b_im;
	FPC_MUL_NOALIAS(b_re, b_im, a_re, a_im, l01_re, l01_im);
	fpr x0 = fpr_add(t0[0], b_re);
	fpr x1 = fpr_add(t0[1], b_im);
	t1[0] = y0;
	t1[1] = y1;

	/* Second recursive invocation, on the split tb0, using
	   the left sub-tree. tb0 is [x0, x1], and the split is
	   trivial since logn = 1. */
	FPR_SQRT_MUL_INV_SIGMA(leaf, d00_re, ss->logn);
	t0[0] = fpr_of32(sampler_next(ss, x0, leaf));
	t0[1] = fpr_of32(sampler_next(ss, x1, leaf));
#endif
	return;
}

#if FNDSA_SYMP_TREE

/* ffsamp_fft_deepest() on a mirrored subtree (Symplectree). */
void
ffsamp_fft_deepest_mirrored(sampler_state *ss, fpr *tmp)
{
	fpr *t0 = tmp;
	fpr *t1 = tmp + 2;
#if FNDSA_SIGN_L10_64
	/* Same band layout as ffsamp_fft_inner, at logn = 1. */
	fpr *g00 = (fpr *)((uint8_t *)tmp + 2 * (sizeof(fpr) << 1)
		+ (sizeof(l10_t) << 1));
	fpr *g11 = (fpr *)((uint8_t *)g00 + sizeof(fpr));
	fpr *g01 = (fpr *)((uint8_t *)g11 + sizeof(fpr));
#else
	fpr *g01 = tmp + 4;
	fpr *g00 = tmp + 6;
	fpr *g11 = tmp + 7;
#endif
	/* Decompose G into LDL. g00 and g11 are self-adjoint,
	   thus only one (real) coefficient each. */
	fpr g00_re = g00[0];
	fpr g01_re = g01[0], g01_im = g01[1];
	fpr g11_re = g11[0];
	fpr_rcp inv_g00;
	fpr mu_re, mu_im;
	FPR_RECIP_X(inv_g00, g00_re);
	FPR_MUL_RECIP(mu_re, g01_re, inv_g00);
	FPR_MUL_RECIP(mu_im, g01_im, inv_g00);
	fpr zo_re = fpr_add(
		fpr_mul(mu_re, g01_re),
		fpr_mul(mu_im, g01_im));
	fpr d00_re = g00_re;
	fpr d11_re = fpr_sub(g11_re, zo_re);


	/* Mirror: leaves are q^2/x, with d00 and d11 swapped. */
    fpr_rcp inv_d11;
    fpr q2c = Q2;
    fpr tree2_d00_re, tree2_d11_re;
    FPR_RECIP_X(inv_d11, d11_re);
    FPR_MUL_RECIP(tree2_d00_re, q2c, inv_d11);
    FPR_MUL_RECIP(tree2_d11_re, q2c, inv_g00);
	/* The mirror also negates l10. */
	fpr l01_re = fpr_neg(mu_re);
	fpr l01_im = mu_im;

	/* No split on d00 and d11, since they are one-coeff each. */

	/* The half-size Gram matrices for the recursive LDL tree
	   exploration are now:
	     - left sub-tree:   d00_re, zero, d00_re
	     - right sub-tree:  d11_re, zero, d11_re
	   t1 split is trivial. */
	fpr w0 = t1[0];
	fpr w1 = t1[1];
	fpr leaf;
	FPR_SQRT_MUL_INV_SIGMA(leaf, tree2_d11_re, ss->logn);
	fpr y0 = fpr_of32(sampler_next(ss, w0, leaf));
	fpr y1 = fpr_of32(sampler_next(ss, w1, leaf));

	/* Merge is trivial, since logn = 1. */

	/* At this point:
	     t0 and t1 are unmodified; t1 is also [w0, w1]
	     l10 is in [l10_re, l10_im]
	     z1 is [y0, y1]
	   Compute tb0 = t0 + (t1 - z1)*l10  (into [x0, x1]).
	   z1 is moved into t1. */
	fpr a_re = fpr_sub(w0, y0);
	fpr a_im = fpr_sub(w1, y1);
	fpr b_re, b_im;
	FPC_MUL_NOALIAS(b_re, b_im, a_re, a_im, l01_re, l01_im);
	fpr x0 = fpr_add(t0[0], b_re);
	fpr x1 = fpr_add(t0[1], b_im);
	t1[0] = y0;
	t1[1] = y1;

	/* Second recursive invocation, on the split tb0, using
	   the left sub-tree. tb0 is [x0, x1], and the split is
	   trivial since logn = 1. */
	FPR_SQRT_MUL_INV_SIGMA(leaf, tree2_d00_re, ss->logn);
	t0[0] = fpr_of32(sampler_next(ss, x0, leaf));
	t0[1] = fpr_of32(sampler_next(ss, x1, leaf));

	return;
}

/* ffsamp_fft_inner() on a mirrored subtree (Symplectree). */
static void
ffsamp_fft_mirror(sampler_state *ss, unsigned logn, fpr *tmp)
{
	/* When logn = 1, arrays have length 2; we unroll the last steps
	   in a dedicated function. */
	if (logn == 1) {
		ffsamp_fft_deepest_mirrored(ss, tmp);
		return;
	}

	/* General case: logn >= 2 */

	/* Same layout as ffsamp_fft_inner, with the diagonals swapped:
	   d00 (for the right subtree) is transient, d11 (for the left
	   subtree) persistent. The mirror is never the root. */
	size_t sz_fpr = sizeof(fpr) << logn;
	size_t sz_l10 = sizeof(l10_t) << logn;
	size_t child_sz_fpr = sz_fpr >> 1;
	size_t child_sz_l10 = sizeof(l10_t) << (logn - 1);
	uint8_t *base = (uint8_t *)tmp;

	fpr *t0 = tmp;
	fpr *t1 = (fpr *)(base + sz_fpr);
	l10_t *l10 = (l10_t *)(base + 2 * sz_fpr);
	/* d11 (from g11) is persistent, d00 (= g00) transient. */
	fpr *d11 = (fpr *)(base + 2 * sz_fpr + sz_l10);
	uint8_t *child = base + 2 * sz_fpr + sz_l10 + sz_fpr / 2;
	fpr *d00 = (fpr *)child;
#if FNDSA_SIGN_L10_64
	fpr *g01 = (fpr *)(child + sz_fpr / 2);
#else
	fpr *g01 = (fpr *)l10;
#endif

	/* Where the child reads its inputs; its two diagonal slots
	   receive the same self-adjoint polynomial. */
	fpr *child_t0 = (fpr *)child;
	fpr *child_t1 = (fpr *)(child + child_sz_fpr);
	fpr *child_diag = (fpr *)(child + 2 * child_sz_fpr + child_sz_l10);
	uint8_t *grandchild = (uint8_t *)child_diag + child_sz_fpr / 2;
	fpr *child_diag_copy = (fpr *)grandchild;
#if FNDSA_SIGN_L10_64
	fpr *child_g01 = (fpr *)(grandchild + child_sz_fpr / 2);
#else
	fpr *child_g01 = (fpr *)(child + 2 * child_sz_fpr);
#endif

	/* Scratch, in the child region while it is still free. */
	fpr *corr = (fpr *)child;             /* (t1 - z1) * l10 */
	fpr *z1 = (fpr *)(child + sz_fpr);

	fpoly_LDL_fft(logn, d00, g01, d11, l10);
	/* The mirror negates l10. */
	fpoly_l10_neg(logn, l10);

	/* Split d00 (split yields right_00 and right_01, right_11 is a
	   copy of right_00). */
	fpoly_split_selfadj_fft(logn, child_diag, child_g01, d00);
	memcpy(child_diag_copy, child_diag, sizeof(fpr) << (logn - 2));

	/* Split t1 and make the first recursive call on the two
	   halves, using the right sub-tree, then merge the result
	   into z1. */
	fpoly_split_fft(logn, child_t0, child_t1, t1);
	ffsamp_fft_mirror(ss, logn - 1, (fpr *)child);
	fpoly_merge_fft(logn, z1, child_t0, child_t1);

	/* Compute tb0 = t0 + (t1 - z1)*l10 (into t0) and move z1 into t1. */
	memcpy(corr, t1, sizeof(fpr) << logn);
	fpoly_sub(logn, corr, z1);
	memcpy(t1, z1, sizeof(fpr) << logn);
	fpoly_mul_l10_fft(logn, corr, l10);
	fpoly_add(logn, t0, corr);

	/* Split d11 to obtain the left sub-tree. */
	fpoly_split_selfadj_fft(logn, child_diag, child_g01, d11);
	memcpy(child_diag_copy, child_diag, sizeof(fpr) << (logn - 2));

	/* Split tb0 and perform the second recursive call on the
	   split output; the final merge produces z0, which we write
	   into t0. */
	fpoly_split_fft(logn, child_t0, child_t1, t0);
	ffsamp_fft_mirror(ss, logn - 1, (fpr *)child);
	fpoly_merge_fft(logn, t0, child_t0, child_t1);
}


#endif


#if FNDSA_FFSAMP_ASM
#define ffsamp_fft_inner   fndsa_ffsamp_fft_inner
void ffsamp_fft_inner(sampler_state *ss, unsigned logn, fpr *tmp);
#else
TARGET_SSE2 TARGET_NEON
static void
ffsamp_fft_inner(sampler_state *ss, unsigned logn, fpr *tmp)
{
	/* When logn = 1, arrays have length 2; we unroll the last steps
	   in a dedicated function. */
	if (logn == 1) {
		ffsamp_fft_deepest(ss, tmp);
		return;
	}

	/* General case: logn >= 2 */

	/* Byte layout of this node: the bands kept across the recursive
	   calls, then the child region.

	      t0     | t1     | l10    | d00      | child region ...
	      sz_fpr   sz_fpr   sz_l10   sz_fpr/2

	   g11 (turned into d11) starts the child region. g01 follows it
	   with L10 in 64 bits below the root; otherwise it sits in the
	   l10 band and the LDL overwrites it. */
	int root = (logn == ss->logn);
	size_t sz_fpr = sizeof(fpr) << logn;       /* full fpr band */
	size_t sz_l10 = root ? sz_fpr : (sizeof(l10_t) << logn);
	size_t child_sz_fpr = sz_fpr >> 1;
	size_t child_sz_l10 = sizeof(l10_t) << (logn - 1);
	uint8_t *base = (uint8_t *)tmp;

	fpr *t0 = tmp;
	fpr *t1 = (fpr *)(base + sz_fpr);
	l10_t *l10 = (l10_t *)(base + 2 * sz_fpr);
	/* d00 (= g00) is persistent, d11 (from g11) transient. */
	fpr *d00 = (fpr *)(base + 2 * sz_fpr + sz_l10);
	uint8_t *child = base + 2 * sz_fpr + sz_l10 + sz_fpr / 2;
	fpr *d11 = (fpr *)child;
#if FNDSA_SIGN_L10_64
	fpr *g01 = root ? (fpr *)l10 : (fpr *)(child + sz_fpr / 2);
#else
	fpr *g01 = (fpr *)l10;
#endif

	/* Where the child reads its inputs (the child is never the
	   root); its two diagonal slots receive the same self-adjoint
	   polynomial. */
	fpr *child_t0 = (fpr *)child;
	fpr *child_t1 = (fpr *)(child + child_sz_fpr);
	fpr *child_diag = (fpr *)(child + 2 * child_sz_fpr + child_sz_l10);
	uint8_t *grandchild = (uint8_t *)child_diag + child_sz_fpr / 2;
	fpr *child_diag_copy = (fpr *)grandchild;
#if FNDSA_SIGN_L10_64
	fpr *child_g01 = (fpr *)(grandchild + child_sz_fpr / 2);
#else
	fpr *child_g01 = (fpr *)(child + 2 * child_sz_fpr);
#endif

	/* Scratch, in the child region while it is still free. */
	fpr *corr = (fpr *)child;             /* (t1 - z1) * l10 */
	fpr *z1 = (fpr *)(child + sz_fpr);

	/* Decompose G into LDL; the decomposed matrix replaces G. */
#if FNDSA_SYMP_ROOT
	/* Symplecroot: d11 = q^2/g00 holds at the root only. */
	fpr *d_right = d11;
	if (root) {
		fpoly_LDL_root_fft(logn, d00, g01, d11);
	} else {
		fpoly_LDL_fft(logn, d00, g01, d11, l10);
	}
#elif FNDSA_SYMP_TREE
	/* The mirrored right subtree is built from d00. */
	fpr *d_right;
	if (root) {
		d_right = d00;
		fpoly_LDL_root_fft(logn, d00, g01);
	} else {
		d_right = d11;
		fpoly_LDL_fft(logn, d00, g01, d11, l10);
	}
#else
	fpoly_LDL_fft(logn, d00, g01, d11, l10);
	fpr *d_right = d11;
#endif

	/* Split the right sub-tree's diagonal (split yields right_00
	   and right_01, right_11 is a copy of right_00). */
	fpoly_split_selfadj_fft(logn, child_diag, child_g01, d_right);
	memcpy(child_diag_copy, child_diag, sizeof(fpr) << (logn - 2));

	/* Split t1 and make the first recursive call on the two
	   halves, using the right sub-tree, then merge the result
	   into z1. */
	fpoly_split_fft(logn, child_t0, child_t1, t1);
#if FNDSA_SYMP_TREE
	if (root) {
		ffsamp_fft_mirror(ss, logn - 1, (fpr *)child);
	} else {
		ffsamp_fft_inner(ss, logn - 1, (fpr *)child);
	}
#else
	ffsamp_fft_inner(ss, logn - 1, (fpr *)child);
#endif
	fpoly_merge_fft(logn, z1, child_t0, child_t1);

	/* Compute tb0 = t0 + (t1 - z1)*l10 (into t0) and move z1 into t1. */
	memcpy(corr, t1, sizeof(fpr) << logn);
	fpoly_sub(logn, corr, z1);
	memcpy(t1, z1, sizeof(fpr) << logn);
	/* At the root l10 is a full fpr polynomial (it has an integer
	   part); below it, an l10_t. */
	if (root) {
		fpoly_mul_fft(logn, corr, (fpr *)l10);
	} else {
		fpoly_mul_l10_fft(logn, corr, l10);
	}
	fpoly_add(logn, t0, corr);

	/* Split d00 to obtain the left-subtree. */
	fpoly_split_selfadj_fft(logn, child_diag, child_g01, d00);
	memcpy(child_diag_copy, child_diag, sizeof(fpr) << (logn - 2));

	/* Split tb0 and perform the second recursive call on the
	   split output; the final merge produces z0, which we write
	   into t0. */
	fpoly_split_fft(logn, child_t0, child_t1, t0);
	ffsamp_fft_inner(ss, logn - 1, (fpr *)child);
	fpoly_merge_fft(logn, t0, child_t0, child_t1);
}
#endif


/* see sign_inner.h */
void
ffsamp_fft(sampler_state *ss, fpr *tmp)
{
	ffsamp_fft_inner(ss, ss->logn, tmp);
}
