/*
 * Floating-point polynomials.
 */

#include "sign_inner.h"
#include <stddef.h>

/* The GM[] twiddle table is defined per representation, in
   fpr/<rep>/fpr_const_*.c. */

/* see sign_inner.h */
TARGET_SSE2 TARGET_NEON
void
fpoly_FFT(unsigned logn, fpr *f)
{
#if FNDSA_SSE2
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	size_t t = hn;
	double *ff = (double *)f;
	/* We separate the last iteration of the loop. With that change,
	   t >= 4 and ht >= 2 in all iteration of the loop. */
	for (unsigned lm = 1; lm < (logn - 1); lm ++) {
		size_t m = (size_t)1 << lm;
		size_t hm = m >> 1;
		size_t ht = t >> 1;
		size_t j0 = 0;
		for (size_t i = 0; i < hm; i ++) {
			__m128d s = _mm_loadu_pd(
				(const double *)GM + ((m + i) << 1));
			__m128d s_re = _mm_shuffle_pd(s, s, 0);
			__m128d s_im = _mm_shuffle_pd(s, s, 3);
			for (size_t j = 0; j < ht; j += 2) {
				size_t j1 = j0 + j;
				size_t j2 = j1 + ht;
				__m128d x_re = _mm_loadu_pd(ff + j1);
				__m128d x_im = _mm_loadu_pd(ff + j1 + hn);
				__m128d y_re = _mm_loadu_pd(ff + j2);
				__m128d y_im = _mm_loadu_pd(ff + j2 + hn);
				__m128d z_re = _mm_sub_pd(
					_mm_mul_pd(y_re, s_re),
					_mm_mul_pd(y_im, s_im));
				__m128d z_im = _mm_add_pd(
					_mm_mul_pd(y_re, s_im),
					_mm_mul_pd(y_im, s_re));
				_mm_storeu_pd(ff + j1,
					_mm_add_pd(x_re, z_re));
				_mm_storeu_pd(ff + j1 + hn,
					_mm_add_pd(x_im, z_im));
				_mm_storeu_pd(ff + j2,
					_mm_sub_pd(x_re, z_re));
				_mm_storeu_pd(ff + j2 + hn,
					_mm_sub_pd(x_im, z_im));
			}
			j0 += t;
		}
		t = ht;
	}

	/* Last iteration: m = n/2, hm = n/4, t = 2, ht = 1 */
	if (logn >= 2) {
		__m128d cz = _mm_castsi128_pd(
			_mm_setr_epi32(0, 0, 0, -0x80000000));
		for (size_t i = 0; i < hn; i += 2) {
			/* s <- re(s):im(s) */
			__m128d s = _mm_loadu_pd((const double *)GM + n + i);
			/* xy_re <- re(x):re(y) */
			__m128d xy_re = _mm_loadu_pd(ff + i);
			/* xy_im <- im(x):im(y) */
			__m128d xy_im = _mm_loadu_pd(ff + i + hn);
			/* y1 <- re(y):im(y) */
			__m128d y1 = _mm_shuffle_pd(xy_re, xy_im, 3);
			/* y2 <- im(y):re(y) */
			__m128d y2 = _mm_shuffle_pd(xy_im, xy_re, 3);
			/* z_re <- re(y)*re(s):im(y)*im(s) */
			__m128d z_re = _mm_mul_pd(y1, s);
			/* z_im <- im(y)*re(s):re(y)*im(s) */
			__m128d z_im = _mm_mul_pd(y2, s);
			/* With u = y*s (complex):
			   u_re <- re(u):-re(u)
			   u_im <- im(u):-im(u)  */
			__m128d u_re = _mm_sub_pd(
				z_re,
				_mm_shuffle_pd(z_re, z_re, 1));
			__m128d u_im = _mm_xor_pd(cz,
				_mm_add_pd(
					z_im,
					_mm_shuffle_pd(z_im, z_im, 1)));
			/* u_re <- re(x+z):re(x-z)
			   u_im <- im(x+z):im(x-z) */
			u_re = _mm_add_pd(u_re,
				_mm_shuffle_pd(xy_re, xy_re, 0));
			u_im = _mm_add_pd(u_im,
				_mm_shuffle_pd(xy_im, xy_im, 0));
			_mm_storeu_pd(ff + i, u_re);
			_mm_storeu_pd(ff + i + hn, u_im);
		}
	}
#elif FNDSA_NEON
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	size_t t = hn;
	float64_t *ff = (float64_t *)f;
	/* We separate the last iteration of the loop. With that change,
	   t >= 4 and ht >= 2 in all iteration of the loop. */
	for (unsigned lm = 1; lm < (logn - 1); lm ++) {
		size_t m = (size_t)1 << lm;
		size_t hm = m >> 1;
		size_t ht = t >> 1;
		size_t j0 = 0;
		for (size_t i = 0; i < hm; i ++) {
			float64x2_t s = vld1q_f64(
				(const float64_t *)GM + ((m + i) << 1));
			float64x2_t s_re = vzip1q_f64(s, s);
			float64x2_t s_im = vzip2q_f64(s, s);
			for (size_t j = 0; j < ht; j += 2) {
				size_t j1 = j0 + j;
				size_t j2 = j1 + ht;
				float64x2_t x_re = vld1q_f64(ff + j1);
				float64x2_t x_im = vld1q_f64(ff + j1 + hn);
				float64x2_t y_re = vld1q_f64(ff + j2);
				float64x2_t y_im = vld1q_f64(ff + j2 + hn);
				float64x2_t z_re = vsubq_f64(
					vmulq_f64(y_re, s_re),
					vmulq_f64(y_im, s_im));
				float64x2_t z_im = vaddq_f64(
					vmulq_f64(y_re, s_im),
					vmulq_f64(y_im, s_re));
				vst1q_f64(ff + j1, vaddq_f64(x_re, z_re));
				vst1q_f64(ff + j1 + hn, vaddq_f64(x_im, z_im));
				vst1q_f64(ff + j2, vsubq_f64(x_re, z_re));
				vst1q_f64(ff + j2 + hn, vsubq_f64(x_im, z_im));
			}
			j0 += t;
		}
		t = ht;
	}

	/* Last iteration: m = n/2, hm = n/4, t = 2, ht = 1 */
	if (logn >= 2) {
		static const union { fpr f[2]; float64x2_t x; }
			cz = { { FPR_ZERO, FPR_NZERO } };
		for (size_t i = 0; i < hn; i += 2) {
			/* s <- re(s):im(s) */
			float64x2_t s = vld1q_f64(
				(const float64_t *)GM + n + i);
			/* xy_re <- re(x):re(y) */
			float64x2_t xy_re = vld1q_f64(ff + i);
			/* xy_im <- im(x):im(y) */
			float64x2_t xy_im = vld1q_f64(ff + i + hn);
			/* y1 <- re(y):im(y) */
			float64x2_t y1 = vzip2q_f64(xy_re, xy_im);
			/* y2 <- im(y):re(y) */
			float64x2_t y2 = vzip2q_f64(xy_im, xy_re);
			/* z_re <- re(y)*re(s):im(y)*im(s) */
			float64x2_t z_re = vmulq_f64(y1, s);
			/* z_im <- im(y)*re(s):re(y)*im(s) */
			float64x2_t z_im = vmulq_f64(y2, s);
			/* With u = y*s (complex):
			   u_re <- re(u):-re(u)
			   u_im <- im(u):-im(u)  */
			float64x2_t u_re = vsubq_f64(
				z_re,
				vextq_f64(z_re, z_re, 1));
			float64x2_t u_im = vreinterpretq_f64_u64(
				veorq_u64(cz.x, vreinterpretq_u64_f64(
					vaddq_f64(
						z_im,
						vextq_f64(z_im, z_im, 1)))));
			/* u_re <- re(x+z):re(x-z)
			   u_im <- im(x+z):im(x-z) */
			u_re = vaddq_f64(u_re, vdupq_laneq_f64(xy_re, 0));
			u_im = vaddq_f64(u_im, vdupq_laneq_f64(xy_im, 0));
			vst1q_f64(ff + i, u_re);
			vst1q_f64(ff + i + hn, u_im);
		}
	}
#elif FNDSA_RV64D
	size_t hn = (size_t)1 << (logn - 1);
	size_t t = hn;
	f64 *ff = (f64 *)f;
	for (unsigned lm = 1; lm < logn; lm ++) {
		size_t m = (size_t)1 << lm;
		size_t hm = m >> 1;
		size_t ht = t >> 1;
		size_t j0 = 0;
		for (size_t i = 0; i < hm; i ++) {
			f64 s_re = ((const f64 *)GM)[((m + i) << 1) + 0];
			f64 s_im = ((const f64 *)GM)[((m + i) << 1) + 1];
			for (size_t j = 0; j < ht; j ++) {
				size_t j1 = j0 + j;
				size_t j2 = j1 + ht;
				f64 x_re = ff[j1];
				f64 x_im = ff[j1 + hn];
				f64 y_re = ff[j2];
				f64 y_im = ff[j2 + hn];
				f64 z_re = f64_sub(
					f64_mul(y_re, s_re),
					f64_mul(y_im, s_im));
				f64 z_im = f64_add(
					f64_mul(y_im, s_re),
					f64_mul(y_re, s_im));
				ff[j1] = f64_add(x_re, z_re);
				ff[j1 + hn] = f64_add(x_im, z_im);
				ff[j2] = f64_sub(x_re, z_re);
				ff[j2 + hn] = f64_sub(x_im, z_im);
			}
			j0 += t;
		}
		t = ht;
	}
#else
	size_t hn = (size_t)1 << (logn - 1);
	size_t t = hn;
	for (unsigned lm = 1; lm < logn; lm ++) {
		size_t m = (size_t)1 << lm;
		size_t hm = m >> 1;
		size_t ht = t >> 1;
		size_t j0 = 0;
		for (size_t i = 0; i < hm; i ++) {
			FPOLY_GM_TW(s, m + i);
			for (size_t j = 0; j < ht; j ++) {
				size_t j1 = j0 + j;
				size_t j2 = j1 + ht;
				fpr x_re = f[j1];
				fpr x_im = f[j1 + hn];
				fpr y_re = f[j2];
				fpr y_im = f[j2 + hn];
				FPOLY_GM_BFLY_FWD(f[j1], f[j1 + hn],
					f[j2], f[j2 + hn],
					x_re, x_im, y_re, y_im, s);
			}
			j0 += t;
		}
		t = ht;
	}
#endif
}

#if FNDSA_FFT_NEW
/* SIMD and RV64D builds do not link the scalar fpr routines, so they
   compute in native double here. */
static inline fpr
fpoly_int_add(fpr x, fpr y)
{
#if FNDSA_SSE2 || FNDSA_NEON || FNDSA_RV64D
	union { fpr f; double d; } ux, uy, uz;
	ux.f = x;
	uy.f = y;
	uz.d = ux.d + uy.d;
	return uz.f;
#else
	return fpr_add(x, y);
#endif
}

static inline fpr
fpoly_int_sub(fpr x, fpr y)
{
#if FNDSA_SSE2 || FNDSA_NEON || FNDSA_RV64D
	union { fpr f; double d; } ux, uy, uz;
	ux.f = x;
	uy.f = y;
	uz.d = ux.d - uy.d;
	return uz.f;
#else
	return fpr_sub(x, y);
#endif
}

static inline fpr
fpoly_int_mul(fpr x, fpr y)
{
#if FNDSA_SSE2 || FNDSA_NEON || FNDSA_RV64D
	union { fpr f; double d; } ux, uy, uz;
	ux.f = x;
	uy.f = y;
	uz.d = ux.d * uy.d;
	return uz.f;
#else
	return fpr_mul(x, y);
#endif
}

/* a = x + y, b = x - y. */
#if FNDSA_SSE2 || FNDSA_NEON || FNDSA_RV64D
#define FPOLY_INT_ADD_SUB(a, b, x, y)   do { \
		(a) = fpoly_int_add((x), (y)); \
		(b) = fpoly_int_sub((x), (y)); \
	} while (0)
#else
#define FPOLY_INT_ADD_SUB(a, b, x, y)   FPR_ADD_SUB((a), (b), (x), (y))
#endif

/* Split an RNJ coefficient into sign, integer part and 64-bit fraction
   of its absolute value. The sign is public. */
#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 \
	&& FNDSA_TWIDDLE_U64 && FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
static inline int
fpoly_rnj_parts(fpr r, uint64_t *frac, uint32_t *integer)
{
	uint64_t u = (uint64_t)r.w[0] | ((uint64_t)r.w[1] << 32);
	uint32_t h = r.w[2];
	int neg = (int)(h >> 31);
	if (neg) {
		u = (uint64_t)-u;
		h = ~h + (uint32_t)(u == 0);
	}
	*frac = u;
	*integer = h;
	return neg;
}
#endif

/* Compact FFT of a self-adjoint polynomial from int32_t input. With
   scale_half, the result is also halved (in the final merge). */
static void
fpoly_FFT_int32_selfadj_inner(unsigned logn, fpr *d, const int32_t *a,
	int scale_half)
{
	if (logn == 1) {
		d[0] = fpr_of32(a[0]);
		if (scale_half) {
			d[0] = fpr_half(d[0]);
		}
		return;
	}

	/* From c-fn-dsa-alt. The splits stay in native integers. */
	size_t hn = (size_t)1 << (logn - 1);
	size_t qn = hn >> 1;

	/* aa is the tail of the first half of d; the source may already
	   be there. */
	int32_t *aa = (int32_t *)(void *)(d + hn) - hn;
	if (a != aa) {
		memmove(aa, a, hn * sizeof *aa);
	}

	/* All splits down to self-adjoint polynomials modulo X^4+1. */
	for (unsigned lk = 0; lk < (logn - 2); lk ++) {
		size_t k = (size_t)1 << lk;
		size_t dk = k << 1;
		for (size_t i = hn - dk; i > 0; i -= dk) {
			for (size_t j = 0; j < k; j ++) {
				aa[i + j + k] += aa[i + j - k];
			}
		}
		for (size_t j = 0; j < k; j ++) {
			aa[j + k] <<= 1;
		}
	}

	/* aa aliases d, so the conversion is a pass of its own. */
	for (size_t j = 0; j < hn; j ++) {
		d[j] = fpr_of32(aa[j]);
	}

	/* The omitted final split and the first merge combine into one
	   multiplication by sqrt(2). */
	#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 \
	&& FNDSA_TWIDDLE_U64 && FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
	/* RNJ[0] = sqrt(2), so the batched RNJ kernel applies. */
	fpr sqrt2 = RNJ[0];
	uint64_t sqrt2_frac;
	uint32_t sqrt2_integer;
	int sqrt2_neg = fpoly_rnj_parts(sqrt2, &sqrt2_frac, &sqrt2_integer);
	FNDSA_FPR_BFLY_RNJ_BATCH(d, d + qn, qn,
		sqrt2_frac, sqrt2_integer, sqrt2_neg);
	#else
	for (size_t j = 0; j < qn; j ++) {
		/* RNJ[0] is sqrt(2) rounded once, not 2*RNJ[1]. */
		fpr x = d[j];
		fpr z = fpoly_int_mul(d[j + qn], RNJ[0]);
		FPOLY_INT_ADD_SUB(d[j], d[j + qn], x, z);
	}
	#endif
	if (scale_half && logn == 2) {
		for (size_t j = 0; j < hn; j ++) {
			d[j] = fpr_half(d[j]);
		}
	}

	/* Real-only recombination. */
	for (unsigned lm = 1; lm < (logn - 1); lm ++) {
		unsigned lk = logn - 2 - lm;
		size_t k = (size_t)1 << lk;
		size_t m = (size_t)1 << lm;
		int final_merge = scale_half && (lm == logn - 2);
		for (size_t j = 0; j < m; j ++) {
			size_t mj = j << (lk + 1);
#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 \
	&& FNDSA_TWIDDLE_U64 && FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
			fpr r = RNJ[j + m];
			uint64_t rfrac;
			uint32_t rinteger;
			int rneg = fpoly_rnj_parts(r, &rfrac, &rinteger);
			if (rinteger == 0) {
				/* The whole group shares r; one assembly call. */
				if (final_merge) {
					FNDSA_FPR_BFLY_GM_BATCH_HALF(d + mj, d + mj + k,
						k, rfrac, rneg);
				} else {
					FNDSA_FPR_BFLY_GM_BATCH(d + mj, d + mj + k,
						k, rfrac, rneg);
				}
			} else {
				/* Same, with r = rinteger + rfrac/2^64. */
				if (final_merge) {
					FNDSA_FPR_BFLY_RNJ_BATCH_HALF(d + mj, d + mj + k,
						k, rfrac, rinteger, rneg);
				} else {
					FNDSA_FPR_BFLY_RNJ_BATCH(d + mj, d + mj + k,
						k, rfrac, rinteger, rneg);
				}
			}
			continue;
#endif
			fpr rr = RNJ[j + m];
			for (size_t i = 0; i < k; i ++) {
				fpr x = d[mj + i];
				fpr z = fpoly_int_mul(d[mj + i + k], rr);
				FPOLY_INT_ADD_SUB(d[mj + i], d[mj + i + k], x, z);
				if (final_merge) {
					d[mj + i] = fpr_half(d[mj + i]);
					d[mj + i + k] = fpr_half(d[mj + i + k]);
				}
			}
		}
	}
}

/* see sign_inner.h */
void
fpoly_FFT_int32_selfadj(unsigned logn, fpr *d, const int32_t *a)
{
	fpoly_FFT_int32_selfadj_inner(logn, d, a, 0);
}

/* see sign_inner.h */
void
fpoly_FFT_int32(unsigned logn, fpr *d, const int32_t *a)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;

	/* The end of d holds n int32_t of scratch. */
	int32_t *cc = (int32_t *)(void *)(d + n) - n;
	if (a != cc) {
		memmove(cc, a, n * sizeof *cc);
	}

	/* c = (a + X^(n/2)b)/2, with a and b both self-adjoint.
	   The halving is done after the transforms, so the inputs stay
	   integral. */
	int32_t *aa = (int32_t *)(void *)(d + hn) - hn;
	aa[0] = cc[0] << 1;
	for (size_t j = 1; j < hn; j ++) {
		aa[j] = cc[j] - cc[n - j];
	}
	fpoly_FFT_int32_selfadj_inner(logn, d, aa, 1);

	int32_t *bb = cc + hn;
	bb[0] <<= 1;
	for (size_t j = 1; j < hn; j ++) {
		bb[j] += cc[hn - j];
	}
	fpoly_FFT_int32_selfadj_inner(logn, d + hn, bb, 1);
}

/* As fpoly_FFT_int32(), with aa and bb already in place. */
static void
fpoly_FFT_int32_finish(unsigned logn, fpr *d, int32_t *aa, int32_t *bb)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;

	fpoly_FFT_int32_selfadj_inner(logn, d, aa, 1);
	fpoly_FFT_int32_selfadj_inner(logn, d + hn, bb, 1);
}

/* see sign_inner.h */
void
fpoly_FFT_int8(unsigned logn, fpr *d, const int8_t *a)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	int32_t *aa = (int32_t *)(void *)(d + hn) - hn;
	int32_t *bb = (int32_t *)(void *)(d + n) - hn;

	/* aa and bb as in fpoly_FFT_int32(), built directly from a. */
	aa[0] = 2 * (int32_t)a[0];
	bb[0] = 2 * (int32_t)a[hn];
	for (size_t j = 1; j < hn; j ++) {
		aa[j] = (int32_t)a[j] - (int32_t)a[n - j];
		bb[j] = (int32_t)a[j + hn] + (int32_t)a[hn - j];
	}
	fpoly_FFT_int32_finish(logn, d, aa, bb);
}

/* see sign_inner.h */
void
fpoly_FFT_uint16(unsigned logn, fpr *d, const uint16_t *a)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	int32_t *aa = (int32_t *)(void *)(d + hn) - hn;
	int32_t *bb = (int32_t *)(void *)(d + n) - hn;

	aa[0] = 2 * (int32_t)a[0];
	bb[0] = 2 * (int32_t)a[hn];
	for (size_t j = 1; j < hn; j ++) {
		aa[j] = (int32_t)a[j] - (int32_t)a[n - j];
		bb[j] = (int32_t)a[j + hn] + (int32_t)a[hn - j];
	}
	fpoly_FFT_int32_finish(logn, d, aa, bb);
}
#endif /* FNDSA_FFT_NEW */

/* see sign_inner.h */
TARGET_SSE2 TARGET_NEON
void
fpoly_iFFT(unsigned logn, fpr *f)
{
#if FNDSA_SSE2
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	double *ff = (double *)f;

	/* Applying the reverse process of fpoly_FFT(), we separate the
	   first iteration (for which t = 1). */
	if (logn >= 2) {
		__m128d cz = _mm_castsi128_pd(
			_mm_setr_epi32(0, 0, 0, -0x80000000));
		for (size_t i = 0; i < hn; i += 2) {
			/* s <- re(s):im(s) */
			__m128d s = _mm_loadu_pd((const double *)GM + n + i);
			/* sc <- re(s):-im(s)  (conjugation) */
			__m128d sc = _mm_xor_pd(s, cz);
			/* xy_re <- re(x):re(y) */
			__m128d xy_re = _mm_loadu_pd(ff + i);
			/* xy_im <- im(x):im(y) */
			__m128d xy_im = _mm_loadu_pd(ff + i + hn);
			/* x <- re(x):im(x)
			   y <- re(y):im(y) */
			__m128d x = _mm_shuffle_pd(xy_re, xy_im, 0);
			__m128d y = _mm_shuffle_pd(xy_re, xy_im, 3);
			/* u <- x + y (complex) */
			__m128d u = _mm_add_pd(x, y);
			/* z <- x - y (complex) */
			__m128d z = _mm_sub_pd(x, y);
			/* v <- z*conj(s) (complex) */
			__m128d v1 = _mm_mul_pd(z, s);
			__m128d v2 = _mm_mul_pd(z, _mm_shuffle_pd(sc, sc, 1));
			__m128d v = _mm_add_pd(
				_mm_shuffle_pd(v1, v2, 2),
				_mm_shuffle_pd(v1, v2, 1));
			/* separate re/im and write */
			_mm_storeu_pd(ff + i, _mm_shuffle_pd(u, v, 0));
			_mm_storeu_pd(ff + i + hn, _mm_shuffle_pd(u, v, 3));
		}
	}

	size_t t = 2;
	for (unsigned lm = 2; lm < logn; lm ++) {
		size_t hm = (size_t)1 << (logn - lm);
		size_t m = hm << 1;
		size_t dt = t << 1;
		size_t j0 = 0;
		for (size_t i = 0; i < hm; i += 2) {
			__m128d s = _mm_loadu_pd((const double *)GM + m + i);
			__m128d s_re = _mm_shuffle_pd(s, s, 0);
			__m128d s_im = _mm_shuffle_pd(s, s, 3);
			for (size_t j = 0; j < t; j += 2) {
				size_t j1 = j0 + j;
				size_t j2 = j1 + t;
				__m128d x_re = _mm_loadu_pd(ff + j1);
				__m128d x_im = _mm_loadu_pd(ff + j1 + hn);
				__m128d y_re = _mm_loadu_pd(ff + j2);
				__m128d y_im = _mm_loadu_pd(ff + j2 + hn);
				_mm_storeu_pd(ff + j1,
					_mm_add_pd(x_re, y_re));
				_mm_storeu_pd(ff + j1 + hn,
					_mm_add_pd(x_im, y_im));
				__m128d u_re = _mm_sub_pd(x_re, y_re);
				__m128d u_im = _mm_sub_pd(x_im, y_im);
				/* Note: we loaded s but we want to
				   multiply with conj(s). */
				__m128d z_re = _mm_add_pd(
					_mm_mul_pd(u_re, s_re),
					_mm_mul_pd(u_im, s_im));
				__m128d z_im = _mm_sub_pd(
					_mm_mul_pd(u_im, s_re),
					_mm_mul_pd(u_re, s_im));
				_mm_storeu_pd(ff + j2, z_re);
				_mm_storeu_pd(ff + j2 + hn, z_im);
			}
			j0 += dt;
		}
		t = dt;
	}

	if (n >= 4) {
		/* MM[i] = 1/2^i */
		static const fpr MM[] = {
			FPR(4503599627370496, -52),
			FPR(4503599627370496, -53),
			FPR(4503599627370496, -54),
			FPR(4503599627370496, -55),
			FPR(4503599627370496, -56),
			FPR(4503599627370496, -57),
			FPR(4503599627370496, -58),
			FPR(4503599627370496, -59),
			FPR(4503599627370496, -60),
			FPR(4503599627370496, -61)
		};

		__m128d e = _mm_load_sd((const double *)MM + (logn - 1));
		e = _mm_shuffle_pd(e, e, 0);
		for (size_t i = 0; i < n; i += 2) {
			__m128d x = _mm_loadu_pd(ff + i);
			x = _mm_mul_pd(x, e);
			_mm_storeu_pd(ff + i, x);
		}
	}
#elif FNDSA_NEON
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	float64_t *ff = (float64_t *)f;

	/* Applying the reverse process of fpoly_FFT(), we separate the
	   first iteration (for which t = 1). */
	if (logn >= 2) {
		static const union { fpr f[2]; float64x2_t x; }
			cz = { { FPR_ZERO, FPR_NZERO } };
		for (size_t i = 0; i < hn; i += 2) {
			/* s <- re(s):im(s) */
			float64x2_t s = vld1q_f64(
				(const float64_t *)GM + n + i);
			/* sc <- re(s):-im(s)  (conjugation) */
			float64x2_t sc = vreinterpretq_f64_u64(
				veorq_u64(cz.x, vreinterpretq_u64_f64(s)));
			/* xy_re <- re(x):re(y) */
			float64x2_t xy_re = vld1q_f64(ff + i);
			/* xy_im <- im(x):im(y) */
			float64x2_t xy_im = vld1q_f64(ff + i + hn);
			/* x <- re(x):im(x)
			   y <- re(y):im(y) */
			float64x2_t x = vzip1q_f64(xy_re, xy_im);
			float64x2_t y = vzip2q_f64(xy_re, xy_im);
			/* u <- x + y (complex) */
			float64x2_t u = vaddq_f64(x, y);
			/* z <- x - y (complex) */
			float64x2_t z = vsubq_f64(x, y);
			/* v <- z*conj(s) (complex) */
			float64x2_t v1 = vmulq_f64(z, s);
			float64x2_t v2 = vmulq_f64(z, vextq_f64(sc, sc, 1));
			float64x2_t v = vpaddq_f64(v1, v2);
			/* separate re/im and write */
			vst1q_f64(ff + i, vzip1q_f64(u, v));
			vst1q_f64(ff + i + hn, vzip2q_f64(u, v));
		}
	}

	size_t t = 2;
	for (unsigned lm = 2; lm < logn; lm ++) {
		size_t hm = (size_t)1 << (logn - lm);
		size_t m = hm << 1;
		size_t dt = t << 1;
		size_t j0 = 0;
		for (size_t i = 0; i < hm; i += 2) {
			float64x2_t s = vld1q_f64(
				(const double *)GM + m + i);
			float64x2_t s_re = vzip1q_f64(s, s);
			float64x2_t s_im = vzip2q_f64(s, s);
			for (size_t j = 0; j < t; j += 2) {
				size_t j1 = j0 + j;
				size_t j2 = j1 + t;
				float64x2_t x_re = vld1q_f64(ff + j1);
				float64x2_t x_im = vld1q_f64(ff + j1 + hn);
				float64x2_t y_re = vld1q_f64(ff + j2);
				float64x2_t y_im = vld1q_f64(ff + j2 + hn);
				vst1q_f64(ff + j1, vaddq_f64(x_re, y_re));
				vst1q_f64(ff + j1 + hn, vaddq_f64(x_im, y_im));
				float64x2_t u_re = vsubq_f64(x_re, y_re);
				float64x2_t u_im = vsubq_f64(x_im, y_im);
				/* Note: we loaded s but we want to
				   multiply with conj(s). */
				float64x2_t z_re = vaddq_f64(
					vmulq_f64(u_re, s_re),
					vmulq_f64(u_im, s_im));
				float64x2_t z_im = vsubq_f64(
					vmulq_f64(u_im, s_re),
					vmulq_f64(u_re, s_im));
				vst1q_f64(ff + j2, z_re);
				vst1q_f64(ff + j2 + hn, z_im);
			}
			j0 += dt;
		}
		t = dt;
	}

	if (n >= 4) {
		/* MM[i] = 1/2^i */
		static const union { fpr f; float64x1_t v; } MM[] = {
			{ FPR(4503599627370496, -52) },
			{ FPR(4503599627370496, -53) },
			{ FPR(4503599627370496, -54) },
			{ FPR(4503599627370496, -55) },
			{ FPR(4503599627370496, -56) },
			{ FPR(4503599627370496, -57) },
			{ FPR(4503599627370496, -58) },
			{ FPR(4503599627370496, -59) },
			{ FPR(4503599627370496, -60) },
			{ FPR(4503599627370496, -61) }
		};

		float64x2_t e = vdupq_lane_f64(MM[logn - 1].v, 0);
		for (size_t i = 0; i < n; i += 2) {
			float64x2_t x = vld1q_f64(ff + i);
			x = vmulq_f64(x, e);
			vst1q_f64(ff + i, x);
		}
	}
#elif FNDSA_RV64D
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	size_t t = 1;
	f64 *ff = (f64 *)f;
	for (unsigned lm = 1; lm < logn; lm ++) {
		size_t hm = (size_t)1 << (logn - lm);
		size_t dt = t << 1;
		size_t j0 = 0;
		for (size_t i = 0; i < (hm >> 1); i ++) {
			f64 s_re = ((const f64 *)GM)[((hm + i) << 1) + 0];
			f64 s_im = ((const f64 *)GM)[((hm + i) << 1) + 1];
			for (size_t j = 0; j < t; j ++) {
				size_t j1 = j0 + j;
				size_t j2 = j1 + t;
				f64 x_re = ff[j1];
				f64 x_im = ff[j1 + hn];
				f64 y_re = ff[j2];
				f64 y_im = ff[j2 + hn];
				ff[j1] = f64_add(x_re, y_re);
				ff[j1 + hn] = f64_add(x_im, y_im);
				x_re = f64_sub(x_re, y_re);
				x_im = f64_sub(x_im, y_im);
				/* Note: multiply with conj(s), not s */
				ff[j2] = f64_add(
					f64_mul(x_re, s_re),
					f64_mul(x_im, s_im));
				ff[j2 + hn] = f64_sub(
					f64_mul(x_im, s_re),
					f64_mul(x_re, s_im));
			}
			j0 += dt;
		}
		t = dt;
	}

	/* MM[i] = 1/2^i */
	static const union { fpr f; f64 v; } MM[] = {
		{ FPR(4503599627370496, -52) },
		{ FPR(4503599627370496, -53) },
		{ FPR(4503599627370496, -54) },
		{ FPR(4503599627370496, -55) },
		{ FPR(4503599627370496, -56) },
		{ FPR(4503599627370496, -57) },
		{ FPR(4503599627370496, -58) },
		{ FPR(4503599627370496, -59) },
		{ FPR(4503599627370496, -60) },
		{ FPR(4503599627370496, -61) }
	};
	f64 z = MM[logn - 1].v;

	for (size_t i = 0; i < n; i ++) {
		ff[i] = f64_mul(ff[i], z);
	}
#else
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	size_t t = 1;
	for (unsigned lm = 1; lm < logn; lm ++) {
		size_t hm = (size_t)1 << (logn - lm);
		size_t dt = t << 1;
		size_t j0 = 0;
		for (size_t i = 0; i < (hm >> 1); i ++) {
			FPOLY_GM_TW_CONJ(s, hm + i);
			for (size_t j = 0; j < t; j ++) {
				size_t j1 = j0 + j;
				size_t j2 = j1 + t;
				fpr x_re = f[j1];
				fpr x_im = f[j1 + hn];
				fpr y_re = f[j2];
				fpr y_im = f[j2 + hn];
				FPOLY_GM_BFLY_INV(f[j1], f[j1 + hn],
					f[j2], f[j2 + hn],
					x_re, x_im, y_re, y_im, s);
			}
			j0 += dt;
		}
		t = dt;
	}

#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264
	if (logn == 9) {
		for (size_t i = 0; i < n; i ++) {
			f[i] = fpr_div2e8_fxp3264(f[i]);
		}
	} else if (logn == 10) {
		for (size_t i = 0; i < n; i ++) {
			f[i] = fpr_div2e9_fxp3264(f[i]);
		}
	} else
#endif
	{
		for (size_t i = 0; i < n; i ++) {
			f[i] = fpr_div2e(f[i], logn - 1);
		}
	}
#endif
}

#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 \
	&& !(FNDSA_SSE2 || FNDSA_NEON || FNDSA_RV64D)
/* see sign_inner.h */
void
fpoly_iFFT_no_div2e(unsigned logn, fpr *f)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	size_t t = 1;

	for (unsigned lm = 1; lm < logn; lm ++) {
		size_t hm = (size_t)1 << (logn - lm);
		size_t dt = t << 1;
		size_t j0 = 0;
		for (size_t i = 0; i < (hm >> 1); i ++) {
			FPOLY_GM_TW_CONJ(s, hm + i);
			for (size_t j = 0; j < t; j ++) {
				size_t j1 = j0 + j;
				size_t j2 = j1 + t;
				fpr x_re = f[j1];
				fpr x_im = f[j1 + hn];
				fpr y_re = f[j2];
				fpr y_im = f[j2 + hn];
				FPOLY_GM_BFLY_INV(f[j1], f[j1 + hn],
					f[j2], f[j2 + hn],
					x_re, x_im, y_re, y_im, s);
			}
			j0 += dt;
		}
		t = dt;
	}
}
#endif

/* see sign_inner.h */
TARGET_SSE2 TARGET_NEON
void
fpoly_set_small(unsigned logn, fpr *d, const int8_t *f)
{
	size_t n = (size_t)1 << logn;
#if FNDSA_SSE2
	size_t i = 0;
	double *dd = (double *)d;
	for (; (i + 4) <= n; i += 4) {
		__m128i x = _mm_setr_epi32(f[i], f[i + 1], f[i + 2], f[i + 3]);
		__m128d y0 = _mm_cvtepi32_pd(x);
		__m128d y1 = _mm_cvtepi32_pd(_mm_bsrli_si128(x, 8));
		_mm_storeu_pd(dd + i, y0);
		_mm_storeu_pd(dd + i + 2, y1);
	}
	__m128d z = _mm_setzero_pd();
	for (; i < n; i ++) {
		_mm_store_sd(dd + i, _mm_cvtsi32_sd(z, f[i]));
	}
#elif FNDSA_NEON
	size_t i = 0;
	float64_t *dd = (float64_t *)d;
	for (; (i + 8) <= n; i += 8) {
		int8x8_t x1 = vld1_s8(f + i);
		int16x8_t x2 = vmovl_s8(x1);
		int32x4_t x30 = vmovl_s16(vget_low_s16(x2));
		int32x4_t x31 = vmovl_high_s16(x2);
		int64x2_t x40 = vmovl_s32(vget_low_s32(x30));
		int64x2_t x41 = vmovl_high_s32(x30);
		int64x2_t x42 = vmovl_s32(vget_low_s32(x31));
		int64x2_t x43 = vmovl_high_s32(x31);

		float64x2_t y0 = vcvtq_f64_s64(x40);
		float64x2_t y1 = vcvtq_f64_s64(x41);
		float64x2_t y2 = vcvtq_f64_s64(x42);
		float64x2_t y3 = vcvtq_f64_s64(x43);
		vst1q_f64(dd + i + 0, y0);
		vst1q_f64(dd + i + 2, y1);
		vst1q_f64(dd + i + 4, y2);
		vst1q_f64(dd + i + 6, y3);
	}
	for (; i < n; i ++) {
		vst1_f64(dd + i, vcvt_f64_s64(vcreate_s64((uint64_t)f[i])));
	}
#elif FNDSA_RV64D
	f64 *dd = (f64 *)d;
	for (size_t i = 0; i < n; i ++) {
		dd[i] = f64_of32(f[i]);
	}
#else
	for (size_t i = 0; i < n; i ++) {
		d[i] = fpr_of32(f[i]);
	}
#endif
}

/* see sign_inner.h */
TARGET_SSE2 TARGET_NEON
void
fpoly_add(unsigned logn, fpr *a, const fpr *b)
{
	size_t n = (size_t)1 << logn;
#if FNDSA_SSE2
	if (n >= 2) {
		for (size_t i = 0; i < n; i += 2) {
			__m128d xa = _mm_loadu_pd((const double *)a + i);
			__m128d xb = _mm_loadu_pd((const double *)b + i);
			_mm_storeu_pd((double *)a + i, _mm_add_pd(xa, xb));
		}
	} else {
		__m128d xa = _mm_load_sd((const double *)a);
		__m128d xb = _mm_load_sd((const double *)b);
		_mm_store_sd((double *)a, _mm_add_sd(xa, xb));
	}
#elif FNDSA_NEON
	if (n >= 2) {
		for (size_t i = 0; i < n; i += 2) {
			float64x2_t xa = vld1q_f64((const float64_t *)a + i);
			float64x2_t xb = vld1q_f64((const float64_t *)b + i);
			vst1q_f64((float64_t *)a + i, vaddq_f64(xa, xb));
		}
	} else {
		float64x1_t xa = vld1_f64((const float64_t *)a);
		float64x1_t xb = vld1_f64((const float64_t *)b);
		vst1_f64((float64_t *)a, vadd_f64(xa, xb));
	}
#elif FNDSA_RV64D
	f64 *aa = (f64 *)a;
	const f64 *bb = (const f64 *)b;
	for (size_t i = 0; i < n; i ++) {
		aa[i] = f64_add(aa[i], bb[i]);
	}
#else
	for (size_t i = 0; i < n; i ++) {
		a[i] = fpr_add(a[i], b[i]);
	}
#endif
}

/* see sign_inner.h */
TARGET_SSE2 TARGET_NEON
void
fpoly_sub(unsigned logn, fpr *a, const fpr *b)
{
	size_t n = (size_t)1 << logn;
#if FNDSA_SSE2
	if (n >= 2) {
		for (size_t i = 0; i < n; i += 2) {
			__m128d xa = _mm_loadu_pd((const double *)a + i);
			__m128d xb = _mm_loadu_pd((const double *)b + i);
			_mm_storeu_pd((double *)a + i, _mm_sub_pd(xa, xb));
		}
	} else {
		__m128d xa = _mm_load_sd((const double *)a);
		__m128d xb = _mm_load_sd((const double *)b);
		_mm_store_sd((double *)a, _mm_sub_sd(xa, xb));
	}
#elif FNDSA_NEON
	if (n >= 2) {
		for (size_t i = 0; i < n; i += 2) {
			float64x2_t xa = vld1q_f64((const float64_t *)a + i);
			float64x2_t xb = vld1q_f64((const float64_t *)b + i);
			vst1q_f64((float64_t *)a + i, vsubq_f64(xa, xb));
		}
	} else {
		float64x1_t xa = vld1_f64((const float64_t *)a);
		float64x1_t xb = vld1_f64((const float64_t *)b);
		vst1_f64((float64_t *)a, vsub_f64(xa, xb));
	}
#elif FNDSA_RV64D
	f64 *aa = (f64 *)a;
	const f64 *bb = (const f64 *)b;
	for (size_t i = 0; i < n; i ++) {
		aa[i] = f64_sub(aa[i], bb[i]);
	}
#else
	for (size_t i = 0; i < n; i ++) {
		a[i] = fpr_sub(a[i], b[i]);
	}
#endif
}

/* see sign_inner.h */
TARGET_SSE2 TARGET_NEON
void
fpoly_neg(unsigned logn, fpr *a)
{
	size_t n = (size_t)1 << logn;
#if FNDSA_SSE2
	__m128d xz = _mm_setzero_pd();
	if (n >= 2) {
		for (size_t i = 0; i < n; i += 2) {
			__m128d xa = _mm_loadu_pd((const double *)a + i);
			_mm_storeu_pd((double *)a + i, _mm_sub_pd(xz, xa));
		}
	} else {
		__m128d xa = _mm_load_sd((const double *)a);
		_mm_store_sd((double *)a, _mm_sub_sd(xz, xa));
	}
#elif FNDSA_NEON
	if (n >= 2) {
		for (size_t i = 0; i < n; i += 2) {
			float64x2_t xa = vld1q_f64((const float64_t *)a + i);
			vst1q_f64((float64_t *)a + i, vnegq_f64(xa));
		}
	} else {
		float64x1_t xa = vld1_f64((const float64_t *)a);
		vst1_f64((float64_t *)a, vneg_f64(xa));
	}
#elif FNDSA_RV64D
	f64 *aa = (f64 *)a;
	for (size_t i = 0; i < n; i ++) {
		aa[i] = f64_neg(aa[i]);
	}
#else
	for (size_t i = 0; i < n; i ++) {
		a[i] = fpr_neg(a[i]);
	}
#endif
}

/* see sign_inner.h */
TARGET_SSE2 TARGET_NEON
void
fpoly_mul_fft(unsigned logn, fpr *a, const fpr *b)
{
	size_t hn = (size_t)1 << (logn - 1);
#if FNDSA_SSE2
	if (hn >= 2) {
		for (size_t i = 0; i < hn; i += 2) {
			__m128d xar = _mm_loadu_pd((const double *)a + i);
			__m128d xai = _mm_loadu_pd((const double *)a + i + hn);
			__m128d xbr = _mm_loadu_pd((const double *)b + i);
			__m128d xbi = _mm_loadu_pd((const double *)b + i + hn);

			__m128d xcr = _mm_sub_pd(
				_mm_mul_pd(xar, xbr),
				_mm_mul_pd(xai, xbi));
			__m128d xci = _mm_add_pd(
				_mm_mul_pd(xar, xbi),
				_mm_mul_pd(xai, xbr));

			_mm_storeu_pd((double *)a + i, xcr);
			_mm_storeu_pd((double *)a + i + hn, xci);
		}
	} else if (hn >= 1) {
		__m128d xa = _mm_loadu_pd((const double *)a);
		__m128d xb = _mm_loadu_pd((const double *)b);
		__m128d xcr = _mm_mul_pd(xa, xb);
		__m128d xci = _mm_mul_pd(xa, _mm_shuffle_pd(xb, xb, 1));
		xcr = _mm_sub_pd(xcr, _mm_shuffle_pd(xcr, xcr, 1));
		xci = _mm_add_pd(xci, _mm_shuffle_pd(xci, xci, 1));
		__m128d xc = _mm_shuffle_pd(xcr, xci, 0);
		_mm_storeu_pd((double *)a, xc);
	}
#elif FNDSA_NEON
	if (hn >= 2) {
		for (size_t i = 0; i < hn; i += 2) {
			float64x2_t xar =
				vld1q_f64((const float64_t *)a + i);
			float64x2_t xai =
				vld1q_f64((const float64_t *)a + i + hn);
			float64x2_t xbr =
				vld1q_f64((const float64_t *)b + i);
			float64x2_t xbi =
				vld1q_f64((const float64_t *)b + i + hn);

			float64x2_t xcr = vsubq_f64(
				vmulq_f64(xar, xbr),
				vmulq_f64(xai, xbi));
			float64x2_t xci = vaddq_f64(
				vmulq_f64(xar, xbi),
				vmulq_f64(xai, xbr));

			vst1q_f64((float64_t *)a + i, xcr);
			vst1q_f64((float64_t *)a + i + hn, xci);
		}
	} else if (hn >= 1) {
		static const union {
			uint64_t u[2];
			float64x2_t x;
		} cz = { { 0, (uint64_t)1 << 63 } };
		float64x2_t xa = vld1q_f64((const float64_t *)a);
		float64x2_t xb = vld1q_f64((const float64_t *)b);
		float64x2_t xcr = vmulq_f64(xa, xb);
		float64x2_t xci = vmulq_f64(xa, vextq_f64(xb, xb, 1));
		xcr = vreinterpretq_f64_u64(
			veorq_u64(vreinterpretq_u64_f64(xcr), cz.x));
		float64x2_t xc = vpaddq_f64(xcr, xci);
		vst1q_f64((float64_t *)a, xc);
	}
#elif FNDSA_RV64D
	f64 *aa = (f64 *)a;
	const f64 *bb = (const f64 *)b;
	for (size_t i = 0; i < hn; i ++) {
		f64 a_re = aa[i];
		f64 a_im = aa[i + hn];
		f64 b_re = bb[i];
		f64 b_im = bb[i + hn];
		aa[i] = f64_sub(f64_mul(a_re, b_re), f64_mul(a_im, b_im));
		aa[i + hn] = f64_add(f64_mul(a_im, b_re), f64_mul(a_re, b_im));
	}
#else
	if (a == b) {
		for (size_t i = 0; i < hn; i ++) {
			FPC_SQR_SELF(a[i], a[i + hn]);
		}
	} else {
		for (size_t i = 0; i < hn; i ++) {
			FPC_MUL_SELF(a[i], a[i + hn], b[i], b[i + hn]);
		}
	}
#endif
}

/* unused
TARGET_SSE2 TARGET_NEON
void
fpoly_muladj_fft(unsigned logn, fpr *a, const fpr *b)
{
	size_t hn = (size_t)1 << (logn - 1);
#if FNDSA_SSE2
	if (hn >= 2) {
		for (size_t i = 0; i < hn; i += 2) {
			__m128d xar = _mm_loadu_pd((const double *)a + i);
			__m128d xai = _mm_loadu_pd((const double *)a + i + hn);
			__m128d xbr = _mm_loadu_pd((const double *)b + i);
			__m128d xbi = _mm_loadu_pd((const double *)b + i + hn);

			__m128d xcr = _mm_add_pd(
				_mm_mul_pd(xar, xbr),
				_mm_mul_pd(xai, xbi));
			__m128d xci = _mm_sub_pd(
				_mm_mul_pd(xai, xbr),
				_mm_mul_pd(xar, xbi));

			_mm_storeu_pd((double *)a + i, xcr);
			_mm_storeu_pd((double *)a + i + hn, xci);
		}
	} else if (hn >= 1) {
		__m128d xa = _mm_loadu_pd((const double *)a);
		__m128d xb = _mm_loadu_pd((const double *)b);
		__m128d xcr = _mm_mul_pd(xa, xb);
		__m128d xci = _mm_mul_pd(xa, _mm_shuffle_pd(xb, xb, 1));
		xcr = _mm_add_pd(xcr, _mm_shuffle_pd(xcr, xcr, 1));
		xci = _mm_sub_pd(xci, _mm_shuffle_pd(xci, xci, 1));
		__m128d xc = _mm_shuffle_pd(xcr, xci, 2);
		_mm_storeu_pd((double *)a, xc);
	}
#elif FNDSA_NEON
	if (hn >= 2) {
		for (size_t i = 0; i < hn; i += 2) {
			float64x2_t xar =
				vld1q_f64((const float64_t *)a + i);
			float64x2_t xai =
				vld1q_f64((const float64_t *)a + i + hn);
			float64x2_t xbr =
				vld1q_f64((const float64_t *)b + i);
			float64x2_t xbi =
				vld1q_f64((const float64_t *)b + i + hn);

			float64x2_t xcr = vaddq_f64(
				vmulq_f64(xar, xbr),
				vmulq_f64(xai, xbi));
			float64x2_t xci = vsubq_f64(
				vmulq_f64(xai, xbr),
				vmulq_f64(xar, xbi));

			vst1q_f64((float64_t *)a + i, xcr);
			vst1q_f64((float64_t *)a + i + hn, xci);
		}
	} else if (hn >= 1) {
		static const union {
			uint64_t u[2];
			float64x2_t x;
		} cz = { { 0, (uint64_t)1 << 63 } };
		float64x2_t xa = vld1q_f64((const float64_t *)a);
		float64x2_t xb = vld1q_f64((const float64_t *)b);
		float64x2_t xcr = vmulq_f64(xa, xb);
		xb = vreinterpretq_f64_u64(
			veorq_u64(vreinterpretq_u64_f64(xb), cz.x));
		float64x2_t xci = vmulq_f64(xa, vextq_f64(xb, xb, 1));
		float64x2_t xc = vpaddq_f64(xcr, xci);
		vst1q_f64((float64_t *)a, xc);
	}
#elif FNDSA_RV64D
	f64 *aa = (f64 *)a;
	const f64 *bb = (const f64 *)b;
	for (size_t i = 0; i < hn; i ++) {
		f64 a_re = aa[i];
		f64 a_im = aa[i + hn];
		f64 b_re = bb[i];
		f64 b_im = bb[i + hn];
		aa[i] = f64_add(f64_mul(a_re, b_re), f64_mul(a_im, b_im));
		aa[i + hn] = f64_sub(f64_mul(a_im, b_re), f64_mul(a_re, b_im));
	}
#else
	for (size_t i = 0; i < hn; i ++) {
		FPC_MUL(a[i], a[i + hn],
			a[i], a[i + hn], b[i], fpr_neg(b[i + hn]));
	}
#endif
}
*/

/* unused
TARGET_SSE2 TARGET_NEON
void
fpoly_mulownadj_fft(unsigned logn, fpr *a)
{
	size_t hn = (size_t)1 << (logn - 1);
#if FNDSA_SSE2
	__m128d xz = _mm_setzero_pd();
	if (hn >= 2) {
		for (size_t i = 0; i < hn; i += 2) {
			__m128d xar = _mm_loadu_pd((const double *)a + i);
			__m128d xai = _mm_loadu_pd((const double *)a + i + hn);
			__m128d xcr = _mm_add_pd(
				_mm_mul_pd(xar, xar),
				_mm_mul_pd(xai, xai));
			_mm_storeu_pd((double *)a + i, xcr);
			_mm_storeu_pd((double *)a + i + hn, xz);
		}
	} else if (hn >= 1) {
		__m128d xa = _mm_loadu_pd((const double *)a);
		__m128d xcr = _mm_mul_pd(xa, xa);
		xcr = _mm_add_pd(xcr, _mm_shuffle_pd(xcr, xcr, 1));
		__m128d xc = _mm_shuffle_pd(xcr, xz, 2);
		_mm_storeu_pd((double *)a, xc);
	}
#elif FNDSA_NEON
	float64x2_t xz = vdupq_lane_f64(vcreate_f64(0), 0);
	if (hn >= 2) {
		for (size_t i = 0; i < hn; i += 2) {
			float64x2_t xar = vld1q_f64(
				(const float64_t *)a + i);
			float64x2_t xai = vld1q_f64(
				(const float64_t *)a + i + hn);
			float64x2_t xcr = vaddq_f64(
				vmulq_f64(xar, xar),
				vmulq_f64(xai, xai));
			vst1q_f64((float64_t *)a + i, xcr);
			vst1q_f64((float64_t *)a + i + hn, xz);
		}
	} else if (hn >= 1) {
		float64x2_t xa = vld1q_f64((const float64_t *)a);
		float64x2_t xcr = vmulq_f64(xa, xa);
		float64x2_t xc = vpaddq_f64(xa, xz);
		vst1q_f64((float64_t *)a, xc);
	}
#elif FNDSA_RV64D
	f64 *aa = (f64 *)a;
	for (size_t i = 0; i < hn; i ++) {
		f64 a_re = aa[i];
		f64 a_im = aa[i + hn];
		aa[i] = f64_add(f64_mul(a_re, a_re), f64_mul(a_im, a_im));
	}
	memset(aa + hn, 0, hn * sizeof(f64));
#else
	for (size_t i = 0; i < hn; i ++) {
		a[i] = fpr_add(fpr_sqr(a[i]), fpr_sqr(a[i + hn]));
		a[i + hn] = FPR_ZERO;
	}
#endif
}
*/

/* see sign_inner.h */
TARGET_SSE2 TARGET_NEON
void
fpoly_mulconst(unsigned logn, fpr *a, fpr x)
{
	size_t n = (size_t)1 << logn;
#if FNDSA_SSE2
	__m128d xx = _mm_load_sd((const double *)&x);
	if (n >= 2) {
		xx = _mm_shuffle_pd(xx, xx, 0);
		for (size_t i = 0; i < n; i += 2) {
			__m128d xa = _mm_loadu_pd((const double *)a + i);
			_mm_storeu_pd((double *)a + i, _mm_mul_pd(xa, xx));
		}
	} else {
		__m128d xa = _mm_load_sd((const double *)a);
		_mm_store_sd((double *)a, _mm_mul_sd(xa, xx));
	}
#elif FNDSA_NEON
	float64x1_t x1 = vld1_f64((const float64_t *)&x);
	if (n >= 2) {
		float64x2_t x2 = vdupq_lane_f64(x1, 0);
		for (size_t i = 0; i < n; i += 2) {
			float64x2_t xa = vld1q_f64((const float64_t *)a + i);
			vst1q_f64((float64_t *)a + i, vmulq_f64(xa, x2));
		}
	} else {
		float64x1_t xa = vld1_f64((const float64_t *)a);
		vst1_f64((float64_t *)a, vmul_f64(xa, x1));
	}
#elif FNDSA_RV64D
	f64 z = f64_from_raw(x);
	f64 *aa = (f64 *)a;
	for (size_t i = 0; i < n; i ++) {
		aa[i] = f64_mul(aa[i], z);
	}
#else
	for (size_t i = 0; i < n; i ++) {
		a[i] = fpr_mul(a[i], x);
	}
#endif
}

#if FNDSA_SIGN_INVQ_WIDE
/* see sign_inner.h */
void
fpoly_mul_inv_q(unsigned logn, fpr *a)
{
	size_t n = (size_t)1 << logn;
	for (size_t i = 0; i < n; i ++) {
		a[i] = FPR_MUL_INV_Q(a[i]);
	}
}
#endif

/* see sign_inner.h */
TARGET_SSE2 TARGET_NEON
void
fpoly_LDL_fft(unsigned logn, const fpr *g00, const fpr *g01, fpr *g11,
	l10_t *l10)
{
	size_t hn = (size_t)1 << (logn - 1);
#if FNDSA_SSE2
	static const union {
		fpr f[2];
		__m128d x;
	} one = { { FPR_ONE, FPR_ONE } };
	__m128d nz = _mm_castsi128_pd(
		_mm_setr_epi32(0, -0x80000000, 0, -0x80000000));
	const double *p00 = (const double *)g00;
	const double *p01 = (const double *)g01;
	double *p11 = (double *)g11;
	double *pl = (double *)l10;
	if (hn >= 2) {
		for (size_t i = 0; i < hn; i += 2) {
			__m128d g00_re = _mm_loadu_pd(p00 + i);
			__m128d g01_re = _mm_loadu_pd(p01 + i);
			__m128d g01_im = _mm_loadu_pd(p01 + i + hn);
			__m128d g11_re = _mm_loadu_pd(p11 + i);
			__m128d inv_g00_re = _mm_div_pd(one.x, g00_re);
			__m128d mu_re = _mm_mul_pd(g01_re, inv_g00_re);
			__m128d mu_im = _mm_mul_pd(g01_im, inv_g00_re);
			__m128d zo_re = _mm_add_pd(
				_mm_mul_pd(mu_re, g01_re),
				_mm_mul_pd(mu_im, g01_im));
			_mm_storeu_pd(p11 + i, _mm_sub_pd(g11_re, zo_re));
			_mm_storeu_pd(pl + i, mu_re);
			_mm_storeu_pd(pl + i + hn, _mm_xor_pd(nz, mu_im));
		}
	} else {
		__m128d g00_re = _mm_load_sd(p00);
		__m128d g01_re = _mm_load_sd(p01);
		__m128d g01_im = _mm_load_sd(p01 + 1);
		__m128d g11_re = _mm_load_sd(p11);
		__m128d inv_g00_re = _mm_div_sd(one.x, g00_re);
		__m128d mu_re = _mm_mul_sd(g01_re, inv_g00_re);
		__m128d mu_im = _mm_mul_sd(g01_im, inv_g00_re);
		__m128d zo_re = _mm_add_sd(
			_mm_mul_sd(mu_re, g01_re),
			_mm_mul_sd(mu_im, g01_im));
		_mm_store_sd(p11, _mm_sub_sd(g11_re, zo_re));
		_mm_store_sd(pl, mu_re);
		_mm_store_sd(pl + 1, _mm_xor_pd(nz, mu_im));
	}
#elif FNDSA_NEON
	static const union {
		fpr f[2];
		float64x1_t s;
		float64x2_t x;
	} one = { { FPR_ONE, FPR_ONE } };
	const float64_t *p00 = (const float64_t *)g00;
	const float64_t *p01 = (const float64_t *)g01;
	float64_t *p11 = (float64_t *)g11;
	float64_t *pl = (float64_t *)l10;
	if (hn >= 2) {
		for (size_t i = 0; i < hn; i += 2) {
			float64x2_t g00_re = vld1q_f64(p00 + i);
			float64x2_t g01_re = vld1q_f64(p01 + i);
			float64x2_t g01_im = vld1q_f64(p01 + i + hn);
			float64x2_t g11_re = vld1q_f64(p11 + i);
			float64x2_t inv_g00_re = vdivq_f64(one.x, g00_re);
			float64x2_t mu_re = vmulq_f64(g01_re, inv_g00_re);
			float64x2_t mu_im = vmulq_f64(g01_im, inv_g00_re);
			float64x2_t zo_re = vaddq_f64(
				vmulq_f64(mu_re, g01_re),
				vmulq_f64(mu_im, g01_im));
			vst1q_f64(p11 + i, vsubq_f64(g11_re, zo_re));
			vst1q_f64(pl + i, mu_re);
			vst1q_f64(pl + i + hn, vnegq_f64(mu_im));
		}
	} else {
		float64x1_t g00_re = vld1_f64(p00);
		float64x1_t g01_re = vld1_f64(p01);
		float64x1_t g01_im = vld1_f64(p01 + 1);
		float64x1_t g11_re = vld1_f64(p11);
		float64x1_t inv_g00_re = vdiv_f64(one.s, g00_re);
		float64x1_t mu_re = vmul_f64(g01_re, inv_g00_re);
		float64x1_t mu_im = vmul_f64(g01_im, inv_g00_re);
		float64x1_t zo_re = vadd_f64(
			vmul_f64(mu_re, g01_re),
			vmul_f64(mu_im, g01_im));
		vst1_f64(p11, vsub_f64(g11_re, zo_re));
		vst1_f64(pl, mu_re);
		vst1_f64(pl + 1, vneg_f64(mu_im));
	}
#elif FNDSA_RV64D
	const f64 *gg00 = (const f64 *)g00;
	const f64 *gg01 = (const f64 *)g01;
	f64 *gg11 = (f64 *)g11;
	f64 *ll = (f64 *)l10;
	for (size_t i = 0; i < hn; i ++) {
		f64 g00_re = gg00[i];
		f64 g01_re = gg01[i], g01_im = gg01[i + hn];
		f64 g11_re = gg11[i];
		f64 inv_g00_re = f64_inv(g00_re);
		f64 mu_re = f64_mul(g01_re, inv_g00_re);
		f64 mu_im = f64_mul(g01_im, inv_g00_re);
		f64 zo_re = f64_add(
			f64_mul(mu_re, g01_re),
			f64_mul(mu_im, g01_im));
		gg11[i] = f64_sub(g11_re, zo_re);
		ll[i] = mu_re;
		ll[i + hn] = f64_neg(mu_im);
	}
#else
	size_t l10_im = L10_IM_OFF(logn);
	for (size_t i = 0; i < hn; i ++) {
		fpr g00_re = g00[i];
		fpr g01_re = g01[i], g01_im = g01[i + hn];
		fpr g11_re = g11[i];
		fpr_rcp inv_g00;
		fpr mu_re, mu_im;
		l10_t s_re, s_im;
		FPR_RECIP_X(inv_g00, g00_re);
		FPR_L10_DIV(s_re, s_im, mu_re, mu_im, g01_re, g01_im, inv_g00);
		fpr zo_re = fpr_add(
			fpr_mul(mu_re, g01_re),
			fpr_mul(mu_im, g01_im));
		g11[i] = fpr_sub(g11_re, zo_re);
		l10[i] = s_re;
		l10[i + l10_im] = s_im;
	}
#endif
}

#if FNDSA_SYMP_ROOT
/* see sign_inner.h
   d11 = q^2/g00 avoids the cancellation in g11 - |g01|^2/g00. The input
   value of g11 is not used. */
void
fpoly_LDL_root_fft(unsigned logn, const fpr *g00, fpr *g01, fpr *g11)
{
	size_t hn = (size_t)1 << (logn - 1);
	fpr q2 = Q2;

	for (size_t i = 0; i < hn; i ++) {
		fpr g00_re = g00[i];
		fpr g01_re = g01[i], g01_im = g01[i + hn];
		fpr_rcp inv_g00;
		fpr mu_re, mu_im;
		FPR_RECIP_X(inv_g00, g00_re);
		FPR_MUL_RECIP(mu_re, g01_re, inv_g00);
		FPR_MUL_RECIP(mu_im, g01_im, inv_g00);
		FPR_MUL_RECIP(g11[i], q2, inv_g00);
		g01[i] = mu_re;
		g01[i + hn] = fpr_neg(mu_im);
	}
}

#elif FNDSA_SYMP_TREE

/* see sign_inner.h */
void
fpoly_LDL_root_fft(unsigned logn, const fpr *g00, fpr *g01)
{
	size_t hn = (size_t)1 << (logn - 1);
	fpr q2 = Q2;

	for (size_t i = 0; i < hn; i ++) {
		fpr g00_re = g00[i];
		fpr g01_re = g01[i], g01_im = g01[i + hn];
		fpr_rcp inv_g00;
		fpr mu_re, mu_im;
		FPR_RECIP_X(inv_g00, g00_re);
		FPR_MUL_RECIP(mu_re, g01_re, inv_g00);
		FPR_MUL_RECIP(mu_im, g01_im, inv_g00);
		g01[i] = mu_re;
		g01[i + hn] = fpr_neg(mu_im);
	}
}

#endif

/* see sign_inner.h */
TARGET_SSE2 TARGET_NEON
void
fpoly_split_fft(unsigned logn, fpr *f0, fpr *f1, const fpr *f)
{
	size_t hn = (size_t)1 << (logn - 1);
	size_t qn = hn >> 1;

	if (logn == 1) {
		f0[0] = f[0];
		f1[0] = f[hn];
	}

#if FNDSA_SSE2
	static union {
		fpr f[2];
		__m128d x;
	} h = { {
		FPR(4503599627370496, -53), FPR(4503599627370496, -53)
	} };
	const double *ff = (const double *)f;
	double *ff0 = (double *)f0;
	double *ff1 = (double *)f1;
	__m128d cz = _mm_castsi128_pd(_mm_setr_epi32(0, 0, 0, -0x80000000));
	for (size_t i = 0; i < qn; i ++) {
		__m128d ab_re = _mm_loadu_pd(ff + (i << 1));
		__m128d ab_im = _mm_loadu_pd(ff + (i << 1) + hn);
		__m128d a = _mm_shuffle_pd(ab_re, ab_im, 0);
		__m128d b = _mm_shuffle_pd(ab_re, ab_im, 3);
		__m128d u = _mm_add_pd(a, b);
		__m128d v = _mm_sub_pd(a, b);
		__m128d s = _mm_loadu_pd((const double *)GM + ((i + hn) << 1));
		__m128d sc = _mm_xor_pd(s, cz);
		/* We compute w = v*conj(s) */
		__m128d w1 = _mm_mul_pd(v, s);
		__m128d w2 = _mm_mul_pd(v, _mm_shuffle_pd(sc, sc, 1));
		__m128d w = _mm_add_pd(
			_mm_shuffle_pd(w1, w2, 0),
			_mm_shuffle_pd(w1, w2, 3));
		u = _mm_mul_pd(u, h.x);
		w = _mm_mul_pd(w, h.x);
		_mm_store_sd(ff0 + i, u);
		_mm_store_sd(ff0 + i + qn, _mm_shuffle_pd(u, u, 1));
		_mm_store_sd(ff1 + i, w);
		_mm_store_sd(ff1 + i + qn, _mm_shuffle_pd(w, w, 1));
	}
#elif FNDSA_NEON
	static union { fpr f[2]; uint64x2_t w; float64x2_t x; }
		h = { {
			FPR(4503599627370496, -53), FPR(4503599627370496, -53)
		} },
		cz = { {
			FPR_ZERO, FPR_NZERO,
		} };
	const float64_t *ff = (const float64_t *)f;
	float64_t *ff0 = (float64_t *)f0;
	float64_t *ff1 = (float64_t *)f1;
	for (size_t i = 0; i < qn; i ++) {
		float64x2_t ab_re = vld1q_f64(ff + (i << 1));
		float64x2_t ab_im = vld1q_f64(ff + (i << 1) + hn);
		float64x2_t a = vzip1q_f64(ab_re, ab_im);
		float64x2_t b = vzip2q_f64(ab_re, ab_im);
		float64x2_t u = vaddq_f64(a, b);
		float64x2_t v = vsubq_f64(a, b);
		float64x2_t s = vld1q_f64(
			(const float64_t *)GM + ((i + hn) << 1));
		float64x2_t sc = vreinterpretq_f64_u64(
			veorq_u64(vreinterpretq_u64_f64(s), cz.w));
		/* We compute w = v*conj(s) */
		float64x2_t w1 = vmulq_f64(v, s);
		float64x2_t w2 = vmulq_f64(v, vextq_f64(sc, sc, 1));
		float64x2_t w = vaddq_f64(
			vzip1q_f64(w1, w2),
			vzip2q_f64(w1, w2));
		u = vmulq_f64(u, h.x);
		w = vmulq_f64(w, h.x);
		vst1_f64(ff0 + i, vget_low_f64(u));
		vst1_f64(ff0 + i + qn, vget_high_f64(u));
		vst1_f64(ff1 + i, vget_low_f64(w));
		vst1_f64(ff1 + i + qn, vget_high_f64(w));
	}
#elif FNDSA_RV64D
	const f64 *ff = (const f64 *)f;
	f64 *ff0 = (f64 *)f0;
	f64 *ff1 = (f64 *)f1;
	for (size_t i = 0; i < qn; i ++) {
		f64 a_re = ff[(i << 1) + 0], a_im = ff[(i << 1) + 0 + hn];
		f64 b_re = ff[(i << 1) + 1], b_im = ff[(i << 1) + 1 + hn];
		f64 t_re, t_im;

		t_re = f64_add(a_re, b_re);
		t_im = f64_add(a_im, b_im);
		ff0[i] = f64_half(t_re);
		ff0[i + qn] = f64_half(t_im);

		t_re = f64_sub(a_re, b_re);
		t_im = f64_sub(a_im, b_im);
		f64 u_re = ((const f64 *)GM)[((i + hn) << 1) + 0];
		f64 u_im = ((const f64 *)GM)[((i + hn) << 1) + 1];
		f64 v_re = f64_add(f64_mul(t_re, u_re), f64_mul(t_im, u_im));
		f64 v_im = f64_sub(f64_mul(t_im, u_re), f64_mul(t_re, u_im));
		ff1[i] = f64_half(v_re);
		ff1[i + qn] = f64_half(v_im);
	}
#else
	for (size_t i = 0; i < qn; i ++) {
		FPOLY_GM_SPLIT(f0, f1, f, i, hn, qn, i + hn);
	}
#endif
}

/* see sign_inner.h */
TARGET_SSE2 TARGET_NEON
void
fpoly_split_selfadj_fft(unsigned logn, fpr *f0, fpr *f1, const fpr *f)
{
	size_t hn = (size_t)1 << (logn - 1);
	size_t qn = hn >> 1;

	if (logn == 1) {
		f0[0] = f[0];
		f1[0] = FPR_ZERO;
	}

#if FNDSA_SSE2
	static union { fpr f[2]; __m128d x; }
		h = { {
			FPR(4503599627370496, -53), FPR(4503599627370496, -53)
		} };
	const double *ff = (const double *)f;
	double *ff0 = (double *)f0;
	double *ff1 = (double *)f1;
	for (size_t i = 0; i < qn; i ++) {
		__m128d ab_re = _mm_loadu_pd(ff + (i << 1));
		__m128d t = _mm_shuffle_pd(ab_re, ab_re, 1);
		__m128d u = _mm_mul_pd(h.x, _mm_add_pd(ab_re, t));
		__m128d v = _mm_mul_pd(h.x, _mm_sub_pd(ab_re, t));
		__m128d s = _mm_loadu_pd((const double *)GM + ((i + hn) << 1));
		__m128d w = _mm_mul_pd(v, s);
		_mm_store_sd(ff0 + i, u);
		_mm_store_sd(ff0 + i + qn, _mm_setzero_pd());
		_mm_store_sd(ff1 + i, w);
		_mm_store_sd(ff1 + i + qn, _mm_shuffle_pd(w, w, 1));
	}
#elif FNDSA_NEON
	static union { fpr f[2]; float64x2_t x; }
		h = { {
			FPR(4503599627370496, -53), FPR(4503599627370496, -53)
		} };
	const float64_t *ff = (const float64_t *)f;
	float64_t *ff0 = (float64_t *)f0;
	float64_t *ff1 = (float64_t *)f1;
	float64x1_t z1 = vcreate_f64(0);
	for (size_t i = 0; i < qn; i ++) {
		float64x2_t ab_re = vld1q_f64(ff + (i << 1));
		float64x2_t t = vextq_f64(ab_re, ab_re, 1);
		float64x2_t u = vmulq_f64(h.x, vaddq_f64(ab_re, t));
		float64x2_t v = vmulq_f64(h.x, vsubq_f64(ab_re, t));
		float64x2_t s = vld1q_f64(
			(const float64_t *)GM + ((i + hn) << 1));
		float64x2_t w = vmulq_f64(v, s);
		vst1_f64(ff0 + i, vget_low_f64(u));
		vst1_f64(ff0 + i + qn, z1);
		vst1_f64(ff1 + i, vget_low_f64(w));
		vst1_f64(ff1 + i + qn, vget_high_f64(w));
	}
#elif FNDSA_RV64D
	const f64 *ff = (const f64 *)f;
	f64 *ff0 = (f64 *)f0;
	f64 *ff1 = (f64 *)f1;
	for (size_t i = 0; i < qn; i ++) {
		f64 a_re = ff[(i << 1) + 0];
		f64 b_re = ff[(i << 1) + 1];
		f64 t_re;

		t_re = f64_add(a_re, b_re);
		ff0[i] = f64_half(t_re);
		ff0[i + qn] = (f64){ 0.0 };

		t_re = f64_half(f64_sub(a_re, b_re));
		ff1[i] = f64_mul(t_re,
			((const f64 *)GM)[((i + hn) << 1) + 0]);
		ff1[i + qn] = f64_mul(t_re,
			f64_neg(((const f64 *)GM)[((i + hn) << 1) + 1]));
	}
#else
	for (size_t i = 0; i < qn; i ++) {
		fpr a_re = f[(i << 1) + 0];
		fpr b_re = f[(i << 1) + 1];
		FPOLY_GM_SPLIT_SELFADJ(f0[i], f0[i + qn],
			f1[i], f1[i + qn], a_re, b_re, i + hn);
	}
#endif
}

/* see sign_inner.h */
TARGET_SSE2 TARGET_NEON
void
fpoly_merge_fft(unsigned logn, fpr *f, const fpr *f0, const fpr *f1)
{
	size_t hn = (size_t)1 << (logn - 1);
	size_t qn = hn >> 1;

	if (logn == 1) {
		f[0] = f0[0];
		f[hn] = f1[0];
	}

#if FNDSA_SSE2
	const double *ff0 = (const double *)f0;
	const double *ff1 = (const double *)f1;
	double *ff = (double *)f;
	__m128d cz = _mm_castsi128_pd(_mm_setr_epi32(0, 0, 0, -0x80000000));
	for (size_t i = 0; i < qn; i ++) {
		__m128d a_re = _mm_load_sd(ff0 + i);
		__m128d a_im = _mm_load_sd(ff0 + i + qn);
		__m128d b_re = _mm_load_sd(ff1 + i);
		__m128d b_im = _mm_load_sd(ff1 + i + qn);

		__m128d s = _mm_loadu_pd((const double *)GM + ((i + hn) << 1));
		__m128d c1 = _mm_mul_pd(s, _mm_shuffle_pd(b_re, b_im, 0));
		__m128d c2 = _mm_mul_pd(s, _mm_shuffle_pd(b_im, b_re, 0));

		/* c_re <- re(b*s):-re(b*s)
		   c_im <- im(b*s):-im(b*s) */
		__m128d c_re = _mm_sub_pd(c1, _mm_shuffle_pd(c1, c1, 1));
		__m128d c_im = _mm_xor_pd(cz,
			_mm_add_pd(c2, _mm_shuffle_pd(c2, c2, 1)));

		_mm_storeu_pd(ff + (i << 1),
			_mm_add_pd(c_re, _mm_shuffle_pd(a_re, a_re, 0)));
		_mm_storeu_pd(ff + (i << 1) + hn,
			_mm_add_pd(c_im, _mm_shuffle_pd(a_im, a_im, 0)));
	}
#elif FNDSA_NEON
	static const union { fpr f[2]; uint64x2_t w; }
		cz = { { FPR_ZERO, FPR_NZERO } };
	const float64_t *ff0 = (const float64_t *)f0;
	const float64_t *ff1 = (const float64_t *)f1;
	float64_t *ff = (float64_t *)f;
	for (size_t i = 0; i < qn; i ++) {
		float64x1_t a_re = vld1_f64(ff0 + i);
		float64x1_t a_im = vld1_f64(ff0 + i + qn);
		float64x1_t b_re = vld1_f64(ff1 + i);
		float64x1_t b_im = vld1_f64(ff1 + i + qn);
		float64x2_t b = vcombine_f64(b_re, b_im);

		float64x2_t s = vld1q_f64(
			(const float64_t *)GM + ((i + hn) << 1));
		float64x2_t c1 = vmulq_f64(s, b);
		float64x2_t c2 = vmulq_f64(s, vextq_f64(b, b, 1));

		/* c_re <- re(b*s):-re(b*s)
		   c_im <- im(b*s):-im(b*s) */
		float64x2_t c_re = vsubq_f64(c1, vextq_f64(c1, c1, 1));
		float64x2_t c_im = vreinterpretq_f64_u64(
			veorq_u64(cz.w, vreinterpretq_u64_f64(
				vaddq_f64(c2, vextq_f64(c2, c2, 1)))));

		vst1q_f64(ff + (i << 1),
			vaddq_f64(c_re, vdupq_lane_f64(a_re, 0)));
		vst1q_f64(ff + (i << 1) + hn,
			vaddq_f64(c_im, vdupq_lane_f64(a_im, 0)));
	}
#elif FNDSA_RV64D
	const f64 *ff0 = (const f64 *)f0;
	const f64 *ff1 = (const f64 *)f1;
	f64 *ff = (f64 *)f;
	for (size_t i = 0; i < qn; i ++) {
		f64 a_re = ff0[i], a_im = ff0[i + qn];
		f64 b_re = ff1[i], b_im = ff1[i + qn];
		f64 s_re = ((const f64 *)GM)[((i + hn) << 1) + 0];
		f64 s_im = ((const f64 *)GM)[((i + hn) << 1) + 1];
		f64 c_re = f64_sub(f64_mul(b_re, s_re), f64_mul(b_im, s_im));
		f64 c_im = f64_add(f64_mul(b_im, s_re), f64_mul(b_re, s_im));
		ff[(i << 1) + 0] = f64_add(a_re, c_re);
		ff[(i << 1) + 0 + hn] = f64_add(a_im, c_im);
		ff[(i << 1) + 1] = f64_sub(a_re, c_re);
		ff[(i << 1) + 1 + hn] = f64_sub(a_im, c_im);
	}
#else
	for (size_t i = 0; i < qn; i ++) {
		FPOLY_GM_MERGE(f, f0, f1, i, hn, qn, i + hn);
	}
#endif
}

/* see sign_inner.h */
TARGET_SSE2 TARGET_NEON
void
fpoly_gram_fft(unsigned logn, fpr *b00, fpr *b01, fpr *b10, const fpr *b11)
{
	size_t hn = (size_t)1 << (logn - 1);
#if FNDSA_SSE2
	double *p00 = (double *)b00;
	double *p01 = (double *)b01;
	double *p10 = (double *)b10;
	const double *p11 = (const double *)b11;
	for (size_t i = 0; i < hn; i += 2) {
		__m128d b00_re = _mm_loadu_pd(p00 + i);
		__m128d b00_im = _mm_loadu_pd(p00 + i + hn);
		__m128d b01_re = _mm_loadu_pd(p01 + i);
		__m128d b01_im = _mm_loadu_pd(p01 + i + hn);
		__m128d b10_re = _mm_loadu_pd(p10 + i);
		__m128d b10_im = _mm_loadu_pd(p10 + i + hn);
		__m128d b11_re = _mm_loadu_pd(p11 + i);
		__m128d b11_im = _mm_loadu_pd(p11 + i + hn);

		/* g00 = b00*adj(b00) + b01*adj(b01) */
		__m128d g00_re = _mm_add_pd(
			_mm_add_pd(
				_mm_mul_pd(b00_re, b00_re),
				_mm_mul_pd(b00_im, b00_im)),
			_mm_add_pd(
				_mm_mul_pd(b01_re, b01_re),
				_mm_mul_pd(b01_im, b01_im)));
		/* g01 = b00*adj(b10) + b01*adj(b11) */
		__m128d u_re = _mm_add_pd(
			_mm_mul_pd(b00_re, b10_re),
			_mm_mul_pd(b00_im, b10_im));
		__m128d u_im = _mm_sub_pd(
			_mm_mul_pd(b00_im, b10_re),
			_mm_mul_pd(b00_re, b10_im));
		__m128d v_re = _mm_add_pd(
			_mm_mul_pd(b01_re, b11_re),
			_mm_mul_pd(b01_im, b11_im));
		__m128d v_im = _mm_sub_pd(
			_mm_mul_pd(b01_im, b11_re),
			_mm_mul_pd(b01_re, b11_im));
		__m128d g01_re = _mm_add_pd(u_re, v_re);
		__m128d g01_im = _mm_add_pd(u_im, v_im);
		/* g11 = b10*adj(b10) + b11*adj(b11) */
		__m128d g11_re = _mm_add_pd(
			_mm_add_pd(
				_mm_mul_pd(b10_re, b10_re),
				_mm_mul_pd(b10_im, b10_im)),
			_mm_add_pd(
				_mm_mul_pd(b11_re, b11_re),
				_mm_mul_pd(b11_im, b11_im)));

		_mm_storeu_pd(p00 + i, g00_re);
		_mm_storeu_pd(p00 + i + hn, _mm_setzero_pd());
		_mm_storeu_pd(p01 + i, g01_re);
		_mm_storeu_pd(p01 + i + hn, g01_im);
		_mm_storeu_pd(p10 + i, g11_re);
		_mm_storeu_pd(p10 + i + hn, _mm_setzero_pd());
	}
#elif FNDSA_NEON
	float64_t *p00 = (float64_t *)b00;
	float64_t *p01 = (float64_t *)b01;
	float64_t *p10 = (float64_t *)b10;
	const float64_t *p11 = (const float64_t *)b11;
	float64x2_t zero = vdupq_lane_f64(vcreate_f64(0), 0);
	for (size_t i = 0; i < hn; i += 2) {
		float64x2_t b00_re = vld1q_f64(p00 + i);
		float64x2_t b00_im = vld1q_f64(p00 + i + hn);
		float64x2_t b01_re = vld1q_f64(p01 + i);
		float64x2_t b01_im = vld1q_f64(p01 + i + hn);
		float64x2_t b10_re = vld1q_f64(p10 + i);
		float64x2_t b10_im = vld1q_f64(p10 + i + hn);
		float64x2_t b11_re = vld1q_f64(p11 + i);
		float64x2_t b11_im = vld1q_f64(p11 + i + hn);

		/* g00 = b00*adj(b00) + b01*adj(b01) */
		float64x2_t g00_re = vaddq_f64(
			vaddq_f64(
				vmulq_f64(b00_re, b00_re),
				vmulq_f64(b00_im, b00_im)),
			vaddq_f64(
				vmulq_f64(b01_re, b01_re),
				vmulq_f64(b01_im, b01_im)));
		/* g01 = b00*adj(b10) + b01*adj(b11) */
		float64x2_t u_re = vaddq_f64(
			vmulq_f64(b00_re, b10_re),
			vmulq_f64(b00_im, b10_im));
		float64x2_t u_im = vsubq_f64(
			vmulq_f64(b00_im, b10_re),
			vmulq_f64(b00_re, b10_im));
		float64x2_t v_re = vaddq_f64(
			vmulq_f64(b01_re, b11_re),
			vmulq_f64(b01_im, b11_im));
		float64x2_t v_im = vsubq_f64(
			vmulq_f64(b01_im, b11_re),
			vmulq_f64(b01_re, b11_im));
		float64x2_t g01_re = vaddq_f64(u_re, v_re);
		float64x2_t g01_im = vaddq_f64(u_im, v_im);
		/* g11 = b10*adj(b10) + b11*adj(b11) */
		float64x2_t g11_re = vaddq_f64(
			vaddq_f64(
				vmulq_f64(b10_re, b10_re),
				vmulq_f64(b10_im, b10_im)),
			vaddq_f64(
				vmulq_f64(b11_re, b11_re),
				vmulq_f64(b11_im, b11_im)));

		vst1q_f64(p00 + i, g00_re);
		vst1q_f64(p00 + i + hn, zero);
		vst1q_f64(p01 + i, g01_re);
		vst1q_f64(p01 + i + hn, g01_im);
		vst1q_f64(p10 + i, g11_re);
		vst1q_f64(p10 + i + hn, zero);
	}
#elif FNDSA_RV64D
	f64 *bb00 = (f64 *)b00;
	f64 *bb01 = (f64 *)b01;
	f64 *bb10 = (f64 *)b10;
	const f64 *bb11 = (const f64 *)b11;
	for (size_t i = 0; i < hn; i ++) {
		f64 b00_re = bb00[i], b00_im = bb00[i + hn];
		f64 b01_re = bb01[i], b01_im = bb01[i + hn];
		f64 b10_re = bb10[i], b10_im = bb10[i + hn];
		f64 b11_re = bb11[i], b11_im = bb11[i + hn];

		/* g00 = b00*adj(b00) + b01*adj(b01) */
		f64 g00_re = f64_add(
			f64_add(f64_sqr(b00_re), f64_sqr(b00_im)),
			f64_add(f64_sqr(b01_re), f64_sqr(b01_im)));
		/* g01 = b00*adj(b10) + b01*adj(b11) */
		f64 u_re = f64_add(
			f64_mul(b00_re, b10_re),
			f64_mul(b00_im, b10_im));
		f64 u_im = f64_sub(
			f64_mul(b00_im, b10_re),
			f64_mul(b00_re, b10_im));
		f64 v_re = f64_add(
			f64_mul(b01_re, b11_re),
			f64_mul(b01_im, b11_im));
		f64 v_im = f64_sub(
			f64_mul(b01_im, b11_re),
			f64_mul(b01_re, b11_im));
		f64 g01_re = f64_add(u_re, v_re);
		f64 g01_im = f64_add(u_im, v_im);
		/* g11 = b10*adj(b10) + b11*adj(b11) */
		f64 g11_re = f64_add(
			f64_add(f64_sqr(b10_re), f64_sqr(b10_im)),
			f64_add(f64_sqr(b11_re), f64_sqr(b11_im)));

		bb00[i] = g00_re;
		bb00[i + hn] = (f64){ 0.0 };
		bb01[i] = g01_re;
		bb01[i + hn] = g01_im;
		bb10[i] = g11_re;
		bb10[i + hn] = (f64){ 0.0 };
	}
#else
	for (size_t i = 0; i < hn; i ++) {
		fpr b00_re = b00[i], b00_im = b00[i + hn];
		fpr b01_re = b01[i], b01_im = b01[i + hn];
		fpr b10_re = b10[i], b10_im = b10[i + hn];
		fpr b11_re = b11[i], b11_im = b11[i + hn];

		/* g00 = b00*adj(b00) + b01*adj(b01) */
		fpr g00_re = fpr_add(
			fpr_add(fpr_sqr(b00_re), fpr_sqr(b00_im)),
			fpr_add(fpr_sqr(b01_re), fpr_sqr(b01_im)));
		/* g01 = b00*adj(b10) + b01*adj(b11) */
		fpr u_re, u_im, v_re, v_im;
		fpr nb10_im = fpr_neg(b10_im);
		fpr nb11_im = fpr_neg(b11_im);
		FPC_MUL_NOALIAS(u_re, u_im,
			b00_re, b00_im, b10_re, nb10_im);
		FPC_MUL_NOALIAS(v_re, v_im,
			b01_re, b01_im, b11_re, nb11_im);
		fpr g01_re = fpr_add(u_re, v_re);
		fpr g01_im = fpr_add(u_im, v_im);
		/* g11 = b10*adj(b10) + b11*adj(b11) */
		fpr g11_re = fpr_add(
			fpr_add(fpr_sqr(b10_re), fpr_sqr(b10_im)),
			fpr_add(fpr_sqr(b11_re), fpr_sqr(b11_im)));

		b00[i] = g00_re;
		b00[i + hn] = FPR_ZERO;
		b01[i] = g01_re;
		b01[i + hn] = g01_im;
		b10[i] = g11_re;
		b10[i + hn] = FPR_ZERO;
	}
#endif
}

/* INV_Q and MINUS_INV_Q (1/q and -1/q) come from the representation. */

/* see sign_inner.h */
TARGET_SSE2 TARGET_NEON
void
fpoly_apply_basis(unsigned logn, fpr *t0, fpr *t1,
	fpr *b01, fpr *b11, const uint16_t *hm)
{
	size_t n = (size_t)1 << logn;
#if FNDSA_FFT_NEW
	/* hm is integral: use the integer-source FFT. */
	fpoly_FFT_uint16(logn, t0, hm);
#else
#if FNDSA_SSE2
	size_t i = 0;
	double *dd = (double *)t0;
	for (; (i + 4) <= n; i += 4) {
		__m128i x = _mm_setr_epi32(
			hm[i], hm[i + 1], hm[i + 2], hm[i + 3]);
		__m128d y0 = _mm_cvtepi32_pd(x);
		__m128d y1 = _mm_cvtepi32_pd(_mm_bsrli_si128(x, 8));
		_mm_storeu_pd(dd + i, y0);
		_mm_storeu_pd(dd + i + 2, y1);
	}
	__m128d z = _mm_setzero_pd();
	for (; i < n; i ++) {
		_mm_store_sd(dd + i, _mm_cvtsi32_sd(z, hm[i]));
	}
#elif FNDSA_NEON
	size_t i = 0;
	float64_t *dd = (float64_t *)t0;
	for (; (i + 4) <= n; i += 4) {
		uint16x4_t m16 = vld1_u16(hm + i);
		uint32x4_t m32 = vmovl_u16(m16);
		float64x2_t y0 = vcvtq_f64_u64(vmovl_u32(vget_low_u32(m32)));
		float64x2_t y1 = vcvtq_f64_u64(vmovl_high_u32(m32));
		vst1q_f64(dd + i, y0);
		vst1q_f64(dd + i + 2, y1);
	}
	for (; i < n; i ++) {
		vst1_f64(dd + i, vcvt_f64_u64(vcreate_u64(hm[i])));
	}
#elif FNDSA_RV64D
	f64 *tt0 = (f64 *)t0;
	for (size_t i = 0; i < n; i ++) {
		tt0[i] = f64_of(hm[i]);
	}
#else
	for (size_t i = 0; i < n; i ++) {
		t0[i] = fpr_of(hm[i]);
	}
#endif
	fpoly_FFT(logn, t0);
#endif
#if FNDSA_FXP3264 && !FNDSA_SIGN_SPLIT
	/* Multiply by 1/q first: FFT(hm)*b_hat would not fit in 32 integer
	   bits. */
#if FNDSA_SIGN_INVQ_WIDE
	fpoly_mul_inv_q(logn, t0);
#else
	fpoly_mulconst(logn, t0, INV_Q);
#endif
	fpoly_mul_fft(logn, b01, t0);
	fpoly_mul_fft(logn, t0, b11);
	if (t1 != b01) {
		memmove(t1, b01, n * sizeof(fpr));
	}
	fpoly_neg(logn, t1);
#else
	fpoly_mul_fft(logn, b01, t0);
	fpoly_mul_fft(logn, t0, b11);
	if (t1 != b01) {
		memmove(t1, b01, n * sizeof(fpr));
	}
	fpoly_mulconst(logn, t1, MINUS_INV_Q);
	fpoly_mulconst(logn, t0, INV_Q);
#endif
}
