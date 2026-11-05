#ifndef FNDSA_FXP_FUSED_H__
#define FNDSA_FXP_FUSED_H__

/*
 * fxp32.64 fused operations: the fused reciprocal (FNDSA_SIGN_RECIP_FUSED),
 * L10 in 64 bits (FNDSA_SIGN_L10_64), the fused sqrt * 1/sigma
 * (FNDSA_SIGN_SQRT_FUSED), the single-rounding complex multiply
 * (FNDSA_FXP_CMUL_FUSED) and the wide 1/q (FNDSA_SIGN_INVQ_WIDE). With
 * RECIP_FUSED, L10_64 or SQRT_FUSED at 0, sign_inner.h supplies the
 * composed version.
 *
 * The cm4 reciprocal and sqrt leaves assume the gamma checks: every LDL
 * diagonal d is at least q/alpha^2 > 4 (alpha <= 55), and d < 2^31.
 */

#include "fpr_fxp.h"

/* ==================================================================== */
/* Reciprocal kernel: 1/x as an unrounded 128-bit mantissa m and shift  */
/* k, for x in [1, 2^31), so k is in [1, 31]. Also used by fpr_inv.     */

#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
/* The cm4 leaf returns k in r12. */
#define FNDSA_FXP96_RECIP_RAW(mm, kk, xw)   do { \
		register uint32_t *frr_m __asm__("r0") = (mm); \
		register const uint32_t *frr_x __asm__("r1") = (xw); \
		register uint32_t frr_k __asm__("r12"); \
		__asm__ volatile ( \
			"@ fndsa_fxp96_recip_fast_raw_x %0 %1\n\t" \
			"bl\tfndsa_fxp96_recip_fast_raw_x" \
			: "+r" (frr_m), "+r" (frr_x), "=r" (frr_k) \
			: \
			: "r2", "r3", "r4", "r5", "r6", "r7", "r8", "r9", \
			  "r10", "r11", "r14", "cc", "memory"); \
		(kk) = frr_k; \
	} while (0)
#else
void fndsa_fxp96_recip_x(uint32_t *m, uint32_t *k, const uint32_t *x);
#define FNDSA_FXP96_RECIP_RAW(mm, kk, xw) \
	fndsa_fxp96_recip_x((mm), &(kk), (xw))
#endif

/* ==================================================================== */
/* Fused reciprocal: FPR_RECIP_X keeps 1/x unrounded, FPR_MUL_RECIP     */
/* multiplies by it and rounds once.                                    */

#if FNDSA_SIGN_RECIP_FUSED

typedef struct {
	uint32_t m[4];
	uint32_t k;
} fpr_rcp;

#define FPR_RECIP_X(rcp, x) \
	FNDSA_FXP96_RECIP_RAW((rcp).m, (rcp).k, (x).w)

#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
#if !FNDSA_KGEN_GAMMA
#error "CM4 fused reciprocal multiplication requires FNDSA_KGEN_GAMMA=1"
#endif
void fndsa_fxp96_mul_recip_x(uint32_t *z, const uint32_t *a,
	const uint32_t *m, uint32_t k);
uint64_t fndsa_fxp96_mul_recip_p63_x(const uint32_t *a,
	const uint32_t *m, uint32_t k);
/* Needs d >= 4, so that k <= 29. */
#define FPR_MUL_RECIP(d, a, rcp)   do { \
		register uint32_t *fmr_d __asm__("r0") = (d).w; \
		register const uint32_t *fmr_a __asm__("r1") = (a).w; \
		register const uint32_t *fmr_m __asm__("r2") = (rcp).m; \
		register uint32_t fmr_k __asm__("r3") = (rcp).k; \
		__asm__ volatile ( \
			"@ fndsa_fxp96_mul_recip_x %0 %1 %2 %3\n\t" \
			"bl\tfndsa_fxp96_mul_recip_x" \
			: "+r" (fmr_d), "+r" (fmr_a), "+r" (fmr_m), "+r" (fmr_k) \
			: \
			: "r4", "r5", "r6", "r7", "r8", "r9", \
			  "r10", "r11", "r12", "r14", \
			  "s9", "s10", "s11", "s12", "s13", "s14", "s15", \
			  "cc", "memory"); \
	} while (0)
/* Same, rounded to a Q0.63 uint64_t. */
#define FPR_MUL_RECIP_P63(d63, a, rcp)   do { \
		register uint32_t fmr_r0 __asm__("r0") = \
			(uint32_t)(uintptr_t)(a).w; \
		register uint32_t fmr_r1 __asm__("r1") = \
			(uint32_t)(uintptr_t)(rcp).m; \
		register uint32_t fmr_r2 __asm__("r2") = (rcp).k; \
		__asm__ volatile ( \
			"@ fndsa_fxp96_mul_recip_p63_x %0 %1 %2\n\t" \
			"bl	fndsa_fxp96_mul_recip_p63_x" \
			: "+r" (fmr_r0), "+r" (fmr_r1), "+r" (fmr_r2) \
			: \
			: "r3", "r4", "r5", "r6", "r7", "r8", "r9", \
			  "r10", "r11", "r12", "r14", \
			  "s10", "s11", "s12", "s15", \
			  "cc", "memory"); \
		(d63) = ((uint64_t)fmr_r1 << 32) | (uint64_t)fmr_r0; \
	} while (0)
#else
void fndsa_fxp96_mul_recip(uint32_t *z, const uint32_t *a,
	const uint32_t *m, uint32_t k);
uint64_t fndsa_fxp96_mul_recip_p63(const uint32_t *a,
	const uint32_t *m, uint32_t k);
#define FPR_MUL_RECIP(d, a, rcp)   do { \
		fndsa_fxp96_mul_recip((d).w, (a).w, (rcp).m, (rcp).k); \
	} while (0)
/* Same, rounded to a Q0.63 uint64_t. */
#define FPR_MUL_RECIP_P63(d63, a, rcp) \
	((d63) = fndsa_fxp96_mul_recip_p63((a).w, (rcp).m, (rcp).k))
#endif /* FNDSA_FPR_IMPL */

#endif /* FNDSA_SIGN_RECIP_FUSED */

/* ==================================================================== */
/* L10 in 64 bits: the non-root l10 is stored as Q0.63.                 */

#if FNDSA_SIGN_L10_64

#if !FNDSA_SIGN_RECIP_FUSED
#error FNDSA_SIGN_L10_64=1 requires FNDSA_SIGN_RECIP_FUSED=1
#endif
/* Without Symplecroot or Symplectree, the root l10 (which has an integer
   part) would be truncated to Q0.63. */
#if !(FNDSA_SYMP_ROOT || FNDSA_SYMP_TREE)
#error FNDSA_SIGN_L10_64=1 requires FNDSA_SYMP_ROOT or FNDSA_SYMP_TREE
#endif

/* Q0.63, loaded as 8 bytes, so the l10 band must be 8-aligned. */
typedef uint64_t l10_t;

/* Index of the imaginary half within the l10 band. */
#define L10_IM_OFF(logn)   ((size_t)1 << ((logn) - 1))

/* Exact widening: Q0.63 lies on the 2^-64 grid. */
static inline fpr
l10_to_fpr(l10_t r)
{
	fpr d;

	d.w[0] = (uint32_t)(r << 1);
	d.w[1] = (uint32_t)(r >> 31);
	d.w[2] = (uint32_t)((int64_t)r >> 63);
	return d;
}

/* s = conj(g/d) in Q0.63, w = g/d widened. */
#define FPR_L10_DIV(s_re, s_im, w_re, w_im, g_re, g_im, rcp)  do { \
		FPR_MUL_RECIP_P63((s_re), (g_re), (rcp)); \
		FPR_MUL_RECIP_P63((s_im), (g_im), (rcp)); \
		(w_re) = l10_to_fpr(s_re); \
		(w_im) = l10_to_fpr(s_im); \
		(s_im) = 0 - (s_im); \
	} while (0)

#define fpoly_mul_l10_fft   fndsa_fpoly_mul_l10_fft
#define fpoly_l10_neg       fndsa_fpoly_l10_neg
void fpoly_mul_l10_fft(unsigned logn, fpr *a, const l10_t *b);
void fpoly_l10_neg(unsigned logn, l10_t *a);

#endif /* FNDSA_SIGN_L10_64 */

/* ==================================================================== */
/* Fused sqrt: sqrt(x) * INV_SIGMA[logn] with one rounding, using the   */
/* wider INV_SIGMA_Q71 copy of the constant.                            */

#if FNDSA_SIGN_SQRT_FUSED

void fndsa_fxp96_sqrt_mul_inv_sigma(uint32_t *d, const uint32_t *x,
	unsigned logn);

/* The cm4 leaf needs the gamma bounds; the C version takes any input. */
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4 && FNDSA_KGEN_GAMMA
void fndsa_fxp96_sqrt_mul_inv_sigma_x(uint32_t *d, const uint32_t *x,
	unsigned logn);
#define FPR_SQRT_MUL_INV_SIGMA(d, x, logn)   do { \
		register uint32_t *fsm_d __asm__("r0") = (d).w; \
		register const uint32_t *fsm_x __asm__("r1") = (x).w; \
		register uint32_t fsm_logn __asm__("r2") = (uint32_t)(logn); \
		__asm__ volatile ( \
			"@ fndsa_fxp96_sqrt_mul_inv_sigma_x %0 %1 %2\n\t" \
			"bl	fndsa_fxp96_sqrt_mul_inv_sigma_x" \
			: "+r" (fsm_d), "+r" (fsm_x), "+r" (fsm_logn) \
			: \
			: "r3", "r4", "r5", "r6", "r7", "r8", "r9", \
			  "r10", "r11", "r12", "r14", \
			  "s9", "s10", "s11", "s12", "s13", "s14", "s15", \
			  "cc", "memory"); \
	} while (0)
#else
#define FPR_SQRT_MUL_INV_SIGMA(d, x, logn) \
	fndsa_fxp96_sqrt_mul_inv_sigma((d).w, (x).w, (logn))
#endif

#endif /* FNDSA_SIGN_SQRT_FUSED */

/* ==================================================================== */
/* Complex multiply: each component is summed exactly and rounded once, */
/* or composed from rounded products with FNDSA_FXP_CMUL_FUSED=0.       */

#if FNDSA_FXP_CMUL_FUSED
/* Portable C version; the cm4 build keeps it as the fprtest reference. */
void fndsa_fxp96_cmul_rne64(
	uint32_t *out, const uint32_t *a, const uint32_t *b);
#endif

#if FNDSA_FXP_CMUL_FUSED && FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4

void fndsa_fxp96_cmul_self_rne64_x(
	uint32_t *a_re, uint32_t *a_im,
	const uint32_t *b_re, const uint32_t *b_im);
void fndsa_fxp96_cmul_noalias_rne64_x(
	uint32_t *d_re, uint32_t *d_im,
	const uint32_t *a_re, const uint32_t *a_im,
	const uint32_t *b_re, const uint32_t *b_im);
#define FPC_MUL(d_re, d_im, a_re, a_im, b_re, b_im)   do { \
		fpr fpct_a_re = (a_re); \
		fpr fpct_a_im = (a_im); \
		fpr fpct_b_re = (b_re); \
		fpr fpct_b_im = (b_im); \
		fpr fpct_out_re, fpct_out_im; \
		FPC_MUL_NOALIAS(fpct_out_re, fpct_out_im, \
			fpct_a_re, fpct_a_im, fpct_b_re, fpct_b_im); \
		(d_re) = fpct_out_re; \
		(d_im) = fpct_out_im; \
	} while (0)
#define FPC_MUL_SELF(d_re, d_im, b_re, b_im)   do { \
		register uint32_t *fpct_re __asm__("r0") = (d_re).w; \
		register uint32_t *fpct_im __asm__("r1") = (d_im).w; \
		register const uint32_t *fpct_bre __asm__("r2") = (b_re).w; \
		register const uint32_t *fpct_bim __asm__("r3") = (b_im).w; \
		__asm__ volatile ( \
			"@ fndsa_fxp96_cmul_self_rne64_x %0 %1 %2 %3\n\t" \
			"bl	fndsa_fxp96_cmul_self_rne64_x" \
			: "+r" (fpct_re), "+r" (fpct_im), \
			  "+r" (fpct_bre), "+r" (fpct_bim) \
			: \
			: "r4", "r5", "r6", "r7", "r8", "r9", \
			  "r10", "r11", "r12", "r14", \
			  "s9", "s11", "s12", \
			  "cc", "memory"); \
	} while (0)
#define FPC_SQR_SELF(d_re, d_im)   do { \
		fpr fpct_sqr_re = (d_re); \
		fpr fpct_sqr_im = (d_im); \
		FPC_MUL_SELF((d_re), (d_im), fpct_sqr_re, fpct_sqr_im); \
	} while (0)
#define FPC_MUL_NOALIAS(d_re, d_im, a_re, a_im, b_re, b_im)   do { \
		register uint32_t *fpct_re __asm__("r0") = (d_re).w; \
		register uint32_t *fpct_im __asm__("r1") = (d_im).w; \
		register const uint32_t *fpct_are __asm__("r2") = (a_re).w; \
		register const uint32_t *fpct_aim __asm__("r3") = (a_im).w; \
		register const uint32_t *fpct_bre __asm__("r4") = (b_re).w; \
		register const uint32_t *fpct_bim __asm__("r5") = (b_im).w; \
		__asm__ volatile ( \
			"@ fndsa_fxp96_cmul_noalias_rne64_x %0 %1 %2 %3 %4 %5\n\t" \
			"bl	fndsa_fxp96_cmul_noalias_rne64_x" \
			: "+r" (fpct_re), "+r" (fpct_im), \
			  "+r" (fpct_are), "+r" (fpct_aim), \
			  "+r" (fpct_bre), "+r" (fpct_bim) \
			: \
			: "r6", "r7", "r8", "r9", "r10", "r11", \
			  "r12", "r14", "s8", "s9", "s10", "s12", "s13", \
			  "cc", "memory"); \
	} while (0)

#elif FNDSA_FXP_CMUL_FUSED

#define FPC_MUL(d_re, d_im, a_re, a_im, b_re, b_im)   do { \
		fxp_fpc fpct_a, fpct_b, fpct_d; \
		fpct_a.z.re = (a_re); \
		fpct_a.z.im = (a_im); \
		fpct_b.z.re = (b_re); \
		fpct_b.z.im = (b_im); \
		fndsa_fxp96_cmul_rne64(fpct_d.w, fpct_a.w, fpct_b.w); \
		(d_re) = fpct_d.z.re; \
		(d_im) = fpct_d.z.im; \
	} while (0)
#define FPC_MUL_SELF(d_re, d_im, b_re, b_im) \
	FPC_MUL((d_re), (d_im), (d_re), (d_im), (b_re), (b_im))
#define FPC_SQR_SELF(d_re, d_im) \
	FPC_MUL((d_re), (d_im), (d_re), (d_im), (d_re), (d_im))
#define FPC_MUL_NOALIAS(d_re, d_im, a_re, a_im, b_re, b_im) \
	FPC_MUL((d_re), (d_im), (a_re), (a_im), (b_re), (b_im))

#else

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

#endif /* FNDSA_FXP_CMUL_FUSED */

/* ==================================================================== */
/* Wide 1/q.                                                            */

/* x/q as x * INV_Q_Q77 / 2^77 with one rounding. INV_Q on the 2^-64
   grid has too few significant bits for the large FFT-domain target.
   Used when FNDSA_SIGN_INVQ_WIDE=1. */
#define INV_Q_Q77   (((uint64_t)INV_Q_Q77_HI << 32) | INV_Q_Q77_LO)
#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
void fndsa_fpr_mul_inv_q_x(void);
static inline fpr
fndsa_fpr_mul_inv_q_inline(fpr x)
{
	register uint32_t r0 __asm__("r0") = x.w[0];
	register uint32_t r1 __asm__("r1") = x.w[1];
	register uint32_t r2 __asm__("r2") = x.w[2];

	__asm__ volatile (
		"@ fndsa_fpr_mul_inv_q_x %0 %1 %2\n\t"
		"bl	fndsa_fpr_mul_inv_q_x"
		: "+r" (r0), "+r" (r1), "+r" (r2)
		:
		: "r3", "r4", "r5", "r6", "r7", "r14", "cc");

	return (fpr){{ r0, r1, r2 }};
}
#define FPR_MUL_INV_Q(x)   fndsa_fpr_mul_inv_q_inline(x)
#else
fpr fndsa_fpr_mul_inv_q(fpr x);
#define FPR_MUL_INV_Q(x)   fndsa_fpr_mul_inv_q(x)
#endif

#endif
