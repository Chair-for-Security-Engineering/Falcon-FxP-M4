#ifndef FNDSA_FXP_FFT_H__
#define FNDSA_FXP_FFT_H__

/*
 * fxp32.64 FFT support: the unsigned 64-bit twiddle table
 * (FNDSA_TWIDDLE_U64), its multiplies, and the calls into the assembly
 * butterflies of sign_fpoly.S (FNDSA_FFT_ASM, which needs the unsigned
 * table and the cm4 implementation).
 */

#include "fpr_fxp.h"

/* Assembly butterflies (sign_fpoly.S), cm4 impl only. */
#define FNDSA_FFT_ASM_CM4 \
	(FNDSA_FFT_ASM && FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4)
#if !FNDSA_TWIDDLE_U64

/* Plain fpr table; sign_inner.h supplies the generic FPOLY_GM_* macros. */
extern const fpr GM[];

#else

/* Unsigned Q0.64 twiddle magnitudes; cos is negative for odd pair
   indices. The FPOLY_GM_* macros below absorb that sign by choosing
   the kernel and swapping ADD_SUB operands. The parity is public. */
extern const uint64_t GM[];

/* ==================================================================== */
/* Twiddle multiplies (and the compact-FFT batch butterflies on cm4).   */

#if FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4
void fndsa_fxp96_cmul_gm_x(uint32_t *out, const uint32_t *a,
	const uint64_t *gm);
void fndsa_fxp96_cmul_gmj_x(uint32_t *out, const uint32_t *a,
	const uint64_t *gm);
void fndsa_fpr_mul_gm_x(void);
/* Multiply by integer + frac/2^64, rounded once. */
void fndsa_fpr_mul_rnj_x(void);
/* Batched Q0.64 compact-FFT butterflies. */
void fndsa_fpr_bfly_gm_batch_x(void);
void fndsa_fpr_bfly_gm_batch_half_x(void);
/* Same with a coefficient integer + frac/2^64. */
void fndsa_fpr_bfly_rnj_batch_x(void);
void fndsa_fpr_bfly_rnj_batch_half_x(void);

static inline fpr
fndsa_fpr_mul_gm_inline(fpr x, uint64_t c)
{
	uint32_t x0 = x.w[0];
	uint32_t x1 = x.w[1];
	uint32_t x2 = x.w[2];
	uint32_t c0 = (uint32_t)c;
	uint32_t c1 = (uint32_t)(c >> 32);
	register uint32_t r0 __asm__("r0") = x0;
	register uint32_t r1 __asm__("r1") = x1;
	register uint32_t r2 __asm__("r2") = x2;
	register uint32_t r3 __asm__("r3") = c0;
	register uint32_t r4 __asm__("r4") = c1;

	__asm__ volatile (
		"@ fndsa_fpr_mul_gm_x %0 %1 %2 %3 %4\n\t"
		"bl	fndsa_fpr_mul_gm_x"
		: "+r" (r0), "+r" (r1), "+r" (r2), "+r" (r3), "+r" (r4)
		:
		: "r5", "r6", "r7", "r14", "cc");

	return (fpr){{ r0, r1, r2 }};
}

#define FNDSA_FPR_MUL_GM(x, c)   fndsa_fpr_mul_gm_inline((x), (c))

static inline void
fndsa_fpr_bfly_gm_batch_inline(fpr *a, fpr *b, size_t num,
	uint64_t frac, int neg)
{
	register fpr *r0 __asm__("r0") = a;
	register fpr *r1 __asm__("r1") = b;
	register uint32_t r2 __asm__("r2") = (uint32_t)frac;
	register uint32_t r3 __asm__("r3") = (uint32_t)(frac >> 32);
	register uint32_t r4 __asm__("r4") = (uint32_t)num;
	register uint32_t r5 __asm__("r5") = (uint32_t)neg;

	__asm__ volatile (
		"@ fndsa_fpr_bfly_gm_batch_x %0 %1 %2 %3 %4 %5\n\t"
		"bl\tfndsa_fpr_bfly_gm_batch_x"
		: "+r" (r0), "+r" (r1), "+r" (r2), "+r" (r3), "+r" (r4),
		  "+r" (r5)
		:
		: "r6", "r7", "r8", "r9", "r10", "r11", "r12", "r14",
		  "s3", "s4", "s5", "cc", "memory");
}

#define FNDSA_FPR_BFLY_GM_BATCH(a, b, num, frac, neg) \
	fndsa_fpr_bfly_gm_batch_inline((a), (b), (num), (frac), (neg))

static inline void
fndsa_fpr_bfly_gm_batch_half_inline(fpr *a, fpr *b, size_t num,
	uint64_t frac, int neg)
{
	register fpr *r0 __asm__("r0") = a;
	register fpr *r1 __asm__("r1") = b;
	register uint32_t r2 __asm__("r2") = (uint32_t)frac;
	register uint32_t r3 __asm__("r3") = (uint32_t)(frac >> 32);
	register uint32_t r4 __asm__("r4") = (uint32_t)num;
	register uint32_t r5 __asm__("r5") = (uint32_t)neg;

	__asm__ volatile (
		"@ fndsa_fpr_bfly_gm_batch_half_x %0 %1 %2 %3 %4 %5\n\t"
		"bl\tfndsa_fpr_bfly_gm_batch_half_x"
		: "+r" (r0), "+r" (r1), "+r" (r2), "+r" (r3), "+r" (r4),
		  "+r" (r5)
		:
		: "r6", "r7", "r8", "r9", "r10", "r11", "r12", "r14",
		  "s3", "s4", "s5", "cc", "memory");
}

#define FNDSA_FPR_BFLY_GM_BATCH_HALF(a, b, num, frac, neg) \
	fndsa_fpr_bfly_gm_batch_half_inline((a), (b), (num), (frac), (neg))

static inline void
fndsa_fpr_bfly_rnj_batch_inline(fpr *a, fpr *b, size_t num,
	uint64_t frac, uint32_t integer, int neg)
{
	register fpr *r0 __asm__("r0") = a;
	register fpr *r1 __asm__("r1") = b;
	register uint32_t r2 __asm__("r2") = (uint32_t)frac;
	register uint32_t r3 __asm__("r3") = (uint32_t)(frac >> 32);
	register uint32_t r4 __asm__("r4") = integer;
	register uint32_t r5 __asm__("r5") = (uint32_t)num;
	register uint32_t r6 __asm__("r6") = (uint32_t)neg;

	__asm__ volatile (
		"@ fndsa_fpr_bfly_rnj_batch_x %0 %1 %2 %3 %4 %5 %6\n\t"
		"bl\tfndsa_fpr_bfly_rnj_batch_x"
		: "+r" (r0), "+r" (r1), "+r" (r2), "+r" (r3), "+r" (r4),
		  "+r" (r5), "+r" (r6)
		:
		: "r7", "r8", "r9", "r10", "r11", "r12", "r14",
		  "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11",
		  "cc", "memory");
}

#define FNDSA_FPR_BFLY_RNJ_BATCH(a, b, num, frac, integer, neg) \
	fndsa_fpr_bfly_rnj_batch_inline((a), (b), (num), (frac), \
		(integer), (neg))

static inline void
fndsa_fpr_bfly_rnj_batch_half_inline(fpr *a, fpr *b, size_t num,
	uint64_t frac, uint32_t integer, int neg)
{
	register fpr *r0 __asm__("r0") = a;
	register fpr *r1 __asm__("r1") = b;
	register uint32_t r2 __asm__("r2") = (uint32_t)frac;
	register uint32_t r3 __asm__("r3") = (uint32_t)(frac >> 32);
	register uint32_t r4 __asm__("r4") = integer;
	register uint32_t r5 __asm__("r5") = (uint32_t)num;
	register uint32_t r6 __asm__("r6") = (uint32_t)neg;

	__asm__ volatile (
		"@ fndsa_fpr_bfly_rnj_batch_half_x %0 %1 %2 %3 %4 %5 %6\n\t"
		"bl\tfndsa_fpr_bfly_rnj_batch_half_x"
		: "+r" (r0), "+r" (r1), "+r" (r2), "+r" (r3), "+r" (r4),
		  "+r" (r5), "+r" (r6)
		:
		: "r7", "r8", "r9", "r10", "r11", "r12", "r14",
		  "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11",
		  "cc", "memory");
}

#define FNDSA_FPR_BFLY_RNJ_BATCH_HALF(a, b, num, frac, integer, neg) \
	fndsa_fpr_bfly_rnj_batch_half_inline((a), (b), (num), (frac), \
		(integer), (neg))

#define FNDSA_FXP96_CMUL_GM(out, a, gm)   do { \
		register uint32_t *fgm_out __asm__("r0") = (out); \
		register const uint32_t *fgm_a __asm__("r1") = (a); \
		register const uint64_t *fgm_gm __asm__("r2") = (gm); \
		__asm__ volatile ( \
			"@ fndsa_fxp96_cmul_gm_x %0 %1 %2\n\t" \
			"bl	fndsa_fxp96_cmul_gm_x" \
			: "+r" (fgm_out), "+r" (fgm_a), "+r" (fgm_gm) \
			: \
			: "r3", "r4", "r5", "r6", "r7", "r8", "r9", \
			  "r10", "r11", "r12", "r14", "s8", "s9", \
			  "cc", "memory"); \
	} while (0)
#define FNDSA_FXP96_CMUL_GMJ(out, a, gm)   do { \
		register uint32_t *fgm_out __asm__("r0") = (out); \
		register const uint32_t *fgm_a __asm__("r1") = (a); \
		register const uint64_t *fgm_gm __asm__("r2") = (gm); \
		__asm__ volatile ( \
			"@ fndsa_fxp96_cmul_gmj_x %0 %1 %2\n\t" \
			"bl	fndsa_fxp96_cmul_gmj_x" \
			: "+r" (fgm_out), "+r" (fgm_a), "+r" (fgm_gm) \
			: \
			: "r3", "r4", "r5", "r6", "r7", "r8", "r9", \
			  "r10", "r11", "r12", "r14", "s8", "s9", \
			  "cc", "memory"); \
	} while (0)
#else
/* Multiply by a GM pair A+iB, by A-iB, and by a single word. out may
   alias a. */
void fndsa_fxp96_cmul_gm(uint32_t *out, const uint32_t *a,
	const uint64_t *gm);
void fndsa_fxp96_cmul_gmj(uint32_t *out, const uint32_t *a,
	const uint64_t *gm);
fpr fndsa_fpr_mul_gm(fpr x, uint64_t c);
fpr fndsa_fpr_mul_rnj(fpr x, uint64_t frac, uint32_t integer);

#define FNDSA_FXP96_CMUL_GM(out, a, gm) \
	fndsa_fxp96_cmul_gm((out), (a), (gm))
#define FNDSA_FXP96_CMUL_GMJ(out, a, gm) \
	fndsa_fxp96_cmul_gmj((out), (a), (gm))
#define FNDSA_FPR_MUL_GM(x, c)   fndsa_fpr_mul_gm((x), (c))
#endif /* FNDSA_FPR_IMPL */

/* ==================================================================== */
/* FFT butterflies. The twiddle state is the GM pair and its parity.    */

#define FPOLY_GM_TW(s, k) \
	const uint64_t *s ## _gm = &GM[(k) << 1]; \
	unsigned s ## _par = (k) & 1

#define FPOLY_GM_TW_CONJ(s, k) FPOLY_GM_TW(s, k)

#if FNDSA_FFT_ASM_CM4
/* Custom ABI: r0-r3 = A, B, O1, O2, r4 = twiddle pair, r5/r6 = byte
   offset to the imaginary half (in/out). */
void fndsa_fxp96_bfly_fwd_even_x(void);
void fndsa_fxp96_bfly_fwd_odd_x(void);
void fndsa_fxp96_bfly_fwd_even_ip_x(void);
void fndsa_fxp96_bfly_fwd_odd_ip_x(void);
void fndsa_fxp96_bfly_inv_even_ip_x(void);
void fndsa_fxp96_bfly_inv_odd_ip_x(void);
void fndsa_fxp96_bfly_split_even_x(void);
void fndsa_fxp96_bfly_split_odd_x(void);

#define FNDSA_FXP96_BFLY_CALL(fn, A, B, O1, O2, gm, off_in, off_out)   do { \
		register const uint32_t *fb_a __asm__("r0") = \
			(const uint32_t *)(A); \
		register const uint32_t *fb_b __asm__("r1") = \
			(const uint32_t *)(B); \
		register uint32_t *fb_o1 __asm__("r2") = (uint32_t *)(O1); \
		register uint32_t *fb_o2 __asm__("r3") = (uint32_t *)(O2); \
		register const uint64_t *fb_gm __asm__("r4") = (gm); \
		register uint32_t fb_oi __asm__("r5") = (uint32_t)(off_in); \
		register uint32_t fb_oo __asm__("r6") = (uint32_t)(off_out); \
		__asm__ volatile ( \
			"@ " #fn " %0 %1 %2 %3 %4 %5 %6\n\t" \
			"bl	" #fn \
			: "+r" (fb_a), "+r" (fb_b), "+r" (fb_o1), "+r" (fb_o2), \
			  "+r" (fb_gm), "+r" (fb_oi), "+r" (fb_oo) \
			: \
			: "r7", "r8", "r9", "r10", "r11", "r12", "r14", \
			  "s5", "s6", "s7", "s8", "s9", "s10", "s11", "s12", \
			  "s13", "s14", "s15", "cc", "memory"); \
	} while (0)

/* In place: r0 = A (x, o1), r1 = B (y, o2), r2 = gm, r3 = byte offset
   to the imaginary halves. */
#define FNDSA_FXP96_BFLY_CALL_IP(fn, A, B, gm, off)   do { \
		register uint32_t *fb_a __asm__("r0") = (uint32_t *)(A); \
		register uint32_t *fb_b __asm__("r1") = (uint32_t *)(B); \
		register const uint64_t *fb_gm __asm__("r2") = (gm); \
		register uint32_t fb_off __asm__("r3") = (uint32_t)(off); \
		__asm__ volatile ( \
			"@ " #fn " %0 %1 %2 %3\n\t" \
			"bl	" #fn \
			: "+r" (fb_a), "+r" (fb_b), "+r" (fb_gm), "+r" (fb_off) \
			: \
			: "r4", "r5", "r6", "r7", "r8", "r9", "r10", "r11", \
			  "r12", "r14", \
			  "s5", "s6", "s7", "s8", "s9", "s10", "s11", "s12", \
			  "s13", "s14", "s15", "cc", "memory"); \
	} while (0)

/* x and y already sit in the o1/o2 slots, so only the addresses are
   passed. */
#define FPOLY_GM_BFLY_FWD(o1_re, o1_im, o2_re, o2_im, x_re, x_im, y_re, y_im, s)  do { \
		uint32_t fbp_off = (uint32_t)((char *)&(o1_im) - (char *)&(o1_re)); \
		(void)(x_re); (void)(x_im); (void)(y_re); (void)(y_im); \
		if (s ## _par) { \
			FNDSA_FXP96_BFLY_CALL_IP(fndsa_fxp96_bfly_fwd_odd_ip_x, \
				&(o1_re), &(o2_re), s ## _gm, fbp_off); \
		} else { \
			FNDSA_FXP96_BFLY_CALL_IP(fndsa_fxp96_bfly_fwd_even_ip_x, \
				&(o1_re), &(o2_re), s ## _gm, fbp_off); \
		} \
	} while (0)
#define FPOLY_GM_BFLY_INV(o1_re, o1_im, o2_re, o2_im, x_re, x_im, y_re, y_im, s)  do { \
		uint32_t fbp_off = (uint32_t)((char *)&(o1_im) - (char *)&(o1_re)); \
		(void)(x_re); (void)(x_im); (void)(y_re); (void)(y_im); \
		if (s ## _par) { \
			FNDSA_FXP96_BFLY_CALL_IP(fndsa_fxp96_bfly_inv_odd_ip_x, \
				&(o1_re), &(o2_re), s ## _gm, fbp_off); \
		} else { \
			FNDSA_FXP96_BFLY_CALL_IP(fndsa_fxp96_bfly_inv_even_ip_x, \
				&(o1_re), &(o2_re), s ## _gm, fbp_off); \
		} \
	} while (0)
#else
/* (o1, o2) = (x + y*s, x - y*s). Even: s = A+iB. Odd: s = -(A-iB), so
   use the conjugate kernel and swap the outputs. */
#define FPOLY_GM_BFLY_FWD(o1_re, o1_im, o2_re, o2_im, x_re, x_im, y_re, y_im, s)  do { \
		fxp_fpc fgm_p; \
		fgm_p.z.re = (y_re); \
		fgm_p.z.im = (y_im); \
		if (s ## _par) { \
			FNDSA_FXP96_CMUL_GMJ(fgm_p.w, fgm_p.w, s ## _gm); \
			FPR_ADD_SUB((o2_re), (o1_re), (x_re), fgm_p.z.re); \
			FPR_ADD_SUB((o2_im), (o1_im), (x_im), fgm_p.z.im); \
		} else { \
			FNDSA_FXP96_CMUL_GM(fgm_p.w, fgm_p.w, s ## _gm); \
			FPR_ADD_SUB((o1_re), (o2_re), (x_re), fgm_p.z.re); \
			FPR_ADD_SUB((o1_im), (o2_im), (x_im), fgm_p.z.im); \
		} \
	} while (0)
/* o1 = x + y, o2 = (x - y)*conj(s). Even: conj(s) = A-iB. Odd:
   conj(s) = -(A+iB), so swap the inputs to negate the difference. */
#define FPOLY_GM_BFLY_INV(o1_re, o1_im, o2_re, o2_im, x_re, x_im, y_re, y_im, s)  do { \
		fxp_fpc fgm_p; \
		if (s ## _par) { \
			FPR_ADD_SUB((o1_re), fgm_p.z.re, (y_re), (x_re)); \
			FPR_ADD_SUB((o1_im), fgm_p.z.im, (y_im), (x_im)); \
			FNDSA_FXP96_CMUL_GM(fgm_p.w, fgm_p.w, s ## _gm); \
		} else { \
			FPR_ADD_SUB((o1_re), fgm_p.z.re, (x_re), (y_re)); \
			FPR_ADD_SUB((o1_im), fgm_p.z.im, (x_im), (y_im)); \
			FNDSA_FXP96_CMUL_GMJ(fgm_p.w, fgm_p.w, s ## _gm); \
		} \
		(o2_re) = fgm_p.z.re; \
		(o2_im) = fgm_p.z.im; \
	} while (0)
#endif /* FNDSA_FFT_ASM_CM4 */

/* ==================================================================== */
/* Merge and split steps.                                               */

/* (o1, o2) = (a + b*GM[k], a - b*GM[k]), parity as above. */
#define FPOLY_GM_MERGE_BFLY(o1_re, o1_im, o2_re, o2_im, a_re, a_im, b_re, b_im, k)  do { \
		fxp_fpc fgm_p; \
		fgm_p.z.re = (b_re); \
		fgm_p.z.im = (b_im); \
		if ((k) & 1) { \
			FNDSA_FXP96_CMUL_GMJ(fgm_p.w, fgm_p.w, &GM[(k) << 1]); \
			FPR_ADD_SUB((o2_re), (o1_re), (a_re), fgm_p.z.re); \
			FPR_ADD_SUB((o2_im), (o1_im), (a_im), fgm_p.z.im); \
		} else { \
			FNDSA_FXP96_CMUL_GM(fgm_p.w, fgm_p.w, &GM[(k) << 1]); \
			FPR_ADD_SUB((o1_re), (o2_re), (a_re), fgm_p.z.re); \
			FPR_ADD_SUB((o1_im), (o2_im), (a_im), fgm_p.z.im); \
		} \
	} while (0)
/* o0 = (a + b)/2, o1 = ((a - b)/2)*conj(GM[k]), parity as in the
   inverse butterfly. */
#define FPOLY_GM_SPLIT_BFLY(o0_re, o0_im, o1_re, o1_im, a_re, a_im, b_re, b_im, k)  do { \
		fxp_fpc fgm_p; \
		fpr fgm_t_re, fgm_t_im; \
		if ((k) & 1) { \
			FPR_ADD_SUB(fgm_t_re, fgm_p.z.re, (b_re), (a_re)); \
			FPR_ADD_SUB(fgm_t_im, fgm_p.z.im, (b_im), (a_im)); \
			FNDSA_FXP96_CMUL_GM(fgm_p.w, fgm_p.w, &GM[(k) << 1]); \
		} else { \
			FPR_ADD_SUB(fgm_t_re, fgm_p.z.re, (a_re), (b_re)); \
			FPR_ADD_SUB(fgm_t_im, fgm_p.z.im, (a_im), (b_im)); \
			FNDSA_FXP96_CMUL_GMJ(fgm_p.w, fgm_p.w, &GM[(k) << 1]); \
		} \
		(o0_re) = fpr_half(fgm_t_re); \
		(o0_im) = fpr_half(fgm_t_im); \
		(o1_re) = fpr_half(fgm_p.z.re); \
		(o1_im) = fpr_half(fgm_p.z.im); \
	} while (0)
/* Same for a self-adjoint input: o0_im is zero and two scalar muls
   suffice. The even arm needs an explicit negation, as -u*B has no
   later add/sub to absorb the sign. */
#define FPOLY_GM_SPLIT_SELFADJ(o0_re, o0_im, o1_re, o1_im, a_re, b_re, k)  do { \
		fpr fgm_t_re, fgm_u_re; \
		if ((k) & 1) { \
			FPR_ADD_SUB(fgm_t_re, fgm_u_re, (b_re), (a_re)); \
			fgm_u_re = fpr_half(fgm_u_re); \
			(o1_re) = FNDSA_FPR_MUL_GM(fgm_u_re, GM[((k) << 1) + 0]); \
			(o1_im) = FNDSA_FPR_MUL_GM(fgm_u_re, GM[((k) << 1) + 1]); \
		} else { \
			FPR_ADD_SUB(fgm_t_re, fgm_u_re, (a_re), (b_re)); \
			fgm_u_re = fpr_half(fgm_u_re); \
			(o1_re) = FNDSA_FPR_MUL_GM(fgm_u_re, GM[((k) << 1) + 0]); \
			(o1_im) = fpr_neg(FNDSA_FPR_MUL_GM(fgm_u_re, \
				GM[((k) << 1) + 1])); \
		} \
		(o0_re) = fpr_half(fgm_t_re); \
		(o0_im) = FPR_ZERO; \
	} while (0)

#if FNDSA_FFT_ASM_CM4
/* Merge/split step i; strides in bytes. */
#define FPOLY_GM_MERGE(f, f0, f1, i, hn, qn, k)   do { \
		uint32_t fbp_oi = (uint32_t)((qn) * sizeof(fpr)); \
		uint32_t fbp_oo = (uint32_t)((hn) * sizeof(fpr)); \
		fpr *fbp_f = (f) + ((i) << 1); \
		const fpr *fbp_f0 = (f0) + (i), *fbp_f1 = (f1) + (i); \
		if ((k) & 1) { \
			FNDSA_FXP96_BFLY_CALL(fndsa_fxp96_bfly_fwd_odd_x, \
				fbp_f0, fbp_f1, fbp_f, fbp_f + 1, &GM[(k) << 1], \
				fbp_oi, fbp_oo); \
		} else { \
			FNDSA_FXP96_BFLY_CALL(fndsa_fxp96_bfly_fwd_even_x, \
				fbp_f0, fbp_f1, fbp_f, fbp_f + 1, &GM[(k) << 1], \
				fbp_oi, fbp_oo); \
		} \
	} while (0)
#define FPOLY_GM_SPLIT(f0, f1, f, i, hn, qn, k)   do { \
		uint32_t fbp_oi = (uint32_t)((hn) * sizeof(fpr)); \
		uint32_t fbp_oo = (uint32_t)((qn) * sizeof(fpr)); \
		const fpr *fbp_f = (f) + ((i) << 1); \
		fpr *fbp_f0 = (f0) + (i), *fbp_f1 = (f1) + (i); \
		if ((k) & 1) { \
			FNDSA_FXP96_BFLY_CALL(fndsa_fxp96_bfly_split_odd_x, \
				fbp_f, fbp_f + 1, fbp_f0, fbp_f1, &GM[(k) << 1], \
				fbp_oi, fbp_oo); \
		} else { \
			FNDSA_FXP96_BFLY_CALL(fndsa_fxp96_bfly_split_even_x, \
				fbp_f, fbp_f + 1, fbp_f0, fbp_f1, &GM[(k) << 1], \
				fbp_oi, fbp_oo); \
		} \
	} while (0)
#endif /* FNDSA_FFT_ASM_CM4 */

#endif /* FNDSA_TWIDDLE_U64 */

#endif
