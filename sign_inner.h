#ifndef FNDSA_SIGN_INNER_H__
#define FNDSA_SIGN_INNER_H__

/* ==================================================================== */
/*
 * This file includes declarations used by the signature generation code
 * only.
 */

#include "inner.h"

/* ==================================================================== */
/*
 * Floating-point values.
 *
 * The 'fpr' type, its constants and all operations on it are provided
 * by a compile-time-selected backend; see fpr/fpr_api.h for the
 * contract and fpr/ieee/ for the IEEE-754 (uint64_t) representation.
 */

#include "fpr/fpr_api.h"

#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264
#define FNDSA_HAS_FPR_DIV2E8_FXP3264   1
#define FNDSA_HAS_FPR_DIV2E9_FXP3264   1
/* x/2^8 and x/2^9, rounded to nearest (ties to even). */
static inline fpr
fpr_div2e8_fxp3264(fpr x)
{
	uint32_t x0 = x.w[0];
	uint32_t x1 = x.w[1];
	uint32_t x2 = x.w[2];
	uint32_t z0 = (x0 >> 8) | (x1 << 24);
	uint32_t z1 = (x1 >> 8) | (x2 << 24);
	uint32_t z2 = (uint32_t)((int32_t)x2 >> 8);
	uint32_t rb = (x0 >> 7) & 1;
	uint32_t inc = rb & ((uint32_t)((x0 & 0x7F) != 0) | (z0 & 1));
	uint64_t cc = (uint64_t)z0 + inc;

	z0 = (uint32_t)cc;
	cc = (cc >> 32) + (uint64_t)z1;
	z1 = (uint32_t)cc;
	z2 += (uint32_t)(cc >> 32);
	return (fpr){{ z0, z1, z2 }};
}

static inline fpr
fpr_div2e9_fxp3264(fpr x)
{
	uint32_t x0 = x.w[0];
	uint32_t x1 = x.w[1];
	uint32_t x2 = x.w[2];
	uint32_t z0 = (x0 >> 9) | (x1 << 23);
	uint32_t z1 = (x1 >> 9) | (x2 << 23);
	uint32_t z2 = (uint32_t)((int32_t)x2 >> 9);
	uint32_t rb = (x0 >> 8) & 1;
	uint32_t inc = rb & ((uint32_t)((x0 & 0xFF) != 0) | (z0 & 1));
	uint64_t cc = (uint64_t)z0 + inc;

	z0 = (uint32_t)cc;
	cc = (cc >> 32) + (uint64_t)z1;
	z1 = (uint32_t)cc;
	z2 += (uint32_t)(cc >> 32);
	return (fpr){{ z0, z1, z2 }};
}
#endif

#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 \
	&& !(FNDSA_SSE2 || FNDSA_NEON || FNDSA_RV64D)
#define FNDSA_HAS_FPR_RINT_DIV2E8_TO_U16   1
#define FNDSA_HAS_FPR_RINT_DIV2E9_TO_U16   1
/* round(x/2^8) and round(x/2^9), truncated to 16 bits. */
static inline uint16_t
fpr_rint_div2e8_to_u16(fpr x)
{
	uint32_t x0 = x.w[0];
	uint32_t x1 = x.w[1];
	uint32_t x2 = x.w[2];
	uint32_t q = (uint32_t)((int32_t)x2 >> 8);
	uint32_t rb = (x2 >> 7) & 1;
	uint32_t st = (uint32_t)((x0 | x1 | (x2 & 0x7F)) != 0);
	uint32_t inc = rb & (st | (q & 1));

	return (uint16_t)(q + inc);
}

static inline uint16_t
fpr_rint_div2e9_to_u16(fpr x)
{
	uint32_t x0 = x.w[0];
	uint32_t x1 = x.w[1];
	uint32_t x2 = x.w[2];
	uint32_t q = (uint32_t)((int32_t)x2 >> 9);
	uint32_t rb = (x2 >> 8) & 1;
	uint32_t st = (uint32_t)((x0 | x1 | (x2 & 0xFF)) != 0);
	uint32_t inc = rb & (st | (q & 1));

	return (uint16_t)(q + inc);
}
#endif

/* Fused reciprocal: FPR_RECIP_X computes 1/x in a wider local format
   (fpr_rcp, never stored), FPR_MUL_RECIP multiplies by it with a single
   rounding. The fallback is a plain multiply by fpr_inv(x). */
#ifndef FPR_RECIP_X
typedef fpr fpr_rcp;
#define FPR_RECIP_X(rcp, x)        ((rcp) = fpr_inv(x))
#define FPR_MUL_RECIP(d, a, rcp)   ((d) = fpr_mul((a), (rcp)))
/* Tells the benchmarks that the seam expands to plain fpr calls. */
#define FNDSA_FPR_RECIP_FALLBACK   1
#endif

/* x/q. With wide 1/q, the backend keeps more significant bits of 1/q
   than an fpr holds. */
#ifndef FPR_MUL_INV_Q
#define FPR_MUL_INV_Q(x)   fpr_mul((x), INV_Q)
#endif

/* The non-root l10 may be stored in a narrower format l10_t. FPR_L10_DIV
   computes mu = g01/g00 and returns both the stored l10 = conj(mu) in
   (s_re, s_im) and mu as fpr in (w_re, w_im). By default l10_t is fpr. */
#if FNDSA_SIGN_L10_64 && !defined(FPR_L10_DIV)
#error "FNDSA_SIGN_L10_64: this fpr representation has no narrow l10"
#endif
#ifndef FPR_L10_DIV
typedef fpr l10_t;
/* Index of the imaginary half within the l10 band. */
#define L10_IM_OFF(logn)   ((size_t)1 << ((logn) - 1))
#define fpoly_mul_l10_fft fpoly_mul_fft
#define fpoly_l10_neg fpoly_neg
#define FPR_L10_DIV(s_re, s_im, w_re, w_im, g_re, g_im, rcp)  do { \
		FPR_MUL_RECIP((w_re), (g_re), (rcp)); \
		FPR_MUL_RECIP((w_im), (g_im), (rcp)); \
		(s_re) = (w_re); \
		(s_im) = fpr_neg(w_im); \
	} while (0)
#endif

/* Special cases of FPC_MUL that a backend may implement more cheaply:
   SQR_SELF is d = d*d, MUL_SELF is d = d*b, NOALIAS needs d distinct
   from a and b. */
#ifndef FPC_SQR_SELF
#define FPC_SQR_SELF(d_re, d_im) \
	FPC_MUL((d_re), (d_im), (d_re), (d_im), (d_re), (d_im))
#endif
#ifndef FPC_MUL_SELF
#define FPC_MUL_SELF(d_re, d_im, b_re, b_im) \
	FPC_MUL((d_re), (d_im), (d_re), (d_im), (b_re), (b_im))
#endif
#ifndef FPC_MUL_NOALIAS
#define FPC_MUL_NOALIAS(d_re, d_im, a_re, a_im, b_re, b_im) \
	FPC_MUL((d_re), (d_im), (a_re), (a_im), (b_re), (b_im))
#endif

/* The FFT code reads the GM table only through these macros, so a
   representation can store the twiddles differently (unsigned twiddle
   table). TW/TW_CONJ load GM[k] or conj(GM[k]) into a local s for the
   BFLY macros (k is evaluated more than once). A backend that overrides
   one of the seven must override all of them. */
#ifndef FPOLY_GM_TW
#define FPOLY_GM_TW(s, k) \
	fpr s ## _re = GM[(k) << 1], s ## _im = GM[((k) << 1) + 1]
#define FPOLY_GM_TW_CONJ(s, k) \
	fpr s ## _re = GM[(k) << 1], s ## _im = fpr_neg(GM[((k) << 1) + 1])
/* Forward butterfly: (o1, o2) = (x + y*s, x - y*s). */
#define FPOLY_GM_BFLY_FWD(o1_re, o1_im, o2_re, o2_im, x_re, x_im, y_re, y_im, s)  do { \
		fpr fgm_re, fgm_im; \
		FPC_MUL_NOALIAS(fgm_re, fgm_im, \
			(y_re), (y_im), s ## _re, s ## _im); \
		FPR_ADD_SUB((o1_re), (o2_re), (x_re), fgm_re); \
		FPR_ADD_SUB((o1_im), (o2_im), (x_im), fgm_im); \
	} while (0)
/* Inverse butterfly: o1 = x + y, o2 = (x - y)*s, with s from TW_CONJ;
   x is consumed as scratch. */
#define FPOLY_GM_BFLY_INV(o1_re, o1_im, o2_re, o2_im, x_re, x_im, y_re, y_im, s)  do { \
		FPR_ADD_SUB((o1_re), (x_re), (x_re), (y_re)); \
		FPR_ADD_SUB((o1_im), (x_im), (x_im), (y_im)); \
		FPC_MUL_SELF((x_re), (x_im), s ## _re, s ## _im); \
		(o2_re) = (x_re); \
		(o2_im) = (x_im); \
	} while (0)
/* Merge step (forward): (o1, o2) = (a + b*GM[k], a - b*GM[k]);
   b is consumed as scratch. */
#define FPOLY_GM_MERGE_BFLY(o1_re, o1_im, o2_re, o2_im, a_re, a_im, b_re, b_im, k)  do { \
		FPC_MUL_SELF((b_re), (b_im), \
			GM[((k) << 1) + 0], GM[((k) << 1) + 1]); \
		FPR_ADD_SUB((o1_re), (o2_re), (a_re), (b_re)); \
		FPR_ADD_SUB((o1_im), (o2_im), (a_im), (b_im)); \
	} while (0)
/* Split step (inverse of the merge step): o0 = (a + b)/2,
   o1 = ((a - b)/2)*conj(GM[k]). */
#define FPOLY_GM_SPLIT_BFLY(o0_re, o0_im, o1_re, o1_im, a_re, a_im, b_re, b_im, k)  do { \
		fpr fgm_t_re, fgm_t_im, fgm_u_re, fgm_u_im, fgm_s_im; \
		FPR_ADD_SUB(fgm_t_re, fgm_u_re, (a_re), (b_re)); \
		FPR_ADD_SUB(fgm_t_im, fgm_u_im, (a_im), (b_im)); \
		(o0_re) = fpr_half(fgm_t_re); \
		(o0_im) = fpr_half(fgm_t_im); \
		fgm_s_im = fpr_neg(GM[((k) << 1) + 1]); \
		FPC_MUL_SELF(fgm_u_re, fgm_u_im, \
			GM[((k) << 1) + 0], fgm_s_im); \
		(o1_re) = fpr_half(fgm_u_re); \
		(o1_im) = fpr_half(fgm_u_im); \
	} while (0)
/* Same for a self-adjoint input (real FFT coefficients): o0_im is
   exactly zero and the twiddle product needs two scalar muls only. */
#define FPOLY_GM_SPLIT_SELFADJ(o0_re, o0_im, o1_re, o1_im, a_re, b_re, k)  do { \
		fpr fgm_t_re, fgm_u_re; \
		FPR_ADD_SUB(fgm_t_re, fgm_u_re, (a_re), (b_re)); \
		(o0_re) = fpr_half(fgm_t_re); \
		(o0_im) = FPR_ZERO; \
		fgm_u_re = fpr_half(fgm_u_re); \
		(o1_re) = fpr_mul(fgm_u_re, GM[((k) << 1) + 0]); \
		(o1_im) = fpr_mul(fgm_u_re, fpr_neg(GM[((k) << 1) + 1])); \
	} while (0)
#endif

/* Merge/split step i between f (n coefficients) and f0, f1, with
   hn = n/2, qn = n/4, k = i + hn. Takes the arrays so that fused
   butterflies can work in place. */
#ifndef FPOLY_GM_SPLIT
#define FPOLY_GM_MERGE(f, f0, f1, i, hn, qn, k)   do { \
		fpr fgp_a_re = (f0)[i], fgp_a_im = (f0)[(i) + (qn)]; \
		fpr fgp_b_re = (f1)[i], fgp_b_im = (f1)[(i) + (qn)]; \
		FPOLY_GM_MERGE_BFLY((f)[((i) << 1) + 0], (f)[((i) << 1) + 0 + (hn)], \
			(f)[((i) << 1) + 1], (f)[((i) << 1) + 1 + (hn)], \
			fgp_a_re, fgp_a_im, fgp_b_re, fgp_b_im, (k)); \
	} while (0)
#define FPOLY_GM_SPLIT(f0, f1, f, i, hn, qn, k)   do { \
		fpr fgp_a_re = (f)[((i) << 1) + 0], fgp_a_im = (f)[((i) << 1) + 0 + (hn)]; \
		fpr fgp_b_re = (f)[((i) << 1) + 1], fgp_b_im = (f)[((i) << 1) + 1 + (hn)]; \
		FPOLY_GM_SPLIT_BFLY((f0)[i], (f0)[(i) + (qn)], (f1)[i], (f1)[(i) + (qn)], \
			fgp_a_re, fgp_a_im, fgp_b_re, fgp_b_im, (k)); \
	} while (0)
#endif

/* Sampler leaves: sqrt(x) * INV_SIGMA[logn]. A backend with the fused
   square root overrides it to round once. */
#ifndef FPR_SQRT_MUL_INV_SIGMA
#define FPR_SQRT_MUL_INV_SIGMA(d, x, logn) \
	((d) = fpr_mul(fpr_sqrt(x), INV_SIGMA[(logn)].f))
#define FNDSA_FPR_SQRT_MUL_INV_SIGMA_FALLBACK   1
#endif

/* ==================================================================== */
/*
 * Signing temporary-buffer layout, in bytes per coefficient. The
 * ffSampling peak is 6 fpr polynomials plus one l10 band; hm (uint16_t)
 * and then G (int8_t) follow it.
 *
 *   stride  8 (ieee):             HM_OFF 56, G_OFF 58, BPC 59
 *   stride 12 (tw, fxp32.64):     HM_OFF 84, G_OFF 86, BPC 87
 *   stride 12, L10 in 64 bits:    HM_OFF 80, G_OFF 82, BPC 83
 */
#define FNDSA_SIGN_HM_OFF   ((size_t)6 * FPR_STRIDE + sizeof(l10_t))
#define FNDSA_SIGN_G_OFF    (FNDSA_SIGN_HM_OFF + 2)
#define FNDSA_SIGN_BPC      (FNDSA_SIGN_G_OFF + 1)

/* The ieee ffSampling assembly (sign_sampler_cm4.s) reads the binary64
   layout; other representations, Symplecroot and Symplectree use C. */
#if FNDSA_ASM_CORTEXM4 && FNDSA_FPR_REP == FNDSA_FPR_REP_IEEE \
	&& !FNDSA_SYMP_ROOT && !FNDSA_SYMP_TREE
#define FNDSA_FFSAMP_ASM   1
#else
#define FNDSA_FFSAMP_ASM   0
#endif

/* Split target t = t_int + t_frac: scalar or SSE2 only. */
#if FNDSA_SIGN_SPLIT && (FNDSA_NEON || FNDSA_RV64D)
#error "FNDSA_SIGN_SPLIT requires the scalar or SSE2 signing path"
#endif

/* fpr_p63(x) = floor(x * 2^63) for x in [0, 1), used by expm_p63().
   The generic form needs 63 integer bits; fxp32.64 and prec override it. */
#ifndef fpr_p63
#define fpr_p63(x)   ((uint64_t)fpr_trunc(fpr_mul2e((x), 63)))
#endif

/* Symplecroot: at the tree root, d11 = q^2/g00, since
   g00*g11 - |g01|^2 = q^2 in every FFT slot. Scalar or SSE2 only. */
#if FNDSA_SYMP_ROOT && (FNDSA_NEON || FNDSA_RV64D)
#error "FNDSA_SYMP_ROOT requires the scalar or SSE2 C signing path"
#endif

/* Symplectree: the right ffLDL subtree is the mirror of the left one
   (internal nodes negated, leaves x replaced by q^2/x). Scalar C only. */
#if FNDSA_SYMP_TREE && (FNDSA_SSE2 || FNDSA_NEON || FNDSA_RV64D)
#error "FNDSA_SYMP_TREE requires the scalar C signing path"
#endif

/* ==================================================================== */
/*
 * Floating-point polynomials.
 *
 * Polynomials are modulo X^n+1 and have real coefficients. When converted
 * to FFT format, they have n/2 complex coefficients, with the real and
 * imaginary parts separated (real in f[i] and imaginary in f[i+n/2]).
 *
 * Unless specified explicitly, the functions below can work with
 * polynomials in both real and FFT representations (if using several
 * operands, then they must all be in the same representation).
 */

/* Convert polynomial f from real to FFT representation (in-place). */
#define fpoly_FFT   fndsa_fpoly_FFT
void fpoly_FFT(unsigned logn, fpr *f);

/* The integer-source FFTs below (FFT-new) need the RNJ table, which
   only these representations provide. */
#define FNDSA_FPR_HAS_RNJ \
	(FNDSA_FPR_REP == FNDSA_FPR_REP_IEEE \
		|| FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 \
		|| FNDSA_FPR_REP == FNDSA_FPR_REP_PREC)
#if FNDSA_FFT_NEW && !FNDSA_FPR_HAS_RNJ
#error "FNDSA_FFT_NEW needs a representation that carries the RNJ table"
#endif

#if FNDSA_FFT_NEW
/* Compact FFT from int32_t input, converted to fpr after the split
   phase. */
#define fpoly_FFT_int32_selfadj   fndsa_fpoly_FFT_int32_selfadj
void fpoly_FFT_int32_selfadj(unsigned logn, fpr *d, const int32_t *a);

#define fpoly_FFT_int32   fndsa_fpoly_FFT_int32
void fpoly_FFT_int32(unsigned logn, fpr *d, const int32_t *a);

/* Same from int8_t or uint16_t input, without the int32_t copy.
   Source and destination must not overlap. */
#define fpoly_FFT_int8   fndsa_fpoly_FFT_int8
void fpoly_FFT_int8(unsigned logn, fpr *d, const int8_t *a);

#define fpoly_FFT_uint16   fndsa_fpoly_FFT_uint16
void fpoly_FFT_uint16(unsigned logn, fpr *d, const uint16_t *a);
#endif

/* Convert polynomial f from FFT to real representation (in-place). */
#define fpoly_iFFT   fndsa_fpoly_iFFT
void fpoly_iFFT(unsigned logn, fpr *f);

#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 \
	&& !(FNDSA_SSE2 || FNDSA_NEON || FNDSA_RV64D)
/* fpoly_iFFT without the final division by n/2. */
#define fpoly_iFFT_no_div2e   fndsa_fpoly_iFFT_no_div2e
void fpoly_iFFT_no_div2e(unsigned logn, fpr *f);
#endif

/* Set polynomial d from small polynomial f with integer coefficients. */
#define fpoly_set_small   fndsa_fpoly_set_small
void fpoly_set_small(unsigned logn, fpr *d, const int8_t *f);

/* Add polynomial b to polynomial a. */
#define fpoly_add   fndsa_fpoly_add
void fpoly_add(unsigned logn, fpr *a, const fpr *b);

/* Subtract polynomial b from polynomial a. */
#define fpoly_sub   fndsa_fpoly_sub
void fpoly_sub(unsigned logn, fpr *a, const fpr *b);

/* Negate polynomial a (in-place). */
#define fpoly_neg   fndsa_fpoly_neg
void fpoly_neg(unsigned logn, fpr *a);

/* Multiply polynomial a with polynomial b (FFT representation only). */
#define fpoly_mul_fft   fndsa_fpoly_mul_fft
void fpoly_mul_fft(unsigned logn, fpr *a, const fpr *b);

/* unused
   Multiply polynomial a with the adjoint of polynomial b (FFT
   representation only).
void fpoly_muladj_fft(unsigned logn, fpr *a, const fpr *b);
*/

/* unused
   Multiply polynomial a with its own adjoint (FFT representation only).
   Coefficients n/2 to n-1 are set to zero.
void fpoly_mulownadj_fft(unsigned logn, fpr *a);
*/

/* Multiply polynomial a with real constant x. */
#define fpoly_mulconst   fndsa_fpoly_mulconst
void fpoly_mulconst(unsigned logn, fpr *a, fpr x);

/* Multiply polynomial a by the wide 1/q (FNDSA_SIGN_INVQ_WIDE only). */
#define fpoly_mul_inv_q   fndsa_fpoly_mul_inv_q
void fpoly_mul_inv_q(unsigned logn, fpr *a);

/* Perform an LDL decomposition of a self-adjoint matrix G. The matrix
   is G = [[g00, g01], [adj(g01), g11]]; g00 and g11 are self-adjoint
   polynomials. The decomposition is G = L*D*adj(L), with:
      D = [[g00, 0], [0, d11]]
      L = [[1, 0], [l10, 1]]
   The output polynomial l10 is written to l10 (which may alias g01:
   each slot is read before it is written), and d11 over g11. g00,
   g11 and d11 are self-adjoint: only their first n/2 coefficients
   are accessed. g00 is unmodified. All polynomials are in FFT
   representation. */
#define fpoly_LDL_fft   fndsa_fpoly_LDL_fft
void fpoly_LDL_fft(unsigned logn, const fpr *g00, const fpr *g01,
	fpr *g11, l10_t *l10);
#if FNDSA_SYMP_ROOT
/* Root LDL with d11 = q^2/g00; l10 is written over g01. */
#define fpoly_LDL_root_fft   fndsa_fpoly_LDL_root_fft
void fpoly_LDL_root_fft(unsigned logn, const fpr *g00, fpr *g01, fpr *g11);
#elif FNDSA_SYMP_TREE
/* Root LDL for the Symplectree: l10 only, written over g01. */
#define fpoly_LDL_root_fft   fndsa_fpoly_LDL_root_fft
void fpoly_LDL_root_fft(unsigned logn, const fpr *g00, fpr *g01);
#endif

/* Split operation on a polynomial: for input polynomial f,
   half-size polynomials f0 and f1 (modulo X^(n/2)+1) are such that
   f = f0(x^2) + x*f1(x^2). All polynomials are in FFT representation. */
#define fpoly_split_fft   fndsa_fpoly_split_fft
void fpoly_split_fft(unsigned logn, fpr *f0, fpr *f1, const fpr *f);

/* Specialized version of fpoly_split_fft() when the source polynomial
   f is self-adjoint. Only the first n/2 coefficients of f are accessed.
   On output, f0 is self-adjoint (all its n/2 coefficients are set), but
   in general f1 is not self-adjoint. */
#define fpoly_split_selfadj_fft   fndsa_fpoly_split_selfadj_fft
void fpoly_split_selfadj_fft(unsigned logn, fpr *f0, fpr *f1, const fpr *f);

/* Merge operation on polynomials: for input half-size polynomials f0
   and f1 (modulo X^(n/2)+1), compute f = f0(x^2) + x*f1(x^2). All
   polynomials are in FFT representation. */
#define fpoly_merge_fft   fndsa_fpoly_merge_fft
void fpoly_merge_fft(unsigned logn, fpr *f, const fpr *f0, const fpr *f1);

/* Given matrix B = [[b00, b01], [b10, b11]], compute the Gram matrix
   G = B*adj(B) = [[g00, g01], [g10, g11]], with:
      g00 = b00*adj(b00) + b01*adj(b01)
      g01 = b00*adj(b10) + b01*adj(b11)
      g10 = b10*adj(b00) + b11*adj(b01)
      g11 = b10*adj(b10) + b11*adj(b11)
   Only g00, g01 and g11 are returned, in b00, b01 and b10, respectively.
   b11 is unmodified. All polynomials are in FFT representation. */
#define fpoly_gram_fft   fndsa_fpoly_gram_fft
void fpoly_gram_fft(unsigned logn,
	fpr *b00, fpr *b01, fpr *b10, const fpr *b11);

/* Given matrix B = [[b00, b01], [b10, b11]], compute the target vector
   [t0,t1] = (1/q)*B*[0,hm], for polynomial hm with coefficients in [0,q-1].
   Only b01 and b11 are needed. b01 and b11 are consumed. t1 may be the
   same pointer as b11, or as b01 (b01 is then multiplied in place and
   the move is skipped). hm is in normal representation; all other
   polynomials are in FFT representation. */
#define fpoly_apply_basis   fndsa_fpoly_apply_basis
void fpoly_apply_basis(unsigned logn, fpr *t0, fpr *t1,
	fpr *b01, fpr *b11, const uint16_t *hm);

/* ==================================================================== */
/*
 * Gaussian sampling.
 *
 * The sampler follows a Gaussian distribution whose centre and standard
 * deviations are not integral, dynamically obtained, and secret. The
 * sampler state includes a PRNG, from which random bytes are obtained
 * to sample value with a given bimodal Gaussian used as source for
 * rejection sampling of the target distribution.
 */

typedef struct {
#if FNDSA_SHAKE256X4
	shake256x4_context pc;
#else
	shake_context pc;
#endif
	unsigned logn;
} sampler_state;

/* Initialize the sampler for a given degree and seed. */
#define sampler_init   fndsa_sampler_init
void sampler_init(sampler_state *ss, unsigned logn,
	const void *seed, size_t seed_len);

/* Sample the next small integer. Parameters are:
      ss       sampler state
      mu       distribution centre
      isigma   inverse of the distribution standard deviation  */
#define sampler_next   fndsa_sampler_next
int32_t sampler_next(sampler_state *ss, fpr mu, fpr isigma);

/* Apply Fast Fourier sampling:
      ss    sampler state (initialized)
      tmp   parameters and free space (total: 7*n fpr slots):
	  t0  (n fpr slots)
	  t1  (n fpr slots)
	  g01 (n fpr slots)
          g00 (hn fpr slots)
	  g11 (hn fpr slots)
	  free space (at least 3*n fpr slots)
   Target vector is [t0, t1]. Gram matrix is [[g00, g01], [adj(g01), g11]]
   with g00 and g11 being both self-adjoint (hence using only hn slots each).
   All polynomials are in FFT representation.
   All values are consumed. On output, the sampled vector [z0,z1] is written
   over [t0,t1]. */
#define ffsamp_fft   fndsa_ffsamp_fft
void ffsamp_fft(sampler_state *ss, fpr *tmp);

/* This function is global on ARM Cortex M4 so that it can be called
   from assembly code. We define its global name here so that test code
   can override it (in test_sampler.c and test_sign.c). */
#if FNDSA_ASM_CORTEXM4
#define ffsamp_fft_deepest   fndsa_ffsamp_fft_deepest
#endif

/* ==================================================================== */
/*
 * Mod-q' (q' = 18433) layer from c-fn-dsa-alt (sign_mq2.c/sign_mq2_cm4.s).
 * On Cortex-M4 the internal representation is relaxed [0,q'].
 */

#define MQ2_MOD         18433
#define MQ2_IQ              3      /*  1/q mod q' */

#define mq2poly_small_to_mod           fndsa_mq2poly_small_to_mod
#define mq2poly_signed_to_mod          fndsa_mq2poly_signed_to_mod
#define mq2poly_unsigned_to_mod        fndsa_mq2poly_unsigned_to_mod
#define mq2poly_mod_to_signed          fndsa_mq2poly_mod_to_signed
#define mq2poly_NTT                    fndsa_mq2poly_NTT
#define mq2poly_iNTT                   fndsa_mq2poly_iNTT
#define mq2poly_mul_ntt                fndsa_mq2poly_mul_ntt
#define mq2poly_muladj_ntt             fndsa_mq2poly_muladj_ntt
#define mq2poly_muladj_x2_ntt          fndsa_mq2poly_muladj_x2_ntt
#define mq2poly_sub                    fndsa_mq2poly_sub
#define mq2poly_add                    fndsa_mq2poly_add
#define mq2poly_mulconst               fndsa_mq2poly_mulconst
#define mq2poly_addmulconst            fndsa_mq2poly_addmulconst
#define mq2_GM                         fndsa_mq2_GM
#define mq2_iGM                        fndsa_mq2_iGM
#define mq2mqpoly_adjust_save          fndsa_mq2mqpoly_adjust_save
#define mq2mqpoly_adjust               fndsa_mq2mqpoly_adjust
#define mq2mqpoly_CRT                  fndsa_mq2mqpoly_CRT

/* Convert a polynomial from small (f) to mod (d). */
void mq2poly_small_to_mod(unsigned logn, const int8_t *f, uint16_t *d);

/* Convert a polynomial from signed to mod (in-place). */
void mq2poly_signed_to_mod(unsigned logn, uint16_t *d);

/* Convert a polynomial from unsigned to mod (in-place). */
#if FNDSA_ASM_CORTEXM4
static inline void
mq2poly_unsigned_to_mod(unsigned logn, uint16_t *d)
{
	/* For integers modulo q', we use the "relaxed" format with
	   values in [0, q'], so the unsigned values in [0, q] are
	   already fine. */
	(void)logn;
	(void)d;
}
#else
void mq2poly_unsigned_to_mod(unsigned logn, uint16_t *d);
#endif

/* Convert a polynomial from mod to signed (in-place). */
void mq2poly_mod_to_signed(unsigned logn, uint16_t *d);

/* Convert a polynomial from mod to ntt (in-place). */
void mq2poly_NTT(unsigned logn, uint16_t *d);

/* Convert a polynomial from ntt to mod (in-place). */
void mq2poly_iNTT(unsigned logn, uint16_t *d);

/* Multiply polynomial a by polynomial b (both in ntt representation). */
void mq2poly_mul_ntt(unsigned logn, uint16_t *a, const uint16_t *b);

/* a <- a*adj(b) (ntt representation.) */
void mq2poly_muladj_ntt(unsigned logn, uint16_t *a, const uint16_t *b);

/* a <- a*adj(a) + b*adj(b) (ntt representation). */
void mq2poly_muladj_x2_ntt(unsigned logn, uint16_t *a, const uint16_t *b);

/* Subtract polynomial b from polynomial a (both must be in int
   representation, or both must be in ntt representation). */
void mq2poly_sub(unsigned logn, uint16_t *a, const uint16_t *b);

/* Add polynomial b to polynomial a (both must be in int
   representation, or both must be in ntt representation). */
void mq2poly_add(unsigned logn, uint16_t *a, const uint16_t *b);

/* For all i = 0 to n-1, replace a[i] with a[i]*m (mod q'). */
void mq2poly_mulconst(unsigned logn, uint16_t *a, unsigned m);

/* For all i = 0 to n-1, replace a[i] with (a[i] + s)*m (mod q'). */
void mq2poly_addmulconst(unsigned logn, uint16_t *a, unsigned m, unsigned s);

/* Given polynomials a and b, representing the same polynomial modulo q'
   and modulo q, respectively (both normalized), adjust the value of a
   so that it matches b. This function assumes that the underlying polynomial
   has coefficients in [-q' + 1, q' - 1]. Adjustment bits are saved into
   adj[] (which must be large enough to contain n bits). */
void mq2mqpoly_adjust_save(unsigned logn,
	int16_t *a, const int16_t *b, uint8_t *adj);

/* Reduced variant of mq2mqpoly_adjust_save(), assuming that adjustment bits
   have already been set in adj[]. */
void mq2mqpoly_adjust(unsigned logn, int16_t *a, const uint8_t *adj);

/* Given a mod q and a mod q' (both normalized), compute a mod q*q' into d
   (normalized) */
void mq2mqpoly_CRT(unsigned logn, int32_t *d,
	const int16_t *aq, const int16_t *aq2);

/* Tables of constants for NTT and inverse NTT modulo q. */
extern const uint16_t mq2_GM[];
extern const uint16_t mq2_iGM[];

/* ==================================================================== */
/*
 * Internal signing function.
 */

/* Internal signing function. The complete signing key (encoded for f,
   g and F, but skipping the leading header byte, and decoded for G) is
   provided, as well as the hashed verifying key, data to sign (context,
   id, hash value), the random seed to work on (optional), the signature
   output buffer, and the temporary area. The signature buffer has been
   verified to be large enough. The temporary area is large enough and
   32-byte aligned.

   Returned value is the signature size (in bytes), or 0 on error. An
   error is possible if seed is NULL and the system RNG fails.

   tmp size: FNDSA_SIGN_G_OFF*n bytes  */
#define sign_core   fndsa_sign_core
size_t sign_core(unsigned logn,
	const uint8_t *sign_key_fgF, const int8_t *G,
	const uint8_t *hashed_vk, const uint8_t *ctx, size_t ctx_len,
	const char *id, const uint8_t *hv, size_t hv_len,
	const uint8_t *seed, size_t seed_len, uint8_t *sig, void *tmp);

/* ==================================================================== */

#endif
