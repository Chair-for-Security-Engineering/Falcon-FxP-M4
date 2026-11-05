#ifndef FNDSA_FPR_PREC_H__
#define FNDSA_FPR_PREC_H__

/* ==================================================================== */
/*
 * 'prec' representation, for precision evaluation (native only).
 *
 * Each fpr holds two values computed in lockstep: ref, an MPFR value at
 * PREC_REF_BITS (default 200) rounded to nearest at every operation,
 * and tgt, the value the embedded target backend computes with its own
 * ops (target_glue.c; PREC_TARGET=none leaves it zero).
 *
 * MPFR memory comes from the Boehm GC and is never freed explicitly.
 * Every op makes a fresh value, so plain copies of fpr values are safe.
 */

#include <mpfr.h>

#if FNDSA_ASM_CORTEXM4 || (defined FNDSA_FPR_IMPL \
	&& FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4)
#error the prec representation is a native-only analysis backend (no cm4)
#endif

/* Reference precision in bits. */
#ifndef PREC_REF_BITS
#define PREC_REF_BITS   200
#endif

/* The target value as opaque bytes, sized for the largest target, so
   FPR_STRIDE is the same for all targets. */
typedef struct {
	mpfr_t ref;
	uint64_t tgt[2];
} prec_fpr;

typedef prec_fpr fpr;

/* sizeof(__mpfr_struct) is 32 on LP64; checked below. */
#define FPR_STRIDE   48
typedef char fndsa_prec_stride_check[(sizeof(prec_fpr) == FPR_STRIDE)
	? 1 : -1];

/* ==================================================================== */
/* Out-of-line operations (ops_mpfr.c). The fndsa_prec_ prefix keeps
   them apart from the target's fndsa_fpr_* symbols.                    */

#define fpr_scaled   fndsa_prec_scaled
#define fpr_of       fndsa_prec_of
#define fpr_of32     fndsa_prec_of32
#define fpr_rint     fndsa_prec_rint
#define fpr_floor    fndsa_prec_floor
#define fpr_trunc    fndsa_prec_trunc
#define fpr_add      fndsa_prec_add
#define fpr_sub      fndsa_prec_sub
#define fpr_neg      fndsa_prec_neg
#define fpr_half     fndsa_prec_half
#define fpr_double   fndsa_prec_double
#define fpr_mul      fndsa_prec_mul
#define fpr_sqr      fndsa_prec_sqr
#define fpr_div      fndsa_prec_div
#define fpr_inv      fndsa_prec_inv
#define fpr_sqrt     fndsa_prec_sqrt
#define fpr_div2e    fndsa_prec_div2e
#define fpr_mul2e    fndsa_prec_mul2e
/* The default fpr_p63 goes through fpr_mul2e(x, 63), which would wrap
   a fixed-point target, so the target's own fpr_p63 is used. */
#define fpr_p63      fndsa_prec_p63

fpr fpr_scaled(int64_t i, int sc);
fpr fpr_of(int64_t i);
fpr fpr_of32(int32_t i);
int64_t fpr_rint(fpr x);
int64_t fpr_floor(fpr x);
int64_t fpr_trunc(fpr x);
uint64_t fpr_p63(fpr x);
fpr fpr_add(fpr x, fpr y);
fpr fpr_sub(fpr x, fpr y);
fpr fpr_neg(fpr x);
fpr fpr_half(fpr x);
fpr fpr_double(fpr x);
fpr fpr_mul(fpr x, fpr y);
fpr fpr_sqr(fpr x);
fpr fpr_div(fpr x, fpr y);
fpr fpr_inv(fpr x);
fpr fpr_sqrt(fpr x);
fpr fpr_div2e(fpr x, unsigned e);
fpr fpr_mul2e(fpr x, unsigned e);

/* Fused reciprocal: the reference is a plain reciprocal and multiply,
   the target runs its own FPR_RECIP_X/FPR_MUL_RECIP. */
typedef struct {
	mpfr_t ref;
	uint64_t tgt[3];
} prec_rcp;
typedef prec_rcp fpr_rcp;
fpr_rcp fndsa_prec_recip_x(fpr x);
fpr fndsa_prec_mul_recip(fpr a, fpr_rcp rcp);
#define FPR_RECIP_X(rcp, x)       ((rcp) = fndsa_prec_recip_x(x))
#define FPR_MUL_RECIP(d, a, rcp)  ((d) = fndsa_prec_mul_recip((a), (rcp)))

/* l10: l10_t stays a prec fpr; only the target value is rounded onto
   the narrow grid. */
#if FNDSA_SIGN_L10_64
fpr fndsa_prec_mul_recip_l10(fpr a, fpr_rcp rcp);
typedef fpr l10_t;
#define L10_IM_OFF(logn)   ((size_t)1 << ((logn) - 1))
#define fpoly_mul_l10_fft  fpoly_mul_fft
#define fpoly_l10_neg      fpoly_neg
#define FPR_L10_DIV(s_re, s_im, w_re, w_im, g_re, g_im, rcp)  do { \
		(w_re) = fndsa_prec_mul_recip_l10((g_re), (rcp)); \
		(w_im) = fndsa_prec_mul_recip_l10((g_im), (rcp)); \
		(s_re) = (w_re); \
		(s_im) = fpr_neg(w_im); \
	} while (0)
#endif

/* Fused sampler leaf: the reference is a plain sqrt and multiply. */
#if defined(FNDSA_SIGN_SQRT_FUSED) && FNDSA_SIGN_SQRT_FUSED
fpr fndsa_prec_sqrt_mul_inv_sigma(fpr x, unsigned logn);
#define FPR_SQRT_MUL_INV_SIGMA(d, x, logn) \
	((d) = fndsa_prec_sqrt_mul_inv_sigma((x), (logn)))
#endif

/* x/q: the reference uses the exact 1/q, the target FPR_MUL_INV_Q. */
fpr fndsa_prec_mul_inv_q(fpr x);
#define FPR_MUL_INV_Q(x)   fndsa_prec_mul_inv_q(x)

/* Sets up the GC and the constant tables. Idempotent; also runs from a
   constructor and from every op. */
void fndsa_prec_init(void);

/* ==================================================================== */
/* Constants (runtime-constructed; see fpr_const_prec.c).               */

/* Constructors for scalar constants. FPR(i, e) = i * 2^e, exact. */
#define FPR(i, e)   fndsa_prec_scaled((int64_t)(i), (e))
fpr fndsa_prec_zero(void);
fpr fndsa_prec_nzero(void);
fpr fndsa_prec_one(void);
fpr fndsa_prec_c_inv_2sqrsigma0(void);
fpr fndsa_prec_c_log2(void);
fpr fndsa_prec_c_inv_log2(void);
fpr fndsa_prec_c_inv_q(void);
fpr fndsa_prec_c_minus_inv_q(void);
fpr fndsa_prec_c_q2(void);

#define FPR_ZERO    (fndsa_prec_zero())
#define FPR_NZERO   (fndsa_prec_nzero())
#define FPR_ONE     (fndsa_prec_one())

/* 1/(2*(1.8205^2)) */
#define INV_2SQRSIGMA0   (fndsa_prec_c_inv_2sqrsigma0())
/* log(2) */
#define LOG2       (fndsa_prec_c_log2())
/* 1/log(2) */
#define INV_LOG2   (fndsa_prec_c_inv_log2())
/* 1/q and -1/q (q = 12289), used to compute the target. */
#define INV_Q         (fndsa_prec_c_inv_q())
#define MINUS_INV_Q   (fndsa_prec_c_minus_inv_q())
/* q^2 = 151019521 (exact), the determinant of the Gram matrix. */
#define Q2            (fndsa_prec_c_q2())

/* ==================================================================== */
/* Constant-time 64-bit shift helpers.                                  */

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
/* Macros required by the contract.                                     */

/* Combined addition/subtraction: a <- x + y, b <- x - y. */
#define FPR_ADD_SUB(a, b, x, y)   do { \
		fpr t_add_sub_x = (x); \
		fpr t_add_sub_y = (y); \
		(a) = fpr_add(t_add_sub_x, t_add_sub_y); \
		(b) = fpr_sub(t_add_sub_x, t_add_sub_y); \
	} while (0)

/* Complex multiplication; the target side uses the target's own
   FPC_MUL, so a fused kernel is measured as such. */
void fndsa_prec_fpc_mul(prec_fpr *d_re, prec_fpr *d_im,
	prec_fpr a_re, prec_fpr a_im, prec_fpr b_re, prec_fpr b_im);
#define FPC_MUL(d_re, d_im, a_re, a_im, b_re, b_im)   do { \
		prec_fpr fpct_d_re, fpct_d_im; \
		fndsa_prec_fpc_mul(&fpct_d_re, &fpct_d_im, \
			(a_re), (a_im), (b_re), (b_im)); \
		(d_re) = fpct_d_re; \
		(d_im) = fpct_d_im; \
	} while (0)

/* ==================================================================== */
/* Data tables, filled at startup in fpr_const_prec.c: the reference
   values by MPFR, the target values from the target's own tables.     */

typedef union {
	fpr f;
} fpr_u;

#define GM          fndsa_prec_GM
#define RNJ         fndsa_prec_RNJ
#define INV_SIGMA   fndsa_prec_INV_SIGMA
#define SIGMA_MIN   fndsa_prec_SIGMA_MIN
extern fpr GM[];
#if FNDSA_FFT_NEW
/* Compact self-adjoint FFT constants. */
extern fpr RNJ[];
#endif
extern fpr_u INV_SIGMA[];
extern fpr_u SIGMA_MIN[];

/* ==================================================================== */

#endif
