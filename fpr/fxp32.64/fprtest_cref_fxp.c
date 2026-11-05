/*
 * C reference fpr operations for the on-board fprtest (fxp3264-cm4).
 * Compiles ops_c.c a second time with the C implementation forced and
 * the symbols renamed, so the C and assembly versions can be compared
 * in one binary.
 */

#undef FNDSA_FPR_IMPL
#define FNDSA_FPR_IMPL   FNDSA_FPR_IMPL_C

#define fndsa_fpr_scaled          fndsa_cref_fpr_scaled
#define fndsa_fpr_mul             fndsa_cref_fpr_mul
#define fndsa_fpr_inv             fndsa_cref_fpr_inv
#define fndsa_fpr_div             fndsa_cref_fpr_div
#define fndsa_fpr_sqrt            fndsa_cref_fpr_sqrt
#define fndsa_fxp96_cmul_rne64    fndsa_cref_fxp96_cmul_rne64
#define fndsa_fxp96_sqrt_mul_inv_sigma fndsa_cref_fxp96_sqrt_mul_inv_sigma
#define fndsa_fxp96_recip_x       fndsa_cref_fxp96_recip_x
#define fndsa_fxp96_mul_recip     fndsa_cref_fxp96_mul_recip
#define fndsa_fxp96_mul_recip_p63 fndsa_cref_fxp96_mul_recip_p63
#define fndsa_fxp96_cmul_gm       fndsa_cref_fxp96_cmul_gm
#define fndsa_fxp96_cmul_gmj      fndsa_cref_fxp96_cmul_gmj
#define fndsa_fpr_mul_gm          fndsa_cref_fpr_mul_gm
#define fndsa_fpr_mul_inv_q       fndsa_cref_fpr_mul_inv_q
#define fndsa_fpoly_mul_l10_fft   fndsa_cref_fpoly_mul_l10_fft
#define fndsa_fpoly_l10_neg       fndsa_cref_fpoly_l10_neg

#include "ops_c.c"

#include "fpr/fprtest.h"

static fpr c_of32(int32_t i)           { return fpr_of32(i); }
static fpr c_scaled(int64_t i, int sc) { return fpr_scaled(i, sc); }
static fpr c_add(fpr x, fpr y)         { return fpr_add(x, y); }
static fpr c_sub(fpr x, fpr y)         { return fpr_sub(x, y); }
static void c_add_sub(fpr *a, fpr *b, fpr x, fpr y) { FPR_ADD_SUB(*a, *b, x, y); }
static fpr c_mul(fpr x, fpr y)         { return fpr_mul(x, y); }
static fpr c_sqr(fpr x)                { return fpr_sqr(x); }
static fpr c_div(fpr x, fpr y)         { return fpr_div(x, y); }
static fpr c_inv(fpr x)                { return fpr_inv(x); }
static fpr c_sqrt(fpr x)               { return fpr_sqrt(x); }
#if FNDSA_SIGN_SQRT_FUSED && FNDSA_KGEN_GAMMA
static void
c_sqrt_inv_sigma(fpr *d, fpr x, unsigned logn)
{
	fndsa_cref_fxp96_sqrt_mul_inv_sigma(d->w, x.w, logn);
}
#endif
static void c_fpc_mul(fpr *d_re, fpr *d_im, fpr a_re, fpr a_im,
	fpr b_re, fpr b_im)
{
	FPC_MUL(*d_re, *d_im, a_re, a_im, b_re, b_im);
}

#ifdef FPR_MUL_RECIP_P63
static fpr
c_p63_to_fpr(uint64_t q)
{
	fpr d;

	d.w[0] = (uint32_t)(q << 1);
	d.w[1] = (uint32_t)(q >> 31);
	d.w[2] = (uint32_t)((int64_t)q >> 63);
	return d;
}
#endif

static fpr
c_mul_recip(fpr a, fpr x)
{
	fpr_rcp rcp;
	fpr d;

	FPR_RECIP_X(rcp, x);
	FPR_MUL_RECIP(d, a, rcp);
	return d;
}

#ifdef FPR_MUL_RECIP_P63
static fpr
c_mul_recip_p63(fpr a, fpr x)
{
	fpr_rcp rcp;
	uint64_t q;

	FPR_RECIP_X(rcp, x);
	FPR_MUL_RECIP_P63(q, a, rcp);
	return c_p63_to_fpr(q);
}
#endif

static fpr
c_mul_inv_q(fpr x)
{
	return FPR_MUL_INV_Q(x);
}

#if FNDSA_TWIDDLE_U64
static fpr
c_fpr_mul_gm(fpr x, uint64_t c)
{
	return fndsa_fpr_mul_gm(x, c);
}
#endif

/* rint/floor/trunc are the same inlines in both builds (0 = skip). */
const fprtest_ops fndsa_fprtest_cref = {
	.op_scaled = c_scaled,
	.op_of32 = c_of32,
	.op_add = c_add,
	.op_sub = c_sub,
	.op_add_sub = c_add_sub,
	.op_mul = c_mul,
	.op_sqr = c_sqr,
	.op_div = c_div,
	.op_inv = c_inv,
	.op_sqrt = c_sqrt,
#if FNDSA_SIGN_SQRT_FUSED && FNDSA_KGEN_GAMMA
	.op_sqrt_inv_sigma = c_sqrt_inv_sigma,
#endif
	.op_rint = 0,
	.op_floor = 0,
	.op_trunc = 0,
	.op_fpc_mul = c_fpc_mul,
	.op_mul_recip = c_mul_recip,
#ifdef FPR_MUL_RECIP_P63
	.op_mul_recip_p63 = c_mul_recip_p63,
#endif
	.op_mul_inv_q = c_mul_inv_q,
#if FNDSA_TWIDDLE_U64
	.op_fxp96_cmul_gm = fndsa_fxp96_cmul_gm,
	.op_fxp96_cmul_gmj = fndsa_fxp96_cmul_gmj,
	.op_fpr_mul_gm = c_fpr_mul_gm,
#endif
};
