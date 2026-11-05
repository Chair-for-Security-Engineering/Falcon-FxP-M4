/*
 * Target glue for the prec backend (see target_glue.h). This file is
 * compiled with the target's representation selected, so 'fpr', the
 * contract ops and the tables are the target's own. PREC_TARGET=none
 * builds it with PREC_TGT_NONE=1 and everything is a no-op.
 */

#include "sign_inner.h"

#include <string.h>
#include "fpr/prec/target_glue.h"

#ifdef PREC_TGT_NONE

const int prec_tgt_enabled = 0;
const int prec_tgt_metric = PREC_TGT_METRIC_FLOAT;
const int prec_tgt_frac_bits = 0;
const char *const prec_tgt_name = "none";

static void
blob_zero(uint64_t *r)
{
	r[0] = 0;
	r[1] = 0;
}

void prec_tgt_scaled(uint64_t *r, int64_t i, int sc)
{ (void)i; (void)sc; blob_zero(r); }
void prec_tgt_recip_x(uint64_t *r, const uint64_t *x)
{ (void)x; r[0] = 0; r[1] = 0; r[2] = 0; }
void prec_tgt_mul_recip(uint64_t *r, const uint64_t *a, const uint64_t *rc)
{ (void)a; (void)rc; blob_zero(r); }
void prec_tgt_mul_recip_l10(uint64_t *r, const uint64_t *a,
	const uint64_t *rc)
{ (void)a; (void)rc; blob_zero(r); }
void prec_tgt_sqrt_mul_inv_sigma(uint64_t *r, const uint64_t *x,
	unsigned logn)
{ (void)x; (void)logn; blob_zero(r); }
void prec_tgt_mul_inv_q(uint64_t *r, const uint64_t *a)
{ (void)a; blob_zero(r); }
void prec_tgt_of(uint64_t *r, int64_t i) { (void)i; blob_zero(r); }
void prec_tgt_of32(uint64_t *r, int32_t i) { (void)i; blob_zero(r); }
void prec_tgt_const_scalar(uint64_t *r, int which)
{ (void)which; blob_zero(r); }
void prec_tgt_gm(uint64_t *r, int i) { (void)i; blob_zero(r); }
void prec_tgt_rnj(uint64_t *r, int i) { (void)i; blob_zero(r); }
void prec_tgt_inv_sigma(uint64_t *r, unsigned logn)
{ (void)logn; blob_zero(r); }
void prec_tgt_sigma_min(uint64_t *r, unsigned logn)
{ (void)logn; blob_zero(r); }
int64_t prec_tgt_rint(const uint64_t *a) { (void)a; return 0; }
int64_t prec_tgt_floor(const uint64_t *a) { (void)a; return 0; }
int64_t prec_tgt_trunc(const uint64_t *a) { (void)a; return 0; }
uint64_t prec_tgt_p63(const uint64_t *a) { (void)a; return 0; }
void prec_tgt_add(uint64_t *r, const uint64_t *a, const uint64_t *b)
{ (void)a; (void)b; blob_zero(r); }
void prec_tgt_cmul(uint64_t *r_re, uint64_t *r_im,
	const uint64_t *a_re, const uint64_t *a_im,
	const uint64_t *b_re, const uint64_t *b_im)
{ (void)a_re; (void)a_im; (void)b_re; (void)b_im;
  blob_zero(r_re); blob_zero(r_im); }
void prec_tgt_sub(uint64_t *r, const uint64_t *a, const uint64_t *b)
{ (void)a; (void)b; blob_zero(r); }
void prec_tgt_mul(uint64_t *r, const uint64_t *a, const uint64_t *b)
{ (void)a; (void)b; blob_zero(r); }
void prec_tgt_div(uint64_t *r, const uint64_t *a, const uint64_t *b)
{ (void)a; (void)b; blob_zero(r); }
void prec_tgt_neg(uint64_t *r, const uint64_t *a) { (void)a; blob_zero(r); }
void prec_tgt_half(uint64_t *r, const uint64_t *a) { (void)a; blob_zero(r); }
void prec_tgt_double(uint64_t *r, const uint64_t *a) { (void)a; blob_zero(r); }
void prec_tgt_sqr(uint64_t *r, const uint64_t *a) { (void)a; blob_zero(r); }
void prec_tgt_inv(uint64_t *r, const uint64_t *a) { (void)a; blob_zero(r); }
void prec_tgt_sqrt(uint64_t *r, const uint64_t *a) { (void)a; blob_zero(r); }
void prec_tgt_div2e(uint64_t *r, const uint64_t *a, unsigned e)
{ (void)a; (void)e; blob_zero(r); }
void prec_tgt_mul2e(uint64_t *r, const uint64_t *a, unsigned e)
{ (void)a; (void)e; blob_zero(r); }
void prec_tgt_to_mpfr(mpfr_t rop, const uint64_t *a)
{ (void)a; mpfr_set_zero(rop, +1); }

#else /* real target below */

#if FNDSA_FPR_REP == FNDSA_FPR_REP_PREC
#error target_glue.c must be compiled with the TARGET representation
#endif

/* The blob must hold the target's fpr. */
typedef char prec_tgt_size_check[(sizeof(fpr) <= 16) ? 1 : -1];

const int prec_tgt_enabled = 1;
#ifdef FXP_FRAC_BITS
const int prec_tgt_metric = PREC_TGT_METRIC_FXP;
const int prec_tgt_frac_bits = FXP_FRAC_BITS;
#else
const int prec_tgt_metric = PREC_TGT_METRIC_FLOAT;
const int prec_tgt_frac_bits = 0;
#endif
#if FNDSA_FPR_REP == FNDSA_FPR_REP_IEEE
const char *const prec_tgt_name = "ieee";
#elif FNDSA_FPR_REP == FNDSA_FPR_REP_TW
const char *const prec_tgt_name = "tw";
#elif FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264
const char *const prec_tgt_name = "fxp32.64";
#else
#error unknown target representation
#endif

static inline fpr
blob_get(const uint64_t *a)
{
	fpr x;

	memcpy(&x, a, sizeof x);
	return x;
}

static inline void
blob_put(uint64_t *r, fpr x)
{
	r[0] = 0;
	r[1] = 0;
	memcpy(r, &x, sizeof x);
}

void
prec_tgt_scaled(uint64_t *r, int64_t i, int sc)
{
	blob_put(r, fpr_scaled(i, sc));
}

void
prec_tgt_of(uint64_t *r, int64_t i)
{
	blob_put(r, fpr_of(i));
}

void
prec_tgt_of32(uint64_t *r, int32_t i)
{
	blob_put(r, fpr_of32(i));
}

void
prec_tgt_const_scalar(uint64_t *r, int which)
{
	switch (which) {
	case PREC_TGT_C_ZERO:            blob_put(r, FPR_ZERO); break;
	case PREC_TGT_C_NZERO:           blob_put(r, FPR_NZERO); break;
	case PREC_TGT_C_ONE:             blob_put(r, FPR_ONE); break;
	case PREC_TGT_C_INV_2SQRSIGMA0:  blob_put(r, INV_2SQRSIGMA0); break;
	case PREC_TGT_C_LOG2:            blob_put(r, LOG2); break;
	case PREC_TGT_C_INV_LOG2:        blob_put(r, INV_LOG2); break;
	case PREC_TGT_C_INV_Q:           blob_put(r, INV_Q); break;
	case PREC_TGT_C_MINUS_INV_Q:     blob_put(r, MINUS_INV_Q); break;
	case PREC_TGT_C_Q2:              blob_put(r, Q2); break;
	default:                         blob_put(r, FPR_ZERO); break;
	}
}

void
prec_tgt_gm(uint64_t *r, int i)
{
	blob_put(r, GM[i]);
}

void
prec_tgt_rnj(uint64_t *r, int i)
{
	/* Only targets with the compact-FFT table. */
#if FNDSA_FFT_NEW && (FNDSA_FPR_REP == FNDSA_FPR_REP_IEEE \
	|| FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264)
	blob_put(r, RNJ[i]);
#else
	(void)i;
	blob_put(r, FPR_ZERO);
#endif
}

void
prec_tgt_inv_sigma(uint64_t *r, unsigned logn)
{
	blob_put(r, INV_SIGMA[logn].f);
}

void
prec_tgt_sigma_min(uint64_t *r, unsigned logn)
{
	blob_put(r, SIGMA_MIN[logn].f);
}

int64_t
prec_tgt_rint(const uint64_t *a)
{
	return fpr_rint(blob_get(a));
}

int64_t
prec_tgt_floor(const uint64_t *a)
{
	return fpr_floor(blob_get(a));
}

int64_t
prec_tgt_trunc(const uint64_t *a)
{
	return fpr_trunc(blob_get(a));
}

uint64_t
prec_tgt_p63(const uint64_t *a)
{
	/* The target's fpr_p63. */
	return fpr_p63(blob_get(a));
}

void
prec_tgt_add(uint64_t *r, const uint64_t *a, const uint64_t *b)
{
	blob_put(r, fpr_add(blob_get(a), blob_get(b)));
}

void
prec_tgt_sub(uint64_t *r, const uint64_t *a, const uint64_t *b)
{
	blob_put(r, fpr_sub(blob_get(a), blob_get(b)));
}

void
prec_tgt_mul(uint64_t *r, const uint64_t *a, const uint64_t *b)
{
	blob_put(r, fpr_mul(blob_get(a), blob_get(b)));
}

void
prec_tgt_div(uint64_t *r, const uint64_t *a, const uint64_t *b)
{
	blob_put(r, fpr_div(blob_get(a), blob_get(b)));
}

void
prec_tgt_cmul(uint64_t *r_re, uint64_t *r_im,
	const uint64_t *a_re, const uint64_t *a_im,
	const uint64_t *b_re, const uint64_t *b_im)
{
	fpr d_re, d_im;

	FPC_MUL(d_re, d_im, blob_get(a_re), blob_get(a_im),
		blob_get(b_re), blob_get(b_im));
	blob_put(r_re, d_re);
	blob_put(r_im, d_im);
}

void
prec_tgt_neg(uint64_t *r, const uint64_t *a)
{
	blob_put(r, fpr_neg(blob_get(a)));
}

void
prec_tgt_half(uint64_t *r, const uint64_t *a)
{
	blob_put(r, fpr_half(blob_get(a)));
}

void
prec_tgt_double(uint64_t *r, const uint64_t *a)
{
	blob_put(r, fpr_double(blob_get(a)));
}

void
prec_tgt_sqr(uint64_t *r, const uint64_t *a)
{
	blob_put(r, fpr_sqr(blob_get(a)));
}

void
prec_tgt_inv(uint64_t *r, const uint64_t *a)
{
	blob_put(r, fpr_inv(blob_get(a)));
}

/* The target's FPR_RECIP_X / FPR_MUL_RECIP; rcp is a 24-byte blob. */
typedef char prec_tgt_rcp_size_check[
	(sizeof(fpr_rcp) <= 3 * sizeof(uint64_t)) ? 1 : -1];

void
prec_tgt_recip_x(uint64_t *r, const uint64_t *x)
{
	fpr xv = blob_get(x);
	fpr_rcp rcp;

	r[0] = 0;
	r[1] = 0;
	r[2] = 0;
	FPR_RECIP_X(rcp, xv);
	memcpy(r, &rcp, sizeof rcp);
}

void
prec_tgt_mul_recip(uint64_t *r, const uint64_t *a, const uint64_t *rc)
{
	fpr av = blob_get(a);
	fpr z;
	fpr_rcp rcp;

	memcpy(&rcp, rc, sizeof rcp);
	FPR_MUL_RECIP(z, av, rcp);
	blob_put(r, z);
}

#if FNDSA_SIGN_L10_64
void
prec_tgt_mul_recip_l10(uint64_t *r, const uint64_t *a, const uint64_t *rc)
{
	fpr av = blob_get(a);
	fpr_rcp rcp;
	l10_t d;

	memcpy(&rcp, rc, sizeof rcp);
	FPR_MUL_RECIP_P63(d, av, rcp);
	blob_put(r, l10_to_fpr(d));
}
#endif

#if FNDSA_SIGN_SQRT_FUSED
void
prec_tgt_sqrt_mul_inv_sigma(uint64_t *r, const uint64_t *x, unsigned logn)
{
	fpr xv = blob_get(x);
	fpr d;

	FPR_SQRT_MUL_INV_SIGMA(d, xv, logn);
	blob_put(r, d);
}
#endif

void
prec_tgt_mul_inv_q(uint64_t *r, const uint64_t *a)
{
	blob_put(r, FPR_MUL_INV_Q(blob_get(a)));
}

void
prec_tgt_sqrt(uint64_t *r, const uint64_t *a)
{
	blob_put(r, fpr_sqrt(blob_get(a)));
}

void
prec_tgt_div2e(uint64_t *r, const uint64_t *a, unsigned e)
{
	blob_put(r, fpr_div2e(blob_get(a), e));
}

void
prec_tgt_mul2e(uint64_t *r, const uint64_t *a, unsigned e)
{
	blob_put(r, fpr_mul2e(blob_get(a), e));
}

/* Exact conversion of a target value to MPFR. */
void
prec_tgt_to_mpfr(mpfr_t rop, const uint64_t *a)
{
#if FNDSA_FPR_REP == FNDSA_FPR_REP_IEEE
	/* fpr is the binary64 bit pattern in a uint64_t. */
	fpr x = blob_get(a);
	double d;
	memcpy(&d, &x, sizeof d);
	mpfr_set_d(rop, d, MPFR_RNDN);
#elif FNDSA_FPR_REP == FNDSA_FPR_REP_TW
	/* Triple-word: the value is the exact sum of three floats.
	   Each add is exact at rop's precision (>= ~160 bits). */
	fpr x = blob_get(a);
	mpfr_set_d(rop, (double)x.x[0], MPFR_RNDN);
	mpfr_add_d(rop, rop, (double)x.x[1], MPFR_RNDN);
	mpfr_add_d(rop, rop, (double)x.x[2], MPFR_RNDN);
#else
	/* Fixed point: value = I + w[1]*2^-32 + w[0]*2^-64, with I the
	   signed integer part. */
	fpr x = blob_get(a);
#if FXP_INT_BITS == 64
	int64_t ip = (int64_t)(((uint64_t)x.w[3] << 32) | x.w[2]);
#else
	int64_t ip = (int64_t)*(int32_t *)&x.w[2];
#endif
	mpfr_set_sj(rop, (intmax_t)ip, MPFR_RNDN);
	mpfr_mul_2si(rop, rop, 32, MPFR_RNDN);
	mpfr_add_ui(rop, rop, (unsigned long)x.w[1], MPFR_RNDN);
	mpfr_mul_2si(rop, rop, 32, MPFR_RNDN);
	mpfr_add_ui(rop, rop, (unsigned long)x.w[0], MPFR_RNDN);
	mpfr_mul_2si(rop, rop, -64, MPFR_RNDN);
#endif
}

#endif /* PREC_TGT_NONE */
