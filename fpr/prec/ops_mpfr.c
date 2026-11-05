/*
 * fpr operations for the 'prec' representation. Each op computes the
 * MPFR reference and, with a target embedded, the target value through
 * the target's own ops (target_glue.c). Errors are measured only at the
 * sampler checkpoint.
 *
 * fpr_rint/floor/trunc compute both integers, count disagreements, and
 * return the target's result (the reference's with PREC_MASTER=ref or
 * no target), so execution follows the target's path.
 *
 * As in TWFalcon's precision test, MPFR/GMP memory comes from the Boehm
 * GC and is never freed explicitly. mpfr limb pointers are interior
 * pointers, so GC_all_interior_pointers must stay enabled.
 */

#include "sign_inner.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gc.h>

#include "target_glue.h"
#include "prec_eval.h"

/* ==================================================================== */
/* Initialization and GC glue.                                          */

static int prec_ready = 0;
static int prec_alloc_uncollectable = 0;
static int prec_master_ref = 1;

static prec_eval_stats prec_stats;
static prec_eval_sink prec_sink = 0;
static void *prec_sink_ctx = 0;

static void *
prec_gc_malloc(size_t n)
{
	if (prec_alloc_uncollectable) {
		return GC_MALLOC_UNCOLLECTABLE(n);
	}
	return GC_MALLOC(n);
}

static void *
prec_gc_realloc(void *p, size_t old_n, size_t new_n)
{
	(void)old_n;
	return GC_REALLOC(p, new_n);
}

static void
prec_gc_free(void *p, size_t n)
{
	/* Deliberately nothing: the GC reclaims unreachable limbs. */
	(void)p;
	(void)n;
}

/* Defined in fpr_const_prec.c. */
void fndsa_prec_init_consts(void);

void
fndsa_prec_init(void)
{
	if (prec_ready) {
		return;
	}
	prec_ready = 1;
	GC_INIT();
	{
		/* Large initial heap so collections stay rare; override
		   with PREC_GC_HEAP_MB. */
		size_t mb = 256;
		const char *e = getenv("PREC_GC_HEAP_MB");
		if (e != NULL && *e != 0) {
			long v = strtol(e, NULL, 10);
			if (v > 0) {
				mb = (size_t)v;
			}
		}
		GC_expand_hp(mb << 20);
	}
	mp_set_memory_functions(&prec_gc_malloc, &prec_gc_realloc,
		&prec_gc_free);
	/* The target drives control flow unless PREC_MASTER=ref or there
	   is no target. */
	prec_master_ref = !prec_tgt_enabled;
	{
		const char *m = getenv("PREC_MASTER");
		if (m != NULL && strcmp(m, "ref") == 0) {
			prec_master_ref = 1;
		}
	}
	/* Constant tables are allocated uncollectable. */
	prec_alloc_uncollectable = 1;
	fndsa_prec_init_consts();
	prec_alloc_uncollectable = 0;
}

__attribute__((constructor))
static void
prec_ctor(void)
{
	fndsa_prec_init();
}

static inline void
prec_ensure(void)
{
	if (!prec_ready) {
		fndsa_prec_init();
	}
}

/* A fresh value with zero tgt; existing values are never changed. */
static inline prec_fpr
prec_new(void)
{
	prec_fpr r;

	prec_ensure();
	mpfr_init2(r.ref, PREC_REF_BITS);
	r.tgt[0] = 0;
	r.tgt[1] = 0;
	return r;
}

/* ==================================================================== */
/* Conversions and rounding.                                            */

prec_fpr
fndsa_prec_scaled(int64_t i, int sc)
{
	prec_fpr r = prec_new();

	mpfr_set_si_2exp(r.ref, i, sc, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_scaled(r.tgt, i, sc);
	}
	return r;
}

prec_fpr
fndsa_prec_of(int64_t i)
{
	prec_fpr r = prec_new();

	mpfr_set_si(r.ref, i, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_of(r.tgt, i);
	}
	return r;
}

prec_fpr
fndsa_prec_of32(int32_t i)
{
	prec_fpr r = prec_new();

	mpfr_set_si(r.ref, i, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_of32(r.tgt, i);
	}
	return r;
}

int64_t
fndsa_prec_rint(prec_fpr x)
{
	prec_ensure();
	/* Round to nearest, ties to even (matches ieee fpr_rint). */
	int64_t ir = mpfr_get_si(x.ref, MPFR_RNDN);
	if (prec_tgt_enabled) {
		int64_t it = prec_tgt_rint(x.tgt);
		prec_stats.rint_calls ++;
		if (it != ir) {
			prec_stats.rint_div ++;
		}
		if (!prec_master_ref) {
			return it;
		}
	}
	return ir;
}

int64_t
fndsa_prec_floor(prec_fpr x)
{
	prec_ensure();
	int64_t ir = mpfr_get_si(x.ref, MPFR_RNDD);
	if (prec_tgt_enabled) {
		int64_t it = prec_tgt_floor(x.tgt);
		prec_stats.floor_calls ++;
		if (it != ir) {
			prec_stats.floor_div ++;
		}
		if (!prec_master_ref) {
			return it;
		}
	}
	return ir;
}

int64_t
fndsa_prec_trunc(prec_fpr x)
{
	prec_ensure();
	int64_t ir = mpfr_get_si(x.ref, MPFR_RNDZ);
	if (prec_tgt_enabled) {
		int64_t it = prec_tgt_trunc(x.tgt);
		prec_stats.trunc_calls ++;
		if (it != ir) {
			prec_stats.trunc_div ++;
		}
		if (!prec_master_ref) {
			return it;
		}
	}
	return ir;
}

uint64_t
fndsa_prec_p63(prec_fpr x)
{
	mpfr_t t;
	uint64_t ir;

	prec_ensure();
	/* floor(x * 2^63), x in [0, 1). */
	mpfr_init2(t, PREC_REF_BITS);
	mpfr_mul_2si(t, x.ref, 63, MPFR_RNDN);
	ir = (uint64_t)mpfr_get_sj(t, MPFR_RNDD);
	if (prec_tgt_enabled) {
		uint64_t it = prec_tgt_p63(x.tgt);
		/* Own counters: ref and target often differ in the last
		   bit here, which would hide real fpr_trunc divergences. */
		prec_stats.p63_calls ++;
		if (it != ir) {
			prec_stats.p63_div ++;
		}
		if (!prec_master_ref) {
			return it;
		}
	}
	return ir;
}

/* ==================================================================== */
/* Arithmetic.                                                          */

prec_fpr
fndsa_prec_add(prec_fpr x, prec_fpr y)
{
	prec_fpr r = prec_new();

	mpfr_add(r.ref, x.ref, y.ref, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_add(r.tgt, x.tgt, y.tgt);
	}
	return r;
}

prec_fpr
fndsa_prec_sub(prec_fpr x, prec_fpr y)
{
	prec_fpr r = prec_new();

	mpfr_sub(r.ref, x.ref, y.ref, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_sub(r.tgt, x.tgt, y.tgt);
	}
	return r;
}

prec_fpr
fndsa_prec_neg(prec_fpr x)
{
	prec_fpr r = prec_new();

	mpfr_neg(r.ref, x.ref, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_neg(r.tgt, x.tgt);
	}
	return r;
}

prec_fpr
fndsa_prec_half(prec_fpr x)
{
	prec_fpr r = prec_new();

	mpfr_div_2si(r.ref, x.ref, 1, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_half(r.tgt, x.tgt);
	}
	return r;
}

prec_fpr
fndsa_prec_double(prec_fpr x)
{
	prec_fpr r = prec_new();

	mpfr_mul_2si(r.ref, x.ref, 1, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_double(r.tgt, x.tgt);
	}
	return r;
}

prec_fpr
fndsa_prec_mul(prec_fpr x, prec_fpr y)
{
	prec_fpr r = prec_new();

	mpfr_mul(r.ref, x.ref, y.ref, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_mul(r.tgt, x.tgt, y.tgt);
	}
	return r;
}

prec_fpr
fndsa_prec_sqr(prec_fpr x)
{
	prec_fpr r = prec_new();

	mpfr_sqr(r.ref, x.ref, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_sqr(r.tgt, x.tgt);
	}
	return r;
}

prec_fpr
fndsa_prec_div(prec_fpr x, prec_fpr y)
{
	prec_fpr r = prec_new();

	mpfr_div(r.ref, x.ref, y.ref, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_div(r.tgt, x.tgt, y.tgt);
	}
	return r;
}

prec_fpr
fndsa_prec_inv(prec_fpr x)
{
	prec_fpr r = prec_new();

	mpfr_si_div(r.ref, 1, x.ref, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_inv(r.tgt, x.tgt);
	}
	return r;
}

/* Complex product (see fpr_prec.h): ref at PREC_REF_BITS, tgt by the
   target's own FPC_MUL. */
void
fndsa_prec_fpc_mul(prec_fpr *d_re, prec_fpr *d_im,
	prec_fpr a_re, prec_fpr a_im, prec_fpr b_re, prec_fpr b_im)
{
	prec_fpr r_re = prec_new();
	prec_fpr r_im = prec_new();
	mpfr_t t;

	mpfr_init2(t, PREC_REF_BITS);
	mpfr_mul(r_re.ref, a_re.ref, b_re.ref, MPFR_RNDN);
	mpfr_mul(t, a_im.ref, b_im.ref, MPFR_RNDN);
	mpfr_sub(r_re.ref, r_re.ref, t, MPFR_RNDN);
	mpfr_mul(r_im.ref, a_re.ref, b_im.ref, MPFR_RNDN);
	mpfr_mul(t, a_im.ref, b_re.ref, MPFR_RNDN);
	mpfr_add(r_im.ref, r_im.ref, t, MPFR_RNDN);
	mpfr_clear(t);
	if (prec_tgt_enabled) {
		prec_tgt_cmul(r_re.tgt, r_im.tgt,
			a_re.tgt, a_im.tgt, b_re.tgt, b_im.tgt);
	}
	*d_re = r_re;
	*d_im = r_im;
}

/* Fused reciprocal seam (see fpr_prec.h). */
fpr_rcp
fndsa_prec_recip_x(prec_fpr x)
{
	fpr_rcp r;

	prec_ensure();
	mpfr_init2(r.ref, PREC_REF_BITS);
	mpfr_si_div(r.ref, 1, x.ref, MPFR_RNDN);
	r.tgt[0] = 0;
	r.tgt[1] = 0;
	r.tgt[2] = 0;
	if (prec_tgt_enabled) {
		prec_tgt_recip_x(r.tgt, x.tgt);
	}
	return r;
}

prec_fpr
fndsa_prec_mul_recip(prec_fpr a, fpr_rcp rcp)
{
	prec_fpr r = prec_new();

	mpfr_mul(r.ref, a.ref, rcp.ref, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_mul_recip(r.tgt, a.tgt, rcp.tgt);
	}
	return r;
}

#if FNDSA_SIGN_L10_64
prec_fpr
fndsa_prec_mul_recip_l10(prec_fpr a, fpr_rcp rcp)
{
	prec_fpr r = prec_new();

	mpfr_mul(r.ref, a.ref, rcp.ref, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_mul_recip_l10(r.tgt, a.tgt, rcp.tgt);
	}
	return r;
}
#endif

prec_fpr
fndsa_prec_sqrt(prec_fpr x)
{
	prec_fpr r = prec_new();

	mpfr_sqrt(r.ref, x.ref, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_sqrt(r.tgt, x.tgt);
	}
	return r;
}

#if defined(FNDSA_SIGN_SQRT_FUSED) && FNDSA_SIGN_SQRT_FUSED
prec_fpr
fndsa_prec_sqrt_mul_inv_sigma(prec_fpr x, unsigned logn)
{
	prec_fpr r = prec_new();

	mpfr_sqrt(r.ref, x.ref, MPFR_RNDN);
	mpfr_mul(r.ref, r.ref, INV_SIGMA[logn].f.ref, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_sqrt_mul_inv_sigma(r.tgt, x.tgt, logn);
	}
	return r;
}
#endif

prec_fpr
fndsa_prec_div2e(prec_fpr x, unsigned e)
{
	prec_fpr r = prec_new();

	mpfr_div_2si(r.ref, x.ref, (long)e, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_div2e(r.tgt, x.tgt, e);
	}
	return r;
}

prec_fpr
fndsa_prec_mul2e(prec_fpr x, unsigned e)
{
	prec_fpr r = prec_new();

	mpfr_mul_2si(r.ref, x.ref, (long)e, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_mul2e(r.tgt, x.tgt, e);
	}
	return r;
}

/* ==================================================================== */
/* Scalar constants (cached; built once at init).                       */

static prec_fpr c_zero, c_nzero, c_one;
static prec_fpr c_inv_2sqrsigma0, c_log2, c_inv_log2, c_inv_q, c_minus_inv_q;
static prec_fpr c_q2;

/* Called from fndsa_prec_init_consts(). */
void
fndsa_prec_init_scalars(void)
{
	mpfr_init2(c_zero.ref, PREC_REF_BITS);
	mpfr_set_zero(c_zero.ref, +1);
	mpfr_init2(c_nzero.ref, PREC_REF_BITS);
	mpfr_set_zero(c_nzero.ref, -1);
	mpfr_init2(c_one.ref, PREC_REF_BITS);
	mpfr_set_si(c_one.ref, 1, MPFR_RNDN);

	/* 1/(2*(1.8205^2)); 1.8205 = 18205/10000 exactly. */
	mpfr_init2(c_inv_2sqrsigma0.ref, PREC_REF_BITS);
	{
		mpfr_t t;
		mpfr_init2(t, PREC_REF_BITS + 64);
		mpfr_set_si(t, 18205, MPFR_RNDN);
		mpfr_div_si(t, t, 10000, MPFR_RNDN);
		mpfr_sqr(t, t, MPFR_RNDN);
		mpfr_mul_2si(t, t, 1, MPFR_RNDN);
		mpfr_si_div(c_inv_2sqrsigma0.ref, 1, t, MPFR_RNDN);
	}

	mpfr_init2(c_log2.ref, PREC_REF_BITS);
	mpfr_free_cache();
	mpfr_const_log2(c_log2.ref, MPFR_RNDN);
	mpfr_init2(c_inv_log2.ref, PREC_REF_BITS);
	{
		mpfr_t t;
		mpfr_init2(t, PREC_REF_BITS + 64);
		mpfr_free_cache();
		mpfr_const_log2(t, MPFR_RNDN);
		mpfr_si_div(c_inv_log2.ref, 1, t, MPFR_RNDN);
	}

	mpfr_init2(c_inv_q.ref, PREC_REF_BITS);
	mpfr_set_si(c_inv_q.ref, 1, MPFR_RNDN);
	mpfr_div_si(c_inv_q.ref, c_inv_q.ref, 12289, MPFR_RNDN);
	mpfr_init2(c_minus_inv_q.ref, PREC_REF_BITS);
	mpfr_neg(c_minus_inv_q.ref, c_inv_q.ref, MPFR_RNDN);

	/* q^2 = 151019521, an integer: exact. */
	mpfr_init2(c_q2.ref, PREC_REF_BITS);
	mpfr_set_ui(c_q2.ref, 12289UL * 12289, MPFR_RNDN);

	/* Target halves from the target's own compile-time macros. */
	prec_tgt_const_scalar(c_zero.tgt, PREC_TGT_C_ZERO);
	prec_tgt_const_scalar(c_nzero.tgt, PREC_TGT_C_NZERO);
	prec_tgt_const_scalar(c_one.tgt, PREC_TGT_C_ONE);
	prec_tgt_const_scalar(c_inv_2sqrsigma0.tgt,
		PREC_TGT_C_INV_2SQRSIGMA0);
	prec_tgt_const_scalar(c_log2.tgt, PREC_TGT_C_LOG2);
	prec_tgt_const_scalar(c_inv_log2.tgt, PREC_TGT_C_INV_LOG2);
	prec_tgt_const_scalar(c_inv_q.tgt, PREC_TGT_C_INV_Q);
	prec_tgt_const_scalar(c_minus_inv_q.tgt, PREC_TGT_C_MINUS_INV_Q);
	prec_tgt_const_scalar(c_q2.tgt, PREC_TGT_C_Q2);
}

prec_fpr
fndsa_prec_zero(void)
{
	prec_ensure();
	return c_zero;
}

prec_fpr
fndsa_prec_nzero(void)
{
	prec_ensure();
	return c_nzero;
}

prec_fpr
fndsa_prec_one(void)
{
	prec_ensure();
	return c_one;
}

prec_fpr
fndsa_prec_c_inv_2sqrsigma0(void)
{
	prec_ensure();
	return c_inv_2sqrsigma0;
}

prec_fpr
fndsa_prec_c_log2(void)
{
	prec_ensure();
	return c_log2;
}

prec_fpr
fndsa_prec_c_inv_log2(void)
{
	prec_ensure();
	return c_inv_log2;
}

prec_fpr
fndsa_prec_c_inv_q(void)
{
	prec_ensure();
	return c_inv_q;
}

prec_fpr
fndsa_prec_c_minus_inv_q(void)
{
	prec_ensure();
	return c_minus_inv_q;
}

prec_fpr
fndsa_prec_c_q2(void)
{
	prec_ensure();
	return c_q2;
}

/* Wide 1/q seam (see fpr_prec.h): exact 1/q on the reference side. */
prec_fpr
fndsa_prec_mul_inv_q(prec_fpr x)
{
	prec_fpr r = prec_new();

	mpfr_mul(r.ref, x.ref, c_inv_q.ref, MPFR_RNDN);
	if (prec_tgt_enabled) {
		prec_tgt_mul_inv_q(r.tgt, x.tgt);
	}
	return r;
}

/* ==================================================================== */
/* Sampler checkpoint, called at the top of sampler_next with mu and
   1/sigma: absolute error of mu, relative error of 1/sigma.            */

/* Records err into h and returns its code (see prec_eval.h). */
static int
prec_record_err(prec_eval_hist *h, mpfr_srcptr err)
{
	if (mpfr_zero_p(err)) {
		h->exact ++;
		return PREC_EVAL_E_NONE;
	}
	/* A double is enough: the largest errors dominate the RMSE. */
	double d = mpfr_get_d(err, MPFR_RNDN);
	h->sum_sq += d * d;
	long e = -(long)mpfr_get_exp(err);
	if (e <= 0) {
		h->huge ++;
		return (e < -9999) ? -9999 : (int)e;
	}
	if (e > PREC_EVAL_E_MAX) {
		h->tiny ++;
		return PREC_EVAL_E_MAX + 1;
	}
	h->hist[e] ++;
	return (int)e;
}

void
fndsa_prec_checkpoint_mu_isigma(prec_fpr mu, prec_fpr isigma)
{
	mpfr_t t_mu, t_is, err_mu, err_sig, err_comb;
	int e_mu, e_sig;

	prec_ensure();
	if (!prec_tgt_enabled) {
		return;
	}
	prec_stats.calls ++;

	mpfr_inits2(PREC_REF_BITS, t_mu, t_is, err_mu, err_sig, err_comb,
		(mpfr_ptr)0);
	prec_tgt_to_mpfr(t_mu, mu.tgt);
	prec_tgt_to_mpfr(t_is, isigma.tgt);

	/* mu: absolute error. */
	mpfr_sub(err_mu, t_mu, mu.ref, MPFR_RNDN);
	mpfr_abs(err_mu, err_mu, MPFR_RNDN);

	/* sigma: relative error for every target. The relative error of
	   1/sigma equals that of sigma to first order. */
	mpfr_sub(err_sig, t_is, isigma.ref, MPFR_RNDN);
	mpfr_abs(err_sig, err_sig, MPFR_RNDN);
	if (!mpfr_zero_p(isigma.ref)) {
		mpfr_t d;
		mpfr_init2(d, PREC_REF_BITS);
		mpfr_abs(d, isigma.ref, MPFR_RNDN);
		mpfr_div(err_sig, err_sig, d, MPFR_RNDN);
		mpfr_clear(d);
	}

	/* combined = delta_mu + delta_sigma as values. */
	mpfr_add(err_comb, err_mu, err_sig, MPFR_RNDN);

	e_mu = prec_record_err(&prec_stats.mu, err_mu);
	e_sig = prec_record_err(&prec_stats.sigma, err_sig);
	/* combined is accumulated but not reported. */
	(void)prec_record_err(&prec_stats.combined, err_comb);

	if (prec_sink != 0) {
		prec_sink(prec_sink_ctx, e_mu, e_sig);
	}

	mpfr_clears(t_mu, t_is, err_mu, err_sig, err_comb, (mpfr_ptr)0);
}

/* ==================================================================== */
/* Analysis API (prec_eval.h).                                          */

const prec_eval_stats *
fndsa_prec_stats(void)
{
	return &prec_stats;
}

void
fndsa_prec_stats_reset(void)
{
	memset(&prec_stats, 0, sizeof prec_stats);
}

void
fndsa_prec_set_sink(prec_eval_sink fn, void *ctx)
{
	prec_sink = fn;
	prec_sink_ctx = ctx;
}

int
fndsa_prec_metric(void)
{
	return prec_tgt_metric;
}

int
fndsa_prec_frac_bits(void)
{
	return prec_tgt_frac_bits;
}

const char *
fndsa_prec_target(void)
{
	return prec_tgt_name;
}

int
fndsa_prec_master_is_ref(void)
{
	prec_ensure();
	return prec_master_ref;
}

int
fndsa_prec_ref_bits(void)
{
	return PREC_REF_BITS;
}
