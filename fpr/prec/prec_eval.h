#ifndef FNDSA_PREC_EVAL_H__
#define FNDSA_PREC_EVAL_H__

/*
 * What the evaluation harness (eval_precision.c) reads from the prec
 * backend.
 *
 * An error 0.xxx * 2^-e is recorded as the exponent e, so larger e means
 * smaller error. Exact matches count in .exact, errors >= 0.5 in .huge,
 * e > PREC_EVAL_E_MAX in .tiny. .sum_sq adds the squared errors of all
 * checkpoints, for the RMSE.
 */

#include <stdint.h>

#define PREC_EVAL_E_MAX   400

/* Sink code for an exact match (err == 0): "no error", never a
   worst case. */
#define PREC_EVAL_E_NONE  INT32_MAX

typedef struct {
	uint64_t exact;                       /* error == 0 */
	uint64_t huge;                        /* error >= 0.5 */
	uint64_t tiny;                        /* e > PREC_EVAL_E_MAX */
	double sum_sq;                        /* sum of squared errors */
	uint64_t hist[PREC_EVAL_E_MAX + 1];   /* count per e (index 0 unused) */
} prec_eval_hist;

typedef struct {
	/* sampler checkpoint (mu / 1/sigma at sampler_next entry) */
	uint64_t calls;
	prec_eval_hist mu, sigma, combined;
	/* rounding-divergence sanity counters (master-path validity) */
	uint64_t rint_calls, rint_div;
	uint64_t floor_calls, floor_div;
	uint64_t trunc_calls, trunc_div;
	uint64_t p63_calls, p63_div;   /* floor(x*2^63) in the sampler */
} prec_eval_stats;

/* Idempotent backend init (also runs from a constructor). */
void fndsa_prec_init(void);

/* Accumulated statistics (live; reset with _reset). */
const prec_eval_stats *fndsa_prec_stats(void);
void fndsa_prec_stats_reset(void);

/* Called once per sampler call with the mu and sigma error codes:
   PREC_EVAL_E_NONE for exact, e <= 0 (clamped at -9999) for huge,
   PREC_EVAL_E_MAX + 1 for tiny. */
typedef void (*prec_eval_sink)(void *ctx, int e_mu, int e_sigma);
void fndsa_prec_set_sink(prec_eval_sink fn, void *ctx);

/* Build/runtime configuration of this binary. */
int fndsa_prec_metric(void);          /* PREC_TGT_METRIC_* (target_glue.h) */
int fndsa_prec_frac_bits(void);       /* fxp targets: FXP_FRAC_BITS */
const char *fndsa_prec_target(void);  /* "none", "ieee", "tw", ... */
int fndsa_prec_master_is_ref(void);   /* 1 if ref drives control flow */
int fndsa_prec_ref_bits(void);        /* PREC_REF_BITS */

#endif
