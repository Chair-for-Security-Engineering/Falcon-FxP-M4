#ifndef FNDSA_PREC_TARGET_GLUE_H__
#define FNDSA_PREC_TARGET_GLUE_H__

/*
 * Interface between the prec backend (ops_mpfr.c) and the embedded
 * target (target_glue.c, compiled with the target's representation, so
 * its own ops and tables are used unchanged). Target values cross as
 * opaque 16-byte blobs, so this header has no fpr type.
 */

#include <stdint.h>
#include <mpfr.h>

/* Whether the target is fixed-point (labels the summary only). */
#define PREC_TGT_METRIC_FLOAT     0   /* ieee, tw */
#define PREC_TGT_METRIC_FXP       1   /* fxp32.64 */

extern const int prec_tgt_enabled;      /* 0 for PREC_TARGET=none */
extern const int prec_tgt_metric;       /* PREC_TGT_METRIC_* */
extern const int prec_tgt_frac_bits;    /* FXP_FRAC_BITS, 0 if float */
extern const char *const prec_tgt_name; /* "ieee", "tw", ... */

/* Constructors. */
void prec_tgt_scaled(uint64_t *r, int64_t i, int sc);
void prec_tgt_of(uint64_t *r, int64_t i);
void prec_tgt_of32(uint64_t *r, int32_t i);

/* Scalar constants (the target header's compile-time macros). */
#define PREC_TGT_C_ZERO             0
#define PREC_TGT_C_NZERO            1
#define PREC_TGT_C_ONE              2
#define PREC_TGT_C_INV_2SQRSIGMA0   3
#define PREC_TGT_C_LOG2             4
#define PREC_TGT_C_INV_LOG2         5
#define PREC_TGT_C_INV_Q            6
#define PREC_TGT_C_MINUS_INV_Q      7
#define PREC_TGT_C_Q2               8
void prec_tgt_const_scalar(uint64_t *r, int which);

/* Table entries (the target's own tables). */
void prec_tgt_gm(uint64_t *r, int i);
void prec_tgt_rnj(uint64_t *r, int i);
void prec_tgt_inv_sigma(uint64_t *r, unsigned logn);
void prec_tgt_sigma_min(uint64_t *r, unsigned logn);

/* Rounding / integer conversions. */
int64_t prec_tgt_rint(const uint64_t *a);
int64_t prec_tgt_floor(const uint64_t *a);
int64_t prec_tgt_trunc(const uint64_t *a);
/* floor(x * 2^63) for x in [0, 1), with the target's own fpr_p63. */
uint64_t prec_tgt_p63(const uint64_t *a);

/* Arithmetic. */
void prec_tgt_add(uint64_t *r, const uint64_t *a, const uint64_t *b);
void prec_tgt_sub(uint64_t *r, const uint64_t *a, const uint64_t *b);
void prec_tgt_mul(uint64_t *r, const uint64_t *a, const uint64_t *b);
void prec_tgt_div(uint64_t *r, const uint64_t *a, const uint64_t *b);
void prec_tgt_neg(uint64_t *r, const uint64_t *a);
void prec_tgt_half(uint64_t *r, const uint64_t *a);
void prec_tgt_double(uint64_t *r, const uint64_t *a);
void prec_tgt_sqr(uint64_t *r, const uint64_t *a);
/* The target's FPC_MUL (fused on fxp32.64 unless FNDSA_FXP_CMUL_FUSED=0). */
void prec_tgt_cmul(uint64_t *r_re, uint64_t *r_im,
	const uint64_t *a_re, const uint64_t *a_im,
	const uint64_t *b_re, const uint64_t *b_im);
void prec_tgt_inv(uint64_t *r, const uint64_t *a);
/* Fused reciprocal seam: r is a 24-byte (3-word) rcp blob. */
void prec_tgt_recip_x(uint64_t *r, const uint64_t *x);
void prec_tgt_mul_recip(uint64_t *r, const uint64_t *a,
	const uint64_t *rc);
/* Same, rounded once onto the target's l10 grid and widened back. */
void prec_tgt_mul_recip_l10(uint64_t *r, const uint64_t *a,
	const uint64_t *rc);
/* sqrt(x)*INV_SIGMA[logn] with a single rounding. */
void prec_tgt_sqrt_mul_inv_sigma(uint64_t *r, const uint64_t *x,
	unsigned logn);
/* x/q through the target's FPR_MUL_INV_Q. */
void prec_tgt_mul_inv_q(uint64_t *r, const uint64_t *a);
void prec_tgt_sqrt(uint64_t *r, const uint64_t *a);
void prec_tgt_div2e(uint64_t *r, const uint64_t *a, unsigned e);
void prec_tgt_mul2e(uint64_t *r, const uint64_t *a, unsigned e);

/* Exact conversion of a target value to MPFR (rop needs about 160
   bits). */
void prec_tgt_to_mpfr(mpfr_t rop, const uint64_t *a);

#endif
