/*
 * Per-operation test table for fxp32.64, checked against the MPFR oracle.
 * Exact ops and the multiplies must match bit for bit (ENT, ENT_MUL);
 * inv, div, sqrt and the fused chains pass within 1 or 2 ulp OR a
 * relative bound (ENT_TOL).
 *
 * dom_int_bits bounds the operands (|value| < 2^dom); results that
 * overflow the format are discarded. dom_min_bits is a lower bound on
 * the raw bit length of a divisor or radicand, which keeps the
 * intermediate reciprocal in range.
 */

#include "fpr/test_fpr.h"

const char fxp_test_backend[] = "fxp32.64";

#define ENT(nm, knd, dom) \
	{ .name = (nm), .kind = (knd), .enabled = 1, .rnd = FXT_RND_N, \
	  .ulp_bound = 0, .dom_int_bits = (dom), .e_max = 32, \
	  .un = 0, .bin = 0, .rel_bits = 0, .dom_min_bits = 0 }

/* Bit-exact multiply, rounded as the FNDSA_FXP_TIES_UP build does. */
#if FNDSA_FXP_TIES_UP
#define ENT_MUL(nm, knd, dom) \
	{ .name = (nm), .kind = (knd), .enabled = 1, .rnd = FXT_RND_NU, \
	  .ulp_bound = 0, .dom_int_bits = (dom), .e_max = 32, \
	  .un = 0, .bin = 0, .rel_bits = 0, .dom_min_bits = 0 }
#else
#define ENT_MUL   ENT
#endif

/* Pass within 'ulp' ulp OR relative 2^-rel, for operands of raw bit
   length in [dmin, dom+64]. */
#define ENT_TOL(nm, knd, dom, ulp, rel, dmin) \
	{ .name = (nm), .kind = (knd), .enabled = 1, .rnd = FXT_RND_N, \
	  .ulp_bound = (ulp), .dom_int_bits = (dom), .e_max = 32, \
	  .un = 0, .bin = 0, .rel_bits = (rel), .dom_min_bits = (dmin) }

#if FNDSA_SIGN_RECIP_FUSED
/* FPR_RECIP_X + FPR_MUL_RECIP: faithful for small results, relative
   error for large ones. */
static fpr
mul_recip_chain(fpr a, fpr x)
{
	fpr_rcp rcp;
	fpr z;

	FPR_RECIP_X(rcp, x);
	FPR_MUL_RECIP(z, a, rcp);
	return z;
}
#define ENT_MUL_RECIP \
	{ .name = "mul_recip", .kind = FXT_MUL_RECIP, .enabled = 1, \
	  .rnd = FXT_RND_N, .ulp_bound = 1, .dom_int_bits = 15, \
	  .e_max = 32, .un = 0, .bin = mul_recip_chain, .rel_bits = 58, \
	  .dom_min_bits = 65 },

/* FPR_MUL_RECIP_P63, widened to an fpr. Rounds at 2^-63, so the bounds
   are one bit coarser than mul_recip's. */
static fpr
mul_recip_p63_chain(fpr a, fpr x)
{
	fpr_rcp rcp;
	uint64_t r;
	fpr z;

	FPR_RECIP_X(rcp, x);
	FPR_MUL_RECIP_P63(r, a, rcp);
	z.w[0] = (uint32_t)(r << 1);
	z.w[1] = (uint32_t)(r >> 31);
	z.w[2] = (uint32_t)((int64_t)r >> 63);
	return z;
}
#define ENT_MUL_RECIP_P63 \
	{ .name = "mul_recip_p63", .kind = FXT_MUL_RECIP_P63, .enabled = 1, \
	  .rnd = FXT_RND_N, .ulp_bound = 2, .dom_int_bits = 15, \
	  .e_max = 32, .un = 0, .bin = mul_recip_p63_chain, .rel_bits = 57, \
	  .dom_min_bits = 65 },
#else
#define ENT_MUL_RECIP
#define ENT_MUL_RECIP_P63
#endif

#if FNDSA_SIGN_SQRT_FUSED && FNDSA_KGEN_GAMMA
#define ENT_SQRT_SIGMA \
	{ .name = "sqrt_inv_sigma", .kind = FXT_SQRT_INV_SIGMA, \
	  .enabled = 1, .rnd = FXT_RND_N, .ulp_bound = 2, \
	  .dom_int_bits = 31, .e_max = 32, .un = 0, .bin = 0, \
	  .rel_bits = 44, .dom_min_bits = 4 },
#else
#define ENT_SQRT_SIGMA
#endif

const fxp_test_entry fxp_test_table[] = {
	/* conversions */
	ENT("fpr_of",       FXT_OF,      31),
	ENT("fpr_of32",     FXT_OF32,    31),
	ENT("fpr_scaled",   FXT_SCALED,  31),
	ENT("fpr_rint",     FXT_RINT,    31),
	ENT("fpr_floor",    FXT_FLOOR,   31),
	ENT("fpr_trunc",    FXT_TRUNC,   31),
	/* floor(x * 2^63) for the sampler */
	ENT("fpr_p63",      FXT_P63,      0),
	/* exact arithmetic */
	ENT("fpr_add",      FXT_ADD,     30),
	ENT("fpr_sub",      FXT_SUB,     30),
	ENT("FPR_ADD_SUB",  FXT_ADD_SUB, 30),
	ENT("fpr_neg",      FXT_NEG,     31),
	ENT("fpr_half",     FXT_HALF,    31),
	ENT("fpr_double",   FXT_DOUBLE,  30),
	ENT("fpr_mul2e",    FXT_MUL2E,   31),
	ENT("fpr_div2e",    FXT_DIV2E,   31),
	/* rounded arithmetic (bit-exact) */
	ENT_MUL("fpr_mul",  FXT_MUL,     15),
	ENT_MUL("fpr_sqr",  FXT_SQR,     15),
	/* x/q via the Q77 word, exact against the product with that word */
	ENT_MUL("mul_inv_q", FXT_MUL_INV_Q, 31),
	/* reciprocal, operand >= 1 */
	ENT_TOL("fpr_inv",  FXT_INV,     31, 1, 60, 65),
	/* x * (1/y), divisor >= 1 */
	ENT_TOL("fpr_div",  FXT_DIV,     15, 2, 46, 65),
	ENT_MUL_RECIP
	ENT_MUL_RECIP_P63
	/* x * rsqrt(x), rounded once */
	ENT_TOL("fpr_sqrt", FXT_SQRT,    31, 2, 44, 4),
	ENT_SQRT_SIGMA
#if FNDSA_FXP_CMUL_FUSED
	ENT_MUL("FPC_MUL",  FXT_FPC_MUL_FUSED, 15),
#else
	ENT_MUL("FPC_MUL",  FXT_FPC_MUL, 15),
#endif
};

const size_t fxp_test_table_len =
	sizeof fxp_test_table / sizeof fxp_test_table[0];
