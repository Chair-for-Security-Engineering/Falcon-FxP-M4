#ifndef FNDSA_TEST_FPR_H__
#define FNDSA_TEST_FPR_H__

/*
 * Per-operation fpr test table (fxp32.64 only). The backend lists its
 * operations in fpr/fxp32.64/fxp_test_table.c; test_fpr.c checks them
 * against MPFR (fpr/fxp_mpfr.h). Disabled entries report SKIP.
 */

#include "sign_inner.h"

/* Operation kinds; each has its own operands, default op and oracle. */
typedef enum {
	FXT_ADD,        /* fpr_add(x, y) */
	FXT_SUB,        /* fpr_sub(x, y) */
	FXT_MUL,        /* fpr_mul(x, y) */
	FXT_SQR,        /* fpr_sqr(x) */
	FXT_DIV,        /* fpr_div(x, y), y != 0 */
	FXT_SQRT,       /* fpr_sqrt(x), x >= 0 */
	FXT_INV,        /* fpr_inv(x), x != 0 */
	FXT_NEG,        /* fpr_neg(x) */
	FXT_HALF,       /* fpr_half(x) */
	FXT_DOUBLE,     /* fpr_double(x) */
	FXT_MUL2E,      /* fpr_mul2e(x, e), 0 <= e <= e_max */
	FXT_DIV2E,      /* fpr_div2e(x, e), 0 <= e <= e_max */
	FXT_OF,         /* fpr_of(i) */
	FXT_OF32,       /* fpr_of32(i) */
	FXT_SCALED,     /* fpr_scaled(i, sc) */
	FXT_RINT,       /* fpr_rint(x), as int64 */
	FXT_FLOOR,      /* fpr_floor(x), as int64 */
	FXT_TRUNC,      /* fpr_trunc(x), as int64 */
	FXT_ADD_SUB,    /* FPR_ADD_SUB(a, b, x, y) */
	FXT_FPC_MUL,    /* FPC_MUL (composition of grid-rounded steps) */
	FXT_FPC_MUL_FUSED, /* FPC_MUL, exact products, one rounding per
	                   component (fxp32.64) */
	FXT_MUL_RECIP,  /* FPR_RECIP_X + FPR_MUL_RECIP: a/b, one rounding */
	FXT_MUL_RECIP_P63, /* FPR_RECIP_X + FPR_MUL_RECIP_P63: a/b rounded
	                   once to Q0.63, |a/b| < 1 - 2^-64 */
	FXT_SQRT_INV_SIGMA, /* FPR_SQRT_MUL_INV_SIGMA: sqrt(x) times
	                   INV_SIGMA[logn], rounded once */
	FXT_MUL_INV_Q,  /* FPR_MUL_INV_Q(x): x times the 1/q word, rounded
	                   once */
	FXT_P63         /* fpr_p63(x) = floor(x * 2^63), x in [0, 1); use
	                   dom_int_bits 0 */
} fxp_test_kind;

/* Rounding the oracle applies onto the grid: FXT_RND_N, or FXT_RND_NU
   for the multiplies of a ties-up build. Integer results ignore it. */
typedef enum {
	FXT_RND_N,   /* to nearest, ties to even (MPFR_RNDN) */
	FXT_RND_Z,   /* toward zero (MPFR_RNDZ) */
	FXT_RND_D,   /* toward -inf (MPFR_RNDD) */
	FXT_RND_U,   /* toward +inf (MPFR_RNDU) */
	FXT_RND_NU   /* to nearest, ties toward +inf (FNDSA_FXP_TIES_UP) */
} fxp_test_rnd;

typedef struct {
	const char *name;    /* display name, e.g. "fpr_mul" or "mul_recip" */
	fxp_test_kind kind;
	int enabled;         /* 0 = not implemented, SKIP */
	fxp_test_rnd rnd;    /* grid rounding the oracle applies */
	unsigned ulp_bound;  /* allowed ulp distance from the oracle; 0 = exact */
	int dom_int_bits;    /* operands drawn with |value| < 2^dom_int_bits */
	int e_max;           /* FXT_MUL2E/DIV2E: max exponent; others: ignored */
	/* Optional replacement for the default op; NULL = default. */
	fpr (*un)(fpr);
	fpr (*bin)(fpr, fpr);
	/* If nonzero, also pass when |got - exp| * 2^rel_bits <= |exp|. */
	unsigned rel_bits;
	/* Minimum raw bit length of a divisor or radicand (0 = none);
	   smaller operands do not occur in signing. */
	int dom_min_bits;
} fxp_test_entry;

/* Supplied by fpr/<rep>/fxp_test_table.c. */
extern const char fxp_test_backend[];
extern const fxp_test_entry fxp_test_table[];
extern const size_t fxp_test_table_len;

#endif
