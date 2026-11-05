/*
 * C reference fpr operations for fprtest, ieee representation. Includes
 * sign_fpr.c again with FNDSA_ASM_CORTEXM4=0 and renamed symbols, so the
 * C versions sit next to the assembly in one binary.
 */

#undef FNDSA_ASM_CORTEXM4
#define FNDSA_ASM_CORTEXM4   0
#undef FNDSA_SHAKE_HUANG_BUF
#define FNDSA_SHAKE_HUANG_BUF   0

#define fndsa_fpr_scaled   fndsa_cref_fpr_scaled
#define fndsa_fpr_add      fndsa_cref_fpr_add
#define fndsa_fpr_mul      fndsa_cref_fpr_mul
#define fndsa_fpr_div      fndsa_cref_fpr_div
#define fndsa_fpr_sqrt     fndsa_cref_fpr_sqrt

#include "sign_fpr.c"

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

/* rint/floor/trunc are the same inline code in both builds. */
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
	.op_rint = 0,
	.op_floor = 0,
	.op_trunc = 0,
};
