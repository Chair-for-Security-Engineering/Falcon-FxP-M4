/*
 * C reference fpr operations for fprtest, tw representation: the
 * constant-time triple_float C routines, which every tw build links.
 * There is no fpr_scaled on the board.
 */

#include "fpr/fprtest.h"

static fpr c_of32(int32_t i)           { return tw_from_int(i); }
static fpr c_add(fpr x, fpr y)         { return tw_sum_ct(x, y); }
static fpr c_sub(fpr x, fpr y)         { return tw_sub_ct(x, y); }
static void c_add_sub(fpr *a, fpr *b, fpr x, fpr y) { tw_add_sub_ct(x, y, a, b); }
static fpr c_mul(fpr x, fpr y)         { return tw_prod_fast_ct(x, y); }
static fpr c_sqr(fpr x)                { return tw_prod_fast_ct(x, x); }
static fpr c_div(fpr x, fpr y)         { return tw_div_fast_ct(x, y); }
static fpr c_inv(fpr x)                { return tw_reci_fast_ct(x); }
static fpr c_sqrt(fpr x)               { return tw_sqrt_fast_ct(x); }
static int64_t c_rint(fpr x)           { return (int64_t)tw_rint_fast(x); }
static int64_t c_floor(fpr x)          { return (int64_t)tw_floor_fast(x); }
static int64_t c_trunc(fpr x)          { return (int64_t)tw_trunc_full(x); }

const fprtest_ops fndsa_fprtest_cref = {
	.op_scaled = 0,
	.op_of32 = c_of32,
	.op_add = c_add,
	.op_sub = c_sub,
	.op_add_sub = c_add_sub,
	.op_mul = c_mul,
	.op_sqr = c_sqr,
	.op_div = c_div,
	.op_inv = c_inv,
	.op_sqrt = c_sqrt,
	.op_rint = c_rint,
	.op_floor = c_floor,
	.op_trunc = c_trunc,
};
