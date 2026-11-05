/*
 * Copyright 2025 NXP
 * SPDX-License-Identifier: MIT
 *
 * Out-of-line fpr operations for the 'tw' representation in C: wrappers
 * over the constant-time triple_float routines, plus a test-only
 * fpr_scaled. From TWFalcon's sign_fpr.c (MIT; see LICENSES/MIT-TWFalcon.txt).
 */

#include "sign_inner.h"

/* --- helpers for fpr_scaled (from sign_fpr.c) --------------------------- */

/* Like make() in fpr/ieee/sign_fpr.c, but sets e to -1076 when m = 0.
   m MUST be 0 or in [2^54, 2^55-1]. */
static inline uint64_t
make_z(uint64_t s, int32_t e, uint64_t m)
{
	uint32_t eu = (uint32_t)(e + 1076) & -(uint32_t)(m >> 54);
	uint64_t cc = (0xC8u >> ((unsigned)m & 7)) & 1;
	return (s << 63) + ((uint64_t)(eu) << 52) + (m >> 2) + cc;
}

/* Count of leading zeros in a 64-bit non-zero word. */
static inline uint32_t
lzcnt64_nonzero(uint64_t x)
{
	uint32_t x0 = (uint32_t)x;
	uint32_t x1 = (uint32_t)(x >> 32);
	uint32_t m = ~tbmask(x1 | -x1);
	x1 |= x0 & m;
	return lzcnt_nonzero(x1) + (m & 32);
}

/* Normalize m into [2^63, 2^64-1] while preserving m*2^e. */
#define NORM64(m, e)   do { \
		uint64_t norm64_m = (uint64_t)(m); \
		uint32_t norm64_c = lzcnt64_nonzero(norm64_m | 1); \
		(m) = fpr_ulsh(norm64_m, norm64_c); \
		(e) -= (int32_t)norm64_c; \
	} while (0)

/* --- out-of-line operations --------------------------------------------- */

/* Given integer i and scale sc, return i*2^sc (tests only), through a
   binary64 bit pattern. */
fpr
fpr_scaled(int64_t i, int sc)
{
	if (i == 0) {
		return FPR_ZERO;
	}

	/* Get sign and absolute value. */
	uint64_t s = (uint64_t)(i >> 63);
	uint64_t m = ((uint64_t)i ^ s) - s;

	/* Normalize m to [2^63, 2^64-1]. */
	NORM64(m, sc);

	/* Divide m by 2^9 to get it in [2^54, 2^55-1], with a sticky low bit. */
	sc += 9;
	m = (m | ((m & 0x1FF) + 0x1FF)) >> 9;

	uint64_t rz = make_z(s, sc, m);
	return create_from_u64(rz);
}

fpr
fpr_add(fpr x, fpr y)
{
	return tw_sum_ct(x, y);
}

fpr
fpr_mul(fpr x, fpr y)
{
	return tw_prod_fast_ct(x, y);
}

fpr
fpr_div(fpr x, fpr y)
{
	return tw_div_fast_ct(x, y);
}

fpr
fpr_sqrt(fpr x)
{
	return tw_sqrt_fast_ct(x);
}
