/*
 * Per-operation fpr tests for fxp32.64. Each enabled entry of
 * fpr/fxp32.64/fxp_test_table.c runs on edge and seeded random vectors; the
 * result must equal the MPFR value rounded onto the fixed-point grid, or
 * lie within the entry's ulp or relative bound. Vectors whose exact result
 * does not fit the format are discarded.
 * Disabled entries report SKIP. The harness first checks itself with
 * GMP integer versions of add/sub/neg/mul/div/sqrt.
 *
 * Usage: ./test_fpr [nvec] [seed]      (default 10000 vectors, per op)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fpr/test_fpr.h"
#include "fpr/fxp_mpfr.h"

static unsigned long g_nvec = 10000;
static const char *g_seed = "test_fpr";

/* ==================================================================== */
/* Random vector generation (SHAKE256, deterministic per entry).        */

static uint32_t
rnd_u32(shake_context *sc)
{
	uint8_t b[4];

	shake_extract(sc, b, sizeof b);
	return (uint32_t)b[0] | ((uint32_t)b[1] << 8)
		| ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

/* Random raw value with |value| < 2^dom_int_bits, log-uniform over bit
   lengths so all magnitudes are covered. Occasionally zero. */
static void
gen_raw(shake_context *sc, mpz_t z, int dom_int_bits)
{
	unsigned span, nbits, nbytes;
	uint8_t buf[64];

	span = (unsigned)dom_int_bits + FXP_FRAC_BITS;
	if ((rnd_u32(sc) & 0x7F) == 0) {
		mpz_set_ui(z, 0);
		return;
	}
	nbits = 1 + rnd_u32(sc) % span;
	nbytes = (nbits + 7) / 8;
	shake_extract(sc, buf, nbytes);
	mpz_import(z, nbytes, 1, 1, 0, 0, buf);
	mpz_fdiv_r_2exp(z, z, nbits);
	mpz_setbit(z, nbits - 1);
	if (rnd_u32(sc) & 1) {
		mpz_neg(z, z);
	}
}

static fpr
gen_fpr(shake_context *sc, mpz_t tmp, int dom_int_bits)
{
	gen_raw(sc, tmp, dom_int_bits);
	return fxp_from_mpz(tmp);
}

/* Random integer with |i| < 2^bits (bits <= 62), log-uniform. */
static int64_t
gen_i64(shake_context *sc, int bits)
{
	unsigned nbits;
	uint64_t v;
	uint8_t b[8];

	if ((rnd_u32(sc) & 0x3F) == 0) {
		return 0;
	}
	nbits = 1 + rnd_u32(sc) % (unsigned)bits;
	shake_extract(sc, b, sizeof b);
	v = 0;
	for (int i = 0; i < 8; i++) {
		v |= (uint64_t)b[i] << (8 * i);
	}
	v &= ((uint64_t)1 << nbits) - 1;
	v |= (uint64_t)1 << (nbits - 1);
	return (rnd_u32(sc) & 1) ? -(int64_t)v : (int64_t)v;
}

/* ==================================================================== */
/* Edge vectors: interesting raw values within the entry's domain.      */

#define MAX_EDGES   16

static mpz_t edge_pool[MAX_EDGES];
static mpfr_t va, vb, vr, vt, vu, vw;   /* shared temporaries */
static mpz_t za, zb, zr, zt;

#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264
/* Input ranges of the fused leaf kernels: reciprocal d in
   [q/55^2, q*55^2], rsqrt d in [q/110^2, q*110^2]. Random vectors use
   whole-number bounds; the exact endpoints are in the edge tables. */
#define RECIP_ALPHA55_LO       5u
#define RECIP_ALPHA55_HI       37174225u
#define RSQRT_ALPHA110_LO      2u
#define RSQRT_ALPHA110_HI      148696900u

/* Rotates through the leaf kernels' normalization shifts. */
static unsigned long g_leaf_bucket_index;

static fpr
fpr_from_q32_64_words(uint32_t w2, uint32_t w1, uint32_t w0)
{
	mpz_set_ui(za, w2);
	mpz_mul_2exp(za, za, FXP_FRAC_BITS);
	mpz_set_ui(zr, w1);
	mpz_mul_2exp(zr, zr, 32);
	mpz_add(za, za, zr);
	mpz_add_ui(za, za, w0);
	return fxp_from_mpz(za);
}

static fpr
fpr_from_integer_with_fraction(shake_context *sc, uint32_t hi)
{
	mpz_set_ui(za, hi);
	mpz_mul_2exp(za, za, FXP_FRAC_BITS);
	mpz_set_ui(zr, rnd_u32(sc));
	mpz_mul_2exp(zr, zr, 32);
	mpz_add(za, za, zr);
	mpz_add_ui(za, za, rnd_u32(sc));
	return fxp_from_mpz(za);
}

static fpr
gen_recip_alpha55(shake_context *sc)
{
	unsigned k = 6 + (unsigned)(g_leaf_bucket_index++ % 24);
	uint32_t lo = (uint32_t)1 << (31 - k);
	uint32_t hi = (uint32_t)1 << (32 - k);

	if (lo < RECIP_ALPHA55_LO) {
		lo = RECIP_ALPHA55_LO;
	}
	if (hi > RECIP_ALPHA55_HI) {
		hi = RECIP_ALPHA55_HI;
	}
	return fpr_from_integer_with_fraction(sc, lo + rnd_u32(sc) % (hi - lo));
}

static fpr
gen_rsqrt_alpha110(shake_context *sc)
{
	unsigned t = 4 + 2 * (unsigned)(g_leaf_bucket_index++ % 14);
	uint32_t lo = (uint32_t)1 << (30 - t);
	uint32_t hi = (uint32_t)1 << (32 - t);

	if (lo < RSQRT_ALPHA110_LO) {
		lo = RSQRT_ALPHA110_LO;
	}
	if (hi > RSQRT_ALPHA110_HI) {
		hi = RSQRT_ALPHA110_HI;
	}
	return fpr_from_integer_with_fraction(sc, lo + rnd_u32(sc) % (hi - lo));
}

static fpr
recip_alpha55_edge(size_t i)
{
	static const uint32_t v[][3] = {
		/* d = 4, the fused reciprocal's smallest divisor (k = 29). */
		{ 4u,        0x00000000, 0x00000000 },
		/* ceil_Q32.64(q / 55^2), the alpha-55 lower endpoint. */
		{ 4u,        0x0FFEA55D, 0x019BA18F },
		/* First whole value strictly above the alpha-55 lower endpoint. */
		{ 5u,        0x00000000, 0x00000000 },
		/* q / 16, the n=512 alpha_hybrid lower endpoint. */
		{ 768u,      0x10000000, 0x00000000 },
		/* 16q, the n=512 upper endpoint. */
		{ 196624u,   0x00000000, 0x00000000 },
		/* q^2 / 2^19, the n=1024 alpha_hybrid lower endpoint (exact). */
		{ 288u,      0x0C002000, 0x00000000 },
		/* 2^19 - 2^-64, just below the n=1024 upper endpoint. */
		{ 524287u,   0xFFFFFFFF, 0xFFFFFFFF },
		/* q * 55^2, the alpha-55 upper endpoint. */
		{ 37174225u, 0x00000000, 0x00000000 }
	};
	const uint32_t *w = v[i % (sizeof v / sizeof v[0])];

	return fpr_from_q32_64_words(w[0], w[1], w[2]);
}

static fpr
rsqrt_alpha110_edge(size_t i)
{
	static const uint32_t v[][3] = {
		/* d = 1, the rsqrt kernel's smallest radicand (t = 30). */
		{ 1u,         0x00000000, 0x00000000 },
		/* ceil_Q32.64(q / 110^2), the alpha-110 lower endpoint. */
		{ 1u,         0x03FFA957, 0x4066E864 },
		/* First whole value strictly above the alpha-110 lower endpoint. */
		{ 2u,         0x00000000, 0x00000000 },
		/* q / 16, the n=512 alpha_hybrid lower endpoint. */
		{ 768u,       0x10000000, 0x00000000 },
		/* 16q, the n=512 upper endpoint. */
		{ 196624u,    0x00000000, 0x00000000 },
		/* q^2 / 2^19, the n=1024 alpha_hybrid lower endpoint (exact). */
		{ 288u,       0x0C002000, 0x00000000 },
		/* 2^19 - 2^-64, just below the n=1024 upper endpoint. */
		{ 524287u,    0xFFFFFFFF, 0xFFFFFFFF },
		/* q * 110^2, the alpha-110 upper endpoint. */
		{ 148696900u, 0x00000000, 0x00000000 }
	};
	const uint32_t *w = v[i % (sizeof v / sizeof v[0])];

	return fpr_from_q32_64_words(w[0], w[1], w[2]);
}
#endif

static size_t
build_edges(int dom_int_bits)
{
	mpz_t lim;
	size_t n = 0, i, m;

	mpz_init(lim);
	mpz_set_ui(lim, 0);
	mpz_setbit(lim, (mp_bitcnt_t)(dom_int_bits + FXP_FRAC_BITS));

	/* candidates, as raw integers */
	mpz_set_ui(edge_pool[n++], 0);                          /* 0 */
	mpz_set_ui(edge_pool[n++], 1);                          /* 1 ulp */
	mpz_set_si(edge_pool[n++], -1);
	mpz_set_ui(edge_pool[n++], 3);                          /* odd */
	mpz_set_ui(edge_pool[n], 0);                            /* 0.5 */
	mpz_setbit(edge_pool[n++], FXP_FRAC_BITS - 1);
	mpz_set_ui(edge_pool[n], 0);                            /* 1.0 */
	mpz_setbit(edge_pool[n++], FXP_FRAC_BITS);
	mpz_neg(edge_pool[n], edge_pool[n - 1]);                /* -1.0 */
	n++;
	mpz_add_ui(edge_pool[n], edge_pool[n - 2], 1);          /* 1.0 + ulp */
	n++;
	mpz_sub_ui(edge_pool[n], lim, 1);                       /* dom max */
	n++;
	mpz_neg(edge_pool[n], edge_pool[n - 1]);                /* dom min */
	n++;
	mpz_fdiv_q_2exp(edge_pool[n], lim, 1);                  /* dom max / 2 */
	mpz_sub_ui(edge_pool[n], edge_pool[n], 1);
	n++;

	/* keep only |raw| < 2^(dom+frac) */
	m = 0;
	for (i = 0; i < n; i++) {
		if (mpz_cmpabs(edge_pool[i], lim) < 0) {
			if (m != i) {
				mpz_set(edge_pool[m], edge_pool[i]);
			}
			m++;
		}
	}
	mpz_clear(lim);
	return m;
}

/* ==================================================================== */
/* Integer reference ops (GMP) for the harness self-check.             */

/* q <- x / 2^s, round to nearest, ties to even (on the real value). */
static void
mpz_rshift_rn(mpz_t q, const mpz_t x, unsigned s)
{
	mpz_t r, half;
	int c;

	mpz_init(r);
	mpz_init(half);
	mpz_fdiv_q_2exp(q, x, s);       /* floor */
	mpz_fdiv_r_2exp(r, x, s);       /* nonnegative remainder */
	mpz_setbit(half, s - 1);
	c = mpz_cmp(r, half);
	if (c > 0 || (c == 0 && mpz_odd_p(q))) {
		mpz_add_ui(q, q, 1);
	}
	mpz_clear(r);
	mpz_clear(half);
}

#if defined FNDSA_HAS_FPR_DIV2E8_FXP3264
typedef fpr (*fpr_unary_fn)(fpr);

static int
check_fpr_div2e_fixed_raw(const mpz_t raw, const char *kind, unsigned idx,
	unsigned e, fpr_unary_fn fn)
{
	fpr x = fxp_from_mpz(raw);
	fpr got = fn(x);
	fpr exp = fpr_div2e(x, e);

	if (fxp_eq(got, exp)) {
		return 0;
	}
	gmp_printf("  fpr_div2e%u %s[%u] FAILED: raw=%Zd\n",
		e, kind, idx, raw);
	printf("    got=%08lx:%08lx:%08lx exp=%08lx:%08lx:%08lx\n",
		(unsigned long)got.w[2], (unsigned long)got.w[1],
		(unsigned long)got.w[0],
		(unsigned long)exp.w[2], (unsigned long)exp.w[1],
		(unsigned long)exp.w[0]);
	return 1;
}

static void
set_u88(mpz_t d, uint32_t hi24, uint32_t mid, uint32_t lo)
{
	mpz_set_ui(d, hi24);
	mpz_mul_2exp(d, d, 32);
	mpz_add_ui(d, d, mid);
	mpz_mul_2exp(d, d, 32);
	mpz_add_ui(d, d, lo);
}

static int
test_fpr_div2e_fixed_fxp3264(
	const char *name, unsigned e, fpr_unary_fn fn)
{
	static const uint32_t qpos[][3] = {
		{ 0x00000000, 0x00000000, 0x00000000 },
		{ 0x00000000, 0x00000000, 0x00000001 },
		{ 0x00000000, 0x00000000, 0xFFFFFFFF },
		{ 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF },
		{ 0x007FFFFF, 0xFFFFFFFF, 0xFFFFFFFF }
	};
	static const int32_t qneg[] = {
		-1, -2, -3, -255, -256, -257
	};
	static const uint32_t rv[] = {
		0, 1, 63, 64, 65, 127, 128, 129,
		254, 255, 256, 257, 510, 511
	};
	shake_context sc;
	mpz_t raw, q, lim;
	unsigned failures = 0;
	unsigned idx = 0;
	unsigned edge_count;

	mpz_inits(raw, q, lim, (mpz_ptr)0);
	mpz_setbit(lim, FXP_TOTAL_BITS - 1);

	for (size_t i = 0; i < sizeof qpos / sizeof qpos[0]; i++) {
		for (size_t j = 0; j < sizeof rv / sizeof rv[0]; j++) {
			if (rv[j] >= ((uint32_t)1 << e)) {
				continue;
			}
			set_u88(q, qpos[i][0], qpos[i][1], qpos[i][2]);
			mpz_mul_2exp(raw, q, e);
			mpz_add_ui(raw, raw, rv[j]);
			if (mpz_cmpabs(raw, lim) >= 0) {
				continue;
			}
			failures += check_fpr_div2e_fixed_raw(
				raw, "edge", idx++, e, fn);
		}
	}
	for (size_t i = 0; i < sizeof qneg / sizeof qneg[0]; i++) {
		for (size_t j = 0; j < sizeof rv / sizeof rv[0]; j++) {
			if (rv[j] >= ((uint32_t)1 << e)) {
				continue;
			}
			mpz_set_si(q, qneg[i]);
			mpz_mul_2exp(raw, q, e);
			mpz_add_ui(raw, raw, rv[j]);
			if (mpz_cmpabs(raw, lim) >= 0) {
				continue;
			}
			failures += check_fpr_div2e_fixed_raw(
				raw, "edge", idx++, e, fn);
		}
	}
	edge_count = idx;

	shake_init(&sc, 256);
	shake_inject(&sc, g_seed, strlen(g_seed));
	shake_inject(&sc, name, strlen(name));
	shake_flip(&sc);
	for (idx = 0; idx < g_nvec; idx++) {
		uint8_t b[12];

		shake_extract(&sc, b, sizeof b);
		mpz_import(raw, sizeof b, -1, 1, 0, 0, b);
		if (mpz_cmp(raw, lim) >= 0) {
			mpz_submul_ui(raw, lim, 2);
		}
		failures += check_fpr_div2e_fixed_raw(
			raw, "rand", idx, e, fn);
	}

	mpz_clears(raw, q, lim, (mpz_ptr)0);
	if (failures != 0) {
		printf("  %s FAILED (%u mismatches)\n", name, failures);
		return -1;
	}
	printf("  %s PASS (%u edge + %lu random)\n",
		name, edge_count, g_nvec);
	return 0;
}
#endif

#if defined FNDSA_HAS_FPR_RINT_DIV2E8_TO_U16
typedef uint16_t (*fpr_to_u16_fn)(fpr);

static uint16_t
ref_rint_div2e_to_u16(fpr x, unsigned e)
{
	mpz_t raw, q, r, half;
	uint16_t y;

	mpz_inits(raw, q, r, half, (mpz_ptr)0);
	fxp_to_mpz(raw, x);
	mpz_fdiv_q_2exp(q, raw, FXP_FRAC_BITS + e);
	mpz_fdiv_r_2exp(r, raw, FXP_FRAC_BITS + e);
	mpz_setbit(half, FXP_FRAC_BITS + e - 1);
	if (mpz_cmp(r, half) > 0
		|| (mpz_cmp(r, half) == 0 && mpz_odd_p(q)))
	{
		mpz_add_ui(q, q, 1);
	}
	y = (uint16_t)mpz_get_si(q);
	mpz_clears(raw, q, r, half, (mpz_ptr)0);
	return y;
}

static int
check_rint_div2e_raw(const mpz_t raw, const char *kind, unsigned idx,
	unsigned e, fpr_to_u16_fn fn)
{
	fpr x = fxp_from_mpz(raw);
	uint16_t got = fn(x);
	uint16_t exp = ref_rint_div2e_to_u16(x, e);

	if (got == exp) {
		return 0;
	}
	gmp_printf("  rint_div2e%u %s[%u] FAILED:"
		" raw=%Zd got=%u exp=%u\n",
		e, kind, idx, raw, (unsigned)got, (unsigned)exp);
	return 1;
}

static int
test_rint_div2e_to_u16(const char *name, unsigned e, fpr_to_u16_fn fn)
{
	static const int32_t qv[] = {
		-32768, -32767, -257, -256, -255, -3, -2, -1,
		0, 1, 2, 3, 255, 256, 257, 32766, 32767
	};
	static const unsigned rv[] = {
		0, 1, 63, 64, 65, 127, 128, 129,
		254, 255, 256, 257, 510, 511
	};
	shake_context sc;
	mpz_t raw, r, lim;
	unsigned failures = 0;
	unsigned idx = 0;
	unsigned edge_count;

	mpz_inits(raw, r, lim, (mpz_ptr)0);

	for (size_t i = 0; i < sizeof qv / sizeof qv[0]; i++) {
		for (size_t j = 0; j < sizeof rv / sizeof rv[0]; j++) {
			if (rv[j] >= ((unsigned)1 << e)) {
				continue;
			}
			mpz_set_si(raw, qv[i]);
			mpz_mul_2exp(raw, raw, FXP_FRAC_BITS + e);
			mpz_set_ui(r, rv[j]);
			mpz_mul_2exp(r, r, FXP_FRAC_BITS);
			if ((rv[j] & (((unsigned)1 << (e - 1)) - 1)) != 0) {
				mpz_add_ui(r, r, 1);
			}
			mpz_add(raw, raw, r);
			failures += check_rint_div2e_raw(
				raw, "edge", idx++, e, fn);
		}
	}
	edge_count = idx;

	/* Deterministic full-width raw inputs, reduced to the signed range. */
	shake_init(&sc, 256);
	shake_inject(&sc, g_seed, strlen(g_seed));
	shake_inject(&sc, name, strlen(name));
	shake_flip(&sc);
	mpz_setbit(lim, FXP_TOTAL_BITS - 1);
	for (idx = 0; idx < g_nvec; idx++) {
		uint8_t b[12];

		shake_extract(&sc, b, sizeof b);
		mpz_import(raw, sizeof b, -1, 1, 0, 0, b);
		if (mpz_cmp(raw, lim) >= 0) {
			mpz_submul_ui(raw, lim, 2);
		}
		failures += check_rint_div2e_raw(raw, "rand", idx, e, fn);
	}

	mpz_clears(raw, r, lim, (mpz_ptr)0);
	if (failures != 0) {
		printf("  %s FAILED (%u mismatches)\n", name, failures);
		return -1;
	}
	printf("  %s PASS (%u edge + %lu random)\n",
		name, edge_count, g_nvec);
	return 0;
}
#endif

static fpr
href_add(fpr x, fpr y)
{
	fxp_to_mpz(za, x);
	fxp_to_mpz(zb, y);
	mpz_add(zr, za, zb);
	return fxp_from_mpz(zr);
}

static fpr
href_sub(fpr x, fpr y)
{
	fxp_to_mpz(za, x);
	fxp_to_mpz(zb, y);
	mpz_sub(zr, za, zb);
	return fxp_from_mpz(zr);
}

static fpr
href_neg(fpr x)
{
	fxp_to_mpz(za, x);
	mpz_neg(za, za);
	return fxp_from_mpz(za);
}

static fpr
href_mul(fpr x, fpr y)
{
	fxp_to_mpz(za, x);
	fxp_to_mpz(zb, y);
	mpz_mul(zt, za, zb);
	mpz_rshift_rn(zr, zt, FXP_FRAC_BITS);
	return fxp_from_mpz(zr);
}

/* Round-to-nearest-ties-even division of raw quantities:
   result = RN((x << FRAC) / y). */
static fpr
href_div(fpr x, fpr y)
{
	mpz_t n, d, q, r, r2;
	int neg, c;

	mpz_init(n);
	mpz_init(d);
	mpz_init(q);
	mpz_init(r);
	mpz_init(r2);
	fxp_to_mpz(za, x);
	fxp_to_mpz(zb, y);
	neg = (mpz_sgn(za) < 0) != (mpz_sgn(zb) < 0);
	mpz_abs(n, za);
	mpz_mul_2exp(n, n, FXP_FRAC_BITS);
	mpz_abs(d, zb);
	mpz_tdiv_qr(q, r, n, d);
	mpz_mul_2exp(r2, r, 1);
	c = mpz_cmp(r2, d);
	if (c > 0 || (c == 0 && mpz_odd_p(q))) {
		mpz_add_ui(q, q, 1);
	}
	if (neg) {
		mpz_neg(q, q);
	}
	{
		fpr res = fxp_from_mpz(q);
		mpz_clear(n);
		mpz_clear(d);
		mpz_clear(q);
		mpz_clear(r);
		mpz_clear(r2);
		return res;
	}
}

/* RN(sqrt(x)): s = isqrt(x << FRAC), rounded up when the remainder
   exceeds s (exact half-way points cannot occur for integer radicands). */
static fpr
href_sqrt(fpr x)
{
	mpz_t n, s, r;

	mpz_init(n);
	mpz_init(s);
	mpz_init(r);
	fxp_to_mpz(za, x);
	mpz_mul_2exp(n, za, FXP_FRAC_BITS);
	mpz_sqrtrem(s, r, n);
	if (mpz_cmp(r, s) > 0) {
		mpz_add_ui(s, s, 1);
	}
	{
		fpr res = fxp_from_mpz(s);
		mpz_clear(n);
		mpz_clear(s);
		mpz_clear(r);
		return res;
	}
}

/* ==================================================================== */
/* Running one entry.                                                   */

static mpfr_rnd_t
map_rnd(fxp_test_rnd r)
{
	switch (r) {
	case FXT_RND_Z: return MPFR_RNDZ;
	case FXT_RND_D: return MPFR_RNDD;
	case FXT_RND_U: return MPFR_RNDU;
	default:        return MPFR_RNDN;
	}
}

/* Oracle value to grid, in the entry's rounding mode (FXT_RND_NU is
   not an MPFR mode: add half an ulp, then floor). */
static int
grid_round(fpr *out, const mpfr_t v, fxp_test_rnd r)
{
	if (r == FXT_RND_NU) {
		return fxp_from_mpfr_half_up(out, v);
	}
	return fxp_from_mpfr(out, v, map_rnd(r));
}

static void
print_fpr_val(const char *tag, fpr x)
{
	fprintf(stderr, "    %s = 0x", tag);
	for (int i = FXP_WORDS - 1; i >= 0; i--) {
		fprintf(stderr, "%08lx", (unsigned long)x.w[i]);
	}
	fxp_to_mpfr(vu, x);
	mpfr_fprintf(stderr, "  (%.40Rg)\n", vu);
}

static int
cmp_fpr(const fxp_test_entry *ent, fpr got, fpr exp)
{
	if (fxp_eq(got, exp)) {
		return 1;
	}
	/* Passes within either ulp_bound or rel_bits; both 0 means exact. */
	if (ent->ulp_bound > 0) {
		fxp_ulp_dist(zt, got, exp);
		if (mpz_cmp_ui(zt, ent->ulp_bound) <= 0) {
			return 1;
		}
	}
	if (ent->rel_bits > 0) {
		/* |got - exp| * 2^rel_bits <= |exp|. exp == 0 is already handled
		   by the exact compare above, so a nonzero got here fails it. */
		fxp_ulp_dist(zt, got, exp);            /* |raw(got) - raw(exp)| */
		fxp_to_mpz(zr, exp);
		mpz_abs(zr, zr);                       /* |raw(exp)| */
		mpz_mul_2exp(zt, zt, ent->rel_bits);
		if (mpz_cmp(zt, zr) <= 0) {
			return 1;
		}
	}
	return 0;
}

/* Result of one vector: 1 = pass, 0 = discard (out of domain), -1 = fail
   (details already printed). */
static int
run_one(const fxp_test_entry *ent, shake_context *sc,
	int use_edges, size_t ei, size_t ej)
{
	fpr a, b, got, exp;
	int64_t ia = 0, igot, iexp;
	int sc_arg = 0;
	unsigned e_arg = 0;
	fxp_test_rnd rnd = ent->rnd;
	int is_int_result = 0, is_two_results = 0;
	fpr got2, exp2;
	fpr c = FPR_ZERO, d = FPR_ZERO;

	memset(&got, 0, sizeof got);
	memset(&exp, 0, sizeof exp);
	got2 = exp2 = FPR_ZERO;

	/* --- operands ------------------------------------------------- */
	if (use_edges) {
		a = fxp_from_mpz(edge_pool[ei]);
		b = fxp_from_mpz(edge_pool[ej]);
	} else {
		a = gen_fpr(sc, zt, ent->dom_int_bits);
		b = gen_fpr(sc, zt, ent->dom_int_bits);
	}
	switch (ent->kind) {
	case FXT_SQRT:
	case FXT_P63:
		fxp_to_mpz(za, a);
		mpz_abs(za, za);
		a = fxp_from_mpz(za);
		break;
	case FXT_SQRT_INV_SIGMA:
		/* Range of the rsqrt leaf; all even shifts t in [4,30]. */
#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264
		a = use_edges ? rsqrt_alpha110_edge(ei) : gen_rsqrt_alpha110(sc);
#endif
		break;
	case FXT_MUL_RECIP:
	case FXT_MUL_RECIP_P63:
		/* Range of the reciprocal leaf; all shifts k in [6,29]. */
#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264
		b = use_edges ? recip_alpha55_edge(ej) : gen_recip_alpha55(sc);
#endif
		break;
	case FXT_DIV:
	case FXT_INV:
		fxp_to_mpz(zb, ent->kind == FXT_INV ? a : b);
		if (mpz_sgn(zb) == 0) {
			return 0;
		}
		break;
	default:
		break;
	}

	/* Discard reciprocated operands below dom_min_bits. */
	if (ent->dom_min_bits > 0
		&& (ent->kind == FXT_DIV || ent->kind == FXT_INV
			|| ent->kind == FXT_SQRT || ent->kind == FXT_MUL_RECIP
			|| ent->kind == FXT_MUL_RECIP_P63
			|| ent->kind == FXT_SQRT_INV_SIGMA))
	{
		fxp_to_mpz(zt, (ent->kind == FXT_DIV
			|| ent->kind == FXT_MUL_RECIP
			|| ent->kind == FXT_MUL_RECIP_P63) ? b : a);
		mpz_abs(zt, zt);
		if (mpz_sgn(zt) == 0
			|| mpz_sizeinbase(zt, 2) < (size_t)ent->dom_min_bits)
		{
			return 0;
		}
	}

	/* A Q0.63 result cannot hold 1, so keep only |a/b| < 1 - 2^-64,
	   i.e. |b| < 2^64 * (|b| - |a|). */
	if (ent->kind == FXT_MUL_RECIP_P63) {
		fxp_to_mpz(zt, a);
		mpz_abs(zt, zt);
		fxp_to_mpz(zb, b);
		mpz_abs(zb, zb);
		mpz_sub(zt, zb, zt);
		if (mpz_sgn(zt) <= 0) {
			return 0;
		}
		mpz_mul_2exp(zt, zt, 64);
		if (mpz_cmp(zb, zt) >= 0) {
			return 0;
		}
	}

	switch (ent->kind) {
	case FXT_OF:
		ia = use_edges ? (ei == 0 ? 0 : (int64_t)1 << (ent->dom_int_bits - 1))
			: gen_i64(sc, ent->dom_int_bits > 62 ? 62 : ent->dom_int_bits);
		break;
	case FXT_OF32:
		ia = use_edges ? (ei == 0 ? 0 : INT32_MAX)
			: (int32_t)gen_i64(sc,
				ent->dom_int_bits > 31 ? 31 : ent->dom_int_bits);
		break;
	case FXT_SCALED:
		ia = gen_i64(sc, ent->dom_int_bits > 62 ? 62 : ent->dom_int_bits);
		sc_arg = 8 - (int)(rnd_u32(sc) % (FXP_FRAC_BITS + 17));
		break;
	case FXT_MUL2E:
	case FXT_DIV2E:
		e_arg = rnd_u32(sc) % (unsigned)(ent->e_max + 1);
		break;
	case FXT_SQRT_INV_SIGMA:
		e_arg = 1 + (rnd_u32(sc) % 10);
		break;
	case FXT_INV:
		/* Positive operands only, as in signing. */
		fxp_to_mpz(zt, a);
		if (mpz_sgn(zt) < 0) {
			a = fpr_neg(a);
		}
		break;
	case FXT_DIV:
	case FXT_MUL_RECIP:
	case FXT_MUL_RECIP_P63:
		/* Same domain for the divisor (div = mul by reciprocal). */
		fxp_to_mpz(zt, b);
		if (mpz_sgn(zt) < 0) {
			b = fpr_neg(b);
		}
		break;
	case FXT_FPC_MUL:
	case FXT_FPC_MUL_FUSED:
		c = gen_fpr(sc, zt, ent->dom_int_bits);
		d = gen_fpr(sc, zt, ent->dom_int_bits);
		break;
	default:
		break;
	}

	/* --- oracle + operation --------------------------------------- */
	fxp_to_mpfr(va, a);
	fxp_to_mpfr(vb, b);
	switch (ent->kind) {
	case FXT_ADD:
		mpfr_add(vr, va, vb, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		got = ent->bin ? ent->bin(a, b) : fpr_add(a, b);
		break;
	case FXT_SUB:
		mpfr_sub(vr, va, vb, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		got = ent->bin ? ent->bin(a, b) : fpr_sub(a, b);
		break;
	case FXT_MUL:
		mpfr_mul(vr, va, vb, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		got = ent->bin ? ent->bin(a, b) : fpr_mul(a, b);
		break;
	case FXT_SQR:
		mpfr_sqr(vr, va, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		got = ent->un ? ent->un(a) : fpr_sqr(a);
		break;
	case FXT_DIV:
		mpfr_div(vr, va, vb, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		got = ent->bin ? ent->bin(a, b) : fpr_div(a, b);
		break;
	case FXT_MUL_RECIP:
	case FXT_MUL_RECIP_P63:
		/* Exact quotient, rounded once. */
		mpfr_div(vr, va, vb, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		got = ent->bin(a, b);
		break;
	case FXT_SQRT:
		mpfr_sqrt(vr, va, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		got = ent->un ? ent->un(a) : fpr_sqrt(a);
		break;
	case FXT_SQRT_INV_SIGMA:
		mpfr_sqrt(vr, va, MPFR_RNDN);
#if FNDSA_SIGN_SQRT_FUSED && defined(INV_SIGMA_Q71)
		mpz_set_ui(zt,
			(unsigned long)(uint32_t)(INV_SIGMA_Q71[e_arg] >> 32));
		mpz_mul_2exp(zt, zt, 32);
		mpz_add_ui(zt, zt,
			(unsigned long)(uint32_t)INV_SIGMA_Q71[e_arg]);
		mpfr_set_z(vt, zt, MPFR_RNDN);
		mpfr_div_2ui(vt, vt, 71, MPFR_RNDN);
#else
		fxp_to_mpfr(vt, INV_SIGMA[e_arg].f);
#endif
		mpfr_mul(vr, vr, vt, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		FPR_SQRT_MUL_INV_SIGMA(got, a, e_arg);
		break;
	case FXT_MUL_INV_Q:
#ifdef INV_Q_Q77
		/* The word must be round(2^77/q), i.e. |c*q - 2^77| <= q/2;
		   the oracle is then the exact product with it, rounded once. */
		mpz_set_ui(zt, (unsigned long)INV_Q_Q77_HI);
		mpz_mul_2exp(zt, zt, 32);
		mpz_add_ui(zt, zt, (unsigned long)INV_Q_Q77_LO);
		mpz_mul_ui(zb, zt, 12289);
		mpz_set_ui(zr, 0);
		mpz_setbit(zr, 77);
		mpz_sub(zb, zb, zr);
		mpz_abs(zb, zb);
		if (mpz_cmp_ui(zb, 12289 / 2) > 0) {
			fprintf(stderr, "FAIL %s: INV_Q_Q77 is not round(2^77/q)\n",
				ent->name);
			return -1;
		}
		mpfr_set_z(vt, zt, MPFR_RNDN);
		mpfr_div_2ui(vt, vt, 77, MPFR_RNDN);
		mpfr_mul(vr, va, vt, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		got = FPR_MUL_INV_Q(a);
		break;
#else
		fprintf(stderr, "%s: no wide 1/q in this representation\n",
			ent->name);
		return -1;
#endif
	case FXT_INV:
		mpfr_ui_div(vr, 1, va, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		got = ent->un ? ent->un(a) : fpr_inv(a);
		break;
	case FXT_NEG:
		mpfr_neg(vr, va, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		got = ent->un ? ent->un(a) : fpr_neg(a);
		break;
	case FXT_HALF:
		mpfr_div_2si(vr, va, 1, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		got = ent->un ? ent->un(a) : fpr_half(a);
		break;
	case FXT_DOUBLE:
		mpfr_mul_2si(vr, va, 1, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		got = ent->un ? ent->un(a) : fpr_double(a);
		break;
	case FXT_MUL2E:
		mpfr_mul_2si(vr, va, (long)e_arg, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		got = fpr_mul2e(a, e_arg);
		break;
	case FXT_DIV2E:
		mpfr_div_2si(vr, va, (long)e_arg, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		got = fpr_div2e(a, e_arg);
		break;
	case FXT_OF:
		mpfr_set_sj(vr, (intmax_t)ia, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		got = fpr_of(ia);
		break;
	case FXT_OF32:
		mpfr_set_sj(vr, (intmax_t)ia, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		got = fpr_of32((int32_t)ia);
		break;
	case FXT_SCALED:
		mpfr_set_sj(vr, (intmax_t)ia, MPFR_RNDN);
		mpfr_mul_2si(vr, vr, sc_arg, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		got = fpr_scaled(ia, sc_arg);
		break;
	case FXT_RINT:
	case FXT_FLOOR:
	case FXT_TRUNC: {
		mpfr_rnd_t irnd = ent->kind == FXT_RINT ? MPFR_RNDN
			: ent->kind == FXT_FLOOR ? MPFR_RNDD : MPFR_RNDZ;
		is_int_result = 1;
		if (!mpfr_fits_intmax_p(va, irnd)) {
			return 0;
		}
		iexp = (int64_t)mpfr_get_sj(va, irnd);
		igot = ent->kind == FXT_RINT ? fpr_rint(a)
			: ent->kind == FXT_FLOOR ? fpr_floor(a) : fpr_trunc(a);
		break;
	}
	case FXT_P63:
		/* Exact integer oracle: with x = raw * 2^-FRAC in [0, 1)
		   (dom_int_bits 0 plus the abs above), floor(x * 2^63) is
		   raw >> (FRAC - 63). */
		is_int_result = 1;
		fxp_to_mpz(za, a);
		mpz_fdiv_q_2exp(za, za, FXP_FRAC_BITS - 63);
		iexp = (int64_t)mpz_get_ui(za);
		igot = (int64_t)fpr_p63(a);
		break;
	case FXT_ADD_SUB:
		is_two_results = 1;
		mpfr_add(vr, va, vb, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		mpfr_sub(vr, va, vb, MPFR_RNDN);
		if (!grid_round(&exp2, vr, rnd)) return 0;
		FPR_ADD_SUB(got, got2, a, b);
		break;
	case FXT_FPC_MUL: {
		/* Each product rounded, then the sum rounded. */
		fpr m1, m2, m3, m4;
		is_two_results = 1;
		fxp_to_mpfr(vt, c);
		fxp_to_mpfr(vu, d);
		/* a=a_re, b=a_im, c=b_re, d=b_im */
		mpfr_mul(vr, va, vt, MPFR_RNDN);              /* a_re*b_re */
		if (!grid_round(&m1, vr, rnd)) return 0;
		mpfr_mul(vr, vb, vu, MPFR_RNDN);              /* a_im*b_im */
		if (!grid_round(&m2, vr, rnd)) return 0;
		mpfr_mul(vr, va, vu, MPFR_RNDN);              /* a_re*b_im */
		if (!grid_round(&m3, vr, rnd)) return 0;
		mpfr_mul(vr, vb, vt, MPFR_RNDN);              /* a_im*b_re */
		if (!grid_round(&m4, vr, rnd)) return 0;
		fxp_to_mpfr(vr, m1);
		fxp_to_mpfr(vu, m2);
		mpfr_sub(vr, vr, vu, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		fxp_to_mpfr(vr, m3);
		fxp_to_mpfr(vu, m4);
		mpfr_add(vr, vr, vu, MPFR_RNDN);
		if (!grid_round(&exp2, vr, rnd)) return 0;
		FPC_MUL(got, got2, a, b, c, d);
		break;
	}
	case FXT_FPC_MUL_FUSED:
		/* Fused: products and sums are exact at 256 bits, so each
		   component is rounded once. */
		is_two_results = 1;
		fxp_to_mpfr(vt, c);
		fxp_to_mpfr(vu, d);
		/* a=a_re, b=a_im, c=b_re, d=b_im */
		mpfr_mul(vr, va, vt, MPFR_RNDN);              /* a_re*b_re */
		mpfr_mul(vw, vb, vu, MPFR_RNDN);              /* a_im*b_im */
		mpfr_sub(vr, vr, vw, MPFR_RNDN);
		if (!grid_round(&exp, vr, rnd)) return 0;
		mpfr_mul(vr, va, vu, MPFR_RNDN);              /* a_re*b_im */
		mpfr_mul(vw, vb, vt, MPFR_RNDN);              /* a_im*b_re */
		mpfr_add(vr, vr, vw, MPFR_RNDN);
		if (!grid_round(&exp2, vr, rnd)) return 0;
		FPC_MUL(got, got2, a, b, c, d);
		break;
	default:
		fprintf(stderr, "unknown kind %d\n", (int)ent->kind);
		return -1;
	}

	/* --- compare --------------------------------------------------- */
	if (is_int_result) {
		if (igot == iexp) {
			return 1;
		}
		fprintf(stderr, "FAIL %s:\n", ent->name);
		print_fpr_val("x", a);
		fprintf(stderr, "    got %lld, expected %lld\n",
			(long long)igot, (long long)iexp);
		return -1;
	}
	if (!cmp_fpr(ent, got, exp)
		|| (is_two_results && !cmp_fpr(ent, got2, exp2))) {
		fprintf(stderr, "FAIL %s:\n", ent->name);
		print_fpr_val("x", a);
		if (ent->kind != FXT_NEG && ent->kind != FXT_HALF
			&& ent->kind != FXT_SQRT
			&& ent->kind != FXT_SQRT_INV_SIGMA
			&& ent->kind != FXT_MUL_INV_Q
			&& ent->kind != FXT_P63)
		{
			print_fpr_val("y", b);
		}
		if (ent->kind == FXT_FPC_MUL
			|| ent->kind == FXT_FPC_MUL_FUSED)
		{
			print_fpr_val("z", c);
			print_fpr_val("w", d);
		}
		if (ent->kind == FXT_MUL2E || ent->kind == FXT_DIV2E) {
			fprintf(stderr, "    e = %u\n", e_arg);
		}
		if (ent->kind == FXT_SQRT_INV_SIGMA) {
			fprintf(stderr, "    logn = %u\n", e_arg);
		}
		if (ent->kind == FXT_OF || ent->kind == FXT_OF32
			|| ent->kind == FXT_SCALED) {
			fprintf(stderr, "    i = %lld, sc = %d\n",
				(long long)ia, sc_arg);
		}
		print_fpr_val("got     ", got);
		print_fpr_val("expected", exp);
		if (is_two_results) {
			print_fpr_val("got2    ", got2);
			print_fpr_val("expected2", exp2);
		}
		return -1;
	}
	return 1;
}

/* Run all vectors for one entry; returns 0 on pass, -1 on failure. */
static int
run_entry(const fxp_test_entry *ent)
{
	shake_context sc;
	unsigned long checked = 0, discarded = 0, attempts = 0;
	size_t nedges, i, j;
	uint8_t kb;
	int binary_edges;

	if (!ent->enabled) {
		printf("  %-22s SKIP (not implemented)\n", ent->name);
		return 0;
	}

	shake_init(&sc, 256);
	shake_inject(&sc, g_seed, strlen(g_seed));
	shake_inject(&sc, ent->name, strlen(ent->name));
	kb = (uint8_t)ent->kind;
	shake_inject(&sc, &kb, 1);
	shake_flip(&sc);
#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264
	g_leaf_bucket_index = 0;
#endif

	/* Edge vectors (only meaningful for fpr-operand kinds). */
	switch (ent->kind) {
	case FXT_OF: case FXT_OF32: case FXT_SCALED:
		nedges = 0;
		break;
	default:
		nedges = build_edges(ent->dom_int_bits);
		break;
	}
	binary_edges = ent->kind == FXT_ADD || ent->kind == FXT_SUB
		|| ent->kind == FXT_MUL || ent->kind == FXT_DIV
		|| ent->kind == FXT_ADD_SUB || ent->kind == FXT_FPC_MUL
		|| ent->kind == FXT_FPC_MUL_FUSED
		|| ent->kind == FXT_MUL_RECIP
		|| ent->kind == FXT_MUL_RECIP_P63;
	for (i = 0; i < nedges; i++) {
		for (j = 0; j < (binary_edges ? nedges : 1); j++) {
			switch (run_one(ent, &sc, 1, i, j)) {
			case 1: checked++; break;
			case 0: discarded++; break;
			default: return -1;
			}
		}
	}

	/* Random vectors. */
	while (checked < g_nvec && attempts < 4 * g_nvec) {
		attempts++;
		switch (run_one(ent, &sc, 0, 0, 0)) {
		case 1: checked++; break;
		case 0: discarded++; break;
		default: return -1;
		}
	}

	printf("  %-22s PASS (%lu vectors, %lu discarded)\n",
		ent->name, checked, discarded);
	return 0;
}

/* ==================================================================== */

/* fndsa_fpr_mul_rnj (the C version of the asm RNJ leaf) must match
   fpr_mul for every table coefficient. */
#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 && FNDSA_TWIDDLE_U64 \
	&& FNDSA_FFT_NEW
static int
test_rnj_decomposition(void)
{
	shake_context sc;
	unsigned long checked = 0;

	shake_init(&sc, 256);
	shake_inject(&sc, "rnj-decomposition", 17);
	shake_flip(&sc);
	for (size_t j = 1; j < 512; j ++) {
		fpr r = RNJ[j];
		if (r.w[2] >> 31) {
			r = fpr_neg(r);
		}
		uint64_t frac = (uint64_t)r.w[0]
			| ((uint64_t)r.w[1] << 32);
		uint32_t integer = r.w[2];
		for (unsigned k = 0; k < 260; k ++) {
			fpr x;
			if (k == 0) {
				x = (fpr){{ 1, 0, 0 }};
			} else if (k == 1) {
				x = (fpr){{ UINT32_MAX, UINT32_MAX, UINT32_MAX }};
			} else if (k == 2) {
				x = fpr_of32(1);
			} else if (k == 3) {
				x = fpr_of32(-1);
			} else {
				x.w[0] = rnd_u32(&sc);
				x.w[1] = rnd_u32(&sc);
				x.w[2] = rnd_u32(&sc) & 0x000FFFFF;
				if (rnd_u32(&sc) & 1) {
					x = fpr_neg(x);
				}
			}
			fpr got = fndsa_fpr_mul_rnj(x, frac, integer);
			fpr exp = fpr_mul(x, r);
			if (memcmp(got.w, exp.w, sizeof got.w) != 0) {
				printf("  RNJ decomposition FAIL (table %zu, vector %u)\n",
					j, k);
				return -1;
			}
			checked ++;
		}
	}

	/* The compact transform's initial merge uses sqrt(2) = RNJ[0]. */
	{
		fpr r = RNJ[0];
		uint64_t frac = (uint64_t)r.w[0]
			| ((uint64_t)r.w[1] << 32);
		for (unsigned k = 0; k < 260; k ++) {
			fpr x = {{ rnd_u32(&sc), rnd_u32(&sc),
				rnd_u32(&sc) & 0x000FFFFF }};
			if (rnd_u32(&sc) & 1) {
				x = fpr_neg(x);
			}
			fpr got = fndsa_fpr_mul_rnj(x, frac, r.w[2]);
			fpr exp = fpr_mul(x, r);
			if (memcmp(got.w, exp.w, sizeof got.w) != 0) {
				printf("  RNJ decomposition FAIL (sqrt2, vector %u)\n", k);
				return -1;
			}
			checked ++;
		}
	}
	printf("  RNJ decomposition      PASS (%lu C-reference vectors)\n",
		checked);
	return 0;
}
#endif

/* ==================================================================== */

int
main(int argc, char *argv[])
{
	static const fxp_test_entry harness[] = {
		{ .name = "harness:add", .kind = FXT_ADD, .enabled = 1,
			.dom_int_bits = FXP_INT_BITS - 2, .bin = href_add },
		{ .name = "harness:sub", .kind = FXT_SUB, .enabled = 1,
			.dom_int_bits = FXP_INT_BITS - 2, .bin = href_sub },
		{ .name = "harness:neg", .kind = FXT_NEG, .enabled = 1,
			.dom_int_bits = FXP_INT_BITS - 1, .un = href_neg },
		{ .name = "harness:mul", .kind = FXT_MUL, .enabled = 1,
			.dom_int_bits = FXP_INT_BITS / 2 - 1, .bin = href_mul },
		{ .name = "harness:div", .kind = FXT_DIV, .enabled = 1,
			.dom_int_bits = FXP_INT_BITS / 2 - 1, .bin = href_div },
		{ .name = "harness:sqrt", .kind = FXT_SQRT, .enabled = 1,
			.dom_int_bits = FXP_INT_BITS - 1, .un = href_sqrt },
	};
	size_t i;
	int failures = 0;

	if (argc >= 2) {
		g_nvec = strtoul(argv[1], 0, 0);
		if (g_nvec == 0) {
			g_nvec = 10000;
		}
	}
	if (argc >= 3) {
		g_seed = argv[2];
	}

	fxp_oracle_init();
	mpfr_inits2(FXP_ORACLE_PREC, va, vb, vr, vt, vu, vw, (mpfr_ptr)0);
	mpz_inits(za, zb, zr, zt, (mpz_ptr)0);
	for (i = 0; i < MAX_EDGES; i++) {
		mpz_init(edge_pool[i]);
	}

	printf("fpr per-op tests, backend %s (%d.%d fixed point), "
		"%lu vectors/op, seed \"%s\"\n",
		fxp_test_backend, FXP_INT_BITS, FXP_FRAC_BITS, g_nvec, g_seed);

	printf("oracle self-test: ");
	if (fxp_oracle_selftest() != 0) {
		printf("FAILED\n");
		return 2;
	}
	printf("ok\n");

#if defined FNDSA_HAS_FPR_DIV2E8_FXP3264
	printf("iFFT fixed div2e helper self-check:\n");
	if (test_fpr_div2e_fixed_fxp3264(
		"fpr_div2e8_fxp3264", 8, fpr_div2e8_fxp3264) != 0)
	{
		return 2;
	}
#if defined FNDSA_HAS_FPR_DIV2E9_FXP3264
	if (test_fpr_div2e_fixed_fxp3264(
		"fpr_div2e9_fxp3264", 9, fpr_div2e9_fxp3264) != 0)
	{
		return 2;
	}
#endif
#endif

#if defined FNDSA_HAS_FPR_RINT_DIV2E8_TO_U16
	printf("split-tint helper self-check:\n");
	if (test_rint_div2e_to_u16(
		"rint_div2e8_to_u16", 8, fpr_rint_div2e8_to_u16) != 0)
	{
		return 2;
	}
#if defined FNDSA_HAS_FPR_RINT_DIV2E9_TO_U16
	if (test_rint_div2e_to_u16(
		"rint_div2e9_to_u16", 9, fpr_rint_div2e9_to_u16) != 0)
	{
		return 2;
	}
#endif
#endif

	printf("harness self-check (exact integer references):\n");
	for (i = 0; i < sizeof harness / sizeof harness[0]; i++) {
		if (run_entry(&harness[i]) != 0) {
			printf("harness self-check FAILED -- "
				"fix the harness before trusting any result\n");
			return 2;
		}
	}

#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 && FNDSA_TWIDDLE_U64 \
	&& FNDSA_FFT_NEW
	printf("RNJ delayed-rounding C-reference check:\n");
	if (test_rnj_decomposition() != 0) {
		return 2;
	}
#endif

	printf("backend operations:\n");
	for (i = 0; i < fxp_test_table_len; i++) {
		if (run_entry(&fxp_test_table[i]) != 0) {
			failures++;
		}
	}

	{
		unsigned nskip = 0, npass = 0;
		for (i = 0; i < fxp_test_table_len; i++) {
			if (!fxp_test_table[i].enabled) {
				nskip++;
			} else {
				npass++;
			}
		}
		if (failures) {
			printf("RESULT: %d operation(s) FAILED\n", failures);
			return 1;
		}
		printf("RESULT: all enabled operations pass "
			"(%u enabled, %u still skipped)\n",
			npass - (unsigned)failures, nskip);
	}
	return 0;
}
