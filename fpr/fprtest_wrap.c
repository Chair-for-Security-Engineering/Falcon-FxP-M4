/*
 * Driver for the on-board asm-vs-C fpr test (see fprtest.h). Results are
 * compared byte for byte. Vectors come from SHAKE256 with a fixed seed
 * and are built through the contract ops, with magnitudes kept inside
 * every representation's range.
 */

#include <stdio.h>
#include <string.h>

#include "fpr/fprtest.h"

/* tw has no fpr_scaled on the board. */
#if FNDSA_FPR_REP == FNDSA_FPR_REP_TW
#define FPRTEST_HAVE_SCALED   0
#else
#define FPRTEST_HAVE_SCALED   1
#endif

/* --- device under test: the contract ops of this build ---------------- */

static fpr d_of32(int32_t i)   { return fpr_of32(i); }
#if FPRTEST_HAVE_SCALED
static fpr d_scaled(int64_t i, int sc) { return fpr_scaled(i, sc); }
#endif
static fpr d_add(fpr x, fpr y) { return fpr_add(x, y); }
static fpr d_sub(fpr x, fpr y) { return fpr_sub(x, y); }
static void d_add_sub(fpr *a, fpr *b, fpr x, fpr y) { FPR_ADD_SUB(*a, *b, x, y); }
static fpr d_mul(fpr x, fpr y) { return fpr_mul(x, y); }
static void d_fpc_mul(fpr *d_re, fpr *d_im, fpr a_re, fpr a_im,
	fpr b_re, fpr b_im)
{
#ifdef FPC_MUL_NOALIAS
	FPC_MUL_NOALIAS(*d_re, *d_im, a_re, a_im, b_re, b_im);
#else
	FPC_MUL(*d_re, *d_im, a_re, a_im, b_re, b_im);
#endif
}
#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 && FNDSA_TWIDDLE_U64
static void
d_fxp96_cmul_gm(uint32_t *out, const uint32_t *a, const uint64_t *gm)
{
	FNDSA_FXP96_CMUL_GM(out, a, gm);
}

static void
d_fxp96_cmul_gmj(uint32_t *out, const uint32_t *a, const uint64_t *gm)
{
	FNDSA_FXP96_CMUL_GMJ(out, a, gm);
}

static fpr
d_fpr_mul_gm(fpr x, uint64_t c)
{
	return FNDSA_FPR_MUL_GM(x, c);
}
#endif
static fpr d_sqr(fpr x)        { return fpr_sqr(x); }
static fpr d_div(fpr x, fpr y) { return fpr_div(x, y); }
static fpr d_inv(fpr x)        { return fpr_inv(x); }
static fpr d_sqrt(fpr x)       { return fpr_sqrt(x); }
#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 \
	&& FNDSA_SIGN_SQRT_FUSED && FNDSA_KGEN_GAMMA
static void
d_sqrt_inv_sigma(fpr *d, fpr x, unsigned logn)
{
	FPR_SQRT_MUL_INV_SIGMA(*d, x, logn);
}
#endif
static fpr d_mul_inv_q(fpr x)  { return FPR_MUL_INV_Q(x); }
static int64_t d_rint(fpr x)   { return fpr_rint(x); }
static int64_t d_floor(fpr x)  { return fpr_floor(x); }
static int64_t d_trunc(fpr x)  { return fpr_trunc(x); }

static fpr
d_mul_recip(fpr a, fpr x)
{
	fpr_rcp rcp;
	fpr d;

	FPR_RECIP_X(rcp, x);
	FPR_MUL_RECIP(d, a, rcp);
	return d;
}

#ifdef FPR_MUL_RECIP_P63
#define FPRTEST_HAVE_MUL_RECIP_P63   1
static fpr
p63_to_fpr(uint64_t q)
{
	fpr d;

	d.w[0] = (uint32_t)(q << 1);
	d.w[1] = (uint32_t)(q >> 31);
	d.w[2] = (uint32_t)((int64_t)q >> 63);
	return d;
}

static fpr
d_mul_recip_p63(fpr a, fpr x)
{
	fpr_rcp rcp;
	uint64_t q;

	FPR_RECIP_X(rcp, x);
	FPR_MUL_RECIP_P63(q, a, rcp);
	return p63_to_fpr(q);
}
#else
#define FPRTEST_HAVE_MUL_RECIP_P63   0
#endif

/* --- deterministic vector generation ----------------------------------- */

static shake_context g_sc;
static fpr g_zero;
static void (*g_emit)(const char *);
static unsigned long g_nvec;
static int g_fail_ops;

static uint32_t
rnd32(void)
{
	uint8_t b[4];

	shake_extract(&g_sc, b, sizeof b);
	return (uint32_t)b[0] | ((uint32_t)b[1] << 8)
		| ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 && FNDSA_KGEN_GAMMA
/* Input ranges of the fused leaf kernels: reciprocal d in
   [q/55^2, q*55^2], rsqrt d in [q/110^2, q*110^2]. Random vectors use
   whole-number bounds; the exact endpoints are in the edge tables. */
#define RECIP_ALPHA55_LO       5u
#define RECIP_ALPHA55_HI       37174225u
#define RSQRT_ALPHA110_LO      2u
#define RSQRT_ALPHA110_HI      148696900u

static fpr
fpr_from_q32_64_words(uint32_t w2, uint32_t w1, uint32_t w0)
{
	fpr x;

	x.w[0] = w0;
	x.w[1] = w1;
	x.w[2] = w2;
	return x;
}

static fpr
gen_recip_alpha55(unsigned long i)
{
	unsigned k = 6 + (unsigned)(i % 24);
	uint32_t lo = (uint32_t)1 << (31 - k);
	uint32_t hi = (uint32_t)1 << (32 - k);
	fpr x;

	if (lo < RECIP_ALPHA55_LO) {
		lo = RECIP_ALPHA55_LO;
	}
	if (hi > RECIP_ALPHA55_HI) {
		hi = RECIP_ALPHA55_HI;
	}
	x.w[0] = rnd32();
	x.w[1] = rnd32();
	x.w[2] = lo + rnd32() % (hi - lo);
	return x;
}

static fpr
gen_rsqrt_alpha110(unsigned long i)
{
	unsigned t = 4 + 2 * (unsigned)(i % 14);
	uint32_t lo = (uint32_t)1 << (30 - t);
	uint32_t hi = (uint32_t)1 << (32 - t);
	fpr x;

	if (lo < RSQRT_ALPHA110_LO) {
		lo = RSQRT_ALPHA110_LO;
	}
	if (hi > RSQRT_ALPHA110_HI) {
		hi = RSQRT_ALPHA110_HI;
	}
	x.w[0] = rnd32();
	x.w[1] = rnd32();
	x.w[2] = lo + rnd32() % (hi - lo);
	return x;
}

static fpr
recip_alpha55_edge(unsigned long i)
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
rsqrt_alpha110_edge(unsigned long i)
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

/* Magnitude in about [2^-emax, 2^31]; half the time with more than 32
   mantissa bits. */
static fpr
gen_val(unsigned emax)
{
	fpr x = fpr_div2e(fpr_of32((int32_t)rnd32()), rnd32() % (emax + 1));
	if (rnd32() & 1) {
		fpr y = fpr_div2e(fpr_of32((int32_t)rnd32()),
			rnd32() % (emax + 1));
		x = fpr_add(x, y);
	}
	return x;
}

static fpr
gen_nonneg(unsigned emax)
{
	fpr x = gen_val(emax);
	if (fpr_floor(x) < 0) {
		x = fpr_neg(x);
	}
	return x;
}

static fpr
gen_nonzero(unsigned emax)
{
	fpr x;

	do {
		x = gen_val(emax);
	} while (memcmp(&x, &g_zero, sizeof x) == 0);
	return x;
}

/* --- reporting ---------------------------------------------------------- */

#define MAX_DETAIL   3   /* mismatch details printed per operation */

static void
hex_fpr(char *dst, fpr x)
{
	const unsigned char *p = (const unsigned char *)&x;

	for (size_t i = 0; i < sizeof x; i++) {
		sprintf(dst + 2 * i, "%02x", p[i]);
	}
}

static void
hex_i64(char *dst, int64_t v)
{
	sprintf(dst, "%08lx%08lx",
		(unsigned long)(uint32_t)((uint64_t)v >> 32),
		(unsigned long)(uint32_t)(uint64_t)v);
}

static void
emit_fpr(const char *tag, fpr v)
{
	char line[96];
	size_t k;

	sprintf(line, "    %s = 0x", tag);
	k = strlen(line);
	hex_fpr(line + k, v);
	g_emit(line);
}

static void
op_done(const char *name, unsigned long bad)
{
	char line[96];

	if (bad == 0) {
		sprintf(line, "  %-12s PASS (%lu vectors)", name, g_nvec);
	} else {
		sprintf(line, "  %-12s FAIL (%lu of %lu vectors differ)",
			name, bad, g_nvec);
		g_fail_ops++;
	}
	g_emit(line);
}

static void
op_skip(const char *name, const char *why)
{
	char line[96];

	sprintf(line, "  %-12s SKIP (%s)", name, why);
	g_emit(line);
}

/* --- per-shape check loops ---------------------------------------------- */

static void
check_un(const char *name, fpr (*dut)(fpr), fpr (*cref)(fpr),
	unsigned emax, int nonneg, int nonzero)
{
	unsigned long bad = 0;

	if (cref == 0) {
		op_skip(name, "no separate C reference for this backend");
		return;
	}
	for (unsigned long i = 0; i < g_nvec; i++) {
		fpr x;

		do {
			x = nonneg ? gen_nonneg(emax)
				: nonzero ? gen_nonzero(emax) : gen_val(emax);
		} while (nonzero && memcmp(&x, &g_zero, sizeof x) == 0);
		fpr g = dut(x), e = cref(x);
		if (memcmp(&g, &e, sizeof g) != 0) {
			if (bad < MAX_DETAIL) {
				char line[64];
				sprintf(line, "  %s mismatch, vector %lu:",
					name, i);
				g_emit(line);
				emit_fpr("x  ", x);
				emit_fpr("asm", g);
				emit_fpr("C  ", e);
			}
			bad++;
		}
	}
	op_done(name, bad);
}

static void
check_bin(const char *name, fpr (*dut)(fpr, fpr), fpr (*cref)(fpr, fpr),
	unsigned emax, int pos_y)
{
	unsigned long bad = 0;

	if (cref == 0) {
		op_skip(name, "no separate C reference for this backend");
		return;
	}
	for (unsigned long i = 0; i < g_nvec; i++) {
		fpr x = gen_val(emax);
		fpr y;

		/* pos_y: positive nonzero (the divisor domain of fpr_div). */
		do {
			y = pos_y ? gen_nonneg(emax) : gen_val(emax);
		} while (pos_y && memcmp(&y, &g_zero, sizeof y) == 0);
		fpr g = dut(x, y), e = cref(x, y);
		if (memcmp(&g, &e, sizeof g) != 0) {
			if (bad < MAX_DETAIL) {
				char line[64];
				sprintf(line, "  %s mismatch, vector %lu:",
					name, i);
				g_emit(line);
				emit_fpr("x  ", x);
				emit_fpr("y  ", y);
				emit_fpr("asm", g);
				emit_fpr("C  ", e);
			}
			bad++;
		}
	}
	op_done(name, bad);
}

static void
check_toint(const char *name, int64_t (*dut)(fpr), int64_t (*cref)(fpr),
	unsigned emax)
{
	unsigned long bad = 0;

	if (cref == 0) {
		op_skip(name, "no separate C reference for this backend");
		return;
	}
	for (unsigned long i = 0; i < g_nvec; i++) {
		fpr x = gen_val(emax);
		int64_t g = dut(x), e = cref(x);
		if (g != e) {
			if (bad < MAX_DETAIL) {
				char line[96];
				size_t k;
				sprintf(line, "  %s mismatch, vector %lu:",
					name, i);
				g_emit(line);
				emit_fpr("x  ", x);
				sprintf(line, "    asm = 0x");
				k = strlen(line);
				hex_i64(line + k, g);
				g_emit(line);
				sprintf(line, "    C   = 0x");
				k = strlen(line);
				hex_i64(line + k, e);
				g_emit(line);
			}
			bad++;
		}
	}
	op_done(name, bad);
}

static void
check_add_sub(void)
{
	unsigned long bad = 0;

	if (fndsa_fprtest_cref.op_add_sub == 0) {
		op_skip("FPR_ADD_SUB",
			"no separate C reference for this backend");
		return;
	}
	for (unsigned long i = 0; i < g_nvec; i++) {
		fpr x = gen_val(34), y = gen_val(34);
		fpr ga, gb, ea, eb;
		d_add_sub(&ga, &gb, x, y);
		fndsa_fprtest_cref.op_add_sub(&ea, &eb, x, y);
		if (memcmp(&ga, &ea, sizeof ga) != 0
			|| memcmp(&gb, &eb, sizeof gb) != 0) {
			if (bad < MAX_DETAIL) {
				char line[64];
				sprintf(line,
					"  FPR_ADD_SUB mismatch, vector %lu:",
					i);
				g_emit(line);
				emit_fpr("x    ", x);
				emit_fpr("y    ", y);
				emit_fpr("asm a", ga);
				emit_fpr("C   a", ea);
				emit_fpr("asm b", gb);
				emit_fpr("C   b", eb);
			}
			bad++;
		}
	}
	op_done("FPR_ADD_SUB", bad);
}

/* FPC_MUL operands up to about 2^15, so products and sums stay in
   range. */
static fpr
gen_fpc(void)
{
	return fpr_div2e(gen_val(34), 16);
}

static void
fpc_vec(const fprtest_ops *cr, const char *tag, unsigned long i,
	fpr a_re, fpr a_im, fpr b_re, fpr b_im, unsigned long *bad)
{
	fpr g_re, g_im, e_re, e_im;

	d_fpc_mul(&g_re, &g_im, a_re, a_im, b_re, b_im);
	if (cr->op_fpc_mul != 0) {
		cr->op_fpc_mul(&e_re, &e_im, a_re, a_im, b_re, b_im);
	} else {
		e_re = cr->op_sub(cr->op_mul(a_re, b_re),
			cr->op_mul(a_im, b_im));
		e_im = cr->op_add(cr->op_mul(a_re, b_im),
			cr->op_mul(a_im, b_re));
	}
	if (memcmp(&g_re, &e_re, sizeof g_re) != 0
		|| memcmp(&g_im, &e_im, sizeof g_im) != 0) {
		if (*bad < MAX_DETAIL) {
			char line[64];
			sprintf(line,
				"  FPC_MUL mismatch, %s vector %lu:", tag, i);
			g_emit(line);
			emit_fpr("a_re ", a_re);
			emit_fpr("a_im ", a_im);
			emit_fpr("b_re ", b_re);
			emit_fpr("b_im ", b_im);
			emit_fpr("asm re", g_re);
			emit_fpr("C   re", e_re);
			emit_fpr("asm im", g_im);
			emit_fpr("C   im", e_im);
		}
		(*bad)++;
	}
}

static void
check_fpc_mul(void)
{
	const fprtest_ops *cr = &fndsa_fprtest_cref;
	unsigned long bad = 0;

	if (cr->op_fpc_mul == 0
		&& (cr->op_mul == 0 || cr->op_add == 0 || cr->op_sub == 0))
	{
		op_skip("FPC_MUL", "no separate C reference for this backend");
		return;
	}

	/* Rounding ties for fxp32.64 (raw values in the comments). Random
	   vectors almost never hit a tie, where the sticky limb p0 and the
	   parity of q0 decide. Vectors 0 and 2 put a nonzero p0 in one
	   product while the combined p0 is zero; vector 1 the reverse. */
#if FPRTEST_HAVE_SCALED
	{
		fpr ulp   = fpr_scaled(1, -64);           /* raw 1 */
		fpr nulp  = fpr_neg(ulp);
		fpr max63 = fpr_scaled(INT64_MAX, -64);   /* raw 2^63-1 */

		/* 0: imag = 1 + (2^63-1), tie, q0 even; ar*bi has p0=1. */
		fpc_vec(cr, "edge", 0, ulp, ulp, max63, ulp, &bad);
		/* 1: imag = 2^65 + (2^63+5): sticky lives only in the
		   combined p0 (both products have p1 low bits zero). */
		fpc_vec(cr, "edge", 1, fpr_scaled(1, -32), ulp,
			fpr_add(fpr_scaled(1, -1), fpr_scaled(5, -64)),
			fpr_scaled(2, -32), &bad);
		/* 2: real = (2^63-1) - (-1), tie, q0 even; the subtracted
		   product supplies the nonzero per-product p0. */
		fpc_vec(cr, "edge", 2, ulp, nulp, max63, ulp, &bad);
		/* 3: imag = 1 + (3*2^63-1), tie with q0 odd: rounds up. */
		fpc_vec(cr, "edge", 3, ulp, ulp,
			fpr_sub(fpr_scaled(3, -1), ulp), ulp, &bad);
		/* 4: imag = 1 + (2^63-2): one ulp below the tie. */
		fpc_vec(cr, "edge", 4, ulp, ulp,
			fpr_sub(fpr_scaled(1, -1), fpr_scaled(2, -64)),
			ulp, &bad);
		/* 5: negative tie: imag = -1 - (2^63-1) = -2^63. */
		fpc_vec(cr, "edge", 5, nulp, nulp, max63, ulp, &bad);
	}
#endif

	for (unsigned long i = 0; i < g_nvec; i++) {
		fpc_vec(cr, "random", i, gen_fpc(), gen_fpc(),
			gen_fpc(), gen_fpc(), &bad);
	}
	op_done("FPC_MUL", bad);
}

#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 && FNDSA_TWIDDLE_U64
static void
gm_cmul_vec(const char *name, const char *tag, unsigned long i,
	void (*dut)(uint32_t *, const uint32_t *,
	const uint64_t *), void (*cref)(uint32_t *, const uint32_t *,
	const uint64_t *), fpr a_re, fpr a_im, unsigned k, int alias,
	unsigned long *bad)
{
	fxp_fpc ga, ea, in;
	const uint64_t *gm = &GM[k << 1];

	in.z.re = a_re;
	in.z.im = a_im;
	if (alias) {
		ga = in;
		ea = in;
		dut(ga.w, ga.w, gm);
		cref(ea.w, ea.w, gm);
	} else {
		dut(ga.w, in.w, gm);
		cref(ea.w, in.w, gm);
	}

	if (memcmp(&ga, &ea, sizeof ga) != 0) {
		if (*bad < MAX_DETAIL) {
			char line[96];
			sprintf(line,
				"  %s mismatch, %s vector %lu k=%u%s:",
				name, tag, i, k, alias ? " alias" : "");
			g_emit(line);
			emit_fpr("a_re  ", a_re);
			emit_fpr("a_im  ", a_im);
			emit_fpr("asm re", ga.z.re);
			emit_fpr("C   re", ea.z.re);
			emit_fpr("asm im", ga.z.im);
			emit_fpr("C   im", ea.z.im);
		}
		(*bad)++;
	}
}

static void
check_gm_cmul(const char *name,
	void (*dut)(uint32_t *, const uint32_t *, const uint64_t *),
	void (*cref)(uint32_t *, const uint32_t *, const uint64_t *))
{
	static const unsigned edges[] = { 1, 2, 3, 4, 5, 511, 512, 1023 };
	unsigned long bad = 0;

	if (cref == 0) {
		op_skip(name, "no separate C reference for this backend");
		return;
	}

	for (unsigned i = 0; i < sizeof edges / sizeof edges[0]; i ++) {
		fpr a_re = gen_fpc();
		fpr a_im = gen_fpc();
		gm_cmul_vec(name, "edge", i, dut, cref,
			a_re, a_im, edges[i], 0, &bad);
		gm_cmul_vec(name, "edge", i, dut, cref,
			a_re, a_im, edges[i], 1, &bad);
	}
	for (unsigned long i = 0; i < g_nvec; i ++) {
		fpr a_re = gen_fpc();
		fpr a_im = gen_fpc();
		unsigned k = (unsigned)(rnd32() % 1023) + 1;
		gm_cmul_vec(name, "random", i, dut, cref,
			a_re, a_im, k, 0, &bad);
		gm_cmul_vec(name, "random", i, dut, cref,
			a_re, a_im, k, 1, &bad);
	}
	op_done(name, bad);
}

static void
check_fpr_mul_gm(void)
{
	static const unsigned edges[] = { 1, 2, 3, 4, 5, 511, 512, 1023 };
	unsigned long bad = 0;
	fpr (*cref)(fpr, uint64_t) = fndsa_fprtest_cref.op_fpr_mul_gm;

	if (cref == 0) {
		op_skip("MUL_GM", "no separate C reference for this backend");
		return;
	}

	for (unsigned i = 0; i < sizeof edges / sizeof edges[0]; i ++) {
		for (unsigned j = 0; j < 2; j ++) {
			fpr x = gen_fpc();
			uint64_t c = GM[(edges[i] << 1) + j];
			fpr g = d_fpr_mul_gm(x, c);
			fpr e = cref(x, c);
			if (memcmp(&g, &e, sizeof g) != 0) {
				if (bad < MAX_DETAIL) {
					char line[96];
					sprintf(line,
						"  MUL_GM mismatch, edge vector %u k=%u.%u:",
						i, edges[i], j);
					g_emit(line);
					emit_fpr("x  ", x);
					emit_fpr("asm", g);
					emit_fpr("C  ", e);
				}
				bad++;
			}
		}
	}
	for (unsigned long i = 0; i < g_nvec; i ++) {
		unsigned k = (unsigned)(rnd32() % 1023) + 1;
		unsigned j = (unsigned)(rnd32() & 1);
		fpr x = gen_fpc();
		uint64_t c = GM[(k << 1) + j];
		fpr g = d_fpr_mul_gm(x, c);
		fpr e = cref(x, c);
		if (memcmp(&g, &e, sizeof g) != 0) {
			if (bad < MAX_DETAIL) {
				char line[96];
				sprintf(line,
					"  MUL_GM mismatch, random vector %lu k=%u.%u:",
					i, k, j);
				g_emit(line);
				emit_fpr("x  ", x);
				emit_fpr("asm", g);
				emit_fpr("C  ", e);
			}
			bad++;
		}
	}
	op_done("MUL_GM", bad);
}

static void
check_gm(void)
{
	check_gm_cmul("CMUL_GM", d_fxp96_cmul_gm,
		fndsa_fprtest_cref.op_fxp96_cmul_gm);
	check_gm_cmul("CMUL_GMJ", d_fxp96_cmul_gmj,
		fndsa_fprtest_cref.op_fxp96_cmul_gmj);
	check_fpr_mul_gm();
}
#endif

static void
check_of32(void)
{
	unsigned long bad = 0;

	if (fndsa_fprtest_cref.op_of32 == 0) {
		op_skip("fpr_of32", "no separate C reference for this backend");
		return;
	}
	for (unsigned long i = 0; i < g_nvec; i++) {
		int32_t v = (int32_t)rnd32();
		fpr g = d_of32(v), e = fndsa_fprtest_cref.op_of32(v);
		if (memcmp(&g, &e, sizeof g) != 0) {
			if (bad < MAX_DETAIL) {
				char line[64];
				sprintf(line,
					"  fpr_of32 mismatch, vector %lu:", i);
				g_emit(line);
				emit_fpr("asm", g);
				emit_fpr("C  ", e);
			}
			bad++;
		}
	}
	op_done("fpr_of32", bad);
}

static void
check_scaled(void)
{
#if FPRTEST_HAVE_SCALED
	unsigned long bad = 0;

	if (fndsa_fprtest_cref.op_scaled == 0) {
		op_skip("fpr_scaled",
			"no separate C reference for this backend");
		return;
	}
	for (unsigned long i = 0; i < g_nvec; i++) {
		int64_t v = (int64_t)(int32_t)rnd32();
		int sc = (int)(rnd32() % 96) - 64;
		fpr g = d_scaled(v, sc);
		fpr e = fndsa_fprtest_cref.op_scaled(v, sc);
		if (memcmp(&g, &e, sizeof g) != 0) {
			if (bad < MAX_DETAIL) {
				char line[64];
				sprintf(line,
					"  fpr_scaled mismatch, vector %lu:",
					i);
				g_emit(line);
				emit_fpr("asm", g);
				emit_fpr("C  ", e);
			}
			bad++;
		}
	}
	op_done("fpr_scaled", bad);
#else
	op_skip("fpr_scaled", "no on-target implementation for this backend");
#endif
}

static void
check_mul_recip(const char *name, fpr (*dut)(fpr, fpr),
	fpr (*cref)(fpr, fpr), int p63)
{
	unsigned long bad = 0;

	if (cref == 0) {
		op_skip(name, "no separate C reference for this backend");
		return;
	}
	for (unsigned long i = 0; i < g_nvec; i++) {
		fpr x;
		fpr a = p63 ? fpr_div2e(gen_val(8), 40) : gen_val(8);
		fpr g, e;

#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 && FNDSA_KGEN_GAMMA
		if (!p63) {
			/* All normalization shifts k = 6..29; edge cases
			   first. */
			x = i < 8 ? recip_alpha55_edge(i) : gen_recip_alpha55(i);
		} else {
			x = fpr_add(fpr_of32(1), gen_nonneg(8));
		}
#else
		x = fpr_add(fpr_of32(1), gen_nonneg(8));
#endif
		g = dut(a, x);
		e = cref(a, x);
		if (memcmp(&g, &e, sizeof g) != 0) {
			if (bad < MAX_DETAIL) {
				char line[64];
				sprintf(line, "  %s mismatch, vector %lu:",
					name, i);
				g_emit(line);
				emit_fpr("a  ", a);
				emit_fpr("x  ", x);
				emit_fpr("asm", g);
				emit_fpr("C  ", e);
			}
			bad++;
		}
	}
	op_done(name, bad);
}

#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 \
	&& FNDSA_SIGN_SQRT_FUSED && FNDSA_KGEN_GAMMA
/* x in [q/110^2, q*110^2]; all even normalization shifts t = 4..30,
   edge cases first. */
static void
check_sqrt_inv_sigma(void)
{
	unsigned long bad = 0;
	void (*cref)(fpr *, fpr, unsigned) =
		fndsa_fprtest_cref.op_sqrt_inv_sigma;

	if (cref == 0) {
		op_skip("sqrt_inv_sigma", "no separate C reference for this backend");
		return;
	}
	for (unsigned long i = 0; i < g_nvec; i++) {
		fpr x, g, e;
		unsigned logn = 9 + (rnd32() & 1);

		x = i < 8 ? rsqrt_alpha110_edge(i) : gen_rsqrt_alpha110(i);
		d_sqrt_inv_sigma(&g, x, logn);
		cref(&e, x, logn);
		if (memcmp(&g, &e, sizeof g) != 0) {
			if (bad < MAX_DETAIL) {
				char line[80];

				sprintf(line,
					"  sqrt_inv_sigma mismatch, vector %lu logn=%u:",
					i, logn);
				g_emit(line);
				emit_fpr("x  ", x);
				emit_fpr("asm", g);
				emit_fpr("C  ", e);
			}
			bad++;
		}
	}
	op_done("sqrt_inv_sigma", bad);
}

#endif

/* --- entry point --------------------------------------------------------- */

int
fndsa_fprtest_run(void (*emit)(const char *), unsigned long nvec)
{
	static const char seed[] = "fndsa-fprtest";
	char line[64];

	g_emit = emit;
	g_nvec = nvec;
	g_fail_ops = 0;
	shake_init(&g_sc, 256);
	shake_inject(&g_sc, seed, sizeof seed - 1);
	shake_flip(&g_sc);
	g_zero = fpr_of32(0);

	sprintf(line, "sizeof(fpr)=%u, %lu vectors/op",
		(unsigned)sizeof(fpr), nvec);
	emit(line);

	/* Conversions first: the generator itself is built on them. */
	check_of32();
	check_scaled();

	check_bin("fpr_add", d_add, fndsa_fprtest_cref.op_add, 34, 0);
	check_bin("fpr_sub", d_sub, fndsa_fprtest_cref.op_sub, 34, 0);
	check_add_sub();
	check_bin("fpr_mul", d_mul, fndsa_fprtest_cref.op_mul, 34, 0);
	check_fpc_mul();
#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 && FNDSA_TWIDDLE_U64
	check_gm();
#endif
	check_un("fpr_sqr", d_sqr, fndsa_fprtest_cref.op_sqr, 34, 0, 0);
	/* Bounded exponents keep quotients in range; divisors are positive,
	   as in signing. */
	check_bin("fpr_div", d_div, fndsa_fprtest_cref.op_div, 8, 1);
	check_un("fpr_inv", d_inv, fndsa_fprtest_cref.op_inv, 8, 1, 1);
	check_mul_recip("MUL_RECIP", d_mul_recip,
		fndsa_fprtest_cref.op_mul_recip, 0);
#if FPRTEST_HAVE_MUL_RECIP_P63
	check_mul_recip("MUL_RECIP_P63", d_mul_recip_p63,
		fndsa_fprtest_cref.op_mul_recip_p63, 1);
#else
	op_skip("MUL_RECIP_P63", "not implemented by this backend");
#endif
	check_un("fpr_sqrt", d_sqrt, fndsa_fprtest_cref.op_sqrt, 8, 1, 0);
	check_un("MUL_INV_Q", d_mul_inv_q, fndsa_fprtest_cref.op_mul_inv_q,
		60, 0, 0);
#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 \
	&& FNDSA_SIGN_SQRT_FUSED && FNDSA_KGEN_GAMMA
	check_sqrt_inv_sigma();
#endif
	check_toint("fpr_rint", d_rint, fndsa_fprtest_cref.op_rint, 34);
	check_toint("fpr_floor", d_floor, fndsa_fprtest_cref.op_floor, 34);
	check_toint("fpr_trunc", d_trunc, fndsa_fprtest_cref.op_trunc, 34);

	return g_fail_ops;
}
