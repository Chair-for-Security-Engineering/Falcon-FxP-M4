/*
 * Generator for the fxp32.64 constant tables (MPFR).
 *
 * Build/run:  make gen-const        (see the top-level Makefile)
 * Manually:   cc -O2 -o gen_const gen_const.c -lmpfr -lgmp
 *             ./gen_const --twiddle-u64 fpr/fxp32.64 32 64
 *             ./gen_const --check-ieee
 *
 * The ieee tables carry binary64 rounding errors, so every constant is
 * recomputed here from its definition at 512 bits and rounded to nearest
 * (ties to even) onto the target grid. Arguments: output directory,
 * integer bits, fraction bits.
 *
 * Output:
 *   <dir>/fpr_const_fxp.c   GM[2048], INV_SIGMA[11], SIGMA_MIN[11]
 *                           and, for fxp32.64, INV_SIGMA_Q71[11] and
 *                           RNJ[512] (compiled only under FFT_NEW)
 *   <dir>/fxp_const_gen.h   INV_2SQRSIGMA0, LOG2, INV_LOG2, INV_Q,
 *                           MINUS_INV_Q, Q2 (compile-time FPR() macros)
 *                           and, for fxp32.64, INV_Q_Q77_HI/_LO
 *
 * With --twiddle-u64 (64 fraction bits only), GM is also emitted as
 * unsigned Q0.64 magnitudes, selected by FNDSA_TWIDDLE_U64.
 *
 * Definitions (as in fpr/ieee/fpr_const_ieee.c and the FN-DSA spec):
 *
 *   GM[2*j+0] = cos(rev(j)*pi/1024)     j = 1..1023, rev = bit reversal
 *   GM[2*j+1] = sin(rev(j)*pi/1024)     over 10 bits; GM[0], GM[1] unused
 *
 *   q = 12289, n = 2^logn                        (logn = 1..10)
 *   gs_norm   = (117/100)*sqrt(q)
 *   bitsec    = max(2, n/4)
 *   eps       = 1/sqrt(bitsec*2^64)
 *   smoothz2n = sqrt(log(4*n*(1 + 1/eps))/pi)/sqrt(2*pi)
 *   sigma     = smoothz2n*gs_norm
 *   SIGMA_MIN[logn] = smoothz2n
 *   INV_SIGMA[logn] = 1/sigma
 *
 *   INV_2SQRSIGMA0 = 1/(2*1.8205^2)   LOG2 = log(2)   INV_LOG2 = 1/log(2)
 *   INV_Q = 1/q                       MINUS_INV_Q = -1/q
 *   Q2 = q^2 (integer, exact in every representation)
 *
 * --check-ieee prints the same values rounded to binary64, to compare
 * with the ieee tables (made independently, with Sage).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <gmp.h>
#include <mpfr.h>

#define PREC   512

static int INT_BITS = 64;
static int FRAC_BITS = 64;
static int TOTAL_BITS = 128;
static int WORDS = 4;
static int TWIDDLE_U64 = 0;

/* Bit reversal over 10 bits. */
static unsigned
rev10(unsigned x)
{
	unsigned r = 0;

	for (int i = 0; i < 10; i ++) {
		r = (r << 1) | ((x >> i) & 1);
	}
	return r;
}

/* Bit reversal over n bits, with n <= 10. */
static unsigned
revn(unsigned x, unsigned n)
{
	unsigned r = 0;

	while (n -- > 0) {
		r = (r << 1) | (x & 1);
		x >>= 1;
	}
	return r;
}

/* Round v onto the grid (to nearest, ties to even) and return its raw
   two's-complement limbs. Exits if the value does not fit. */
static void
to_limbs(uint32_t *w, mpfr_srcptr v, const char *what)
{
	mpfr_t t;
	mpz_t z, lim;

	mpfr_init2(t, PREC);
	mpz_init(z);
	mpz_init(lim);

	/* Scaling is exact, so the rint is the only rounding. */
	mpfr_mul_2ui(t, v, (unsigned long)FRAC_BITS, MPFR_RNDN);
	mpfr_rint(t, t, MPFR_RNDN);
	mpfr_get_z(z, t, MPFR_RNDN);

	mpz_setbit(lim, (mp_bitcnt_t)(TOTAL_BITS - 1));
	if (mpz_cmpabs(z, lim) > 0
		|| (mpz_cmpabs(z, lim) == 0 && mpz_sgn(z) > 0))
	{
		fprintf(stderr, "constant out of range: %s\n", what);
		exit(1);
	}

	/* Two's complement on TOTAL_BITS bits, split into 32-bit limbs
	   (least significant first). */
	mpz_fdiv_r_2exp(z, z, (mp_bitcnt_t)TOTAL_BITS);
	for (int i = 0; i < WORDS; i ++) {
		mpz_t t2;
		mpz_init(t2);
		mpz_fdiv_q_2exp(t2, z, (mp_bitcnt_t)(32 * i));
		mpz_fdiv_r_2exp(t2, t2, 32);
		w[i] = (uint32_t)mpz_get_ui(t2);
		mpz_clear(t2);
	}

	mpfr_clear(t);
	mpz_clear(z);
	mpz_clear(lim);
}

static void
print_fpr(FILE *f, mpfr_srcptr v, const char *what)
{
	uint32_t w[8];

	to_limbs(w, v, what);
	fprintf(f, "FPR(");
	for (int i = 0; i < WORDS; i ++) {
		fprintf(f, "%s0x%08lX", i == 0 ? "" : ", ", (unsigned long)w[i]);
	}
	fprintf(f, ")");
}

/* Round a nonnegative v onto the Q0.64 grid as one 64-bit word; equals
   the low two limbs from to_limbs. Exits if it does not fit. */
static uint64_t
to_u64(mpfr_srcptr v, const char *what)
{
	uint32_t w[8];

	to_limbs(w, v, what);
	for (int i = 2; i < WORDS; i ++) {
		if (w[i] != 0) {
			fprintf(stderr, "constant does not fit u64: %s\n",
				what);
			exit(1);
		}
	}
	return (uint64_t)w[0] | ((uint64_t)w[1] << 32);
}

/* round(v * 2^sh) as one unsigned 64-bit word, for the normalized
   constants (INV_SIGMA_Q71, INV_Q_Q77). Exits if it does not fit. */
static uint64_t
to_u64_qn(mpfr_srcptr v, unsigned sh, const char *what)
{
	mpfr_t t;
	mpz_t z, limb;
	uint32_t lo, hi;

	mpfr_init2(t, PREC);
	mpz_init(z);
	mpz_init(limb);

	mpfr_mul_2ui(t, v, sh, MPFR_RNDN);
	mpfr_rint(t, t, MPFR_RNDN);
	mpfr_get_z(z, t, MPFR_RNDN);
	if (mpz_sgn(z) < 0 || mpz_sizeinbase(z, 2) > 64) {
		fprintf(stderr, "constant out of uint64 range: %s\n", what);
		exit(1);
	}

	mpz_fdiv_r_2exp(limb, z, 32);
	lo = (uint32_t)mpz_get_ui(limb);
	mpz_fdiv_q_2exp(limb, z, 32);
	hi = (uint32_t)mpz_get_ui(limb);

	mpfr_clear(t);
	mpz_clear(z);
	mpz_clear(limb);
	return ((uint64_t)hi << 32) | (uint64_t)lo;
}

/* --- the constants ---------------------------------------------------- */

/* v <- cos or sin of rev(j)*pi/1024 (im = 0 for cos, 1 for sin). */
static void
gm_entry(mpfr_ptr v, unsigned j, int im)
{
	mpfr_t a, pi;

	if (j == 0) {
		/* GM[0] and GM[1] are unused (zero, as in the ieee table). */
		mpfr_set_ui(v, 0, MPFR_RNDN);
		return;
	}
	mpfr_init2(a, PREC);
	mpfr_init2(pi, PREC);
	mpfr_const_pi(pi, MPFR_RNDN);
	mpfr_mul_ui(a, pi, (unsigned long)rev10(j), MPFR_RNDN);
	mpfr_div_ui(a, a, 1024, MPFR_RNDN);
	if (im) {
		mpfr_sin(v, a, MPFR_RNDN);
	} else {
		mpfr_cos(v, a, MPFR_RNDN);
	}
	mpfr_clear(a);
	mpfr_clear(pi);
}

/* RNJ[m + j] = 1/(2*cos((4*rev(j)+1)*pi/2^(log2(m)+2))), the
   recombination constants of the compact self-adjoint FFT; RNJ[0] is
   sqrt(2) for the first merge. */
static void
rnj_entry(mpfr_ptr v, unsigned j)
{
	unsigned lm = 0, m = 1;
	mpfr_t a, pi;

	if (j == 0) {
		mpfr_sqrt_ui(v, 2, MPFR_RNDN);
		return;
	}
	while ((m << 1) <= j) {
		m <<= 1;
		lm ++;
	}
	mpfr_init2(a, PREC);
	mpfr_init2(pi, PREC);
	mpfr_const_pi(pi, MPFR_RNDN);
	mpfr_mul_ui(a, pi, 4 * revn(j - m, lm) + 1, MPFR_RNDN);
	mpfr_div_2ui(a, a, lm + 2, MPFR_RNDN);
	mpfr_cos(a, a, MPFR_RNDN);
	mpfr_mul_ui(a, a, 2, MPFR_RNDN);
	mpfr_ui_div(v, 1, a, MPFR_RNDN);
	mpfr_clear(a);
	mpfr_clear(pi);
}

/* smoothz2n(logn) = sqrt(log(4*n*(1 + 1/eps))/pi)/sqrt(2*pi), with
   eps = 1/sqrt(max(2, n/4)*2^64). This is SIGMA_MIN[logn]. */
static void
smoothz2n(mpfr_ptr v, int logn)
{
	mpfr_t t, u, pi, bitsec;
	unsigned long n = 1UL << logn;

	mpfr_init2(t, PREC);
	mpfr_init2(u, PREC);
	mpfr_init2(pi, PREC);
	mpfr_init2(bitsec, PREC);
	mpfr_const_pi(pi, MPFR_RNDN);

	/* bitsec = max(2, n/4) */
	mpfr_set_ui(bitsec, n, MPFR_RNDN);
	mpfr_div_ui(bitsec, bitsec, 4, MPFR_RNDN);
	if (mpfr_cmp_ui(bitsec, 2) < 0) {
		mpfr_set_ui(bitsec, 2, MPFR_RNDN);
	}

	/* 1/eps = sqrt(bitsec*2^64) */
	mpfr_mul_2ui(t, bitsec, 64, MPFR_RNDN);
	mpfr_sqrt(t, t, MPFR_RNDN);

	/* u = 4*n*(1 + 1/eps) */
	mpfr_add_ui(u, t, 1, MPFR_RNDN);
	mpfr_mul_ui(u, u, 4 * n, MPFR_RNDN);

	/* v = sqrt(log(u)/pi)/sqrt(2*pi) */
	mpfr_log(u, u, MPFR_RNDN);
	mpfr_div(u, u, pi, MPFR_RNDN);
	mpfr_sqrt(u, u, MPFR_RNDN);
	mpfr_mul_ui(t, pi, 2, MPFR_RNDN);
	mpfr_sqrt(t, t, MPFR_RNDN);
	mpfr_div(v, u, t, MPFR_RNDN);

	mpfr_clear(t);
	mpfr_clear(u);
	mpfr_clear(pi);
	mpfr_clear(bitsec);
}

/* 1/sigma, sigma = smoothz2n(logn) * (117/100)*sqrt(12289). */
static void
inv_sigma(mpfr_ptr v, int logn)
{
	mpfr_t t;

	mpfr_init2(t, PREC);
	smoothz2n(v, logn);
	mpfr_set_ui(t, 12289, MPFR_RNDN);
	mpfr_sqrt(t, t, MPFR_RNDN);
	mpfr_mul_ui(t, t, 117, MPFR_RNDN);
	mpfr_div_ui(t, t, 100, MPFR_RNDN);
	mpfr_mul(v, v, t, MPFR_RNDN);
	mpfr_ui_div(v, 1, v, MPFR_RNDN);
	mpfr_clear(t);
}

static void
scalar(mpfr_ptr v, const char *name)
{
	if (strcmp(name, "INV_2SQRSIGMA0") == 0) {
		/* 1/(2*sigma0^2), sigma0 = 1.8205 (exactly, as a decimal) */
		mpfr_set_ui(v, 18205, MPFR_RNDN);
		mpfr_div_ui(v, v, 10000, MPFR_RNDN);
		mpfr_sqr(v, v, MPFR_RNDN);
		mpfr_mul_ui(v, v, 2, MPFR_RNDN);
		mpfr_ui_div(v, 1, v, MPFR_RNDN);
	} else if (strcmp(name, "LOG2") == 0) {
		mpfr_const_log2(v, MPFR_RNDN);
	} else if (strcmp(name, "INV_LOG2") == 0) {
		mpfr_const_log2(v, MPFR_RNDN);
		mpfr_ui_div(v, 1, v, MPFR_RNDN);
	} else if (strcmp(name, "INV_Q") == 0) {
		mpfr_set_ui(v, 1, MPFR_RNDN);
		mpfr_div_ui(v, v, 12289, MPFR_RNDN);
	} else if (strcmp(name, "MINUS_INV_Q") == 0) {
		mpfr_set_si(v, -1, MPFR_RNDN);
		mpfr_div_ui(v, v, 12289, MPFR_RNDN);
	} else if (strcmp(name, "Q2") == 0) {
		mpfr_set_ui(v, 12289, MPFR_RNDN);
		mpfr_sqr(v, v, MPFR_RNDN);
	} else {
		fprintf(stderr, "unknown scalar %s\n", name);
		exit(1);
	}
}

static const char *SCALARS[] = {
	"INV_2SQRSIGMA0", "LOG2", "INV_LOG2", "INV_Q", "MINUS_INV_Q", "Q2", NULL
};

/* --- emit ------------------------------------------------------------- */

static void
emit_tables(const char *dir, const char *rep)
{
	char path[512];
	FILE *f;
	mpfr_t v;

	mpfr_init2(v, PREC);

	snprintf(path, sizeof path, "%s/fpr_const_fxp.c", dir);
	f = fopen(path, "w");
	if (f == NULL) {
		perror(path);
		exit(1);
	}
	fprintf(f,
"/*\n"
" * Constant tables for the '%s' representation (%d integer bits,\n"
" * %d fraction bits; value = raw * 2^-%d).\n"
" *\n"
" * Generated by fpr/gen/gen_const.c; do not edit, run 'make gen-const'.\n"
" * Values are computed with MPFR at 512 bits and rounded to nearest\n"
" * (ties to even); limbs are least significant first. The symbol\n"
" * renames are in fpr/%s/fpr_fxp.h.\n"
" */\n"
"\n"
"#include \"sign_inner.h\"\n"
"\n"
"/* FFT twiddle factors: with rev() the bit reversal over 10 bits,\n"
"     GM[2*j+0] = cos(rev(j)*pi/1024)\n"
"     GM[2*j+1] = sin(rev(j)*pi/1024)\n"
"   for j = 1 to 1023. GM[0] and GM[1] are not used (zero here, as in the\n"
"   ieee table). */\n", rep, INT_BITS, FRAC_BITS, FRAC_BITS, rep);
	if (TWIDDLE_U64) {
		fprintf(f,
"/* With FNDSA_TWIDDLE_U64 the twiddles are unsigned Q0.64 magnitudes\n"
"   (value = word * 2^-64). sin is always positive; cos is negative for\n"
"   odd pair indices, which the FFT knows from its loop counter. Pairs 0\n"
"   and 1 are unused. */\n"
"#if FNDSA_TWIDDLE_U64\n"
"const uint64_t GM[2048] = {\n");
		for (unsigned j = 0; j < 1024; j ++) {
			uint64_t c, s;

			if (j < 2) {
				fprintf(f, "\t0, 0,\n");
				continue;
			}
			gm_entry(v, j, 0);
			if ((mpfr_sgn(v) < 0) != (int)(j & 1)) {
				fprintf(stderr,
					"GM cos sign != index parity at %u\n",
					j);
				exit(1);
			}
			mpfr_abs(v, v, MPFR_RNDN);
			c = to_u64(v, "GM cos");
			gm_entry(v, j, 1);
			if (mpfr_sgn(v) <= 0) {
				fprintf(stderr,
					"GM sin not positive at %u\n", j);
				exit(1);
			}
			s = to_u64(v, "GM sin");
			fprintf(f, "\t0x%016llX, 0x%016llX,\n",
				(unsigned long long)c, (unsigned long long)s);
		}
		fprintf(f, "};\n#else\n");
	}
	fprintf(f, "const fpr GM[2048] = {\n");
	for (unsigned j = 0; j < 1024; j ++) {
		fprintf(f, "\t");
		gm_entry(v, j, 0);
		print_fpr(f, v, "GM cos");
		fprintf(f, ",\n\t");
		gm_entry(v, j, 1);
		print_fpr(f, v, "GM sin");
		fprintf(f, ",\n");
	}
	fprintf(f, "};\n");
	if (TWIDDLE_U64) {
		fprintf(f, "#endif\n");
	}
	fprintf(f, "\n");

	/* Only fxp32.64 has the compact-FFT table. */
	if (INT_BITS == 32 && FRAC_BITS == 64) {
		fprintf(f,
"#if FNDSA_FFT_NEW\n"
"/* Compact self-adjoint FFT recombination constants.  For m = 2^k:\n"
"     RNJ[m+j] = 1/(2*cos((4*rev_k(j)+1)*pi/2^(k+2))).\n"
"   RNJ[0] holds sqrt(2), the first merge's constant. */\n"
"const fpr RNJ[512] = {\n");
		for (unsigned j = 0; j < 512; j ++) {
			fprintf(f, "\t");
			rnj_entry(v, j);
			print_fpr(f, v, "RNJ");
			fprintf(f, ",\n");
		}
		fprintf(f, "};\n#endif\n\n");
	}
	fprintf(f,
"/* Sampler constants, indexed by logn (index 0 unused):\n"
"     SIGMA_MIN[logn] = smoothz2n = sqrt(log(4*n*(1 + 1/eps))/pi)/sqrt(2*pi)\n"
"     INV_SIGMA[logn] = 1/(smoothz2n*(117/100)*sqrt(12289))\n"
"   with n = 2^logn, eps = 1/sqrt(max(2, n/4)*2^64). */\n"
"const fpr_u INV_SIGMA[11] = {\n\t{ FPR_ZERO },   /* unused */\n");
	for (int logn = 1; logn <= 10; logn ++) {
		inv_sigma(v, logn);
		fprintf(f, "\t{ ");
		print_fpr(f, v, "INV_SIGMA");
		mpfr_fprintf(f, " },   /* %.25Rg */\n", v);
	}
	fprintf(f, "};\n\n");
	if (INT_BITS == 32 && FRAC_BITS == 64) {
		fprintf(f,
	"/* INV_SIGMA as unsigned Q71 words, round(value*2^71), for the fused\n"
	"   square root; seven more significant bits than the Q64 table. */\n"
	"const uint64_t INV_SIGMA_Q71[11] = {\n"
	"\t0,   /* unused */\n");
		for (int logn = 1; logn <= 10; logn ++) {
			uint64_t u;

			inv_sigma(v, logn);
			u = to_u64_qn(v, 71, "INV_SIGMA_Q71");
			fprintf(f, "\t0x%016llX,   /* %d */\n",
				(unsigned long long)u, logn);
		}
		fprintf(f, "};\n\n");
	}
	fprintf(f, "const fpr_u SIGMA_MIN[11] = {\n\t{ FPR_ZERO },   /* unused */\n");
	for (int logn = 1; logn <= 10; logn ++) {
		smoothz2n(v, logn);
		fprintf(f, "\t{ ");
		print_fpr(f, v, "SIGMA_MIN");
		mpfr_fprintf(f, " },   /* %.25Rg */\n", v);
	}
	fprintf(f, "};\n");
	fclose(f);
	printf("wrote %s\n", path);

	snprintf(path, sizeof path, "%s/fxp_const_gen.h", dir);
	f = fopen(path, "w");
	if (f == NULL) {
		perror(path);
		exit(1);
	}
	fprintf(f,
"#ifndef FNDSA_FXP_CONST_GEN_H__\n"
"#define FNDSA_FXP_CONST_GEN_H__\n"
"\n"
"/*\n"
" * Scalar constants for the '%s' representation (%d.%d fixed point).\n"
" *\n"
" * Generated by fpr/gen/gen_const.c; do not edit, run 'make gen-const'.\n"
" * Included by fpr/%s/fpr_fxp.h, which defines FPR().\n"
" */\n"
"\n", rep, INT_BITS, FRAC_BITS, rep);
	for (int i = 0; SCALARS[i] != NULL; i ++) {
		scalar(v, SCALARS[i]);
		mpfr_fprintf(f, "/* %.25Rg */\n#define %-16s ", v, SCALARS[i]);
		print_fpr(f, v, SCALARS[i]);
		fprintf(f, "\n");
	}
	if (INT_BITS == 32 && FRAC_BITS == 64) {
		uint64_t u;

		scalar(v, "INV_Q");
		u = to_u64_qn(v, 77, "INV_Q_Q77");
		fprintf(f,
"\n/* round(2^77/q) as an unsigned Q77 word in two halves: 64 significant\n"
"   bits where INV_Q keeps 51. */\n"
"#define INV_Q_Q77_HI     0x%08lX\n"
"#define INV_Q_Q77_LO     0x%08lX\n",
			(unsigned long)(u >> 32), (unsigned long)(u & 0xFFFFFFFF));
	}
	fprintf(f, "\n#endif\n");
	fclose(f);
	printf("wrote %s\n", path);

	mpfr_clear(v);
}

/* The same values rounded to binary64; matching the independently made
   ieee tables confirms the definitions, notably the GM index mapping. */
static void
check_ieee(void)
{
	mpfr_t v;
	union { double d; uint64_t u; } c;

	mpfr_init2(v, PREC);
	for (unsigned j = 0; j < 1024; j ++) {
		for (int im = 0; im < 2; im ++) {
			gm_entry(v, j, im);
			c.d = mpfr_get_d(v, MPFR_RNDN);
			printf("GM[%u] %016llx\n", 2 * j + (unsigned)im,
				(unsigned long long)c.u);
		}
	}
	for (int logn = 1; logn <= 10; logn ++) {
		inv_sigma(v, logn);
		c.d = mpfr_get_d(v, MPFR_RNDN);
		printf("INV_SIGMA[%d] %016llx\n", logn,
			(unsigned long long)c.u);
	}
	for (int logn = 1; logn <= 10; logn ++) {
		smoothz2n(v, logn);
		c.d = mpfr_get_d(v, MPFR_RNDN);
		printf("SIGMA_MIN[%d] %016llx\n", logn,
			(unsigned long long)c.u);
	}
	for (int i = 0; SCALARS[i] != NULL; i ++) {
		scalar(v, SCALARS[i]);
		c.d = mpfr_get_d(v, MPFR_RNDN);
		printf("%s %016llx\n", SCALARS[i], (unsigned long long)c.u);
	}
	mpfr_clear(v);
}

int
main(int argc, char *argv[])
{
	const char *dir, *rep;

	if (argc == 2 && strcmp(argv[1], "--check-ieee") == 0) {
		check_ieee();
		return 0;
	}
	if (argc >= 2 && strcmp(argv[1], "--twiddle-u64") == 0) {
		TWIDDLE_U64 = 1;
		argv ++;
		argc --;
	}
	if (argc != 4) {
		fprintf(stderr,
			"usage: %s [--twiddle-u64] <dir> <int_bits> <frac_bits>\n"
			"       %s --check-ieee\n", argv[0], argv[0]);
		return 1;
	}
	dir = argv[1];
	INT_BITS = atoi(argv[2]);
	FRAC_BITS = atoi(argv[3]);
	TOTAL_BITS = INT_BITS + FRAC_BITS;
	if (TOTAL_BITS % 32 != 0 || TOTAL_BITS <= 0 || TOTAL_BITS > 256) {
		fprintf(stderr, "bad widths %d.%d\n", INT_BITS, FRAC_BITS);
		return 1;
	}
	if (TWIDDLE_U64 && FRAC_BITS != 64) {
		fprintf(stderr, "--twiddle-u64 requires 64 fraction bits\n");
		return 1;
	}
	WORDS = TOTAL_BITS / 32;

	rep = strrchr(dir, '/');
	rep = rep == NULL ? dir : rep + 1;
	emit_tables(dir, rep);
	return 0;
}
