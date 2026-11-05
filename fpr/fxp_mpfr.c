/*
 * MPFR/GMP oracle glue for the fxp representation (see fxp_mpfr.h).
 */

#include <stdio.h>
#include <string.h>

#include "fxp_mpfr.h"

/* Cached constants: 2^TOTAL, 2^(TOTAL-1)-1 (max raw), -2^(TOTAL-1)
   (min raw). */
static mpz_t z_two_pow_total, z_raw_max, z_raw_min;
static int inited = 0;

void
fxp_oracle_init(void)
{
	if (inited) {
		return;
	}
	mpfr_set_default_prec(FXP_ORACLE_PREC);
	mpz_init(z_two_pow_total);
	mpz_init(z_raw_max);
	mpz_init(z_raw_min);
	mpz_setbit(z_two_pow_total, FXP_TOTAL_BITS);
	mpz_setbit(z_raw_max, FXP_TOTAL_BITS - 1);
	mpz_neg(z_raw_min, z_raw_max);
	mpz_sub_ui(z_raw_max, z_raw_max, 1);
	inited = 1;
}

void
fxp_to_mpz(mpz_t z, fpr x)
{
	mpz_import(z, FXP_WORDS, -1, sizeof(uint32_t), 0, 0, x.w);
	if (mpz_tstbit(z, FXP_TOTAL_BITS - 1)) {
		mpz_sub(z, z, z_two_pow_total);
	}
}

fpr
fxp_from_mpz(const mpz_t z)
{
	mpz_t t;
	fpr r;
	size_t count;

	mpz_init(t);
	/* fdiv_r yields the representative in [0, 2^TOTAL): exactly the
	   two's-complement wrap. */
	mpz_fdiv_r_2exp(t, z, FXP_TOTAL_BITS);
	memset(&r, 0, sizeof r);
	mpz_export(r.w, &count, -1, sizeof(uint32_t), 0, 0, t);
	mpz_clear(t);
	return r;
}

void
fxp_to_mpfr(mpfr_t v, fpr x)
{
	mpz_t z;

	mpz_init(z);
	fxp_to_mpz(z, x);
	mpfr_set_z(v, z, MPFR_RNDN);   /* exact: prec >= TOTAL bits */
	mpfr_mul_2si(v, v, -(long)FXP_FRAC_BITS, MPFR_RNDN);   /* exact */
	mpz_clear(z);
}

/* Common tail: v*2^FRAC (exact), plus 1/2 if half_up, rounded to an
   integer with rnd, range-checked and stored. */
static int
from_mpfr_scaled(fpr *out, const mpfr_t v, int half_up, mpfr_rnd_t rnd)
{
	mpfr_t t;
	mpz_t z;
	int ok;

	if (!mpfr_number_p(v)) {
		return 0;
	}
	mpfr_init2(t, FXP_ORACLE_PREC);
	mpfr_mul_2si(t, v, (long)FXP_FRAC_BITS, MPFR_RNDN);   /* exact */
	if (half_up) {
		mpfr_add_d(t, t, 0.5, MPFR_RNDN);   /* exact at oracle prec */
	}
	mpfr_rint(t, t, rnd);
	mpz_init(z);
	mpfr_get_z(z, t, MPFR_RNDN);   /* t is an integer: exact */
	ok = mpz_cmp(z, z_raw_min) >= 0 && mpz_cmp(z, z_raw_max) <= 0;
	if (ok) {
		*out = fxp_from_mpz(z);
	}
	mpz_clear(z);
	mpfr_clear(t);
	return ok;
}

int
fxp_from_mpfr(fpr *out, const mpfr_t v, mpfr_rnd_t rnd)
{
	return from_mpfr_scaled(out, v, 0, rnd);
}

int
fxp_from_mpfr_half_up(fpr *out, const mpfr_t v)
{
	return from_mpfr_scaled(out, v, 1, MPFR_RNDD);
}

int
fxp_eq(fpr a, fpr b)
{
	return memcmp(a.w, b.w, sizeof a.w) == 0;
}

void
fxp_ulp_dist(mpz_t d, fpr a, fpr b)
{
	mpz_t zb;

	mpz_init(zb);
	fxp_to_mpz(d, a);
	fxp_to_mpz(zb, b);
	mpz_sub(d, d, zb);
	mpz_abs(d, d);
	mpz_clear(zb);
}

/* ==================================================================== */
/* Self-test.                                                           */

static int st_failures;

static void
st_check(int cond, const char *what)
{
	if (!cond) {
		fprintf(stderr, "oracle self-test FAILED: %s\n", what);
		st_failures++;
	}
}

/* Build the real value raw_num * 2^-FRAC (raw_num given as a small
   integer times 2^shift so large values are expressible). */
static void
st_val(mpfr_t v, long num, long shift)
{
	mpfr_set_si(v, num, MPFR_RNDN);
	mpfr_mul_2si(v, v, shift - (long)FXP_FRAC_BITS, MPFR_RNDN);
}

/* Grid-round v with mode rnd and check the raw result equals num*2^shift
   (expect_ok = 1), or that the conversion reports overflow. */
static void
st_round_case(const char *what, mpfr_t v, mpfr_rnd_t rnd,
	int expect_ok, long num, long shift)
{
	fpr r;
	mpz_t z, e;
	int ok;

	ok = fxp_from_mpfr(&r, v, rnd);
	if (!expect_ok) {
		st_check(!ok, what);
		return;
	}
	if (!ok) {
		st_check(0, what);
		return;
	}
	mpz_init(z);
	mpz_init(e);
	fxp_to_mpz(z, r);
	mpz_set_si(e, num);
	if (shift >= 0) {
		mpz_mul_2exp(e, e, (mp_bitcnt_t)shift);
	} else {
		mpz_fdiv_q_2exp(e, e, (mp_bitcnt_t)(-shift));
	}
	st_check(mpz_cmp(z, e) == 0, what);
	mpz_clear(z);
	mpz_clear(e);
}

int
fxp_oracle_selftest(void)
{
	mpfr_t v;
	mpz_t z;
	fpr r;
	int i;

	fxp_oracle_init();
	st_failures = 0;
	mpfr_init2(v, FXP_ORACLE_PREC);
	mpz_init(z);

	/* Zero and one, both directions. */
	fxp_to_mpfr(v, FPR_ZERO);
	st_check(mpfr_zero_p(v), "FPR_ZERO decodes to 0");
	fxp_to_mpfr(v, FPR_ONE);
	st_check(mpfr_cmp_ui(v, 1) == 0, "FPR_ONE decodes to 1");
	mpfr_set_ui(v, 1, MPFR_RNDN);
	st_check(fxp_from_mpfr(&r, v, MPFR_RNDN) && fxp_eq(r, FPR_ONE),
		"1.0 encodes to FPR_ONE");

	/* -1: raw is -2^FRAC; check the two's-complement round trip. */
	mpfr_set_si(v, -1, MPFR_RNDN);
	st_check(fxp_from_mpfr(&r, v, MPFR_RNDN), "-1.0 encodes");
	fxp_to_mpz(z, r);
	mpz_set_si(z, 0);   /* recompute expected below */
	mpz_setbit(z, FXP_FRAC_BITS);
	mpz_neg(z, z);
	{
		mpz_t got;
		mpz_init(got);
		fxp_to_mpz(got, r);
		st_check(mpz_cmp(got, z) == 0, "-1.0 raw is -2^FRAC");
		mpz_clear(got);
	}
	fxp_to_mpfr(v, r);
	st_check(mpfr_cmp_si(v, -1) == 0, "-1.0 round-trips");

	/* Ties-to-even at the grid: n+0.5 ulp. */
	st_val(v, 5, -1);   /* 2.5 ulp */
	st_round_case("RNDN 2.5ulp -> 2 (even)", v, MPFR_RNDN, 1, 2, 0);
	st_val(v, 3, -1);   /* 1.5 ulp */
	st_round_case("RNDN 1.5ulp -> 2 (even)", v, MPFR_RNDN, 1, 2, 0);
	st_val(v, 1, -1);   /* 0.5 ulp */
	st_round_case("RNDN 0.5ulp -> 0 (even)", v, MPFR_RNDN, 1, 0, 0);
	st_val(v, -5, -1);  /* -2.5 ulp */
	st_round_case("RNDN -2.5ulp -> -2 (even)", v, MPFR_RNDN, 1, -2, 0);

	/* Directed modes. */
	st_val(v, 5, -1);
	st_round_case("RNDZ 2.5ulp -> 2", v, MPFR_RNDZ, 1, 2, 0);
	st_val(v, -5, -1);
	st_round_case("RNDZ -2.5ulp -> -2", v, MPFR_RNDZ, 1, -2, 0);
	st_round_case("RNDD -2.5ulp -> -3", v, MPFR_RNDD, 1, -3, 0);
	st_val(v, 5, -1);
	st_round_case("RNDU 2.5ulp -> 3", v, MPFR_RNDU, 1, 3, 0);

	/* Range boundaries. Max raw = 2^(TOTAL-1)-1; min raw = -2^(TOTAL-1). */
	mpfr_set_z(v, z_raw_max, MPFR_RNDN);
	mpfr_mul_2si(v, v, -(long)FXP_FRAC_BITS, MPFR_RNDN);
	st_check(fxp_from_mpfr(&r, v, MPFR_RNDN), "max value encodes");
	fxp_to_mpfr(v, r);
	mpfr_mul_2si(v, v, (long)FXP_FRAC_BITS, MPFR_RNDN);
	mpfr_get_z(z, v, MPFR_RNDN);
	st_check(mpz_cmp(z, z_raw_max) == 0, "max value round-trips");
	mpfr_set_z(v, z_raw_min, MPFR_RNDN);
	mpfr_mul_2si(v, v, -(long)FXP_FRAC_BITS, MPFR_RNDN);
	st_check(fxp_from_mpfr(&r, v, MPFR_RNDN), "min value encodes");
	mpfr_set_si(v, 1, MPFR_RNDN);
	mpfr_mul_2si(v, v, FXP_INT_BITS - 1, MPFR_RNDN);
	st_round_case("2^(INT-1) overflows", v, MPFR_RNDN, 0, 0, 0);
	mpfr_neg(v, v, MPFR_RNDN);
	st_check(fxp_from_mpfr(&r, v, MPFR_RNDN),
		"-2^(INT-1) is representable (min)");

	/* Raw round trips through fpr limbs for assorted bit patterns. */
	for (i = 0; i < 64; i++) {
		mpz_t got;
		mpz_set_ui(z, 0);
		mpz_setbit(z, (mp_bitcnt_t)((i * 37) % (FXP_TOTAL_BITS - 1)));
		mpz_setbit(z, (mp_bitcnt_t)((i * 11 + 5) % (FXP_TOTAL_BITS - 1)));
		if (i & 1) {
			mpz_neg(z, z);
		}
		r = fxp_from_mpz(z);
		mpz_init(got);
		fxp_to_mpz(got, r);
		st_check(mpz_cmp(got, z) == 0, "raw limb round trip");
		mpz_clear(got);
	}

	mpz_clear(z);
	mpfr_clear(v);
	return st_failures;
}
