#ifndef FNDSA_FXP_MPFR_H__
#define FNDSA_FXP_MPFR_H__

/*
 * MPFR/GMP oracle glue for the fxp representation, used only by the
 * native test_fpr harness. An fpr is a two's-complement integer "raw" of
 * FXP_TOTAL_BITS bits with value raw / 2^FXP_FRAC_BITS.
 *
 * At FXP_ORACLE_PREC bits, sums and products of grid values are exact,
 * and quotients and square roots are close enough that rounding them to
 * the grid gives the correctly rounded result.
 */

#include <stdint.h>
#define MPFR_USE_INTMAX_T   1
#include <gmp.h>
#include <mpfr.h>

#include "sign_inner.h"

#if FNDSA_FPR_REP != FNDSA_FPR_REP_FXP3264
#error the MPFR oracle glue only applies to the fxp representation
#endif

#define FXP_WORDS         (FPR_STRIDE / 4)
#define FXP_TOTAL_BITS    (FXP_INT_BITS + FXP_FRAC_BITS)
#define FXP_ORACLE_PREC   (2 * FXP_TOTAL_BITS + 64)

/* Initialize cached bignum constants and the default mpfr precision.
   Call once before anything else. */
void fxp_oracle_init(void);

/* fpr to signed raw integer (value * 2^FXP_FRAC_BITS). */
void fxp_to_mpz(mpz_t z, fpr x);

/* Signed raw integer to fpr, reduced mod 2^FXP_TOTAL_BITS. */
fpr fxp_from_mpz(const mpz_t z);

/* fpr to exact real value (v needs at least FXP_TOTAL_BITS). */
void fxp_to_mpfr(mpfr_t v, fpr x);

/* Real value to fpr, rounded to the grid with rnd. Returns 0 if v is
   out of range or not finite. */
int fxp_from_mpfr(fpr *out, const mpfr_t v, mpfr_rnd_t rnd);

/* Same, rounding to nearest with ties toward +inf. */
int fxp_from_mpfr_half_up(fpr *out, const mpfr_t v);

/* Bit compare of two fpr values. */
int fxp_eq(fpr a, fpr b);

/* d <- |raw(a) - raw(b)|, i.e. their distance in ulp (2^-FRAC units). */
void fxp_ulp_dist(mpz_t d, fpr a, fpr b);

/* Self-test against hand-worked values; returns 0 if all pass. */
int fxp_oracle_selftest(void);

#endif
