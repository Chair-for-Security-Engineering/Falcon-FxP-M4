/*
 * Constant tables for the 'prec' representation, built at startup.
 * The reference values are computed with MPFR from their definitions
 * (see fpr/ieee/fpr_const_ieee.c), not converted from the ieee doubles:
 *
 *   GM[2k+0] = cos(rev10(k)*pi/1024), GM[2k+1] = sin(rev10(k)*pi/1024)
 *   smoothz2n = sqrt(log(4n(1 + 1/eps))/pi)/sqrt(2*pi),
 *     eps = 1/sqrt(max(2, n/4)*2^64)
 *   SIGMA_MIN[logn] = smoothz2n
 *   INV_SIGMA[logn] = 1/(smoothz2n*1.17*sqrt(12289))
 *
 * A startup check compares them, rounded to binary64, with the ieee
 * tables. The target values come from the target glue.
 */

#include "sign_inner.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "target_glue.h"

/* The data tables (names expand to fndsa_prec_* via fpr_prec.h). */
fpr GM[2048];
#if FNDSA_FFT_NEW
fpr RNJ[512];
#endif
fpr_u INV_SIGMA[11];
fpr_u SIGMA_MIN[11];

/* Working precision for computing the constants. */
#define PREC_CONST_WORK_BITS   1000

static inline int
rev10(int k)
{
	int r = 0;

	for (int i = 0; i < 10; i ++) {
		r |= ((k >> i) & 1) << (9 - i);
	}
	return r;
}

/* smoothz2n (see file header) at working precision into rop. */
static void
compute_smoothz2n(mpfr_t rop, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t bitsec = (n / 4) < 2 ? 2 : (n / 4);
	mpfr_t eps, pi, sq2pi, t;

	mpfr_inits2(PREC_CONST_WORK_BITS, eps, pi, sq2pi, t, (mpfr_ptr)0);

	/* eps = 1/sqrt(bitsec*2^64) */
	mpfr_set_ui(eps, (unsigned long)bitsec, MPFR_RNDN);
	mpfr_mul_2si(eps, eps, 64, MPFR_RNDN);
	mpfr_rec_sqrt(eps, eps, MPFR_RNDN);

	mpfr_free_cache();
	mpfr_const_pi(pi, MPFR_RNDN);
	/* sq2pi = sqrt(2*pi) */
	mpfr_mul_si(sq2pi, pi, 2, MPFR_RNDN);
	mpfr_sqrt(sq2pi, sq2pi, MPFR_RNDN);

	/* t = sqrt(log(4n*(1 + 1/eps))/pi) */
	mpfr_si_div(t, 1, eps, MPFR_RNDN);
	mpfr_add_si(t, t, 1, MPFR_RNDN);
	mpfr_mul_si(t, t, (long)(4 * n), MPFR_RNDN);
	mpfr_log(t, t, MPFR_RNDN);
	mpfr_div(t, t, pi, MPFR_RNDN);
	mpfr_sqrt(t, t, MPFR_RNDN);

	mpfr_div(rop, t, sq2pi, MPFR_RNDN);

	mpfr_clears(eps, pi, sq2pi, t, (mpfr_ptr)0);
}

/* ==================================================================== */
/* Startup check: ref rounded to binary64 must equal the ieee table.   */

static void
check_const(const char *what, unsigned idx, mpfr_srcptr ref,
	int64_t m, int e)
{
	double expect = ldexp((double)m, e);
	double got = mpfr_get_d(ref, MPFR_RNDN);

	if (got != expect) {
		fprintf(stderr, "prec const self-check FAILED: %s[%u]: "
			"ref->double %.17g != ieee %.17g\n",
			what, idx, got, expect);
		abort();
	}
}

/* ieee INV_SIGMA and SIGMA_MIN mantissas (exponents -60 and -52, as in
   fpr/ieee/fpr_const_ieee.c), indexed by logn (entry 0 unused). */
static const int64_t INV_SIGMA_M[11] = {
	0,
	7961475618707097, 7851656902127320, 7746260754658859,
	7595833604889141, 7453842886538220, 7319528409832599,
	7192222552237877, 7071336252758509, 6956347512113097,
	6846791885593314
};
static const int64_t SIGMA_MIN_M[11] = {
	0,
	5028307297130123, 5098636688852518, 5168009084304506,
	5270355833453349, 5370752584786614, 5469306724145091,
	5566116128735780, 5661270305715104, 5754851361258101,
	5846934829975396
};

/* ==================================================================== */

/* Defined in ops_mpfr.c. */
void fndsa_prec_init_scalars(void);

/* Called from fndsa_prec_init() while allocations are uncollectable. */
void
fndsa_prec_init_consts(void)
{
	mpfr_t ang, work;

	fndsa_prec_init_scalars();

	mpfr_inits2(PREC_CONST_WORK_BITS, ang, work, (mpfr_ptr)0);

	/* GM twiddles. */
	mpfr_free_cache();
	for (int i = 0; i < 2048; i ++) {
		mpfr_init2(GM[i].ref, PREC_REF_BITS);
		prec_tgt_gm(GM[i].tgt, i);
		if (i < 2) {
			/* unused entries */
			mpfr_set_zero(GM[i].ref, +1);
			continue;
		}
		int k = rev10(i >> 1);
		mpfr_set_si(ang, k, MPFR_RNDN);
		mpfr_div_si(ang, ang, 1024, MPFR_RNDN);
		if ((i & 1) == 0) {
			mpfr_cospi(work, ang, MPFR_RNDN);
		} else {
			mpfr_sinpi(work, ang, MPFR_RNDN);
		}
		mpfr_set(GM[i].ref, work, MPFR_RNDN);
	}
	/* GM[3] = 1, GM[4] = GM[5] = sqrt(2)/2, GM[8]/GM[9] =
	   cos/sin(pi/8); a wrong bit reversal fails these. */
	check_const("GM", 3, GM[3].ref, 4503599627370496, -52);
	check_const("GM", 4, GM[4].ref, 6369051672525773, -53);
	check_const("GM", 5, GM[5].ref, 6369051672525773, -53);
	check_const("GM", 8, GM[8].ref, 8321567036706118, -53);
	check_const("GM", 9, GM[9].ref, 6893811853601123, -54);

#if FNDSA_FFT_NEW
	/* Compact-FFT recombination constants, for m = 2^k and j < m:
	   RNJ[m+j] = 1/(2*cos((4*rev_k(j)+1)*pi/2^(k+2))); RNJ[0] = sqrt(2). */
	mpfr_init2(RNJ[0].ref, PREC_REF_BITS);
	mpfr_sqrt_ui(RNJ[0].ref, 2, MPFR_RNDN);
	prec_tgt_rnj(RNJ[0].tgt, 0);
	for (int k = 0; k < 9; k ++) {
		int m = 1 << k;
		for (int j = 0; j < m; j ++) {
			int rj = 0;
			for (int b = 0; b < k; b ++) {
				rj |= ((j >> b) & 1) << (k - 1 - b);
			}
			mpfr_init2(RNJ[m + j].ref, PREC_REF_BITS);
			prec_tgt_rnj(RNJ[m + j].tgt, m + j);
			mpfr_set_si(ang, 4 * rj + 1, MPFR_RNDN);
			mpfr_div_2si(ang, ang, k + 2, MPFR_RNDN);
			mpfr_cospi(work, ang, MPFR_RNDN);
			mpfr_mul_2si(work, work, 1, MPFR_RNDN);
			mpfr_si_div(RNJ[m + j].ref, 1, work, MPFR_RNDN);
		}
	}
	/* Anchors against the ieee table: RNJ[0] = sqrt(2), RNJ[1] =
	   sqrt(2)/2, RNJ[2] = 1/(2*cos(pi/8)), RNJ[3] = 1/(2*cos(5*pi/8)). */
	check_const("RNJ", 0, RNJ[0].ref, 6369051672525773, -52);
	check_const("RNJ", 1, RNJ[1].ref, 6369051672525773, -53);
	check_const("RNJ", 2, RNJ[2].ref, 4874661109905557, -53);
	check_const("RNJ", 3, RNJ[3].ref, -5884236481753340, -52);
#endif

	/* INV_SIGMA / SIGMA_MIN. */
	mpfr_init2(INV_SIGMA[0].f.ref, PREC_REF_BITS);
	mpfr_set_zero(INV_SIGMA[0].f.ref, +1);
	INV_SIGMA[0].f.tgt[0] = 0;
	INV_SIGMA[0].f.tgt[1] = 0;
	mpfr_init2(SIGMA_MIN[0].f.ref, PREC_REF_BITS);
	mpfr_set_zero(SIGMA_MIN[0].f.ref, +1);
	SIGMA_MIN[0].f.tgt[0] = 0;
	SIGMA_MIN[0].f.tgt[1] = 0;
	for (unsigned logn = 1; logn <= 10; logn ++) {
		mpfr_t gs_norm;

		compute_smoothz2n(work, logn);

		mpfr_init2(SIGMA_MIN[logn].f.ref, PREC_REF_BITS);
		prec_tgt_sigma_min(SIGMA_MIN[logn].f.tgt, logn);
		mpfr_set(SIGMA_MIN[logn].f.ref, work, MPFR_RNDN);
		check_const("SIGMA_MIN", logn, SIGMA_MIN[logn].f.ref,
			SIGMA_MIN_M[logn], -52);

		/* gs_norm = (117/100)*sqrt(12289);
		   INV_SIGMA = 1/(smoothz2n*gs_norm) */
		mpfr_init2(gs_norm, PREC_CONST_WORK_BITS);
		mpfr_set_si(gs_norm, 12289, MPFR_RNDN);
		mpfr_sqrt(gs_norm, gs_norm, MPFR_RNDN);
		mpfr_mul_si(gs_norm, gs_norm, 117, MPFR_RNDN);
		mpfr_div_si(gs_norm, gs_norm, 100, MPFR_RNDN);
		mpfr_mul(work, work, gs_norm, MPFR_RNDN);

		mpfr_init2(INV_SIGMA[logn].f.ref, PREC_REF_BITS);
		prec_tgt_inv_sigma(INV_SIGMA[logn].f.tgt, logn);
		mpfr_si_div(INV_SIGMA[logn].f.ref, 1, work, MPFR_RNDN);
		check_const("INV_SIGMA", logn, INV_SIGMA[logn].f.ref,
			INV_SIGMA_M[logn], -60);

		mpfr_clear(gs_norm);
	}

	mpfr_clears(ang, work, (mpfr_ptr)0);
}
