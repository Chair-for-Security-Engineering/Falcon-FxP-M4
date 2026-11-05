/*
 * Core signature generation.
 */

#include "sign_inner.h"
#include <inner.h>
#include <stdint.h>

/* There is no wide -1/q, so the c*F residues are negated instead. */
#if FNDSA_SIGN_INVQ_WIDE
#define INVQ_RESIDUE_SIGN   (-1)
#else
#define INVQ_RESIDUE_SIGN   1
#endif

#if !FNDSA_SIGN_INT_GRAM
/* Given f, g, F and G, return the basis [[g, -f], [G, -F]] in FFT
   format (b00, b01, b10 and b11, in that order, are written in the
   destination). */
static void
basis_to_FFT(unsigned logn,
	const int8_t *f, const int8_t *g, const int8_t *F, const int8_t *G,
	fpr *dst)
{
	size_t n = (size_t)1 << logn;
	fpr *b00 = dst;
	fpr *b01 = b00 + n;
	fpr *b10 = b01 + n;
	fpr *b11 = b10 + n;
#if FNDSA_FFT_NEW
	fpoly_FFT_int8(logn, b01, f);
	fpoly_FFT_int8(logn, b00, g);
	fpoly_FFT_int8(logn, b11, F);
	fpoly_FFT_int8(logn, b10, G);
#else
	fpoly_set_small(logn, b01, f);
	fpoly_set_small(logn, b00, g);
	fpoly_set_small(logn, b11, F);
	fpoly_set_small(logn, b10, G);
	fpoly_FFT(logn, b01);
	fpoly_FFT(logn, b00);
	fpoly_FFT(logn, b11);
	fpoly_FFT(logn, b10);
#endif
	fpoly_neg(logn, b01);
	fpoly_neg(logn, b11);
}
#endif

/* see sign_inner.h */
TARGET_SSE2 TARGET_NEON
size_t
sign_core(unsigned logn,
	const uint8_t *sign_key_fgF, const int8_t *G,
	const uint8_t *hashed_vk, const uint8_t *ctx, size_t ctx_len,
	const char *id, const uint8_t *hv, size_t hv_len,
	const uint8_t *seed, size_t seed_len, uint8_t *sig, void *tmp)
{
	/* Output value is 0 on error, or the signature length on success. */
	size_t ret = 0;

#if FNDSA_SSE2
	/* We must ensure that the rounding mode is appropriate (we need
	   the roundTiesToEven policy, which is normally the default, but
	   could have been changed by the calling application. */
	unsigned round_mode = _MM_GET_ROUNDING_MODE();
	_MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
#endif
	/* Note: on aarch64, we use only NEON intrinsics, and NEON always
	   uses "round to nearest" mode, which is the IEEE-754
	   roundTiesToEven (NEON also always flushes denormals to zero,
	   which is not IEEE-754 compliant, but we do not care because
	   we never get any denormal value anywhere in signature
	   generation). */

	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	unsigned nbits;
	switch (logn) {
	case 2: case 3: case 4: case 5:
		nbits = 8;
		break;
	case 6: case 7:
		nbits = 7;
		break;
	case 8: case 9:
		nbits = 6;
		break;
	default:
		nbits = 5;
		break;
	}
	/* Storage format of F is already an array of int8_t. */
	size_t flen = (nbits << logn) >> 3;
	int8_t *F = (int8_t *)(sign_key_fgF + (flen << 1));

	/* The special hash identifier consisting of a byte of value 0xFF
	   followed by a zero denotes the original Falcon algorithm, which
	   we support for now for test vector reproducibility.
	   TODO: remove when final test vectors are available. */
	int orig_falcon = (*(const uint8_t *)id == 0xFF && id[1] == 0);

	/* Buffer for randomness: we need a 40-byte nonce and a 56-byte
	   sub-seed for sampling. */
	uint8_t rndbuf[40 + 56];
	uint8_t *nonce = rndbuf;
	uint8_t *subseed = rndbuf + 40;

	/* TODO: add "hedging" by hashing together the private key, the
	   message and the seed, to get a new seed value (protection against
	   bad RNG). We don't need to use the complete private key; we can
	   use F. */

	for (uint32_t counter = 0;; counter ++) {
		/* Generate the nonce and the sub-seed. In the original
		   Falcon, the nonce was not regenerated in case of
		   restart, but regenerating it makes the algorithm security
		   easier to analyze and prove.

		   When working with an explicit seed: normally, we hash
		   together the seed and a loop counter. For test purposes,
		   if we are using the original Falcon behaviour, this is
		   the first iteration, and the seed length is exactly 96
		   bytes, then we use the seed directly. */
		uint8_t *rndp;
		size_t rndlen;
		if (counter == 0 || !orig_falcon) {
			rndp = rndbuf;
			rndlen = sizeof rndbuf;
		} else {
			rndp = rndbuf + 40;
			rndlen = (sizeof rndbuf) - 40;
		}
		if (seed == NULL) {
			if (!sysrng(rndp, rndlen)) {
				goto sign_exit;
			}
		} else if (orig_falcon && counter == 0 && seed_len == rndlen) {
			memcpy(rndp, seed, rndlen);
		} else {
			/* We can use the tmp buffer for the SHAKE context.
			   It even works at n = 4 (logn = 2) because there
			   are FNDSA_SIGN_G_OFF*4 bytes free in tmp[] at this
			   point, and we need only 208. The larger Huang
			   context uses a local copy when it does not fit. */
			shake_context *sc = (shake_context *)tmp;
	#if FNDSA_SHAKE_HUANG_BUF
			shake_context sc_buf;
			if ((FNDSA_SIGN_G_OFF << logn) < sizeof *sc) {
				sc = &sc_buf;
			}
#endif
			shake_init(sc, 256);
			shake_inject(sc, seed, seed_len);
			uint8_t cbuf[4];
			cbuf[0] = (uint8_t)counter;
			cbuf[1] = (uint8_t)(counter >> 8);
			cbuf[2] = (uint8_t)(counter >> 16);
			cbuf[3] = (uint8_t)(counter >> 24);
			shake_inject(sc, cbuf, 4);
			shake_flip(sc);
			shake_extract(sc, rndp, rndlen);
		}

		/* The nonce (salt) is published with the signature. */
		FNDSA_CT_PUBLIC(nonce, 40);

		/* Hash the message into a polynomial. */
		uint16_t *hm = (uint16_t *)((uint8_t *)tmp + FNDSA_SIGN_HM_OFF * n);
		hash_to_point(logn, nonce, hashed_vk,
			ctx, ctx_len, id, hv, hv_len, hm);

		/* Initialize a sampler state. */
		sampler_state ss;
		sampler_init(&ss, logn, subseed, 56);

		/* Compute the lattice basis B = [[g, -f], [G, -F]] in FFT
		   representation, then compute the Gram matrix G = B*adj(B):
		      g00 = b00*adj(b00) + b01*adj(b01)
		      g01 = b00*adj(b10) + b01*adj(b11)
		      g10 = b10*adj(b00) + b11*adj(b01)
		      g11 = b10*adj(b10) + b11*adj(b11)
		   Note that g10 = adj(g01). For historical reasons, this
		   implementation keeps g01, i.e. the "upper triangle",
		   omitting g10.
		   We want the following layout:
		      free space (2n)
		      g00 (n)
		      g01 (n)
		      g11 (n)
		      b11 (n)
		      b01 (n)  */
		int8_t *f = (int8_t *)tmp;
		int8_t *g = f + n;
		(void)trim_i8_decode(logn, sign_key_fgF, f, nbits);
		(void)trim_i8_decode(logn, sign_key_fgF + flen, g, nbits);
		fpr *t0 = (fpr *)tmp;
		fpr *t1 = t0 + n;
#if FNDSA_SIGN_INT_GRAM
		/* Exact integer Gram, with NTTs mod q and q':
		      a = g00 = f*adj(f) + g*adj(g)          (int16)
		      b = g01 = f*adj(F) + g*adj(G)          (int32)
		   then FFT into the layout below. Exact since keygen
		   ensures ||(f,g)||^2 < 16823, so |g00| < q' and
		   |g01| < q*q'/2. */
		fpr *g01 = t1 + n;
		fpr *g00 = g01 + n;
		/* Free here: scratch for the Gram, then the split target. */
		fpr *b10 = t1 + 3 * n;
		/* Integer scratch scoped to the Gram computation. */
		{
		uint16_t *w0 = (uint16_t *)b10;      /* f mod q', then b mod q' */
		uint16_t *w1 = w0 + n;               /* g mod q', then g*adj(G) */
		uint16_t *w2 = w1 + n;               /* f mod q,  then b mod q  */
		uint16_t *w3 = w2 + n;               /* g mod q */
		uint16_t *v0 = (uint16_t *)t0 + n;   /* a mod q' / F NTT       */
		uint16_t *v1 = v0 + n;               /* a mod q  / G NTT       */
		uint8_t *adj = (uint8_t *)(v1 + n);
		int32_t *bi = (int32_t *)v0;         /* b after CRT (n i32) */

		/* f and g mod q' and mod q (NTT), shared by a and b. */
		mq2poly_small_to_mod(logn, f, w0);
		mq2poly_NTT(logn, w0);
		mq2poly_small_to_mod(logn, g, w1);
		mq2poly_NTT(logn, w1);
		mqpoly_small_to_int(logn, f, w2);
		mqpoly_int_to_ntt(logn, w2);
		mqpoly_small_to_int(logn, g, w3);
		mqpoly_int_to_ntt(logn, w3);

		/* a = f*adj(f) + g*adj(g), mod q' and mod q, on copies:
		   the NTTs of f and g are needed again for g01. */
		memcpy(v0, w0, n * sizeof(uint16_t));
		mq2poly_muladj_x2_ntt(logn, v0, w1);
		mq2poly_iNTT(logn, v0);
		mq2poly_mod_to_signed(logn, v0);
		memcpy(v1, w2, n * sizeof(uint16_t));
		mqpoly_muladj_x2_ntt(logn, v1, w3);
		mqpoly_ntt_to_int(logn, v1);
		mqpoly_int_to_signed(logn, v1);

		/* Adjust the mod-q' value in v0 to the true int16 value,
		   on all n coefficients. */
		mq2mqpoly_adjust_save(logn,
			(int16_t *)v0, (int16_t *)v1, adj);

#if FNDSA_FFT_NEW
		/* g00 is self-adjoint: half-size integer FFT, using the
		   g01 band as scratch. */
		{
			int32_t *a32 = (int32_t *)(void *)g01;
			for (size_t i = 0; i < hn; i ++) {
				a32[i] = *(int16_t *)&v0[i];
			}
			fpoly_FFT_int32_selfadj(logn, g00, a32);
		}
#else
		for (size_t i = 0; i < n; i ++) {
			g01[i] = fpr_of(*(int16_t *)&v0[i]);
		}
		fpoly_FFT(logn, g01);
		memcpy(g00, g01, hn * sizeof(fpr));
#endif

		/* g01 = f*adj(F) + g*adj(G), via CRT into int32. Mod q': */
		mq2poly_small_to_mod(logn, F, v0);
		mq2poly_NTT(logn, v0);
		/* w0 <- f*adj(F) mod q' (NTT) */
		mq2poly_muladj_ntt(logn, w0, v0);
		mq2poly_small_to_mod(logn, G, v0);
		mq2poly_NTT(logn, v0);
		/* w1 <- g*adj(G) mod q' (NTT) */
		mq2poly_muladj_ntt(logn, w1, v0);
		mq2poly_add(logn, w0, w1);
		mq2poly_iNTT(logn, w0);
		mq2poly_mod_to_signed(logn, w0);

		/* Mod q. */
		mqpoly_small_to_int(logn, F, v0);
		mqpoly_int_to_ntt(logn, v0);
		mqpoly_small_to_int(logn, G, v1);
		mqpoly_int_to_ntt(logn, v1);
		/* w2 <- f*adj(F) + g*adj(G) mod q (normalized) */
		mqpoly_muladj_add_muladj(logn, w2, v0, w3, v1);
		mqpoly_ntt_to_int(logn, w2);
		mqpoly_int_to_signed(logn, w2);

		/* CRT to the true int32 value, then FFT. */
		mq2mqpoly_CRT(logn, bi, (int16_t *)w2, (int16_t *)w0);
#if FNDSA_FFT_NEW
		fpoly_FFT_int32(logn, g01, bi);
#else
		for (size_t i = 0; i < n; i ++) {
			g01[i] = fpr_of(bi[i]);
		}
		fpoly_FFT(logn, g01);
#endif
		}

#if !FNDSA_SYMP_ROOT && !FNDSA_SYMP_TREE
		/* The generic LDL needs g11 = (q^2 + |g01|^2)/g00, since
		   the Gram determinant is q^2. */
		{
			fpr *g11 = g00 + hn;
			fpr q2 = Q2;
			for (size_t i = 0; i < hn; i ++) {
				fpr g01_re = g01[i], g01_im = g01[i + hn];
				fpr_rcp inv_g00;
				fpr mu_re, mu_im, d11;
				FPR_RECIP_X(inv_g00, g00[i]);
				FPR_MUL_RECIP(mu_re, g01_re, inv_g00);
				FPR_MUL_RECIP(mu_im, g01_im, inv_g00);
				FPR_MUL_RECIP(d11, q2, inv_g00);
				g11[i] = fpr_add(d11, fpr_add(
					fpr_mul(mu_re, g01_re),
					fpr_mul(mu_im, g01_im)));
			}
		}
#endif

#else
		basis_to_FFT(logn, f, g, F, G, t1 + n);
		fpr *b00 = t1 + n;
		fpr *b01 = b00 + n;
		fpr *b10 = b01 + n;
		fpr *b11 = b10 + n;
#if !FNDSA_SIGN_SPLIT
		/* The Gram computation consumes b01; keep a copy in t1 for
		   the target. */
		memcpy(t1, b01, n * sizeof(fpr));
#endif
		fpoly_gram_fft(logn, b00, b01, b10, b11);

		/* We now move things a bit to get the following (taking
		   into account that g00 and g11 are self-adjoint, hence
		   half-size):
		      free space (2n)
		      g01 (n)
		      g00 (hn)
		      g11 (hn)
		      free space (n)
		      b11 (n)
		      b01 (n)  */
		fpr *g01 = b00;
		fpr *g00 = b01;
		fpr *g11 = b01 + hn;
		/* g00 is saved in the upper half of t0 (f and g are in
		   its first 2n bytes, b01 may be in t1). */
		fpr *sav = t0 + hn;
		memcpy(sav, b00, hn * sizeof(fpr));
		memcpy(g01, b01, n * sizeof(fpr));
		memcpy(g00, sav, hn * sizeof(fpr));
		memcpy(g11, b10, hn * sizeof(fpr));
#endif

#if FNDSA_SIGN_SPLIT
		/* Split target: q*t = (-c*F, c*f). Its centered mod-q
		   residues, divided by q, give t_frac. t0 is written last
		   since f and g are still in it. */
		uint16_t *uc = (uint16_t *)b10;
		uint16_t *qt0 = uc + n;
		uint16_t *qt1 = qt0 + n;
		memcpy(uc, hm, n * sizeof(uint16_t));
		mqpoly_ext_to_int(logn, uc);
		mqpoly_int_to_ntt(logn, uc);
		mqpoly_small_to_int(logn, F, qt0);
		mqpoly_small_to_int(logn, f, qt1);
		mqpoly_int_to_ntt(logn, qt0);
		mqpoly_int_to_ntt(logn, qt1);
		mqpoly_mul_ntt(logn, qt0, uc);
		mqpoly_mul_ntt(logn, qt1, uc);
		mqpoly_ntt_to_int(logn, qt0);
		mqpoly_ntt_to_int(logn, qt1);
		(void)mqpoly_sqnorm_int_to_signed(logn, qt0);
		(void)mqpoly_sqnorm_int_to_signed(logn, qt1);
		for (size_t i = 0; i < n; i ++) {
			t0[i] = fpr_of(INVQ_RESIDUE_SIGN * *(int16_t *)&qt0[i]);
			t1[i] = fpr_of(*(int16_t *)&qt1[i]);
		}
		fpoly_FFT(logn, t0);
		fpoly_FFT(logn, t1);
#if FNDSA_SIGN_INVQ_WIDE
		fpoly_mul_inv_q(logn, t0);
		fpoly_mul_inv_q(logn, t1);
#else
		fpoly_mulconst(logn, t0, MINUS_INV_Q);
		fpoly_mulconst(logn, t1, INV_Q);
#endif
#else
		/* We now set the target [t0,t1] to [hm,0], then apply the
		   lattice basis to obtain the real target vector (after
		   normalization with regard to the modulus q).
		   b11 is unchanged; b01 is the copy in t1. */
		fpoly_apply_basis(logn, t0, t1, t1, b11, hm);
#endif

		/* Current layout:
		      t0  (n)
		      t1  (n)
		      g01 (n)
		      g00 (hn)
		      g11 (hn)
		   We now do the Fast Fourier sampling, which uses
		   up to 3*n slots beyond t1 (hence 7*n total usage
		   in tmp[]). */
		ffsamp_fft(&ss, tmp);

		/*
		 * At this point, [t0,t1] are the FFT representation of
		 * the sampled vector; in normal (non-FFT) representation,
		 * [t0,t1] is integral. We want to apply the lattice basis
		 * Compute the lattice basis B = [[g, -f], [G, -F]] to that
		 * vector, and subtract the result from [hm,0] to get the
		 * signature value. These computations can be done either
		 * in the FFT domain, or in with integers; and integer
		 * computations can be done modulo q = 12289 since the
		 * signature verification also works modulo q. Using
		 * integers is faster than staying in FFT representation
		 * when floating-point operations are emulated.
		 */
#if FNDSA_SIGN_SPLIT

		/* [t0, t1] hold [z0, z1] in FFT form; convert them to int16. */
		uint16_t *z0 = (uint16_t *)(t0);
		uint16_t *z1 = (uint16_t *)(t1);
		uint16_t *u0 = z0 + n;
		uint16_t *u1 = u0 + n;
		uint16_t *u2 = u1 + n;

		/* fxp32.64: the division by n/2 is folded into the rounding. */
#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264 \
	&& !(FNDSA_SSE2 || FNDSA_NEON || FNDSA_RV64D)
		if (logn == 9) {
			fpoly_iFFT_no_div2e(logn, t0);
			fpoly_iFFT_no_div2e(logn, t1);
			for (size_t i = 0; i < n; i ++) {
				z0[i] = fpr_rint_div2e8_to_u16(t0[i]);
				z1[i] = fpr_rint_div2e8_to_u16(t1[i]);
			}
		} else if (logn == 10) {
			fpoly_iFFT_no_div2e(logn, t0);
			fpoly_iFFT_no_div2e(logn, t1);
			for (size_t i = 0; i < n; i ++) {
				z0[i] = fpr_rint_div2e9_to_u16(t0[i]);
				z1[i] = fpr_rint_div2e9_to_u16(t1[i]);
			}
		} else
#endif
		{
			fpoly_iFFT(logn, t0);
			fpoly_iFFT(logn, t1);
			for (size_t i = 0; i < n; i ++) {
				z0[i] = (uint16_t)fpr_rint(t0[i]);
				z1[i] = (uint16_t)fpr_rint(t1[i]);
			}
		}

		/* Build signature: s2 = f*(z0 + r0) + F*(z1 + r1)
		   with:
		      c*f = (c*f mod q) + q*r1
		      -c*F = (-c*F mod q) + q*r0
		   This simplifies to:
		      s2 = f*(z0 + (c*F mod q)/q) + F*(z1 - (c*f mod q)/q)
		   Since s2 is integral, we can do these computations
		   modulo q'. (c-fn-dsa-alt) */

		/* u2  <- c mod q  (NTT mod q) */
		memcpy(u2, hm, n * sizeof(uint16_t));
		mqpoly_ext_to_int(logn, u2);
		mqpoly_int_to_ntt(logn, u2);

		/* nf  <- f mod q  (NTT mod q)
		   (f[] and g[] decoded) */
		f = (int8_t *)(z1 + n);
		g = f + n;
		(void)trim_i8_decode(logn, sign_key_fgF, f, nbits);
		(void)trim_i8_decode(logn, sign_key_fgF + flen, g, nbits);
		uint16_t *nf = (uint16_t *)(g + n);
		uint16_t *qf2 = nf + n;
		uint16_t *qF2 = qf2 + n;
		uint16_t *th = qF2 + n;
		uint16_t *ts = th + n;
		mqpoly_small_to_int(logn, f, nf);
		mqpoly_int_to_ntt(logn, nf);

		/* u0  <- F mod q  (NTT mod q) */
		mqpoly_small_to_int(logn, F, u0);
		mqpoly_int_to_ntt(logn, u0);

		/* u1  <- c*f mod q  (normalized) */
		memcpy(u1, u2, n * sizeof(uint16_t));
		mqpoly_mul_ntt(logn, u1, nf);
		mqpoly_ntt_to_int(logn, u1);
		(void)mqpoly_sqnorm_int_to_signed(logn, u1);

		/* u0  <- c*F mod q  (normalized) */
		mqpoly_mul_ntt(logn, u0, u2);
		mqpoly_ntt_to_int(logn, u0);
		(void)mqpoly_sqnorm_int_to_signed(logn, u0);

		/* u0  <- z0 + (c*F mod q)/q  (NTT mod q') */
		mq2poly_signed_to_mod(logn, u0);
		mq2poly_mulconst(logn, u0, MQ2_IQ);
		mq2poly_signed_to_mod(logn, z0);
		mq2poly_add(logn, u0, z0);
		mq2poly_NTT(logn, u0);

		/* z1  <- z1 - (c*f mod q)/q  (NTT mod q') */
		mq2poly_signed_to_mod(logn, u1);
		mq2poly_mulconst(logn, u1, MQ2_IQ);
		mq2poly_signed_to_mod(logn, z1);
		mq2poly_sub(logn, z1, u1);
		mq2poly_NTT(logn, z1);

		/* qf2 <- f mod q'  (NTT mod q') */
		mq2poly_small_to_mod(logn, f, qf2);
		mq2poly_NTT(logn, qf2);

		/* qF2 <- F mod q'  (NTT mod q') */
		mq2poly_small_to_mod(logn, F, qF2);
		mq2poly_NTT(logn, qF2);

		/* u0  <- s2 = f*(z0 + (c*F mod q)/q)
		               + F*(z1 - (c*f mod q)/q)  (normalized) */
		mq2poly_mul_ntt(logn, u0, qf2);
		mq2poly_mul_ntt(logn, z1, qF2);
		mq2poly_add(logn, u0, z1);
		mq2poly_iNTT(logn, u0);
		mq2poly_mod_to_signed(logn, u0);

		/* At this point:
		   u2    c mod q  (NTT mod q)
		   nf    f mod q  (NTT mod q)
		   u0    s2  (normalized) */

		/* th  <- g mod q  (NTT mod q) */
		mqpoly_small_to_int(logn, g, th);
		mqpoly_int_to_ntt(logn, th);

		/* th  <- h = g/f mod q  (NTT mod q) */
		mqpoly_div_ntt(logn, th, nf);

		/* ts  <- s2 mod q  (NTT mod q) */
		memcpy(ts, u0, n * sizeof(uint16_t));
		mqpoly_signed_to_int(logn, ts);
		mqpoly_int_to_ntt(logn, ts);

		/* u2  <- s1 = c - s2*h mod q  (internal representation) */
		mqpoly_mul_ntt(logn, ts, th);
		mqpoly_sub(logn, u2, ts);
		mqpoly_ntt_to_int(logn, u2);

		uint32_t sqn1 = mqpoly_sqnorm_int(logn, u2);
		uint32_t sqn2 = mqpoly_sqnorm_signed(logn, u0);

		/* If either squared norm saturated, or the sum (i.e.
		   the total squared norm of [s1,s2]) is too high, then
		   we loop. */
		uint32_t sqn = sqn1 + sqn2;
		sqn1 |= sqn2;
		sqn |= (uint32_t)(*(int32_t *)&sqn1 >> 31);
		if (!mqpoly_sqnorm_is_acceptable(logn, sqn)) {
			continue;
		}
		int16_t *s2 = (int16_t *)u0;

#elif !(FNDSA_SSE2 || FNDSA_NEON || FNDSA_RV64D)

		/* Convert [t0, t1] to integers modulo q. */
		uint16_t *ut0 = (uint16_t *)(t1 + n);
		uint16_t *ut1 = ut0 + n;
		uint16_t *ut2 = ut1 + n;
		uint16_t *ut3 = ut2 + n;
#if FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264
		if (logn == 9) {
			fpoly_iFFT_no_div2e(logn, t0);
			fpoly_iFFT_no_div2e(logn, t1);
			for (size_t i = 0; i < n; i ++) {
				ut0[i] = fpr_rint_div2e8_to_u16(t0[i]);
				ut1[i] = fpr_rint_div2e8_to_u16(t1[i]);
			}
		} else if (logn == 10) {
			fpoly_iFFT_no_div2e(logn, t0);
			fpoly_iFFT_no_div2e(logn, t1);
			for (size_t i = 0; i < n; i ++) {
				ut0[i] = fpr_rint_div2e9_to_u16(t0[i]);
				ut1[i] = fpr_rint_div2e9_to_u16(t1[i]);
			}
		} else
#endif
		{
			fpoly_iFFT(logn, t0);
			fpoly_iFFT(logn, t1);
#if FNDSA_SSE2
			/* We inline an fpr_rint() implementation, using SSE2
			   intrinsics (_mm_cvtpd_epi32() for rounding to neareast
			   with roundTiesToEven). */
			for (size_t i = 0; i < n; i += 2) {
				__m128d xt = _mm_loadu_pd((const double *)t0 + i);
				__m128i zt = _mm_cvtpd_epi32(xt);
				ut0[i + 0] = (uint16_t)_mm_cvtsi128_si32(zt);
				ut0[i + 1] = (uint16_t)_mm_cvtsi128_si32(
					_mm_bsrli_si128(zt, 4));
			}
			for (size_t i = 0; i < n; i += 2) {
				__m128d xt = _mm_loadu_pd((const double *)t1 + i);
				__m128i zt = _mm_cvtpd_epi32(xt);
				ut1[i + 0] = (uint16_t)_mm_cvtsi128_si32(zt);
				ut1[i + 1] = (uint16_t)_mm_cvtsi128_si32(
					_mm_bsrli_si128(zt, 4));
			}
#elif FNDSA_NEON
			/* We inline an fpr_rint() implementation, using NEON
			   intrinsics (vcvtnq_s64_f64() for rounding to neareast
			   with roundTiesToEven). */
			for (size_t i = 0; i < n; i += 2) {
				float64x2_t xt = vld1q_f64((const float64_t *)t0 + i);
				int64x2_t zt = vcvtnq_s64_f64(xt);
				ut0[i + 0] = (uint16_t)vgetq_lane_s64(zt, 0);
				ut0[i + 1] = (uint16_t)vgetq_lane_s64(zt, 1);
			}
			for (size_t i = 0; i < n; i += 2) {
				float64x2_t xt = vld1q_f64((const float64_t *)t1 + i);
				int64x2_t zt = vcvtnq_s64_f64(xt);
				ut1[i + 0] = (uint16_t)vgetq_lane_s64(zt, 0);
				ut1[i + 1] = (uint16_t)vgetq_lane_s64(zt, 1);
			}
#elif FNDSA_RV64D
			const f64 *tt0 = (const f64 *)t0;
			const f64 *tt1 = (const f64 *)t1;
			for (size_t i = 0; i < n; i ++) {
				ut0[i] = (uint64_t)f64_rint(tt0[i]);
			}
			for (size_t i = 0; i < n; i ++) {
				ut1[i] = (uint64_t)f64_rint(tt1[i]);
			}
#else
			for (size_t i = 0; i < n; i ++) {
				ut0[i] = (uint16_t)fpr_rint(t0[i]);
				ut1[i] = (uint16_t)fpr_rint(t1[i]);
			}
#endif
		}
		mqpoly_signed_to_int(logn, ut0);
		mqpoly_signed_to_int(logn, ut1);

		/* Convert [t0,t1] to NTT. */
		mqpoly_int_to_ntt(logn, ut0);
		mqpoly_int_to_ntt(logn, ut1);

		/* Decode (f,g) from the signing key. */
		f = (int8_t *)(ut3 + n);
		g = f + n;
		(void)trim_i8_decode(logn, sign_key_fgF, f, nbits);
		(void)trim_i8_decode(logn, sign_key_fgF + flen, g, nbits);

		/* s1 = hm - (g*t0 + G*t1).
		   We compute s1 into ut3; we do not need to keep s1,
		   only its squared norm. */
		mqpoly_small_to_int(logn, g, ut2);
		mqpoly_small_to_int(logn, G, ut3);
		mqpoly_int_to_ntt(logn, ut2);
		mqpoly_int_to_ntt(logn, ut3);
		mqpoly_mul_ntt(logn, ut2, ut0);
		mqpoly_mul_ntt(logn, ut3, ut1);
		mqpoly_add(logn, ut2, ut3);
		mqpoly_ntt_to_int(logn, ut2);
		memcpy(ut3, hm, n * sizeof(uint16_t));
		mqpoly_ext_to_int(logn, ut3);
		mqpoly_sub(logn, ut3, ut2);
		uint32_t sqn1 = mqpoly_sqnorm_int(logn, ut3);

		/* s2 = -(-f*t0 - F*t1) = f*t0 + F*t1
		   We compute s2 into ut3. */
		mqpoly_small_to_int(logn, f, ut2);
		mqpoly_small_to_int(logn, F, ut3);
		mqpoly_int_to_ntt(logn, ut2);
		mqpoly_int_to_ntt(logn, ut3);
		mqpoly_mul_ntt(logn, ut2, ut0);
		mqpoly_mul_ntt(logn, ut3, ut1);
		mqpoly_add(logn, ut3, ut2);
		mqpoly_ntt_to_int(logn, ut3);
		uint32_t sqn2 = mqpoly_sqnorm_int_to_signed(logn, ut3);

		/* If either squared norm saturated, or the sum (i.e.
		   the total squared norm of [s1,s2]) is too high, then
		   we loop. */
		uint32_t sqn = sqn1 + sqn2;
		sqn1 |= sqn2;
		sqn |= (uint32_t)(*(int32_t *)&sqn1 >> 31);
		if (!mqpoly_sqnorm_is_acceptable(logn, sqn)) {
			continue;
		}
		int16_t *s2 = (int16_t *)ut3;

#else
		/*
		 * We stay here in the floating-point domain for the
		 * application of the basis. This happens to be faster
		 * on our test platforms with a hardware FPU.
		 */

		/* Get the lattice point corresponding to the sampled
		   vector. This means computing:
		      [v0,v1] = [t0,t1] * [[g, -f], [G, -F]]
		   hence:
		      v0 = t0*g + t1*G
		      v1 = -t0*f - t1*F  */
		fpr *w0 = t1 + n;
		fpr *w1 = w0 + n;
		f = (int8_t *)(w1 + n);
		g = f + n;
		(void)trim_i8_decode(logn, sign_key_fgF, f, nbits);
		(void)trim_i8_decode(logn, sign_key_fgF + flen, g, nbits);
		fpoly_set_small(logn, w0, g);
		fpoly_set_small(logn, w1, f);
		fpoly_FFT(logn, w0);
		fpoly_FFT(logn, w1);
		fpoly_mul_fft(logn, w1, t0);
		fpoly_mul_fft(logn, t0, w0);
		fpoly_set_small(logn, w0, G);
		fpoly_FFT(logn, w0);
		fpoly_mul_fft(logn, w0, t1);
		fpoly_add(logn, t0, w0);
		fpoly_set_small(logn, w0, F);
		fpoly_FFT(logn, w0);
		fpoly_mul_fft(logn, t1, w0);
		fpoly_add(logn, t1, w1);
		fpoly_neg(logn, t1);
		fpoly_iFFT(logn, t0);
		fpoly_iFFT(logn, t1);

		/* We compute s1, then s2 into buffer s2 (s1 is not
		   retained). We accumulate their squared norm in sqn,
		   with an "overflow" flag in ng. */
		uint32_t sqn = 0;
		uint32_t ng = 0;
		int16_t *s2 = (int16_t *)w0;
#if FNDSA_SSE2
		/* We inline an fpr_rint() implementation, using SSE2
		   intrinsics (_mm_cvtpd_epi32() for rounding to neareast
		   with roundTiesToEven). */
		for (size_t i = 0; i < n; i += 2) {
			__m128d xt = _mm_loadu_pd((const double *)t0 + i);
			__m128i zt = _mm_cvtpd_epi32(xt);
			uint16_t zu0 = hm[i + 0]
				- (uint16_t)_mm_cvtsi128_si32(zt);
			uint16_t zu1 = hm[i + 1]
				- (uint16_t)_mm_cvtsi128_si32(
					_mm_bsrli_si128(zt, 4));
			int32_t z0 = (int32_t)*(int16_t *)&zu0;
			int32_t z1 = (int32_t)*(int16_t *)&zu1;
			sqn += (uint32_t)(z0 * z0);
			ng |= sqn;
			sqn += (uint32_t)(z1 * z1);
			ng |= sqn;
		}
		for (size_t i = 0; i < n; i += 2) {
			__m128d xt = _mm_loadu_pd((const double *)t1 + i);
			__m128i zt = _mm_cvtpd_epi32(xt);
			uint16_t zu0 = -(uint16_t)_mm_cvtsi128_si32(zt);
			uint16_t zu1 = -(uint16_t)_mm_cvtsi128_si32(
				_mm_bsrli_si128(zt, 4));
			int32_t z0 = (int32_t)*(int16_t *)&zu0;
			int32_t z1 = (int32_t)*(int16_t *)&zu1;
			sqn += (uint32_t)(z0 * z0);
			ng |= sqn;
			sqn += (uint32_t)(z1 * z1);
			ng |= sqn;
			s2[i + 0] = (int16_t)z0;
			s2[i + 1] = (int16_t)z1;
		}
#elif FNDSA_NEON
		/* We inline an fpr_rint() implementation, using NEON
		   intrinsics (vcvtnq_s64_f64() for rounding to neareast
		   with roundTiesToEven). */
		for (size_t i = 0; i < n; i += 2) {
			float64x2_t xt = vld1q_f64((const float64_t *)t0 + i);
			int64x2_t zt = vcvtnq_s64_f64(xt);
			uint16_t zu0 = hm[i + 0]
				- (uint16_t)vgetq_lane_s64(zt, 0);
			uint16_t zu1 = hm[i + 1]
				- (uint16_t)vgetq_lane_s64(zt, 1);
			int32_t z0 = (int32_t)*(int16_t *)&zu0;
			int32_t z1 = (int32_t)*(int16_t *)&zu1;
			sqn += (uint32_t)(z0 * z0);
			ng |= sqn;
			sqn += (uint32_t)(z1 * z1);
			ng |= sqn;
		}
		for (size_t i = 0; i < n; i += 2) {
			float64x2_t xt = vld1q_f64((const float64_t *)t1 + i);
			int64x2_t zt = vcvtnq_s64_f64(xt);
			uint16_t zu0 = -(uint16_t)vgetq_lane_s64(zt, 0);
			uint16_t zu1 = -(uint16_t)vgetq_lane_s64(zt, 1);
			int32_t z0 = (int32_t)*(int16_t *)&zu0;
			int32_t z1 = (int32_t)*(int16_t *)&zu1;
			sqn += (uint32_t)(z0 * z0);
			ng |= sqn;
			sqn += (uint32_t)(z1 * z1);
			ng |= sqn;
			s2[i + 0] = (int16_t)z0;
			s2[i + 1] = (int16_t)z1;
		}
#elif FNDSA_RV64D
		const f64 *tt0 = (const f64 *)t0;
		const f64 *tt1 = (const f64 *)t1;
		for (size_t i = 0; i < n; i ++) {
			uint16_t zu = hm[i] - (uint16_t)f64_rint(tt0[i]);
			int32_t z = *(int16_t *)&zu;
			sqn += (uint32_t)(z * z);
			ng |= sqn;
		}
		for (size_t i = 0; i < n; i ++) {
			uint16_t zu = -(uint16_t)f64_rint(tt1[i]);
			int32_t z = *(int16_t *)&zu;
			sqn += (uint32_t)(z * z);
			ng |= sqn;
			s2[i] = (int16_t)z;
		}
#else
		for (size_t i = 0; i < n; i ++) {
			uint16_t zu = hm[i] - (uint16_t)fpr_rint(t0[i]);
			int32_t z = *(int16_t *)&zu;
			sqn += (uint32_t)(z * z);
			ng |= sqn;
		}
		for (size_t i = 0; i < n; i ++) {
			uint16_t zu = -(uint16_t)fpr_rint(t1[i]);
			int32_t z = *(int16_t *)&zu;
			sqn += (uint32_t)(z * z);
			ng |= sqn;
			s2[i] = (int16_t)z;
		}
#endif

		/* If the squared norm exceeds 2^31-1, then at some point
		   the high bit of ng was set, which we use to saturate
		   the squared norm to 2^32-1. If the squared norm is
		   unacceptable, then we loop. */
		sqn |= (uint32_t)(*(int32_t *)&ng >> 31);
		if (!mqpoly_sqnorm_is_acceptable(logn, sqn)) {
			continue;
		}
#endif

		/* We have a candidate signature; we must encode it. This
		   may fail, if the signature cannot be encoded in the
		   target size. */
		size_t sig_len = FNDSA_SIGNATURE_SIZE(logn);
		/* s2 is public: it is the signature, or it is discarded. */
		FNDSA_CT_PUBLIC(s2, n * sizeof *s2);
		if (comp_encode(logn, s2, sig + 41, sig_len - 41)) {
			/* Success! */
			sig[0] = 0x30 + logn;
			memcpy(sig + 1, nonce, 40);
			ret = sig_len;
			goto sign_exit;
		}
	}

sign_exit:
#if FNDSA_SSE2
	_MM_SET_ROUNDING_MODE(round_mode);
#endif
	return ret;
}
