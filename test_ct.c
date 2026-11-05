/*
 * Constant-time check of signing, TIMECOP-style: built with FNDSA_CTGRIND=1
 * and run under Valgrind memcheck. The signing key and seeds are marked
 * secret, so memcheck reports any branch or memory index that depends on
 * them. Key generation is not checked.
 *
 *   ./test_ct [keys] [sigs] [logn...]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "inner.h"
#include "kgen_inner.h"
#include "sign_inner.h"

static void *
xmalloc(size_t len)
{
	void *p = malloc(len);
	if (p == NULL) {
		fprintf(stderr, "out of memory\n");
		exit(EXIT_FAILURE);
	}
	return p;
}

/* Deterministic test bytes: SHAKE256(label || logn || i || j). */
static void
test_bytes(void *dst, size_t len, const char *label,
	unsigned logn, unsigned i, unsigned j)
{
	shake_context sc;
	uint8_t buf[12];
	uint32_t v[3] = { logn, i, j };

	for (int k = 0; k < 12; k ++) {
		buf[k] = (uint8_t)(v[k >> 2] >> ((k & 3) << 3));
	}
	shake_init(&sc, 256);
	shake_inject(&sc, label, strlen(label));
	shake_inject(&sc, buf, sizeof buf);
	shake_flip(&sc);
	shake_extract(&sc, dst, len);
}

static int
test_degree(unsigned logn, unsigned keys, unsigned sigs)
{
	size_t n = (size_t)1 << logn;
	size_t sk_len = FNDSA_SIGN_KEY_SIZE(logn);
	size_t vk_len = FNDSA_VRFY_KEY_SIZE(logn);
	size_t ktmp_len = ((size_t)KGEN_TMP_MUL << logn) + 31;
	size_t stmp_len = ((size_t)FNDSA_SIGN_BPC << logn) + 31;
	uint8_t *sk = xmalloc(sk_len);
	uint8_t *vk = xmalloc(vk_len);
	uint8_t *sig = xmalloc(FNDSA_SIGNATURE_SIZE(logn));
	void *ktmp = xmalloc(ktmp_len);
	void *stmp = xmalloc(stmp_len);
	int ok = 1;

	(void)n;
	for (unsigned i = 0; i < keys && ok; i ++) {
		uint8_t seed[48];

		test_bytes(seed, sizeof seed, "ct-keygen", logn, i, 0);
		if (!fndsa_keygen_seeded_temp(logn, seed, sizeof seed,
			sk, vk, ktmp, ktmp_len))
		{
			printf("logn=%u key %u: keygen failed\n", logn, i);
			ok = 0;
			break;
		}
		/* The signing key is secret apart from its header byte
		   (the degree). */
		FNDSA_CT_SECRET(sk + 1, sk_len - 1);

		for (unsigned j = 0; j < sigs; j ++) {
			uint8_t msg[32], sseed[48];
			size_t sig_len;

			test_bytes(msg, sizeof msg, "ct-msg", logn, i, j);
			test_bytes(sseed, sizeof sseed, "ct-sign", logn, i, j);
			FNDSA_CT_SECRET(sseed, sizeof sseed);
			sig_len = fndsa_sign_seeded_temp(sk, sk_len, NULL, 0,
				FNDSA_HASH_ID_RAW, msg, sizeof msg,
				sseed, sizeof sseed,
				sig, FNDSA_SIGNATURE_SIZE(logn), stmp, stmp_len);
			if (sig_len == 0 || !fndsa_verify(sig, sig_len,
				vk, vk_len, NULL, 0, FNDSA_HASH_ID_RAW,
				msg, sizeof msg))
			{
				printf("logn=%u key %u sig %u: FAILED\n",
					logn, i, j);
				ok = 0;
				break;
			}
		}
	}
	if (ok) {
		printf("logn=%u: %u keys x %u signatures OK\n",
			logn, keys, sigs);
	}
	free(sk);
	free(vk);
	free(sig);
	free(ktmp);
	free(stmp);
	return ok;
}

int
main(int argc, char *argv[])
{
	unsigned keys = argc > 1 ? (unsigned)atoi(argv[1]) : 4;
	unsigned sigs = argc > 2 ? (unsigned)atoi(argv[2]) : 2;
	int ok = 1;

	if (argc > 3) {
		for (int k = 3; k < argc; k ++) {
			ok &= test_degree((unsigned)atoi(argv[k]), keys, sigs);
		}
	} else {
		ok &= test_degree(9, keys, sigs);
		ok &= test_degree(10, keys, sigs);
	}
	return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
