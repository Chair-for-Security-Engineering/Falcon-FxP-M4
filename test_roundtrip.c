/*
 * Keygen, sign and verify round trip through the public API, for every
 * fpr representation. The verifier is integer-only, so a valid signature
 * must verify whatever the representation; a tampered one must not.
 */

#include <stdio.h>
#include <string.h>
#include "fndsa.h"

static int
test_degree(unsigned logn, int iters)
{
	uint8_t sk[FNDSA_SIGN_KEY_SIZE(10)];
	uint8_t vk[FNDSA_VRFY_KEY_SIZE(10)];
	uint8_t sig[FNDSA_SIGNATURE_SIZE(10)];
	size_t sklen = FNDSA_SIGN_KEY_SIZE(logn);
	size_t vklen = FNDSA_VRFY_KEY_SIZE(logn);

	for (int it = 0; it < iters; it++) {
		if (fndsa_keygen(logn, sk, vk) != 1) {
			printf("  logn=%u: keygen FAILED\n", logn);
			return 0;
		}
		char msg[64];
		snprintf(msg, sizeof msg, "roundtrip message %d for logn %u",
			it, logn);
		size_t mlen = strlen(msg);

		/* Raw-message signing (FNDSA_HASH_ID_RAW), empty context. */
		size_t siglen = fndsa_sign(sk, sklen, NULL, 0,
			FNDSA_HASH_ID_RAW, msg, mlen, sig, sizeof sig);
		if (siglen == 0) {
			printf("  logn=%u it=%d: sign FAILED\n", logn, it);
			return 0;
		}
		if (!fndsa_verify(sig, siglen, vk, vklen, NULL, 0,
			FNDSA_HASH_ID_RAW, msg, mlen)) {
			printf("  logn=%u it=%d: VERIFY FAILED (len %zu)\n",
				logn, it, siglen);
			return 0;
		}

		/* Negative control: a flipped byte must not verify. */
		sig[siglen / 2] ^= 0x01;
		if (fndsa_verify(sig, siglen, vk, vklen, NULL, 0,
			FNDSA_HASH_ID_RAW, msg, mlen)) {
			printf("  logn=%u it=%d: tampered sig verified (BAD)\n",
				logn, it);
			return 0;
		}
	}
	printf("  logn=%u (%s): %d sign/verify roundtrips OK, "
		"tamper rejected\n", logn, logn == 9 ? "512" : "1024", iters);
	return 1;
}

int
main(void)
{
	printf("keygen/sign/verify roundtrip (standard degrees)\n");
	int ok = 1;
	ok &= test_degree(9, 20);
	ok &= test_degree(10, 10);
	printf("%s\n", ok ? "ALL OK" : "FAILURE");
	return ok ? 0 : 1;
}
