/* NIST/PQClean signature API for FN-DSA-1024 (logn = 10). Messages are
   signed raw (no pre-hash, empty context); sm = signature || message. */

#include <string.h>

#include "api.h"

#define LOGN   10

int
crypto_sign_keypair(uint8_t *pk, uint8_t *sk)
{
	return fndsa_keygen(LOGN, sk, pk) ? 0 : -1;
}

int
crypto_sign_signature(uint8_t *sig, size_t *siglen,
	const uint8_t *m, size_t mlen, const uint8_t *sk)
{
	size_t r = fndsa_sign(sk, CRYPTO_SECRETKEYBYTES, NULL, 0,
		FNDSA_HASH_ID_RAW, m, mlen, sig, CRYPTO_BYTES);
	if (r != CRYPTO_BYTES) {
		return -1;
	}
	*siglen = r;
	return 0;
}

int
crypto_sign_verify(const uint8_t *sig, size_t siglen,
	const uint8_t *m, size_t mlen, const uint8_t *pk)
{
	return fndsa_verify(sig, siglen, pk, CRYPTO_PUBLICKEYBYTES,
		NULL, 0, FNDSA_HASH_ID_RAW, m, mlen) ? 0 : -1;
}

int
crypto_sign(uint8_t *sm, size_t *smlen,
	const uint8_t *m, size_t mlen, const uint8_t *sk)
{
	size_t siglen;

	memmove(sm + CRYPTO_BYTES, m, mlen);
	if (crypto_sign_signature(sm, &siglen, sm + CRYPTO_BYTES,
		mlen, sk) != 0)
	{
		return -1;
	}
	*smlen = siglen + mlen;
	return 0;
}

int
crypto_sign_open(uint8_t *m, size_t *mlen,
	const uint8_t *sm, size_t smlen, const uint8_t *pk)
{
	if (smlen < CRYPTO_BYTES || crypto_sign_verify(sm, CRYPTO_BYTES,
		sm + CRYPTO_BYTES, smlen - CRYPTO_BYTES, pk) != 0)
	{
		*mlen = (size_t)-1;
		return -1;
	}
	*mlen = smlen - CRYPTO_BYTES;
	memmove(m, sm + CRYPTO_BYTES, *mlen);
	return 0;
}
