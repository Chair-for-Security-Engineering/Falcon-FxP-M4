#ifndef PQM4_FNDSA_API_H__
#define PQM4_FNDSA_API_H__

/* NIST/PQClean signature API for FN-DSA-512 (logn = 9), on top of
   fndsa.h. */

#include <stddef.h>
#include <stdint.h>

#include "fndsa.h"

#define CRYPTO_SECRETKEYBYTES   FNDSA_SIGN_KEY_SIZE(9)
#define CRYPTO_PUBLICKEYBYTES   FNDSA_VRFY_KEY_SIZE(9)
#define CRYPTO_BYTES            FNDSA_SIGNATURE_SIZE(9)
#define CRYPTO_ALGNAME          "FN-DSA-512"

int crypto_sign_keypair(uint8_t *pk, uint8_t *sk);

int crypto_sign_signature(uint8_t *sig, size_t *siglen,
	const uint8_t *m, size_t mlen, const uint8_t *sk);

int crypto_sign_verify(const uint8_t *sig, size_t siglen,
	const uint8_t *m, size_t mlen, const uint8_t *pk);

int crypto_sign(uint8_t *sm, size_t *smlen,
	const uint8_t *m, size_t mlen, const uint8_t *sk);

int crypto_sign_open(uint8_t *m, size_t *mlen,
	const uint8_t *sm, size_t smlen, const uint8_t *pk);

#endif
