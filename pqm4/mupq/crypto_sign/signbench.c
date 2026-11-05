/*
 * Sign-only cycle benchmark, a local addition to pqm4.
 * pqm4's speed binary makes a new key for every signature, and keygen
 * dominates its run time. This one makes SIGNBENCH_KEYS keys and signs
 * MUPQ_ITERATIONS times with each. Output has the same format as speed.c.
 * Built by 'make pqm4-signbench'.
 */

#include "api.h"
#include "hal.h"
#include "sendfn.h"
#include "randombytes.h"

#include <stdio.h>
#include <stdint.h>
#include <string.h>

#define MLEN 59

/* Keys to generate; MUPQ_ITERATIONS signatures are measured under each. */
#ifndef SIGNBENCH_KEYS
#define SIGNBENCH_KEYS 2
#endif

// https://stackoverflow.com/a/1489985/1711232
#define PASTER(x, y) x##y
#define EVALUATOR(x, y) PASTER(x, y)
#define NAMESPACE(fun) EVALUATOR(MUPQ_NAMESPACE, fun)

// use different names so we can have empty namespaces
#define MUPQ_CRYPTO_PUBLICKEYBYTES NAMESPACE(CRYPTO_PUBLICKEYBYTES)
#define MUPQ_CRYPTO_SECRETKEYBYTES NAMESPACE(CRYPTO_SECRETKEYBYTES)
#define MUPQ_CRYPTO_BYTES          NAMESPACE(CRYPTO_BYTES)
#define MUPQ_CRYPTO_ALGNAME        NAMESPACE(CRYPTO_ALGNAME)

#define MUPQ_crypto_sign_keypair NAMESPACE(crypto_sign_keypair)
#define MUPQ_crypto_sign NAMESPACE(crypto_sign)
#define MUPQ_crypto_sign_open NAMESPACE(crypto_sign_open)

#define printcycles(S, U) send_unsignedll((S), (U))

int main(void)
{
  unsigned char sk[MUPQ_CRYPTO_SECRETKEYBYTES];
  unsigned char pk[MUPQ_CRYPTO_PUBLICKEYBYTES];
  unsigned char sm[MLEN+MUPQ_CRYPTO_BYTES];
  size_t smlen;
  unsigned int rc;
  unsigned long long t0, t1;
  int k, i;

  hal_setup(CLOCK_BENCHMARK);

  hal_send_str("==========================");

  for(k=0;k<SIGNBENCH_KEYS;k++)
  {
    t0 = hal_get_time();
    MUPQ_crypto_sign_keypair(pk, sk);
    t1 = hal_get_time();
    printcycles("keypair cycles:", t1-t0);

    for(i=0;i<MUPQ_ITERATIONS;i++)
    {
      randombytes(sm, MLEN);
      t0 = hal_get_time();
      MUPQ_crypto_sign(sm, &smlen, sm, MLEN, sk);
      t1 = hal_get_time();
      printcycles("sign cycles:", t1-t0);

      /* Verify every time, so a broken build cannot pass unnoticed. */
      t0 = hal_get_time();
      rc = MUPQ_crypto_sign_open(sm, &smlen, sm, smlen, pk);
      t1 = hal_get_time();
      printcycles("verify cycles:", t1-t0);

      if (rc) {
        hal_send_str("ERROR Signature did not verify correctly!\n");
      }
      hal_send_str("+");
    }
  }
  hal_send_str("#");
  return 0;
}
