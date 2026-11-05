/*
 * main() of the on-board asm-vs-C fpr test (see fpr/fprtest.h), a
 * local addition to pqm4. Built by
 * 'make pqm4-fprtest[-run]'; only meaningful for the *-m4f and *-cm4
 * implementations.
 */

#include "api.h"
#include "hal.h"

extern int fndsa_fprtest_run(void (*emit)(const char *),
	unsigned long nvec);

#define FPRTEST_NVEC   2000

int
main(void)
{
	int fail_ops;

	hal_setup(CLOCK_FAST);

	hal_send_str("==========================");
	hal_send_str("fpr asm-vs-C differential test");
	hal_send_str(CRYPTO_ALGNAME);

	fail_ops = fndsa_fprtest_run(hal_send_str, FPRTEST_NVEC);

	hal_send_str(fail_ops
		? "FPRTEST RESULT: FAIL (see mismatches above)"
		: "FPRTEST RESULT: all operations bit-identical");
	hal_send_str("#");
	hal_send_str("==========================");
	while (1) {
		/* done */
	}
	return 0;
}
