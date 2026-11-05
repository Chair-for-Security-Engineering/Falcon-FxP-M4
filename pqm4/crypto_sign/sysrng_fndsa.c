/*
 * Replaces sysrng.c inside pqm4: randomness comes from the framework's
 * randombytes() (deterministic in the testvectors binaries).
 */

#include "inner.h"
#include "randombytes.h"

/* see inner.h */
int
sysrng(void *dst, size_t len)
{
	if (len == 0) {
		return 1;
	}
	return randombytes(dst, len) == 0;
}
