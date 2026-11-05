#ifndef FNDSA_FPRTEST_H__
#define FNDSA_FPRTEST_H__

/*
 * On-board test that the Cortex-M4 assembly fpr operations are
 * bit-identical to the portable C of the same representation. One binary
 * links both: the contract ops resolve to the assembly, and
 * fprtest_cref_<rep>.c supplies the C versions under other names (a NULL
 * entry means there is no separate C version; reported as SKIP).
 * Build and run with 'make pqm4-fprtest[-run]'.
 */

#include "../sign_inner.h"

typedef struct {
	fpr (*op_scaled)(int64_t i, int sc);
	fpr (*op_of32)(int32_t i);
	fpr (*op_add)(fpr x, fpr y);
	fpr (*op_sub)(fpr x, fpr y);
	void (*op_add_sub)(fpr *a, fpr *b, fpr x, fpr y);
	fpr (*op_mul)(fpr x, fpr y);
	fpr (*op_sqr)(fpr x);
	fpr (*op_div)(fpr x, fpr y);
	fpr (*op_inv)(fpr x);
	fpr (*op_sqrt)(fpr x);
	void (*op_sqrt_inv_sigma)(fpr *d, fpr x, unsigned logn);
	int64_t (*op_rint)(fpr x);
	int64_t (*op_floor)(fpr x);
	int64_t (*op_trunc)(fpr x);
	/* Only for a fused FPC_MUL (fxp32.64); NULL makes the driver
	   compose it from mul, add and sub. */
	void (*op_fpc_mul)(fpr *d_re, fpr *d_im, fpr a_re, fpr a_im,
		fpr b_re, fpr b_im);
	void (*op_fxp96_cmul_gm)(uint32_t *out, const uint32_t *a,
		const uint64_t *gm);
	void (*op_fxp96_cmul_gmj)(uint32_t *out, const uint32_t *a,
		const uint64_t *gm);
	fpr (*op_fpr_mul_gm)(fpr x, uint64_t c);
	fpr (*op_mul_recip)(fpr a, fpr x);
	fpr (*op_mul_recip_p63)(fpr a, fpr x);
	fpr (*op_mul_inv_q)(fpr x);
} fprtest_ops;

/* Supplied by the representation's fprtest_cref_<rep>.c. */
extern const fprtest_ops fndsa_fprtest_cref;

/* Runs the comparison, one line per operation through 'emit'. Returns
   the number of operations with a mismatch. */
int fndsa_fprtest_run(void (*emit)(const char *), unsigned long nvec);

#endif
