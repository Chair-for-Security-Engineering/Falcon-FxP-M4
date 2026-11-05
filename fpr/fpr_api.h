#ifndef FNDSA_FPR_API_H__
#define FNDSA_FPR_API_H__

/* ==================================================================== */
/*
 * fpr backend contract.
 *
 * The signing code uses real numbers only through the 'fpr' type and the
 * operations of the representation header selected below. A backend is a
 * representation (FNDSA_FPR_REP: type, constants, inline ops) plus an
 * implementation (FNDSA_FPR_IMPL: the out-of-line ops, in C or assembly).
 * Selection is at compile time.
 *
 * Every representation header must provide:
 *   - the 'fpr' type, a plain value type (only the native 'prec'
 *     analysis backend uses heap memory)
 *   - FPR_STRIDE: sizeof(fpr) as used in polynomial buffer layouts
 *   - constants: FPR_ZERO, FPR_NZERO, FPR_ONE, an FPR(...) initializer
 *     whose arguments differ per representation, and the named constants
 *     and tables listed in ieee/fpr_ieee.h (INV_Q, Q2, GM, INV_SIGMA, ...)
 *   - inline ops: fpr_ursh, fpr_irsh, fpr_ulsh, fpr_of, fpr_of32,
 *     fpr_rint, fpr_floor, fpr_trunc, fpr_sub, fpr_neg, fpr_half,
 *     fpr_double, fpr_sqr, fpr_inv, fpr_div2e, fpr_mul2e
 *   - out-of-line ops: fpr_scaled, fpr_add, fpr_mul, fpr_div, fpr_sqrt
 *   - macros: FPR_ADD_SUB(a, b, x, y), FPC_MUL(d_re, d_im, ...)
 * A representation may also override the FPOLY_GM_* macros (defaults in
 * sign_inner.h) to change how the twiddle table is stored.
 */

#include "../inner.h"

/* ==================================================================== */
/* Backend checks (selection and defaults are in knobs.h).              */

/* The SIMD and RISC-V hardware-FP paths treat fpr values as native
   doubles, so they work only with the ieee representation. */
#if FNDSA_FPR_REP != FNDSA_FPR_REP_IEEE
#if FNDSA_AVX2 || FNDSA_SSE2 || FNDSA_NEON || FNDSA_RV64D
#error non-ieee fpr representations are incompatible with FNDSA_AVX2/FNDSA_SSE2/FNDSA_NEON/FNDSA_RV64D
#endif
#endif

/* The ieee cm4 implementation is keyed off FNDSA_ASM_CORTEXM4. */
#if FNDSA_FPR_REP == FNDSA_FPR_REP_IEEE \
	&& FNDSA_FPR_IMPL == FNDSA_FPR_IMPL_CM4 && !FNDSA_ASM_CORTEXM4
#error the (ieee, cm4) backend requires FNDSA_ASM_CORTEXM4=1
#endif

/* ==================================================================== */
/* Representation header.                                                */

#if FNDSA_FPR_REP == FNDSA_FPR_REP_IEEE
#include "ieee/fpr_ieee.h"
#elif FNDSA_FPR_REP == FNDSA_FPR_REP_TW
#include "tw/fpr_tw.h"
#elif FNDSA_FPR_REP == FNDSA_FPR_REP_FXP3264
#include "fxp32.64/fpr_fxp.h"
#include "fxp32.64/fxp_fused.h"
#include "fxp32.64/fxp_fft.h"
#elif FNDSA_FPR_REP == FNDSA_FPR_REP_PREC
#include "prec/fpr_prec.h"
#else
#error unknown or not yet integrated FNDSA_FPR_REP value
#endif

/* ==================================================================== */

#endif
