/* Rounding macros for the fxp32.64 assembly: ties to even, or ties up
 * (round bit alone) with FNDSA_FXP_TIES_UP=1. */
#ifndef __ASSEMBLER__
#error "fxp_round_cm4.h is assembly-only"
#endif
#include "knobs.h"
@ Rounding increment into C for a product cut at 2^64: round bit =
@ bit 31 of \p1, sticky = \p0 | \p1[30:0], quotient low limb \q0,
@ scratch \t (may equal \p0).
.macro fxp_round_carry t, p0, p1, q0
#if FNDSA_FXP_TIES_UP
    lsls    \t, \p1, #1          @ C = round bit
#else
    orr     \t, \p0, \p1, lsl #1 @ sticky bits: p0 | p1[30:0]
    orr     \t, \t, \q0, lsl #31 @ | (q0 & 1) in a spare bit position
    and     \t, \t, \p1, asr #31 @ sticky|odd word, kept iff round bit
    cmp.w   \t, #1               @ C = RNE increment
#endif
.endm

@ fpr_mul has no register left for p0, so ties to even parks (p0 != 0)
@ in the carry flag until the rounding.
.macro fxp_mul_park_sticky p0
#if !FNDSA_FXP_TIES_UP
    cmp     \p0, #1              @ C = (p0 != 0)
#endif
.endm

@ ...and its increment, with C = (p0 != 0) still parked (ties to even).
.macro fxp_mul_round_carry t, p1, q0
#if FNDSA_FXP_TIES_UP
    lsls    \t, \p1, #1          @ C = round bit p1[31]
#else
    adc     \t, \p1, \p1         @ p1<<1 | (p0!=0): nonzero iff sticky
    orr     \t, \t, \q0, lsl #31 @ fold in q0 bit 0 (ties-to-even)
    and     \t, \t, \p1, asr #31 @ keep sticky|odd only if round bit
    cmp     \t, #1               @ C = increment
#endif
.endm
