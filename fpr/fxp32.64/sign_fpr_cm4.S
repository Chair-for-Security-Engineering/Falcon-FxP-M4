	.syntax	unified
	.cpu	cortex-m4
	.file	"sign_fpr_cm4.S"
	.text

/* fxp32.64 arithmetic for the Cortex-M4: multiply, complex and twiddle
   products, 1/q, reciprocal and square root kernels. A value is a signed
   96-bit integer in three limbs, x = x0 + x1*2^32 + x2*2^64, scaled by
   2^-64. All code is constant time. Rounding is to nearest, ties to even,
   or ties up with FNDSA_FXP_TIES_UP (see fxp_round_cm4.h). Most routines
   use a custom register convention, given above each one. */
#include "fpr/fxp32.64/fxp_round_cm4.h"
#include "fpr/fxp32.64/fxp_const_gen.h"

@ =======================================================================
@ fpr fndsa_fpr_mul(fpr x, fpr y)
@
@ Input: x = r0:r1:r2, y = r3:r4:r5. Output: z = r0:r1:r2.
@ Clobbers: r6-r9, r12, flags.
@ =======================================================================

    .align 2
    .global fndsa_fpr_mul
    .thumb_func
    .type fndsa_fpr_mul, %function
fndsa_fpr_mul:
    @ r6 = p0 / scratch, r7 = p1 / scratch, r8:r9:r12 = q0:q1:q2.
    @ Unsigned product first; p0 and p1 only feed the rounding.
    umull   r6, r7, r0, r3      @ r6 = p0, r7 = high(x0*y0)
    umull   r8, r9, r0, r4      @ r8 = low(x0*y1), r9 = high(x0*y1)
    umaal   r7, r8, r1, r3      @ r7 = p1, r8 = carry into column 2

    @ umull/umaal leave the flags alone, so C survives until rounding.
    fxp_mul_park_sticky r6

    @ column 2: carry + x0*y2 + x1*y1 + x2*y0
    umaal   r8, r9, r0, r5
    umull   r6, r12, r1, r4
    umaal   r8, r6, r2, r3

    @ column 3: carry + x1*y2 + x2*y1
    umaal   r9, r12, r1, r5
    umaal   r9, r6, r2, r4

    @ column 4: carry + x2*y2 (p5 is dropped)
    umaal   r12, r6, r2, r5

    @ Round q2:q1:q0 = p4:p3:p2. The signed correction only touches
    @ q1:q2, so it can come after.
    fxp_mul_round_carry r6, r7, r8
    adcs    r8, r8, #0
    adcs    r9, r9, #0
    adcs    r12, r12, #0

    @ Signed correction: subtract (x < 0 ? y : 0) + (y < 0 ? x : 0) at
    @ the 2^96 column; only y0:y1 and x0:x1 reach q1:q2.
    and.w   r6, r3, r2, asr #31 @ y0 if x < 0
    subs    r9, r9, r6
    and.w   r6, r4, r2, asr #31 @ y1 if x < 0
    sbcs    r12, r12, r6

    and.w   r6, r0, r5, asr #31 @ x0 if y < 0
    and.w   r7, r1, r5, asr #31 @ x1 if y < 0
    subs    r1, r9, r6
    sbcs    r2, r12, r7
    mov     r0, r8
    bx      lr

    .size fndsa_fpr_mul, .-fndsa_fpr_mul
@ =======================================================================
@ fndsa_fxp96_cmul_direct_body areg, aim_s, breg, bim_s, outim_s
@
@ Shared body of the two complex multiplies below:
@   imag = round((ar*bi + ai*br) / 2^64), stored to [outim_s]
@   real = round((ar*br - ai*bi) / 2^64), left in r11:r3:r12
@ areg/breg are core registers holding the ar/br pointers, aim_s/bim_s/
@ outim_s VFP registers holding the ai/bi/imag-output pointers.
@ Clobbers: r0, r3-r12, r14, flags, memory; areg and breg are kept.
@ =======================================================================
.macro fndsa_fxp96_cmul_direct_body areg, aim_s, breg, bim_s, outim_s
    @ imag, product 1: x = ar in r4:r5:r6, y = bi in r7:r8:r9
    vmov    r3, \bim_s
    ldr     r4, [\areg, #0]
    ldr     r5, [\areg, #4]
    ldr     r6, [\areg, #8]
    ldr     r7, [r3, #0]
    ldr     r8, [r3, #4]
    ldr     r9, [r3, #8]

    @ unsigned product p0..p4 (p5 dropped)
    umull   r0, r10, r4, r7      @ r0 = p0
    umull   r11, r3, r4, r8
    umaal   r10, r11, r5, r7      @ r10 = p1
    @ column 2: carry + x0*y2 + x1*y1 + x2*y0
    umaal   r11, r3, r4, r9
    umull   lr, r12, r5, r8
    umaal   r11, lr, r6, r7      @ r11 = p2 = q0
    @ column 3: carry + x1*y2 + x2*y1
    umaal   r3, r12, r5, r9
    umaal   r3, lr, r6, r8      @ r3 = p3 = q1
    @ column 4: carry + x2*y2
    umaal   r12, lr, r6, r9     @ r12 = p4 = q2

    @ signed correction, as in fndsa_fpr_mul
    and     r7, r7, r6, asr #31
    subs.w  r3, r3, r7
    and     r7, r8, r6, asr #31
    sbcs    r12, r12, r7
    and     r7, r4, r9, asr #31
    subs.w  r3, r3, r7
    and     r7, r5, r9, asr #31
    sbcs    r12, r12, r7

    @ imag, product 2, accumulated: x = ai in r4:r5:r6, y = br in r7:r8:r9
    vmov    lr, \aim_s
    ldr     r4, [lr, #0]
    ldr     r5, [lr, #4]
    ldr     r6, [lr, #8]
    ldr     r7, [\breg, #0]
    ldr     r8, [\breg, #4]
    ldr     r9, [\breg, #8]

    @ x0y0, x0y1, x0y2, x1y2, x2y2
    mov lr, #0
    umaal   r0, lr, r4, r7
    umaal   r10, lr, r4, r8
    umaal   r11, lr, r4, r9
    umaal   r3, lr, r5, r9
    umaal   r12, lr, r6, r9

    @ x1y0, x1y1, x2y1
    mov     lr, #0
    umaal   r10, lr, r5, r7
    umaal   r11, lr, r5, r8
    umaal   r3, lr, r6, r8
    add.w   r12, r12, lr

    @ x2y0
    mov     lr, #0
    umaal   r11, lr, r6, r7        @ q0 final; lr = col-3 carry, folded below

    @ p0, p1, q0 are final, so round now; the lr carry goes in the same
    @ chain.
    fxp_round_carry r0, r0, r10, r11
    adcs    r11, r11, #0
    adcs    r3, r3, lr
    adc     r12, r12, #0

    @ signed correction
    and     r7, r7, r6, asr #31
    subs.w  r3, r3, r7
    and     r7, r8, r6, asr #31
    sbcs    r12, r12, r7
    and     r7, r4, r9, asr #31
    subs.w  r3, r3, r7
    and     r7, r5, r9, asr #31
    sbcs    r12, r12, r7

    @ Store imag and negate ai for the real component.
    vmov    r7, \outim_s
    mov     lr, #0
    str     r11, [r7, #0]
    rsbs    r4, r4, #0
    str     r3, [r7, #4]
    sbcs    r5, lr, r5
    str     r12, [r7, #8]
    sbcs    r6, lr, r6

    @ real, product 2: x = -ai in r4:r5:r6, y = bi in r7:r8:r9
    vmov    r3, \bim_s
    ldr     r7, [r3, #0]
    ldr     r8, [r3, #4]
    ldr     r9, [r3, #8]

    umull   r0, r10, r4, r7
    umull   r11, r3, r4, r8
    umaal   r10, r11, r5, r7
    umaal   r11, r3, r4, r9
    umull   lr, r12, r5, r8
    umaal   r11, lr, r6, r7
    umaal   r3, r12, r5, r9
    umaal   r3, lr, r6, r8
    umaal   r12, lr, r6, r9

    @ signed correction
    and     r7, r7, r6, asr #31
    subs.w  r3, r3, r7
    and     r7, r8, r6, asr #31
    sbcs    r12, r12, r7
    and     r7, r4, r9, asr #31
    subs.w  r3, r3, r7
    and     r7, r5, r9, asr #31
    sbcs    r12, r12, r7

    @ real, product 1, accumulated: x = ar in r4:r5:r6, y = br in r7:r8:r9
    ldr     r4, [\areg, #0]
    ldr     r5, [\areg, #4]
    ldr     r6, [\areg, #8]
    ldr     r7, [\breg, #0]
    ldr     r8, [\breg, #4]
    ldr     r9, [\breg, #8]

    @ x0y0, x0y1, x0y2, x1y2, x2y2
    mov     lr, #0
    umaal   r0, lr, r4, r7
    umaal   r10, lr, r4, r8
    umaal   r11, lr, r4, r9
    umaal   r3, lr, r5, r9
    umaal   r12, lr, r6, r9

    @ x1y0, x1y1, x2y1
    mov     lr, #0
    umaal   r10, lr, r5, r7
    umaal   r11, lr, r5, r8
    umaal   r3, lr, r6, r8
    add.w   r12, r12, lr

    @ x2y0
    mov     lr, #0
    umaal   r11, lr, r6, r7         @ q0 final; lr = col-3 carry, folded below

    @ round
    fxp_round_carry r0, r0, r10, r11
    adcs    r11, r11, #0
    adcs    r3, r3, lr
    adc     r12, r12, #0

    @ signed correction
    and     r7, r7, r6, asr #31
    subs.w  r3, r3, r7
    and     r7, r8, r6, asr #31
    sbcs    r12, r12, r7
    and     r7, r4, r9, asr #31
    subs.w  r3, r3, r7
    and     r7, r5, r9, asr #31
    sbcs    r12, r12, r7
.endm

#if FNDSA_FXP_CMUL_FUSED
@ =======================================================================
@ void fndsa_fxp96_cmul_self_rne64_x(uint32_t a_re[3],
@                                   uint32_t a_im[3],
@                                   const uint32_t b_re[3],
@                                   const uint32_t b_im[3])
@
@ a <- a * b. Input: r0 = a_re, r1 = a_im, r2 = b_re, r3 = b_im; b must
@ not overlap a.
@ Clobbers: r0-r12, r14, s9, s11, s12, flags, memory.
@ =======================================================================

.align 2
.global fndsa_fxp96_cmul_self_rne64_x
.thumb_func
.type fndsa_fxp96_cmul_self_rne64_x, %function
fndsa_fxp96_cmul_self_rne64_x:
    vmov    s9, r1              @ out_im
    vmov    s11, r3             @ b_im
    vmov    s12, lr
    mov     r1, r0              @ out_re/a_re; r2 remains b_re

    fndsa_fxp96_cmul_direct_body r1, s9, r2, s11, s9

    @ store real
    str     r11, [r1, #0]
    str     r3, [r1, #4]
    str     r12, [r1, #8]

    vmov    lr, s12
    bx      lr

.size fndsa_fxp96_cmul_self_rne64_x, .-fndsa_fxp96_cmul_self_rne64_x
@ =======================================================================
@ void fndsa_fxp96_cmul_noalias_rne64_x(uint32_t d_re[3],
@                                      uint32_t d_im[3],
@                                      const uint32_t a_re[3],
@                                      const uint32_t a_im[3],
@                                      const uint32_t b_re[3],
@                                      const uint32_t b_im[3])
@
@ Input: r0 = d_re, r1 = d_im, r2 = a_re, r3 = a_im, r4 = b_re, r5 = b_im.
@ d must not overlap a or b.
@ Clobbers: r0-r12, r14, s8, s9, s10, s12, s13, flags, memory.
@ =======================================================================

.align 2
.global fndsa_fxp96_cmul_noalias_rne64_x
.thumb_func
.type fndsa_fxp96_cmul_noalias_rne64_x, %function
fndsa_fxp96_cmul_noalias_rne64_x:
    vmov    s8, s9, r0, r1      @ out_re, out_im
    vmov    s10, r3             @ a_im
    vmov    s12, r5             @ b_im
    vmov    s13, lr
    mov     r1, r4              @ b_re; r2 remains a_re

    fndsa_fxp96_cmul_direct_body r2, s10, r1, s12, s9

    @ store real
    vmov    r0, s8
    str     r11, [r0, #0]
    str     r3, [r0, #4]
    str     r12, [r0, #8]

    vmov    lr, s13
    bx      lr

.size fndsa_fxp96_cmul_noalias_rne64_x, .-fndsa_fxp96_cmul_noalias_rne64_x
#endif /* FNDSA_FXP_CMUL_FUSED */
@ =======================================================================
@ void fndsa_fxp96_cmul_gm_x(uint32_t out[6],
@                              const uint32_t a[6],
@                              const uint64_t gm[2])
@
@ out = a * (A + iB) for an unsigned Q0.64 GM pair gm = (A, B). As the
@ complex multiply above, but with 3x2-limb products and only a's sign
@ correction. Input: r0 = out, r1 = a, r2 = gm; out may alias a.
@ Clobbers: r0-r12, r14, s8, s9, flags, memory.
@ =======================================================================

@ Body shared with cmul_gmj. gm computes imag then real, gmj real then
@ imag; the offsets pick the operands and output slots.
.macro fndsa_fxp96_cmul_gm_body first_aoff, second_aoff, first_store, second_store
    vmov    s8, s9, r0, lr

    @ first component (gm: ar*B + ai*A, gmj: ai*B + ar*A)
    @ product 1: x = ar (gm) or ai (gmj) in r4:r5:r6, y = B in r7:r8
    .if \first_aoff == 0
    ldr     r4, [r1, #0]
    ldr     r5, [r1, #4]
    ldr     r6, [r1, #8]
    .else
    ldr     r4, [r1, #12]
    ldr     r5, [r1, #16]
    ldr.w   r6, [r1, #20]
    .endif
    ldr     r7, [r2, #8]   @ B0, resident
    ldr     r8, [r2, #12]
    ldr     lr, [r2, #0]   @ A0, resident
    ldr     r11, [r2, #4]

    umull   r0, r10, r4, r7      @ r0 = p0
    umull   r2, r3, r4, r8
    umaal   r10, r2, r5, r7      @ r10 = p1
    @ column 2: carry + x1*y1 + x2*y0
    umaal   r2, r3, r5, r8
    mov     r9, #0
    umaal   r2, r9, r6, r7       @ r2 = p2 = q0
    @ column 3: carry + x2*y1
    umaal   r3, r9, r6, r8       @ r3 = p3 = q1

    @ signed correction (GM words are unsigned)
    and     r12, r7, r6, asr #31
    subs    r3, r3, r12
    and     r12, r8, r6, asr #31
    sbcs    r9, r9, r12

    @ product 2, accumulated: x = ai (gm) or ar (gmj), y = A in lr:r11
    .if \second_aoff == 0
    ldr     r4, [r1, #0]
    ldr     r5, [r1, #4]
    ldr     r6, [r1, #8]
    .else
    ldr     r4, [r1, #12]
    ldr     r5, [r1, #16]
    ldr.w   r6, [r1, #20]
    .endif

    @ x0y0, x0y1, x1y1, x2y1
    mov     r12, #0
    umaal   r0, r12, r4, lr
    umaal   r10, r12, r4, r11
    umaal   r2, r12, r5, r11
    umaal   r3, r12, r6, r11
    add.w   r9, r9, r12

    @ x1y0, x2y0
    mov     r12, #0
    umaal   r10, r12, r5, lr
    umaal   r2, r12, r6, lr

    @ round, with the r12 carry in the same chain
    fxp_round_carry r0, r0, r10, r2
    adcs    r2, r2, #0
    adcs    r3, r3, r12
    adc     r9, r9, #0

    @ signed correction
    and     r12, lr, r6, asr #31
    subs    r3, r3, r12
    and     r12, r11, r6, asr #31
    sbcs    r9, r9, r12

    @ store, then negate x for the second component
    vmov    r12, s8
    mov     r10, #0
    str     r2, [r12, #\first_store]
    rsbs.w  r4, r4, #0
    str     r3, [r12, #(\first_store + 4)]
    sbcs    r5, r10, r5
    str     r9, [r12, #(\first_store + 8)]
    sbcs    r6, r10, r6

    @ second component (gm: ar*A - ai*B, gmj: ai*A - ar*B)
    @ product 1: -x times B
    umull   r0, r10, r4, r7
    umull   r2, r3, r4, r8
    umaal   r10, r2, r5, r7
    umaal   r2, r3, r5, r8
    mov     r9, #0
    umaal   r2, r9, r6, r7
    umaal   r3, r9, r6, r8

    and     r7, r7, r6, asr #31
    subs.w  r3, r3, r7
    and     r7, r8, r6, asr #31
    sbcs    r9, r9, r7

    @ product 2, accumulated: x = ar (gm) or ai (gmj), y = A in lr:r11
    .if \first_aoff == 0
    ldr     r4, [r1, #0]
    ldr     r5, [r1, #4]
    ldr     r6, [r1, #8]
    .else
    ldr     r4, [r1, #12]
    ldr     r5, [r1, #16]
    ldr.w   r6, [r1, #20]
    .endif

    @ x0y0, x0y1, x1y1, x2y1
    mov     r12, #0
    umaal   r0, r12, r4, lr
    umaal   r10, r12, r4, r11
    umaal   r2, r12, r5, r11
    umaal   r3, r12, r6, r11
    add.w   r9, r9, r12

    @ x1y0, x2y0
    mov     r12, #0
    umaal   r10, r12, r5, lr
    umaal   r2, r12, r6, lr

    @ round
    fxp_round_carry r0, r0, r10, r2
    adcs    r2, r2, #0
    adcs    r3, r3, r12
    adc     r9, r9, #0

    @ signed correction
    and     r12, lr, r6, asr #31
    subs    r3, r3, r12
    and     r12, r11, r6, asr #31
    sbcs    r9, r9, r12

    @ store
    vmov    r12, s8
    str     r2, [r12, #\second_store]
    str     r3, [r12, #(\second_store + 4)]
    str     r9, [r12, #(\second_store + 8)]
    vmov    lr, s9
    bx      lr
.endm

@ Own section, so that it is dropped when only the asm butterflies are
@ used.
.section .text.fndsa_fxp96_cmul_gm_x,"ax",%progbits
.align 2
.global fndsa_fxp96_cmul_gm_x
.thumb_func
.type fndsa_fxp96_cmul_gm_x, %function
fndsa_fxp96_cmul_gm_x:
    fndsa_fxp96_cmul_gm_body 0, 12, 12, 0

.size fndsa_fxp96_cmul_gm_x, .-fndsa_fxp96_cmul_gm_x
@ =======================================================================
@ void fndsa_fxp96_cmul_gmj_x(uint32_t out[6],
@                               const uint32_t a[6],
@                               const uint64_t gm[2])
@
@ out = a * (A - iB). Same convention as fndsa_fxp96_cmul_gm_x.
@ =======================================================================
.section .text.fndsa_fxp96_cmul_gmj_x,"ax",%progbits
.align 2
.global fndsa_fxp96_cmul_gmj_x
.thumb_func
.type fndsa_fxp96_cmul_gmj_x, %function
fndsa_fxp96_cmul_gmj_x:
    fndsa_fxp96_cmul_gm_body 12, 0, 0, 12

.size fndsa_fxp96_cmul_gmj_x, .-fndsa_fxp96_cmul_gmj_x
.text

@ =======================================================================
@ fpr fndsa_fpr_mul_gm_x(fpr x, uint64_t c)
@
@ x * c for one unsigned Q0.64 GM word c, rounded once.
@ Input: x = r0:r1:r2, c = r3:r4. Output: z = r0:r1:r2.
@ Clobbers: r3-r7, flags.
@ =======================================================================

.align 2
.global fndsa_fpr_mul_gm_x
.thumb_func
.type fndsa_fpr_mul_gm_x, %function
fndsa_fpr_mul_gm_x:
    @ p0 = r5, p1 = r6, q0:q1:q2 = r7:r0:r1 (x registers are reused)
    umull   r5, r6, r0, r3      @ x0y0: r5 = p0, r6 = col1
    umull   r7, r0, r0, r4      @ x0y1: r7 col1, r0 = col2 (x0 slot reused)
    umaal   r6, r7, r1, r3      @ x1y0: r6 = p1, r7 = col2
    umaal   r7, r0, r1, r4      @ x1y1: r7 col2, r0 = col3 (x1 slot free)
    mov     r1, #0
    umaal   r7, r1, r2, r3      @ x2y0: r7 = q0 final, r1 = col3
    umaal   r0, r1, r2, r4      @ x2y1: r0 = q1 final, r1 = q2

    @ round, then correct q1:q2
    fxp_round_carry r5, r5, r6, r7
    adcs    r7, r7, #0
    adcs    r0, r0, #0
    adcs    r1, r1, #0

    @ signed correction: subtract c at 2^96 if x < 0
    and     r5, r3, r2, asr #31
    and     r6, r4, r2, asr #31
    subs.w  r0, r0, r5
    sbcs    r2, r1, r6          @ q2 straight into its output slot
    mov     r1, r0
    mov     r0, r7
    bx      lr

.size fndsa_fpr_mul_gm_x, .-fndsa_fpr_mul_gm_x

@ =======================================================================
@ fpr fndsa_fpr_mul_inv_q_x(fpr x)
@
@ x/q as x * INV_Q_Q77 / 2^77, rounded once. The product is that of
@ fndsa_fpr_mul_gm_x, rounded at bit 77 instead of 64.
@ Input: x = r0:r1:r2. Output: z = r0:r1:r2. Clobbers: r3-r7, flags.
@ =======================================================================

.section .text.fndsa_fpr_mul_inv_q_x,"ax",%progbits
.align 2
.global fndsa_fpr_mul_inv_q_x
.thumb_func
.type fndsa_fpr_mul_inv_q_x, %function
fndsa_fpr_mul_inv_q_x:
    movw    r3, #(INV_Q_Q77_LO & 0xFFFF)
    movt    r3, #(INV_Q_Q77_LO >> 16)
    movw    r4, #(INV_Q_Q77_HI & 0xFFFF)
    movt    r4, #(INV_Q_Q77_HI >> 16)
    umull   r5, r6, r0, r3      @ x0c0: r5 = p0, r6 = col1
    umull   r7, r0, r0, r4      @ x0c1: r7 col1, r0 = col2
    umaal   r6, r7, r1, r3      @ x1c0: r6 = p1, r7 = col2
    umaal   r7, r0, r1, r4      @ x1c1: r7 col2, r0 = col3
    mov     r1, #0
    umaal   r7, r1, r2, r3      @ x2c0: r7 = p2, r1 = col3
    umaal   r0, r1, r2, r4      @ x2c1: r0 = p3, r1 = p4

    @ Signed correction at the 2^96 column (c is unsigned).
    and     r3, r3, r2, asr #31
    and     r4, r4, r2, asr #31
    subs    r0, r0, r3
    sbc     r1, r1, r4

    @ Move the round bit (bit 76) to bit 31 of r6, sticky bits 75..64
    @ below it; p0|p1 are the remaining sticky words. Then z = p >> 77.
    orr     r5, r5, r6
    lsl     r6, r7, #19
    asr     r2, r1, #13
    lsl     r1, r1, #19
    orr     r1, r1, r0, lsr #13
    lsl     r0, r0, #19
    orr     r0, r0, r7, lsr #13
    fxp_round_carry r5, r5, r6, r0
    adcs    r0, r0, #0
    adcs    r1, r1, #0
    adc     r2, r2, #0
    bx      lr

.size fndsa_fpr_mul_inv_q_x, .-fndsa_fpr_mul_inv_q_x
.text
#if FNDSA_FFT_NEW
@ =======================================================================
@ fpr fndsa_fpr_mul_rnj_x(fpr x, uint64_t frac, uint32_t integer)
@
@ x * (integer + frac/2^64), rounded once after integer*x is added.
@ Input: x = r0:r1:r2, frac = r3:r4, integer = r5. Output: z = r0:r1:r2.
@ Clobbers: r3-r7, r12, s8-s15, flags.
@ =======================================================================

.align 2
.global fndsa_fpr_mul_rnj_x
.thumb_func
.type fndsa_fpr_mul_rnj_x, %function
fndsa_fpr_mul_rnj_x:
    vmov    s8, s9, r0, r1
    vmov    s10, r2
    vmov    s11, r5
    vmov    s14, s15, r3, r4

    @ x * frac: r5:r6 = p0:p1, r7:r0:r1 = q0:q1:q2
    umull   r5, r6, r0, r3
    umull   r7, r0, r0, r4
    umaal   r6, r7, r1, r3
    umaal   r7, r0, r1, r4
    mov     r1, #0
    umaal   r7, r1, r2, r3
    umaal   r0, r1, r2, r4

    @ low 96 bits of integer*x
    vmov    r2, s11
    vmov    r3, s8
    umull   r3, r4, r3, r2          @ a0:a1 = x0*integer
    vmov    s12, r3
    vmov    s13, r4
    vmov    r3, s9
    umull   r3, r4, r3, r2          @ r3:r4 = x1*integer
    vmov    r12, s13
    adds    r12, r12, r3            @ a1, carry into a2
    adc     r4, r4, #0              @ a2 = high(x1*integer)+carry
    vmov    r3, s10
    umull   r3, r2, r3, r2          @ low(x2*integer) is enough
    adds    r4, r4, r3              @ a2; overflow is outside Q32.64

    @ add integer*x before rounding, so that q0 has the right parity
    vmov    r3, s12
    adds    r7, r7, r3
    adcs    r0, r0, r12
    adcs    r1, r1, r4

    fxp_round_carry r5, r5, r6, r7
    adcs    r7, r7, #0
    adcs    r0, r0, #0
    adcs    r1, r1, #0

    @ signed correction of x*frac, on q1:q2 only
    vmov    r3, r4, s14, s15
    vmov    r2, s10
    and     r5, r3, r2, asr #31
    subs.w  r0, r0, r5
    and     r5, r4, r2, asr #31
    sbcs    r1, r1, r5

    mov     r2, r1
    mov     r1, r0
    mov     r0, r7
    bx      lr

.size fndsa_fpr_mul_rnj_x, .-fndsa_fpr_mul_rnj_x
#endif /* FNDSA_FFT_NEW */

#if FNDSA_FFT_NEW
@ =======================================================================
@ void fndsa_fpr_bfly_gm_batch_x(fpr *a, fpr *b, size_t num,
@                                  uint64_t frac, int neg)
@
@ num compact-FFT butterflies with one Q0.64 coefficient c = frac: with
@ x = a[i], z = b[i]*c, (a[i], b[i]) = (x+z, x-z), or swapped if neg != 0.
@ neg is public. Input: r0 = a, r1 = b, r2:r3 = frac, r4 = num, r5 = neg.
@ Clobbers: r0-r12, r14, s3-s5, flags, memory.
@ =======================================================================

@ lo:mid:hi <- round(lo:mid:hi / 2), as fpr_half.

.macro fndsa_fpr_half96 lo, mid, hi, tmp
    and     \tmp, \lo, \lo, lsr #1
    and     \tmp, \tmp, #1
    asrs    \hi, \hi, #1
    rrxs    \mid, \mid
    rrx     \lo, \lo
    adds    \lo, \lo, \tmp
    adcs    \mid, \mid, #0
    adc     \hi, \hi, #0
.endm

@ Butterfly loop: coefficient in s3:s4, pointers in r11/r12, count in lr.
@ neg and half are assembly-time parameters.
.macro fndsa_fpr_bfly_gm_batch_loop label, neg, half
\label:
    ldrd    r8, r9, [r11, #0]
    ldr     r10, [r11, #8]
    ldrd    r0, r1, [r12, #0]
    ldr     r2, [r12, #8]
    vmov    r3, r4, s3, s4

    umull   r5, r6, r0, r3
    umull   r7, r0, r0, r4
    umaal   r6, r7, r1, r3
    umaal   r7, r0, r1, r4
    mov     r1, #0
    umaal   r7, r1, r2, r3
    umaal   r0, r1, r2, r4
    fxp_round_carry r5, r5, r6, r7
    adcs    r7, r7, #0
    adcs    r0, r0, #0
    adcs    r1, r1, #0
    and     r5, r3, r2, asr #31
    subs.w  r0, r0, r5
    and     r5, r4, r2, asr #31
    sbcs    r1, r1, r5

    @ r3:r4:r5 = x-z and r6:r7:r0 = x+z.
    subs    r3, r8, r7
    sbcs    r4, r9, r0
    sbcs    r5, r10, r1
    adds    r6, r8, r7
    adcs    r7, r9, r0
    adcs    r0, r10, r1
.if \half
    @ final merge: halve the outputs
    fndsa_fpr_half96 r3, r4, r5, r1
    fndsa_fpr_half96 r6, r7, r0, r1
.endif
.if \neg
    strd    r3, r4, [r11, #0]
    str     r5, [r11, #8]
    strd    r6, r7, [r12, #0]
    str     r0, [r12, #8]
.else
    strd    r6, r7, [r11, #0]
    str     r0, [r11, #8]
    strd    r3, r4, [r12, #0]
    str     r5, [r12, #8]
.endif
    add     r11, r11, #12
    add     r12, r12, #12
    subs    lr, lr, #1
    bne     \label
.endm

.align 2
.global fndsa_fpr_bfly_gm_batch_x
.thumb_func
.type fndsa_fpr_bfly_gm_batch_x, %function
fndsa_fpr_bfly_gm_batch_x:
    vmov    s5, lr
    vmov    s3, s4, r2, r3
    mov     r11, r0
    mov     r12, r1
    mov     lr, r4
    cmp     r5, #0
    beq     .Lfndsa_fpr_bfly_gm_batch_pos

    fndsa_fpr_bfly_gm_batch_loop .Lfndsa_fpr_bfly_gm_batch_neg, 1, 0
    b       .Lfndsa_fpr_bfly_gm_batch_done

.Lfndsa_fpr_bfly_gm_batch_pos:
    fndsa_fpr_bfly_gm_batch_loop .Lfndsa_fpr_bfly_gm_batch_pos_loop, 0, 0

.Lfndsa_fpr_bfly_gm_batch_done:
    vmov    lr, s5
    bx      lr

.size fndsa_fpr_bfly_gm_batch_x, .-fndsa_fpr_bfly_gm_batch_x

@ void fndsa_fpr_bfly_gm_batch_half_x(fpr *a, fpr *b, size_t num,
@                                       uint64_t frac, int neg)
@
@ As fndsa_fpr_bfly_gm_batch_x, with halved outputs (final merge).
.align 2
.global fndsa_fpr_bfly_gm_batch_half_x
.thumb_func
.type fndsa_fpr_bfly_gm_batch_half_x, %function
fndsa_fpr_bfly_gm_batch_half_x:
    vmov    s5, lr
    vmov    s3, s4, r2, r3
    mov     r11, r0
    mov     r12, r1
    mov     lr, r4
    cmp     r5, #0
    beq     .Lfndsa_fpr_bfly_gm_batch_half_pos

    fndsa_fpr_bfly_gm_batch_loop .Lfndsa_fpr_bfly_gm_batch_half_neg, 1, 1
    b       .Lfndsa_fpr_bfly_gm_batch_half_done

.Lfndsa_fpr_bfly_gm_batch_half_pos:
    fndsa_fpr_bfly_gm_batch_loop .Lfndsa_fpr_bfly_gm_batch_half_pos_loop, 0, 1

.Lfndsa_fpr_bfly_gm_batch_half_done:
    vmov    lr, s5
    bx      lr

.size fndsa_fpr_bfly_gm_batch_half_x, .-fndsa_fpr_bfly_gm_batch_half_x
#endif /* FNDSA_FFT_NEW */

#if FNDSA_FFT_NEW
@ =======================================================================
@ void fndsa_fpr_bfly_rnj_batch_x(fpr *a, fpr *b, size_t num,
@                                   uint64_t frac, uint32_t integer,
@                                   int neg)
@
@ As fndsa_fpr_bfly_gm_batch_x, with coefficient c = integer + frac/2^64
@ (small public integer); z = b[i]*c is computed as in fndsa_fpr_mul_rnj_x.
@ Input: r0 = a, r1 = b, r2:r3 = frac, r4 = integer, r5 = num, r6 = neg.
@ Clobbers: r0-r12, r14, s3-s11, flags, memory.
@ =======================================================================

.macro fndsa_fpr_bfly_rnj_batch_loop label, neg, half
\label:
    @ x = a[i], y = b[i]
    ldrd    r8, r9, [r11, #0]
    ldr     r10, [r11, #8]
    ldrd    r0, r1, [r12, #0]
    ldr     r2, [r12, #8]
    vmov    r3, r4, s3, s4

    @ y * frac: p0:p1 = r5:r6, quotient r7:r0:r1
    umull   r5, r6, r0, r3
    umull   r7, r0, r0, r4
    umaal   r6, r7, r1, r3
    umaal   r7, r0, r1, r4
    mov     r1, #0
    umaal   r7, r1, r2, r3
    umaal   r0, r1, r2, r4

    @ add low96(integer*y) before rounding (it can change the parity)
    vmov    s7, s8, r5, r6
    ldrd    r2, r3, [r12, #0]
    ldr     r4, [r12, #8]
    vmov    s10, s11, r3, r4
    vmov    r5, s6
    umull   r6, r3, r2, r5          @ r6:r3 = y0*integer
    vmov    r2, s10
    umull   r2, r4, r2, r5          @ r2:r4 = y1*integer
    adds    r3, r3, r2
    adc     r4, r4, #0
    vmov    r2, s11
    mul     r2, r2, r5              @ low(y2*integer) only
    adds    r4, r4, r2
    adds    r7, r7, r6
    adcs    r0, r0, r3
    adcs    r1, r1, r4

    @ round once
    vmov    r5, r6, s7, s8
    fxp_round_carry r5, r5, r6, r7
    adcs    r7, r7, #0
    adcs    r0, r0, #0
    adcs    r1, r1, #0

    @ signed correction of y*frac (integer*y needs none)
    ldr     r2, [r12, #8]
    vmov    r3, r4, s3, s4
    and     r5, r3, r2, asr #31
    subs.w  r0, r0, r5
    and     r5, r4, r2, asr #31
    sbcs    r1, r1, r5

    @ x - z and x + z, stores swapped for neg
    subs    r3, r8, r7
    sbcs    r4, r9, r0
    sbcs    r5, r10, r1
    adds    r6, r8, r7
    adcs    r7, r9, r0
    adcs    r0, r10, r1
.if \half
    fndsa_fpr_half96 r3, r4, r5, r1
    fndsa_fpr_half96 r6, r7, r0, r1
.endif
.if \neg
    strd    r3, r4, [r11, #0]
    str     r5, [r11, #8]
    strd    r6, r7, [r12, #0]
    str     r0, [r12, #8]
.else
    strd    r6, r7, [r11, #0]
    str     r0, [r11, #8]
    strd    r3, r4, [r12, #0]
    str     r5, [r12, #8]
.endif
    add     r11, r11, #12
    add     r12, r12, #12
    subs    lr, lr, #1
    bne     \label
.endm

.align 2
.global fndsa_fpr_bfly_rnj_batch_x
.thumb_func
.type fndsa_fpr_bfly_rnj_batch_x, %function
fndsa_fpr_bfly_rnj_batch_x:
    vmov    s5, lr
    vmov    s3, s4, r2, r3
    vmov    s6, r4
    mov     r11, r0
    mov     r12, r1
    mov     lr, r5
    cmp     r6, #0
    beq     .Lfndsa_fpr_bfly_rnj_batch_pos

    fndsa_fpr_bfly_rnj_batch_loop .Lfndsa_fpr_bfly_rnj_batch_neg, 1, 0
    b       .Lfndsa_fpr_bfly_rnj_batch_done

.Lfndsa_fpr_bfly_rnj_batch_pos:
    fndsa_fpr_bfly_rnj_batch_loop .Lfndsa_fpr_bfly_rnj_batch_pos_loop, 0, 0

.Lfndsa_fpr_bfly_rnj_batch_done:
    vmov    lr, s5
    bx      lr

.size fndsa_fpr_bfly_rnj_batch_x, .-fndsa_fpr_bfly_rnj_batch_x

@ void fndsa_fpr_bfly_rnj_batch_half_x(fpr *a, fpr *b, size_t num,
@                                        uint64_t frac, uint32_t integer,
@                                        int neg)
@
@ As fndsa_fpr_bfly_rnj_batch_x, with halved outputs (final merge).
.align 2
.global fndsa_fpr_bfly_rnj_batch_half_x
.thumb_func
.type fndsa_fpr_bfly_rnj_batch_half_x, %function
fndsa_fpr_bfly_rnj_batch_half_x:
    vmov    s5, lr
    vmov    s3, s4, r2, r3
    vmov    s6, r4
    mov     r11, r0
    mov     r12, r1
    mov     lr, r5
    cmp     r6, #0
    beq     .Lfndsa_fpr_bfly_rnj_batch_half_pos

    fndsa_fpr_bfly_rnj_batch_loop .Lfndsa_fpr_bfly_rnj_batch_half_neg, 1, 1
    b       .Lfndsa_fpr_bfly_rnj_batch_half_done

.Lfndsa_fpr_bfly_rnj_batch_half_pos:
    fndsa_fpr_bfly_rnj_batch_loop .Lfndsa_fpr_bfly_rnj_batch_half_pos_loop, 0, 1

.Lfndsa_fpr_bfly_rnj_batch_half_done:
    vmov    lr, s5
    bx      lr

.size fndsa_fpr_bfly_rnj_batch_half_x, .-fndsa_fpr_bfly_rnj_batch_half_x
#endif /* FNDSA_FFT_NEW */

@ =======================================================================
@ void fndsa_fxp96_mul_recip_x(uint32_t z[3],
@                                  const uint32_t a[3],
@                                  const uint32_t rc[4],
@                                  uint32_t k)
@
@ z = round(a * rc * 2^(k-158)), as fxp_mul_recip in ops_c.c. The divisor
@ is a gamma-checked LDL diagonal d >= 4, so k <= 29 and b = 30-k is in
@ [1, 29]: the result always comes from P[4..7] of the 96x128 product.
@ Input: r0 = z, r1 = a, r2 = rc, r3 = k.
@ Clobbers: r0-r12, r14, s10, s14-s15, flags, memory.
@ =======================================================================

.align 2
.global fndsa_fxp96_mul_recip_x
.thumb_func
.type fndsa_fxp96_mul_recip_x, %function
fndsa_fxp96_mul_recip_x:
#if FNDSA_FXP_TIES_UP
	@ Ties up: no sticky needed, so P[0] and P[1] are dropped early.
	vmov    s15, lr
	rsb     r3, r3, #30

	ldrd    r4, r5, [r1]
	ldr     r6, [r1, #8]
	ldrd    r7, r8, [r2]
	ldrd    r9, r10, [r2, #8]
	umull   r11, r12, r4, r7    @ P[0] dumped in r11, carry to P[1]
	umull   r11, r1, r4, r8
	umaal   r12, r11, r5, r7    @ P[1] final in r12 (dead), r11 = carry
	umaal   r11, r1, r4, r9
	umull   lr, r2, r5, r8
	umaal   r11, lr, r6, r7     @ P[2]
	umaal   r1, r2, r5, r9
	umaal   r1, lr, r6, r8      @ P[3]
	umaal   r2, lr, r6, r9      @ P[4], P[5] in lr
	movs    r12, #0
	umaal   r1, r12, r4, r10
	umaal   r2, r12, r5, r10
	umaal   lr, r12, r6, r10    @ P[6]

	@ signed correction: subtract rc at P[3..6] if a < 0, borrow to P[7]
	and     r5, r7, r6, asr #31
	subs    r1, r1, r5
	and     r5, r8, r6, asr #31
	sbcs    r2, r2, r5
	and     r5, r9, r6, asr #31
	sbcs    lr, lr, r5
	and     r5, r10, r6, asr #31
	sbcs    r12, r12, r5
	sbc     r5, r5, r5          @ P[7] from the final borrow

	@ P >> (128+b); the first shift leaves the round bit in C
	rsb     r6, r3, #32
	lsrs    r8, r2, r3          @ C = round bit P[4][b-1]
	lsl     r7, lr, r6
	orr     r8, r8, r7          @ q0
	lsr     r9, lr, r3
	lsl     r7, r12, r6
	orr     r9, r9, r7          @ q1
	lsr     r10, r12, r3
	lsl     r7, r5, r6
	orr     r10, r10, r7        @ q2
	adcs    r8, r8, #0
	adcs    r9, r9, #0
	adc     r10, r10, #0

	str     r8, [r0]
	strd    r9, r10, [r0, #4]
	vmov    lr, s15
	bx      lr
#else
	@ s15 = return address, s10 = b = 30-k, s14 = z
	vmov    s15, lr
	rsb     r3, r3, #30
	vmov    s10, r3
	vmov    s14, r0             @ output pointer

	@ a = r4:r5:r6, rc = r7:r8:r9:r10, P[0..6] = r11:r12:r0:r1:r2:lr:r3
	ldrd    r4, r5, [r1]
	ldr     r6, [r1, #8]
	ldrd    r7, r8, [r2]
	ldrd    r9, r10, [r2, #8]
	umull   r11, r12, r4, r7
	umull   r0, r1, r4, r8
	umaal   r12, r0, r5, r7
	umaal   r0, r1, r4, r9
	umull   lr, r2, r5, r8
	umaal   r0, lr, r6, r7
	umaal   r1, r2, r5, r9
	umaal   r1, lr, r6, r8
	umaal   r2, lr, r6, r9
	movs    r3, #0
	umaal   r1, r3, r4, r10
	umaal   r2, r3, r5, r10
	umaal   lr, r3, r6, r10

	@ signed correction
	and     r5, r7, r6, asr #31
	subs    r1, r1, r5
	and     r5, r8, r6, asr #31
	sbcs    r2, r2, r5
	and     r5, r9, r6, asr #31
	sbcs    lr, lr, r5
	and     r5, r10, r6, asr #31
	sbcs    r3, r3, r5
	sbc     r5, r5, r5          @ P[7] from the final borrow

	@ sticky from P[0..3] and the round bit, then P >> (128+b)
	vmov    r4, s10
	rsb     r6, r4, #32
	orr     r11, r11, r12
	orr     r11, r11, r0
	orr     r11, r11, r1
	lsl     r0, r2, r6          @ round bit P[4][b-1] in bit 31
	orr     r11, r11, r0, lsl #1
	lsr     r8, r2, r4
	lsl     r7, lr, r6
	orr     r8, r8, r7          @ q0
	lsr     r9, lr, r4
	lsl     r7, r3, r6
	orr     r9, r9, r7          @ q1
	lsr     r10, r3, r4
	lsl     r7, r5, r6
	orr     r10, r10, r7        @ q2
	orr     r11, r11, r8, lsl #31
	and     r11, r11, r0, asr #31
	cmp     r11, #1             @ C = round_mask & (sticky | q0&1)
	adcs    r8, r8, #0
	adcs    r9, r9, #0
	adc     r10, r10, #0

	vmov    r0, s14
	str     r8, [r0]
	strd    r9, r10, [r0, #4]
	vmov    lr, s15
	bx      lr
#endif
.size fndsa_fxp96_mul_recip_x, .-fndsa_fxp96_mul_recip_x
@ =======================================================================
@ uint64_t fndsa_fxp96_mul_recip_p63_x(const uint32_t a[3],
@                                      const uint32_t rc[4],
@                                      uint32_t k)
@
@ round(a * rc * 2^(k-159)) as a Q0.63 value, for the l10 band.
@ Input: r0 = a, r1 = rc, r2 = k. Output: r0 (low), r1 (high).
@ Clobbers: r0-r12, r14, s10-s12, s15, flags.
@ =======================================================================

.align 2
.global fndsa_fxp96_mul_recip_p63_x
.thumb_func
.type fndsa_fxp96_mul_recip_p63_x, %function
fndsa_fxp96_mul_recip_p63_x:
#if FNDSA_FXP_TIES_UP
    @ Ties up: no sticky, so the sign mask can stay in r12.
    vmov    s15, lr             @ return address; lr is a scratch below

    mov     r3, r0
    ldrd    r5, r9, [r3]        @ a0, a1
    ldr     r0, [r0, #8]        @ a2
    ldrd    r6, r4, [r1]        @ rc0, rc1
    asr.w   r12, r0, #31        @ sign mask, kept
    rsb     r10, r2, #32
    lsl.w   r7, r5, r2          @ am0
    lsr.w   r11, r5, r10
    lsl.w   lr, r9, r2
    ldrd    r5, r8, [r1, #8]    @ rc2, rc3
    lsr.w   r1, r0, r10
    lsls    r0, r2
    lsl.w   r2, r12, r2
    orrs    r2, r1              @ am3
    orr.w   lr, lr, r11         @ am1
    lsr.w   r9, r9, r10
    orr.w   r0, r0, r9          @ am2
    mov     r3, r2              @ am3

    @ am * rc, p0..p6 (p7 dropped)
    umull   r10, r2, r7, r6     @ p0 (dead), high(am0*rc0)
    umull   r9, r11, r7, r4     @ low/high(am0*rc1)
    umaal   r2, r9, lr, r6      @ p1, carry into column 2

    umaal   r9, r11, r7, r5     @ partial p2, carry
    umull   r2, r10, lr, r4     @ low/high(am1*rc1)
    umaal   r9, r2, r0, r6      @ p2, carry into column 3

    umull   r1, r9, r7, r8      @ low/high(am0*rc3)
    umaal   r11, r1, lr, r5     @ partial p3, carry
    umaal   r11, r10, r0, r4    @ partial p3, carry
    umaal   r11, r2, r3, r6     @ p3, carry into column 4

    umaal   r9, r1, lr, r8      @ partial p4, carry
    umaal   r9, r10, r0, r5     @ partial p4, carry
    umaal   r9, r2, r3, r4      @ p4, carry into column 5

    umaal   r1, r10, r0, r8     @ partial p5, carry
    umaal   r1, r2, r3, r5      @ p5, carry into column 6
    umaal   r10, r2, r3, r8     @ p6, discarded p7

    @ signed correction at 2^128, then round bit P4[30] into C
    and.w   r7, r12, r6
    subs.w  r9, r9, r7
    and.w   r7, r12, r4
    sbcs.w  r1, r1, r7
    and.w   r7, r12, r5
    sbcs.w  r10, r10, r7

    lsr.w   r6, r1, #31
    orr.w   r6, r6, r10, lsl #1 @ q1 base
    mov.w   r0, r9, lsr #31
    orr.w   r0, r0, r1, lsl #1  @ q0 base
    lsls    r4, r9, #2          @ C = round bit P4[30]
    adcs.w  r0, r0, #0
    adc.w   r1, r6, #0

    vmov    lr, s15
    bx      lr
#else
    vmov    s15, lr             @ return address; lr is a scratch below

    mov     r3, r0
    ldrd    r5, r9, [r3]        @ a0, a1
    ldr     r0, [r0, #8]        @ a2
    ldrd    r6, r4, [r1]        @ rc0, rc1
    asr.w   r12, r0, #31        @ sign mask
    vmov    s10, r12            @ save sign mask for late correction
    rsb     r10, r2, #32
    lsl.w   r7, r5, r2          @ am0
    lsr.w   r11, r5, r10
    lsl.w   lr, r9, r2
    ldrd    r5, r8, [r1, #8]    @ rc2, rc3
    lsr.w   r1, r0, r10
    lsls    r0, r2
    lsl.w   r2, r12, r2
    orrs    r2, r1              @ am3
    orr.w   lr, lr, r11         @ am1
    lsr.w   r9, r9, r10
    orr.w   r0, r0, r9          @ am2
    mov     r3, r2              @ am3

    @ am * rc, p0..p6 (p7 dropped); r10 collects p0|p1|p2|p3 (sticky)
    umull   r10, r2, r7, r6     @ p0, high(am0*rc0)
    umull   r9, r11, r7, r4     @ low/high(am0*rc1)
    umaal   r2, r9, lr, r6      @ p1, carry into column 2
    orr.w   r10, r10, r2

    umaal   r9, r11, r7, r5     @ partial p2, carry
    umull   r2, r12, lr, r4     @ low/high(am1*rc1)
    umaal   r9, r2, r0, r6      @ p2, carry into column 3
    orr.w   r10, r10, r9

    umull   r1, r9, r7, r8      @ low/high(am0*rc3)
    umaal   r11, r1, lr, r5     @ partial p3, carry
    umaal   r11, r12, r0, r4    @ partial p3, carry
    umaal   r11, r2, r3, r6     @ p3, carry into column 4
    orr.w   r10, r10, r11

    umaal   r9, r1, lr, r8      @ partial p4, carry
    umaal   r9, r12, r0, r5     @ partial p4, carry
    umaal   r9, r2, r3, r4      @ p4, carry into column 5

    umaal   r1, r12, r0, r8     @ partial p5, carry
    umaal   r1, r2, r3, r5      @ p5, carry into column 6
    umaal   r12, r2, r3, r8     @ p6, discarded p7

    @ signed correction at 2^128 (and.w keeps the borrow)
    vmov    r3, s10             @ sign mask
    and.w   r7, r3, r6
    subs.w  r9, r9, r7
    and.w   r7, r3, r4
    sbcs.w  r1, r1, r7
    and.w   r7, r3, r5
    sbcs.w  r12, r12, r7

    @ C = round bit & (sticky | q0 & 1)
    lsl.w   r4, r9, #1          @ P4[30] -> bit 31 (round mask)
    orr.w   r10, r10, r4, lsl #1 @ P0..P3 | P4[29:0] (sticky)
    lsrs    r6, r1, #31
    orr.w   r6, r6, r12, lsl #1 @ q1 base
    mov.w   r0, r9, lsr #31
    orr.w   r0, r0, r1, lsl #1  @ q0 base
    orr.w   r10, r10, r0, lsl #31 @ sticky | quotient LSB
    and.w   r10, r10, r4, asr #31 @ keep it iff round bit is set
    cmp.w   r10, #1             @ C = RNE increment
    adcs.w  r0, r0, #0
    adc.w   r1, r6, #0

    vmov    lr, s15
    bx      lr
#endif
	.size fndsa_fxp96_mul_recip_p63_x, .-fndsa_fxp96_mul_recip_p63_x
	@ =======================================================================
	@ void fndsa_fxp96_recip_fast_raw_x(uint32_t m[4],
	@                                   const uint32_t x[3])
	@
	@ m = unrounded 1/x in Q126 for x normalized by k, as fxp_recip_fast in
	@ ops_c.c. Requires 1 <= x < 2^31, so k = clz(x[2]) is in [1, 31].
	@ Input: r0 = m, r1 = x. Output: r12 = k.
	@ Clobbers: r0-r11, flags, memory.
	@ =======================================================================

	.align 2
	.global fndsa_fxp96_recip_fast_raw_x
	.thumb_func
	.type fndsa_fxp96_recip_fast_raw_x, %function
fndsa_fxp96_recip_fast_raw_x:
	@ d = x << k in r6:r5:r4
	ldmia   r1!, {r4, r5, r6}
	clz     r12, r6
	rsb     r1, r12, #32
	lsr     r3, r4, r1
	lsr     r2, r5, r1
	lsl     r4, r4, r12
	lsl     r5, r5, r12
	orr     r5, r5, r3
	lsl     r6, r6, r12
	orr     r6, r6, r2

	@ Seed and two steps on d32 = r6, then one Q62 step on d64. y < 2,
	@ kept in Q31 in r9. Seed constants as in ops_c.c.
	mov.w   r1, #0xF0F0F0F0
	mov.w   r11, #0x69696969
	umull   r10, r9, r6, r1             @ r9 = high32(r1 * d32)
	sub     r9, r11, r9                 @ z0 = 48/17 - that (Q31)

	@ y = high32(y*(2 - d*y)): y0 Q31, y1 Q30, y2 Q28, shifted to Q31
	umull   r1, r2, r6, r9              @ t = d*y0
	rsb     r2, r2, #0                  @ e = 2 - t
	umull   r1, r3, r9, r2              @ r3 = y1 (Q30)
	umull   r1, r2, r6, r3              @ t = d*y1
	rsb     r2, r2, #0x80000000         @ e = 2 - t
	umull   r1, r9, r3, r2              @ r9 = y2 (Q28)
	lsl     r9, r9, #3                  @ y2 in Q31

	@ Q62 step: t = (d64*z)>>32 (Q63), e = 2 - t, y3 = (z*e)>>32 in r8:r7
	movs    r3, #0
	movs    r8, #0
	umull   r1, r2, r5, r9
	umaal   r2, r3, r6, r9              @ t = r3:r2
	rsbs    r2, r2, #0
	sbc     r11, r8, r3                 @ e = r11:r2
	umull   r1, r7, r9, r2
	umaal   r7, r8, r9, r11             @ y3 = r8:r7

	@ P = d*y in r11:r10:r9:r3:r2, then E = (2^159 - P) >> 32
	umull   r2, r3, r4, r7
	umull   r9, r10, r4, r8
	umaal   r3, r9, r5, r7
	umaal   r9, r10, r5, r8
	umull   r1, r11, r6, r7
	adds    r9, r9, r1
	adc     r10, r10, #0
	umaal   r10, r11, r6, r8
	rsbs    r2, r2, #0
	mov.w   r1, #0
	sbcs    r4, r1, r3
	sbcs    r5, r1, r9
	sbcs    r6, r1, r10
	mvn     r9, r11
	adc     r9, r9, #0x80000000         @ 2^31 - P4 - borrow

	@ m = (y*E) >> 62, with E in r9:r6:r5:r4
	umull   r1, r2, r4, r7
	umull   r3, r10, r4, r8
	umaal   r2, r3, r5, r7
	mov     r4, r2
	umaal   r3, r10, r5, r8
	movs    r2, #0
	umaal   r3, r2, r6, r7
	mov     r5, r3
	umaal   r10, r2, r6, r8
	movs    r3, #0
	umaal   r10, r3, r9, r7
	mov     r6, r10
	umaal   r2, r3, r9, r8
	lsrs    r1, r4, #30
	orr     r1, r1, r5, lsl #2
	str     r1, [r0]
	lsrs    r1, r5, #30
	orr     r1, r1, r6, lsl #2
	str     r1, [r0, #4]
	lsrs    r1, r6, #30
	orr     r1, r1, r2, lsl #2
	str     r1, [r0, #8]
	lsrs    r1, r2, #30
	orr     r1, r1, r3, lsl #2
	str     r1, [r0, #12]
	bx      lr
	.size fndsa_fxp96_recip_fast_raw_x, .-fndsa_fxp96_recip_fast_raw_x

	@ =======================================================================
	@ fndsa_rsqrt_fast_g13_body
	@
	@ out = fxp_rsqrt_fast(m) (ops_c.c), the unrounded 1/sqrt(m/2^96) in
	@ Q94, for m in [2^94, 2^96). Unsigned and truncating throughout.
	@ Input: r0 = out, r1 = m. Clobbers: r1-r12, lr, flags.
	@ =======================================================================
	.macro fndsa_rsqrt_fast_g13_body

	@ Seed and three steps on m32 = r5, one Q62 step on m64 = r5:r4.
	@ y < 2, kept in Q31 in r10. fp = 0xC0000000 = 3 in Q30.
	ldrd    r4, r5, [r1, #4]
	ldr     r1, [r1]
	mov.w   r6, #0xAAAAAAAA             @ slope 4/3 in Q31
	mov.w   r8, #0x1B1B1B1B             @ intercept 2.2118 in Q31 mod 2^32
	mov     fp, #-1073741824
	umull   r9, r10, r5, r6             @ r10 = high32(r6 * m32)
	sub     r10, r8, r10                @ y0 (Q31)

	@ a = m*y, b = a*y, y = y*(3 - b): y1 Q30, y2 Q27 << 4, y3 Q30 << 1
	umull   r6, ip, r5, r10
	umull   r6, r8, ip, r10
	rsb     r8, r8, fp
	umull   r6, r9, r10, r8             @ r9 = y1 (Q30)
	umull   r6, ip, r5, r9
	umull   r6, r8, ip, r9
	rsb     r8, r8, #0x30000000
	umull   r6, r10, r9, r8             @ r10 = y2 (Q27)
	lsl     r10, r10, #4                @ Q31
	umull   r6, ip, r5, r10
	umull   r6, r8, ip, r10
	rsb     r8, r8, fp
	umull   r6, r9, r10, r8             @ r9 = y3 (Q30)
	lsl     r10, r9, #1                 @ Q31

	@ Q62 step: a = (m64*z)>>32, b = (a*z)>>32, e = 3 - b, y = (z*e)>>32
	@ in ip:r2
	movs    lr, #0
	movs    r9, #0
	umull   r6, ip, r4, r10
	umaal   ip, lr, r5, r10             @ a = lr:ip
	umull   r6, r8, ip, r10
	umaal   r8, r9, lr, r10             @ b = r9:r8
	rsbs    r8, r8, #0
	sbc     lr, fp, r9                  @ e = lr:r8
	umull   r6, r2, r10, r8
	movs    ip, #0
	umaal   r2, ip, r10, lr             @ y = ip:r2
	movs    r3, #0

	@ full 96-bit step: a = (m*y) >> 64, b = (a*y) >> 62,
	@ out = ((3 - b) * y) >> 63

	@ a = (m * y) >> 64; a = r9:lr:r8.
	umull   r6, r7, r1, r2
	umull   r8, lr, r1, ip
	umaal   r7, r8, r4, r2
	umaal   r8, lr, r4, ip
	mov     r9, #0
	umaal   r8, r9, r5, r2
	umaal   lr, r9, r5, ip

	@ b = (a * y) >> 62 = lr:r9:r4
	umull   r4, r5, r8, r2
	umull   r6, r7, r8, ip
	umaal   r5, r6, lr, r2
	umaal   r6, r7, lr, ip
	mov     r10, #0
	umaal   r6, r10, r9, r2
	umaal   r7, r10, r9, ip
	lsrs    r4, r5, #30
	orr     r4, r4, r6, lsl #2
	lsrs    r9, r6, #30
	orr     r9, r9, r7, lsl #2
	lsrs    lr, r7, #30
	orr     lr, lr, r10, lsl #2

	@ e = THREE_Q94 - b.
	rsbs    r4, r4, #0
	sbcs    r9, r3, r9
	sbc     lr, fp, lr

	@ out = (e * y) >> 63.
	umull   r5, r6, r4, r2
	umull   r7, r8, r4, ip
	umaal   r6, r7, r9, r2
	umaal   r7, r8, r9, ip
	mov     r10, #0
	umaal   r7, r10, lr, r2
	umaal   r8, r10, lr, ip
	lsrs    r6, r6, #31
	orr     r6, r6, r7, lsl #1
	lsrs    r2, r7, #31
	orr     r2, r2, r8, lsl #1
	lsrs    r4, r8, #31
	orr     r4, r4, r10, lsl #1
	strd    r6, r2, [r0]
	str     r4, [r0, #8]
	.endm

	@ =======================================================================
	@ void fndsa_fxp96_sqrt_mul_inv_sigma_x(uint32_t d[3],
	@                                      const uint32_t x[3],
	@                                      unsigned logn)
	@
	@ d = round(x * rsqrt(m) * INV_SIGMA_Q71[logn] >> (181 - t/2)), with
	@ m = x << t, as fndsa_fxp96_sqrt_mul_inv_sigma in ops_c.c. Requires
	@ 1 <= x < 2^31 (the gamma checks give x > 4).
	@ Input: r0 = d, r1 = x, r2 = logn.
	@ Clobbers: r0-r12, r14, s9-s15, flags, memory.
	@
	@ Frame (68 bytes): [4] high Q71 limb, [24] d, [36] xy p3 (ties to
	@ even only), [40..48] m, [52..60] y.
	@ =======================================================================

	.align 2
	.global fndsa_fxp96_sqrt_mul_inv_sigma_x
	.thumb_func
	.type fndsa_fxp96_sqrt_mul_inv_sigma_x, %function
fndsa_fxp96_sqrt_mul_inv_sigma_x:
#if FNDSA_FXP_TIES_UP
	@ t = clz(x[2]) & ~1 is in [0, 30]; m = x << t.
	@ r4 = logn, r5 = t, r6:r9:r10 = x[2]:x[1]:x[0], [sp,#40..48] = m.
	vmov    s15, lr
	sub     sp, sp, #68
	mov     r4, r2
	str     r0, [sp, #24]
	ldrd    r10, r9, [r1]
	ldr     r6, [r1, #8]
	clz     r3, r6
	bic     r5, r3, #1
	rsb     r3, r5, #32
	lsl     r0, r10, r5
	lsr     r1, r10, r3
	lsl     r2, r9, r5
	orr     r1, r1, r2
	lsr     r2, r9, r3
	lsl     r3, r6, r5
	orr     r2, r2, r3
	str     r0, [sp, #40]
	str     r1, [sp, #44]
	str     r2, [sp, #48]
	add     r0, sp, #52
	add     r1, sp, #40

	@ rsqrt, with s9:s10 = logn:t, s11:s12:s13 = x[2]:x[1]:x[0]
	vmov    s9, s10, r4, r5
	vmov    s11, s12, r6, r9
	vmov    s13, r10
	fndsa_rsqrt_fast_g13_body
	vmov    r4, r5, s9, s10
	vmov    r6, r9, s11, s12
	vmov    r10, s13

	@ sh = 181 - t/2 is in [166, 181], so the result is P[5..7] >> b with
	@ b = sh - 160 in [6, 21].
	lsrs    r3, r5, #1
	rsb     r3, r3, #21
	ldr     r5, .Lsqrt_sigma_invq71
	add     r5, r5, r4, lsl #3
	ldrd    r1, r2, [r5]
	@ Ties up: no sticky word, so xy p3 stays in r6.
	str     r2, [sp, #4]
	mov     r5, r3              @ b
	ldr     r2, [sp, #52]

	@ xy = x * y: p0 = r4, p1 = r7, p2 = r3, p3 = r6, p4 = r0, p5 = r9
	umull   r4, r7, r10, r2
	ldr     r11, [sp, #56]
	umull   r8, lr, r10, r11
	umaal   r7, r8, r9, r2
	ldr     r12, [sp, #60]
	umaal   r8, lr, r10, r12
	umull   r0, r3, r9, r11
	umaal   r8, r0, r6, r2
	umaal   lr, r3, r9, r12
	umaal   lr, r0, r6, r11
	umaal   r3, r0, r6, r12
	mov     r9, r0
	mov     r0, r3
	mov     r3, r8
	mov     r6, lr
	ldr     r2, [sp, #4]

	@ xy * INV_SIGMA_Q71; p0..p3 only for their carries
	umull   r12, r8, r4, r1
	umull   r10, r11, r4, r2
	umaal   r8, r10, r7, r1
	umaal   r10, r11, r7, r2
	mov     r12, #0
	umaal   r10, r12, r3, r1
	umaal   r11, r12, r3, r2
	mov     lr, #0
	umaal   r11, lr, r6, r1
	umaal   r12, lr, r6, r2
	mov     fp, #0
	umaal   r12, fp, r0, r1
	umaal   lr, fp, r0, r2
	umlal   lr, fp, r9, r1
	mov     r10, #0
	umaal   fp, r10, r9, r2

	@ P[5..7] = lr:fp:r10 >> b; the first shift leaves the round bit in C
	rsb     r8, r5, #32
	lsrs    r9, lr, r5          @ C = round bit P5[b-1]
	lsl     r4, fp, r8
	orr     r9, r9, r4
	lsr     r4, fp, r5
	lsl     r0, r10, r8
	orr     r11, r4, r0
	lsr     r12, r10, r5
	adcs    r9, r9, #0
	adc     r11, r11, #0
	adc     r12, r12, #0
	ldr     r3, [sp, #24]

	@ store
	strd    r9, r11, [r3]
	str     r12, [r3, #8]
	add     sp, sp, #68
	vmov    lr, s15
	bx      lr
	.align 2
.Lsqrt_sigma_invq71:
	.word   fndsa_fpr_INV_SIGMA_Q71
#else
	@ t = clz(x[2]) & ~1 is in [0, 30]; m = x << t.
	@ r4 = logn, r5 = t, r6:r9:r10 = x[2]:x[1]:x[0], [sp,#40..48] = m.
	vmov    s15, lr
	sub     sp, sp, #68
	mov     r4, r2
	str     r0, [sp, #24]
	ldrd    r10, r9, [r1]
	ldr     r6, [r1, #8]
	clz     r3, r6
	bic     r5, r3, #1
	rsb     r3, r5, #32
	lsl     r0, r10, r5
	lsr     r1, r10, r3
	lsl     r2, r9, r5
	orr     r1, r1, r2
	lsr     r2, r9, r3
	lsl     r3, r6, r5
	orr     r2, r2, r3
	str     r0, [sp, #40]
	str     r1, [sp, #44]
	str     r2, [sp, #48]
	add     r0, sp, #52
	add     r1, sp, #40

	@ rsqrt, with s9:s10 = logn:t, s11:s12:s13 = x[2]:x[1]:x[0]
	vmov    s9, s10, r4, r5
	vmov    s11, s12, r6, r9
	vmov    s13, r10
	fndsa_rsqrt_fast_g13_body
	vmov    r4, r5, s9, s10
	vmov    r6, r9, s11, s12
	vmov    r10, s13

	@ sh = 181 - t/2 is in [166, 181], so the result is P[5..7] >> b with
	@ b = sh - 160 in [6, 21].
	lsrs    r3, r5, #1
	rsb     r3, r3, #21
	ldr     r5, .Lsqrt_sigma_invq71
	add     r5, r5, r4, lsl #3
	ldrd    r1, r2, [r5]
	str     r2, [sp, #4]
	mov     r5, r3              @ b
	ldr     r2, [sp, #52]

	@ xy = x * y: p0 = r4, p1 = r7, p2 = r3, p3 = [sp,#36], p4 = r0,
	@ p5 = r9
	umull   r4, r7, r10, r2
	ldr     r11, [sp, #56]
	umull   r8, lr, r10, r11
	umaal   r7, r8, r9, r2
	ldr     r12, [sp, #60]
	umaal   r8, lr, r10, r12
	umull   r0, r3, r9, r11
	umaal   r8, r0, r6, r2
	umaal   lr, r3, r9, r12
	umaal   lr, r0, r6, r11
	umaal   r3, r0, r6, r12
	str     lr, [sp, #36]
	mov     r9, r0
	mov     r0, r3
	mov     r3, r8
	ldr     r2, [sp, #4]

	@ xy * INV_SIGMA_Q71: r6 = P0|P1|P2|P3 (sticky), r12:lr:fp:r10 = P4..P7.
	@ p5 * Q71 limb 0 is umlal, as lr:fp already holds P5 and the P6 carry.
	umull   r6, r8, r4, r1
	umull   r10, r11, r4, r2
	umaal   r8, r10, r7, r1
	orr     r6, r6, r8
	umaal   r10, r11, r7, r2
	mov     r12, #0
	umaal   r10, r12, r3, r1
	orr     r6, r6, r10
	umaal   r11, r12, r3, r2
	ldr     r4, [sp, #36]
	mov     lr, #0
	umaal   r11, lr, r4, r1
	orr     r6, r6, r11
	umaal   r12, lr, r4, r2
	mov     fp, #0
	umaal   r12, fp, r0, r1
	umaal   lr, fp, r0, r2
	umlal   lr, fp, r9, r1
	mov     r10, #0
	umaal   fp, r10, r9, r2

	@ P[5..7] >> b, C = round bit & (sticky | q0 & 1)
	orr     r6, r6, r12
	rsb     r8, r5, #32
	lsl     r4, lr, r8          @ round bit P5[b-1] in bit 31
	lsr     r9, lr, r5
	lsl     r0, fp, r8
	orr     r9, r9, r0          @ q0
	lsr     r0, fp, r5
	lsl     r2, r10, r8
	orr     r11, r0, r2         @ q1
	lsr     r12, r10, r5        @ q2
	orr     r6, r6, r4, lsl #1  @ P0..P4 | P5 below the round bit
	orr     r6, r6, r9, lsl #31 @ sticky | retained low bit (tie-to-even)
	and     r6, r6, r4, asr #31 @ keep it iff P5 bit b-1 is set
	cmp     r6, #1               @ C = RNE increment
	adcs    r9, r9, #0
	adc     r11, r11, #0
	adc     r12, r12, #0
	ldr     r3, [sp, #24]

	@ store
	strd    r9, r11, [r3]
	str     r12, [r3, #8]
	add     sp, sp, #68
	vmov    lr, s15
	bx      lr
	.align 2
.Lsqrt_sigma_invq71:
	.word   fndsa_fpr_INV_SIGMA_Q71
#endif
	.size fndsa_fxp96_sqrt_mul_inv_sigma_x, .-fndsa_fxp96_sqrt_mul_inv_sigma_x
