	.syntax	unified
	.cpu	cortex-m4
	.file	"sign_fpoly.S"
	.text

#include "fpr/fxp32.64/fxp_round_cm4.h"

@ Fused FFT butterflies for fxp32.64 (FNDSA_TWIDDLE_U64=1, FNDSA_FFT_ASM=1):
@ one call loads x and y, multiplies by the twiddle, does both add/subs
@ and stores o1 and o2. They compute the same results as the FPOLY_GM_*
@ macros in fxp_fft.h.
@
@ Custom ABI, out of place (merge/split, FNDSA_FXP96_BFLY_CALL):
@   r0 = A: x_re, x_im at +r5     r2 = O1: o1_re, o1_im at +r6
@   r1 = B: y_re, y_im at +r5     r3 = O2: o2_re, o2_im at +r6
@   r4 = &GM[k << 1]: cos magnitude at +0, sin at +8, unsigned Q0.64
@ In place (FFT/iFFT, the _ip entries, FNDSA_FXP96_BFLY_CALL_IP):
@   r0 = A (x and o1), r1 = B (y and o2), r2 = gm, r3 = offset to the
@   imaginary halves. The _ip entries just widen that to the form above
@   (O1 = A, O2 = B, one offset) and run the same bodies; every slot is
@   read before it is written.
@ Clobbers r0-r12, r14, s5-s15, flags; returns with bx lr.
@ Twiddle w = A + iB for even k, -(A - iB) for odd k.
@
@ A complex product is two components, P1*B + P2*A then -P2*B + P1*A,
@ with P1 the half loaded first. Each is one call to bfly_kern: operand
@ r4:r5:r6, B in r7:r8, A in lr:r11, accumulator r0 r10 r2:r3:r9,
@ scratch r12. bl clobbers lr, so A0 is kept in s5 and the pointer to
@ the other half in s6.

@ Fresh product x*B into r0 r10 r2:r3:r9, x's sign correction.
.macro bfly_mul_b
    umull   r0, r10, r4, r7
    umull   r2, r3, r4, r8
    umaal   r10, r2, r5, r7
    umaal   r2, r3, r5, r8
    mov     r9, #0
    umaal   r2, r9, r6, r7
    umaal   r3, r9, r6, r8
    and     r12, r7, r6, asr #31
    subs    r3, r3, r12
    and     r12, r8, r6, asr #31
    sbcs    r9, r9, r12
.endm

@ Accumulate x*A; the column-3 carry stays in r12 for bfly_round.
.macro bfly_mul_a_acc
    mov     r12, #0
    umaal   r0, r12, r4, lr
    umaal   r10, r12, r4, r11
    umaal   r2, r12, r5, r11
    umaal   r3, r12, r6, r11
    add     r9, r9, r12
    mov     r12, #0
    umaal   r10, r12, r5, lr
    umaal   r2, r12, r6, lr
.endm

@ Rounding at 2^64 (fxp_round_carry), then x's sign correction against A.
.macro bfly_round
    fxp_round_carry r0, r0, r10, r2
    adcs    r2, r2, #0
    adcs    r3, r3, r12
    adc     r9, r9, #0
    and     r12, lr, r6, asr #31
    subs    r3, r3, r12
    and     r12, r11, r6, asr #31
    sbcs    r9, r9, r12
.endm

@ r4:r5:r6 <- [\p] (\p may be r4).
.macro bfly_ld3 p
    ldr     r6, [\p, #8]
    ldr     r5, [\p, #4]
    ldr     r4, [\p, #0]
.endm

@ r4:r5:r6 <- -r4:r5:r6.
.macro bfly_neg
    mov     r10, #0
    rsbs    r4, r4, #0
    sbcs    r5, r10, r5
    sbcs    r6, r10, r6
.endm

@ fpr_half in place on \a:\b:\c (low to high), temp \t.
.macro bfly_half a, b, c, t
    and     \t, \a, \a, lsr #1
    and     \t, \t, #1
    asrs    \c, \c, #1
    rrxs    \b, \b
    rrx     \a, \a
    adds    \a, \a, \t
    adcs    \b, \b, #0
    adc     \c, \c, #0
.endm

@ Forward add/sub of one half: x at [r0], z in r2:r3:r9;
@ x - z to [\dif], x + z to [\sum] (s-registers holding the pointers).
.macro bfly_addsub sum, dif
    ldr     r6, [r0, #8]
    ldr     r5, [r0, #4]
    ldr     r0, [r0, #0]
    vmov    r12, \dif
    subs    r10, r0, r2
    str     r10, [r12, #0]
    sbcs    r10, r5, r3
    str     r10, [r12, #4]
    sbcs    r10, r6, r9
    str     r10, [r12, #8]
    vmov    r12, \sum
    adds    r0, r0, r2
    str     r0, [r12, #0]
    adcs    r5, r5, r3
    str     r5, [r12, #4]
    adcs    r6, r6, r9
    str     r6, [r12, #8]
.endm

@ One product component, shared by every entry: x (r4:r5:r6) * B (r7:r8)
@ + [s6] * A (s5:r11), rounded once at 2^64 into r2:r3:r9; the other
@ half is left in r4:r5:r6. Keeps r1, r7, r8, r11; clobbers r0, r10,
@ r12, lr, s7. lr is data inside, so the return address goes via s7.
    .align 2
    .thumb_func
    .type fndsa_fxp96_bfly_kern, %function
fndsa_fxp96_bfly_kern:
    vmov    s7, lr
    vmov    lr, s5
    bfly_mul_b
    vmov    r4, s6
    bfly_ld3 r4
    bfly_mul_a_acc
    bfly_round
    vmov    lr, s7
    bx      lr
.size fndsa_fxp96_bfly_kern, .-fndsa_fxp96_bfly_kern

@ Forward butterfly, o1 = x + y*w, o2 = x - y*w. Even: imaginary
@ component first. Odd: real first and the sum/difference targets
@ exchanged, since w = -(A - iB).
@ Parked: s8/s9 x halves for components 2/1, s10/s11 sum targets 2/1,
@ s12/s13 difference targets 2/1, s15 lr.
.macro bfly_fwd odd
    add     r7, r0, r5
    add     r8, r2, r6
    add     r9, r3, r6
  .if \odd
    vmov    s8, s9, r7, r0
    vmov    s10, s11, r9, r3
    vmov    s12, s13, r8, r2
    vmov    s6, r1
    add     r1, r1, r5
  .else
    add     r10, r1, r5
    vmov    s8, s9, r0, r7
    vmov    s10, s11, r2, r8
    vmov    s12, s13, r3, r9
    vmov    s6, r10
  .endif
    vmov    s15, lr

    ldr     r7, [r4, #8]
    ldr     r8, [r4, #12]
    ldr     r11, [r4, #4]
    ldr     r4, [r4, #0]
    vmov    s5, r4
    bfly_ld3 r1
    bl      fndsa_fxp96_bfly_kern

    @ -P2 survives the add/sub: limb 0 in r4, limbs 1-2 in s9/s14.
    bfly_neg
    vmov    r0, s9
    vmov    s9, r5
    vmov    s14, r6
    bfly_addsub s11, s13
    vmov    r5, s9
    vmov    r6, s14

    vmov    s6, r1
    bl      fndsa_fxp96_bfly_kern
    vmov    r0, s8
    bfly_addsub s10, s12
    vmov    lr, s15
    bx      lr
.endm

@ Inverse add/sub of the P2 half: x at [r0], y at [r1], sum to [r2]
@ (halved if \halve), difference (x - y, or y - x if \odd) to [r3]; with
@ \adv every address is offset by r5 (inputs) / r6 (outputs), r4 as the
@ address temp. Leaves the difference slot in s6. Limbs r7-r12, temp lr.
.macro bfly_inv_p2 odd, halve, adv
  .if \adv
    add     r4, r0, r5
    ldr     r7, [r4, #0]
    ldr     r8, [r4, #4]
    ldr     r9, [r4, #8]
    add     r4, r1, r5
    ldr     r10, [r4, #0]
    ldr     r11, [r4, #4]
    ldr     r12, [r4, #8]
    add     r4, r3, r6
  .else
    ldr     r7, [r0, #0]
    ldr     r8, [r0, #4]
    ldr     r9, [r0, #8]
    ldr     r10, [r1, #0]
    ldr     r11, [r1, #4]
    ldr     r12, [r1, #8]
    mov     r4, r3
  .endif
  .if \odd
    subs    lr, r10, r7
    str     lr, [r4, #0]
    sbcs    lr, r11, r8
    str     lr, [r4, #4]
    sbcs    lr, r12, r9
    str     lr, [r4, #8]
  .else
    subs    lr, r7, r10
    str     lr, [r4, #0]
    sbcs    lr, r8, r11
    str     lr, [r4, #4]
    sbcs    lr, r9, r12
    str     lr, [r4, #8]
  .endif
    vmov    s6, r4
  .if \adv
    adds    r4, r2, r6
  .else
    mov     r4, r2
  .endif
    adds    r7, r7, r10
    adcs    r8, r8, r11
    adcs    r9, r9, r12
  .if \halve
    bfly_half r7, r8, r9, lr
  .endif
    str     r7, [r4, #0]
    str     r8, [r4, #4]
    str     r9, [r4, #8]
.endm

@ The P1 half, base addresses: sum to [r2], difference to r4:r6 and [r3]
@ (the kernel multiplies it from the registers and reloads it once).
.macro bfly_inv_p1 odd, halve
    ldr     r7, [r0, #0]
    ldr     r8, [r0, #4]
    ldr     r9, [r0, #8]
    ldr     r10, [r1, #0]
    ldr     r11, [r1, #4]
    ldr     r12, [r1, #8]
  .if \odd
    subs    r4, r10, r7
    sbcs    r5, r11, r8
    sbcs    r6, r12, r9
  .else
    subs    r4, r7, r10
    sbcs    r5, r8, r11
    sbcs    r6, r9, r12
  .endif
    adds    r7, r7, r10
    str     r4, [r3, #0]
    adcs    r8, r8, r11
    str     r5, [r3, #4]
    adcs    r9, r9, r12
    str     r6, [r3, #8]
  .if \halve
    bfly_half r7, r8, r9, r10
  .endif
    str     r7, [r2, #0]
    str     r8, [r2, #4]
    str     r9, [r2, #8]
.endm

@ Inverse butterfly, o1 = x + y, o2 = (x - y)*conj(w). Even: gmj
@ kernel (conj(w) = A - iB), P1 = imaginary half. Odd: y - x and the
@ gm kernel (conj(w) = -(A + iB)), P1 = real half. \halve: split step,
@ every output through fpr_half, the product after its rounding. The
@ P2 half goes first and its difference to the o2 slot; the P1 half
@ leaves its difference in r4:r6 for the kernel (stored too, for the
@ second component's reload).
.macro bfly_inv odd, halve
    vmov    s14, s15, r4, lr
  .if \odd
    bfly_inv_p2 1, \halve, 1
    vmov    lr, s14
    bfly_inv_p1 1, \halve
  .else
    bfly_inv_p2 0, \halve, 0
    vmov    lr, s14
    add     r0, r0, r5
    add     r1, r1, r5
    add     r2, r2, r6
    add     r3, r3, r6
    bfly_inv_p1 0, \halve
  .endif
    mov     r1, r3
    ldr     r7, [lr, #8]
    ldr     r8, [lr, #12]
    ldr     r11, [lr, #4]
    ldr     lr, [lr, #0]
    vmov    s5, lr
    bl      fndsa_fxp96_bfly_kern

    bfly_neg
  .if \halve
    bfly_half r2, r3, r9, r0
  .endif
    vmov    r12, s6
    str     r2, [r12, #0]
    str     r3, [r12, #4]
    str     r9, [r12, #8]

    vmov    s6, r1
    bl      fndsa_fxp96_bfly_kern
  .if \halve
    bfly_half r2, r3, r9, r0
  .endif
    str     r2, [r1, #0]
    str     r3, [r1, #4]
    str     r9, [r1, #8]
    vmov    lr, s15
    bx      lr
.endm

@ In-place entries take r0 = A, r1 = B, r2 = gm, r3 = off and widen
@ that to the 7-register form (O1 = A, O2 = B, off_in = off_out).
.macro bfly_ip_args
    mov     r4, r2
    mov     r5, r3
    mov     r6, r3
    mov     r2, r0
    mov     r3, r1
.endm

@ Entry \ip_name (in place) widens the arguments and falls through into
@ entry \name (7-register), whose body is \body with arguments \args.
.macro bfly_entry name, ip_name, body, args:vararg
    .align 2
    .global \ip_name
    .thumb_func
    .type \ip_name, %function
\ip_name:
    bfly_ip_args
.size \ip_name, .-\ip_name
    .global \name
    .thumb_func
    .type \name, %function
\name:
    \body \args
.size \name, .-\name
.endm

bfly_entry fndsa_fxp96_bfly_fwd_even_x, fndsa_fxp96_bfly_fwd_even_ip_x, bfly_fwd, 0
bfly_entry fndsa_fxp96_bfly_fwd_odd_x,  fndsa_fxp96_bfly_fwd_odd_ip_x,  bfly_fwd, 1
bfly_entry fndsa_fxp96_bfly_inv_even_x, fndsa_fxp96_bfly_inv_even_ip_x, bfly_inv, 0, 0
bfly_entry fndsa_fxp96_bfly_inv_odd_x,  fndsa_fxp96_bfly_inv_odd_ip_x,  bfly_inv, 1, 0

    .align 2
    .global fndsa_fxp96_bfly_split_even_x
    .thumb_func
    .type fndsa_fxp96_bfly_split_even_x, %function
fndsa_fxp96_bfly_split_even_x:
    bfly_inv 0, 1
.size fndsa_fxp96_bfly_split_even_x, .-fndsa_fxp96_bfly_split_even_x

    .align 2
    .global fndsa_fxp96_bfly_split_odd_x
    .thumb_func
    .type fndsa_fxp96_bfly_split_odd_x, %function
fndsa_fxp96_bfly_split_odd_x:
    bfly_inv 1, 1
.size fndsa_fxp96_bfly_split_odd_x, .-fndsa_fxp96_bfly_split_odd_x
