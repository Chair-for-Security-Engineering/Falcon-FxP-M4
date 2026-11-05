@ Copyright 2025 NXP
@ SPDX-License-Identifier: MIT

.syntax	unified
.cpu	cortex-m4
.file	"triple_float_ct.s"
.text

	.align	2
	.global tw_sum_ct_asm	
	.thumb
	.thumb_func
	.type	tw_sum_ct_asm, %function
tw_sum_ct_asm:
         push	{r4}
         sub	sp, #64	@ 0x40
         vmov r2, s0
        vmov s0, lr
         vmov r1, s3
         vmov r3, s2
         vmov r4, s5
         bic.w	r0, r1, #2147483648	@ 0x80000000
         bic.w	lr, r2, #2147483648	@ 0x80000000
         sub.w	lr, lr, r0
         bic.w	ip, r3, #2147483648	@ 0x80000000
         bic.w	r0, r4, #2147483648	@ 0x80000000
         sub.w	ip, ip, r0
         eor.w	r0, r2, r1
         and.w	r0, r0, lr, asr #31
         eors	r2, r0
         vmov	s10, r2
         eor.w	r2, r3, r4
         and.w	r2, r2, ip, asr #31
         eors	r1, r0
         eors	r3, r2
         eors	r2, r4
         vmov	s8, r2
         bic.w	r0, r3, #2147483648	@ 0x80000000
         bic.w	r2, r1, #2147483648	@ 0x80000000
         subs	r2, r0, r2
         eor.w	r4, r1, r3
        nop
         vmov r0, s4
         and.w	r4, r4, r2, asr #31
         vmov r2, s1
         eors	r1, r4
         eors	r3, r4
         bic.w	ip, r2, #2147483648	@ 0x80000000
         bic.w	r4, r0, #2147483648	@ 0x80000000
         sub.w	ip, ip, r4
         eor.w	r4, r2, r0
         and.w	r4, r4, ip, asr #31
         eors	r0, r4
         eors	r2, r4
         bic.w	ip, r1, #2147483648	@ 0x80000000
         bic.w	r4, r0, #2147483648	@ 0x80000000
         sub.w	ip, ip, r4
         eor.w	r4, r1, r0
         and.w	r4, r4, ip, asr #31
         eors	r0, r4
         eors	r1, r4
         vmov	s7, r0
         vmov	s4, r1
         bic.w	r0, r2, #2147483648	@ 0x80000000
         bic.w	r1, r3, #2147483648	@ 0x80000000
         subs	r0, r0, r1
         eor.w	r1, r3, r2
         vadd.f32	s5, s8, s7
         and.w	r1, r1, r0, asr #31
         eors	r3, r1
         vadd.f32	s15, s5, s4
         vmov	s9, r3
         eor.w	r3, r2, r1
         vadd.f32	s14, s15, s9
         vmov	s11, r3
         vadd.f32	s13, s14, s11
         vsub.f32	s3, s14, s15
         vadd.f32	s6, s13, s10
         vsub.f32	s2, s13, s14
         vsub.f32	s1, s6, s13
         vsub.f32	s12, s13, s2
         vsub.f32	s10, s10, s1
         vsub.f32	s1, s6, s1
        vsub.f32	s11, s11, s2
         vsub.f32	s13, s13, s1
        vsub.f32	s12, s14, s12
         vadd.f32	s13, s13, s10
        vadd.f32	s12, s12, s11
         vadd.f32	s10, s6, s13
         vsub.f32	s6, s10, s6
         vmov	r4, s10
        vsub.f32	s13, s13, s6
        vsub.f32	s14, s14, s3
         vmov	r3, s13
         negs	r1, r3
         orrs	r1, r3
         eors	r3, r4
         and.w	r3, r3, r1, asr #31
         eors	r3, r4
         vmov	s13, r3
         vadd.f32	s11, s13, s12
         vsub.f32	s14, s15, s14
         vsub.f32	s13, s11, s13
         vmov	r0, s11
         vsub.f32	s13, s12, s13
         vsub.f32	s9, s9, s3
         vmov	r3, s13
         negs	r2, r3
         orrs	r2, r3
         eors	r3, r0
         and.w	r3, r3, r2, asr #31
         eors	r3, r0
         vadd.f32	s14, s14, s9
         vmov	s13, r3
         vadd.f32	s12, s13, s14
         and.w	r4, r4, r1, asr #31
         vsub.f32	s13, s12, s13
         add.w	ip, sp, #64	@ 0x40
         vsub.f32	s14, s14, s13
         lsrs	r1, r1, #31
         and.w	r0, r0, r2, asr #31
         movs	r3, #0
         add.w	r2, r1, r2, lsr #31
         add.w	r1, ip, r1, lsl #2
         str	r3, [sp, #48]	@ 0x30
         vsub.f32	s11, s15, s5
         vmov	r3, s14
         str	r4, [sp, #40]	@ 0x28
         str.w	r0, [r1, #-24]
         vmov	r1, s12
         vsub.f32	s15, s15, s11
         negs	r0, r3
         orrs	r0, r3
         eors	r3, r1
         vsub.f32	s10, s4, s11
         vsub.f32	s15, s5, s15
         and.w	r3, r3, r0, asr #31
         eors	r3, r1
         vadd.f32	s15, s15, s10
         vmov	s14, r3
         vadd.f32	s13, s14, s15
         and.w	r4, r1, r0, asr #31
         vsub.f32	s14, s13, s14
         add.w	ip, ip, r2, lsl #2
         vsub.f32	s15, s15, s14
         add.w	r2, r2, r0, lsr #31
         vmov	r3, s15
         vmov	r0, s13
         negs	r1, r3
         orrs	r1, r3
         eors	r3, r0
         vsub.f32	s11, s5, s7
         and.w	r3, r3, r1, asr #31
         eors	r3, r0
         vsub.f32	s11, s8, s11
         vmov	s15, r3
         vadd.f32	s14, s15, s11
         and.w	r0, r0, r1, asr #31
         vsub.f32	s15, s14, s15
         add.w	r1, r2, r1, lsr #31
         vsub.f32	s11, s11, s15
         add.w	r2, sp, r2, lsl #2
         add.w	r3, sp, r1, lsl #2
         # vmov r4, s0
        vmov lr, s0
         str.w	r4, [ip, #-24]
         str	r0, [r2, #40]	@ 0x28
         vstr	s14, [r3, #40]	@ 0x28
         vstr	s11, [r3, #44]	@ 0x2c
         vldr	s1, [sp, #44]	@ 0x2c
         vldr	s2, [sp, #48]	@ 0x30
         vldr	s0, [sp, #40]	@ 0x28
         add	sp, #64	@ 0x40
         pop {r4}
         bx lr
         # pop	{pc}

.align	2
.global tw_sum_f_ct_asm
.thumb
.thumb_func
.type	tw_sum_f_ct_asm, %function
tw_sum_f_ct_asm:
         sub	sp, #60	@ 0x3c
         vmov r2, s0
         vmov r3, s3
         vmov s0, r4
         bic.w	r1, r3, #2147483648	@ 0x80000000
         bic.w	r4, r2, #2147483648	@ 0x80000000
         subs	r4, r4, r1
         eor.w	r0, r2, r3
         vmov r1, s1
         and.w	r0, r0, r4, asr #31
         eors	r3, r0
         eors	r2, r0
         vmov	s10, r2
         bic.w	r4, r1, #2147483648	@ 0x80000000
         bic.w	r2, r3, #2147483648	@ 0x80000000
         subs	r4, r4, r2
         eor.w	r0, r3, r1
         vmov r2, s2
         and.w	r0, r0, r4, asr #31
         eors	r3, r0
         eors	r1, r0
         vmov	s12, r1
         bic.w	r0, r2, #2147483648	@ 0x80000000
         bic.w	r1, r3, #2147483648	@ 0x80000000
         subs	r0, r0, r1
         eor.w	r1, r3, r2
         and.w	r1, r1, r0, asr #31
         eors	r3, r1
         vmov	s15, r3
         eor.w	r3, r2, r1
         vmov	s11, r3
         vadd.f32	s9, s11, s15
         movs	r4, #0
         vadd.f32	s14, s9, s12
         vsub.f32	s7, s9, s15
         vadd.f32	s8, s14, s10
         vsub.f32	s6, s14, s9
         vsub.f32	s5, s8, s14
         vsub.f32	s13, s14, s6
         vsub.f32	s4, s8, s5
         vsub.f32	s10, s10, s5
         vsub.f32	s14, s14, s4
         vsub.f32	s13, s9, s13
         vadd.f32	s14, s14, s10
         vsub.f32	s10, s12, s6
         vadd.f32	s12, s8, s14
         vadd.f32	s13, s13, s10
         vsub.f32	s8, s12, s8
         vmov	r0, s12
         vsub.f32	s14, s14, s8
         vsub.f32	s9, s9, s7
         vmov	r3, s14
         rsb	ip, r3, #0
         orr.w	ip, ip, r3
         eors	r3, r0
         and.w	r3, r3, ip, asr #31
         eors	r3, r0
         vmov	s14, r3
         vadd.f32	s12, s13, s14
         vsub.f32	s15, s15, s9
         vsub.f32	s14, s12, s14
         vmov	r1, s12
         vsub.f32	s14, s13, s14
         vsub.f32	s11, s11, s7
         vmov	r3, s14
         negs	r2, r3
         orrs	r2, r3
         eors	r3, r1
         and.w	r3, r3, r2, asr #31
         eors	r3, r1
         vadd.f32	s15, s15, s11
         vmov	s14, r3
         vadd.f32	s13, s15, s14
         mov.w	r3, ip, lsr #31
         and.w	r0, r0, ip, asr #31
         vsub.f32	s14, s13, s14
         add.w	ip, sp, #56	@ 0x38
         and.w	r1, r1, r2, asr #31
         add.w	r2, r3, r2, lsr #31
         add.w	r3, ip, r3, lsl #2
         vsub.f32	s15, s15, s14
         str	r4, [sp, #48]	@ 0x30
         str	r0, [sp, #40]	@ 0x28
         str.w	r1, [r3, #-16]
         add.w	r3, ip, r2, lsl #2
         vmov r4, s0
         vstr	s13, [r3, #-16]
         vstr	s15, [r3, #-12]
         vldr	s1, [sp, #44]	@ 0x2c
         vldr	s2, [sp, #48]	@ 0x30
         vldr	s0, [sp, #40]	@ 0x28
         add	sp, #60	@ 0x3c
         bx	lr


	.align	2
	.global tw_prod_fast_ct_asm	
	.thumb
	.thumb_func
	.type	tw_prod_fast_ct_asm, %function
tw_prod_fast_ct_asm:
        vmul.f32        s13, s3, s1
        vmul.f32        s15, s4, s0
        vmul.f32        s9, s3, s0
        vadd.f32        s14, s13, s15
        vmov.f32        s12, s9
        vfnms.f32       s12, s0, s3
        vsub.f32        s8, s14, s13
        vadd.f32        s11, s14, s12
        vsub.f32        s6, s14, s8
        vsub.f32        s10, s11, s14
        vsub.f32        s8, s15, s8
        vsub.f32        s7, s11, s10
        vsub.f32        s10, s12, s10
        vsub.f32        s12, s13, s6
        vsub.f32        s14, s14, s7
        vadd.f32        s12, s12, s8
        vfnms.f32       s15, s0, s4
        vfnms.f32       s13, s1, s3
        vfma.f32        s15, s2, s3
        vfma.f32        s13, s0, s5
        vadd.f32        s15, s15, s13
        vadd.f32        s13, s14, s10
        vmov.f32        s14, s12
        vfma.f32        s14, s1, s4
        vadd.f32        s15, s15, s14
        vadd.f32        s14, s15, s13
        vadd.f32        s12, s14, s11
        vsub.f32        s6, s14, s13
        vadd.f32        s0, s12, s9
        vsub.f32        s11, s12, s11
        vsub.f32        s9, s0, s9
        vsub.f32        s14, s14, s11
        vsub.f32        s12, s12, s9
        vsub.f32        s15, s15, s6
        vadd.f32        s6, s14, s12
        vsub.f32        s12, s6, s12
        vadd.f32        s11, s6, s15
        vsub.f32        s14, s14, s12
        vsub.f32        s12, s11, s6
        vmov            r1, s14 @ int
        vsub.f32        s10, s15, s12
        vadd.f32        s14, s14, s15
        rsbs            r0, r1, #0
        vmov            r2, s14
        vmov            ip, s11 @ int
        vmov            r3, s10
        orrs            r0, r0, r1
        vmov            r1, s6 @ int
        rsbs            r0, r0, #0
        eor             r1, ip, r1
        eors            r3, r3, r2
        and             r1, r1, r0, lsr #31
        and             r3, r3, r0, lsr #31
        eor             r1, r1, ip
        eors            r3, r3, r2
        vmov            s1, r1  @ int
        vmov            s2, r3  @ int
        bx              lr

	.align	2
	.global tw_prod_fast_f_ct_asm	
	.thumb
	.thumb_func
	.type	tw_prod_fast_f_ct_asm, %function
tw_prod_fast_f_ct_asm:
        vmul.f32        s15, s2, s0
        vmul.f32        s9, s1, s0
        vmov.f32        s14, s9
        vfnms.f32 s14, s0, s1
        vadd.f32        s10, s15, s14
        vsub.f32        s13, s10, s15
        vsub.f32        s12, s10, s13
        vsub.f32        s14, s14, s13
        vsub.f32        s12, s15, s12
        vfnms.f32 s15, s0, s2
        vfma.f32        s15, s3, s0
        vadd.f32        s12, s12, s14
        vadd.f32        s14, s15, s12
        vadd.f32        s13, s14, s10
        vsub.f32        s11, s14, s12
        vadd.f32        s0, s13, s9
        vsub.f32        s10, s13, s10
        vsub.f32        s9, s0, s9
        vsub.f32        s14, s14, s10
        vsub.f32        s12, s13, s9
        vsub.f32        s15, s15, s11
        vadd.f32        s13, s14, s12
        vsub.f32        s12, s13, s12
        vadd.f32        s11, s13, s15
        vsub.f32        s14, s14, s12
        vsub.f32        s12, s11, s13
        vmov    r1, s14 @ int
        vadd.f32        s14, s14, s15
        vsub.f32        s15, s15, s12
        rsbs    r0, r1, #0
        vmov    r2, s14
        vmov    ip, s11 @ int
        vmov    r3, s15
        orrs    r0, r0, r1
        vmov    r1, s13 @ int
        rsbs    r0, r0, #0
        eor     r1, ip, r1
        eors    r3, r3, r2
        and     r1, r1, r0, lsr #31
        and     r3, r3, r0, lsr #31
        eor     r1, r1, ip
        eors    r3, r3, r2
        vmov    s1, r1  @ int
        vmov    s2, r3  @ int
        bx      lr


	.align	2
	.global tw_reci_fast_ct_asm	
	.thumb
	.thumb_func
	.type	tw_reci_fast_ct_asm, %function
tw_reci_fast_ct_asm:
        vldr.32 s15, =1065353217
        vldr.32 s11, =1065353214
        vdiv.f32        s12, s15, s0
        vldr.32 s15, =-1082130431
        vfma.f32 s15, s12, s0
        vmul.f32        s14, s12, s11
        vneg.f32        s10, s12
        vmov.f32        s13, s14
        vneg.f32        s15, s15
        vfnms.f32 s13, s12, s11
        vfma.f32 s15, s10, s1
        vfma.f32 s13, s12, s15
        vadd.f32        s12, s13, s14
        vsub.f32        s14, s12, s14
        vmul.f32        s15, s1, s12
        vsub.f32        s13, s13, s14
        vmul.f32        s6, s12, s0
        vmul.f32        s11, s13, s0
        vmov.f32        s8, s6
        vadd.f32        s14, s15, s11
        vfnms.f32       s8, s12, s0
        vsub.f32        s7, s14, s11
        vadd.f32        s9, s8, s14
        vsub.f32        s10, s14, s7
        vsub.f32        s5, s9, s14
        vsub.f32        s7, s15, s7
        vsub.f32        s10, s11, s10
        vsub.f32        s4, s9, s5
        vadd.f32        s10, s10, s7
        vfnms.f32       s15, s12, s1
        vfnms.f32       s11, s13, s0
        vfma.f32        s11, s12, s2
        vsub.f32        s14, s14, s4
        vadd.f32        s15, s15, s11
        vmov.f32        s11, s10
        vfma.f32        s11, s13, s1
        vsub.f32        s8, s8, s5
        vadd.f32        s15, s15, s11
        vadd.f32        s11, s14, s8
        vneg.f32        s15, s15
        vsub.f32        s14, s15, s11
        vsub.f32        s10, s14, s9
        vadd.f32        s11, s11, s14
        vsub.f32        s8, s10, s6
        vadd.f32        s9, s9, s10
        vadd.f32        s8, s8, s6
        vsub.f32        s14, s14, s9
        vsub.f32        s10, s10, s8
        vsub.f32        s15, s15, s11 @
        vadd.f32        s8, s14, s10
        vsub.f32        s10, s8, s10
        vadd.f32        s9, s15, s8
        vsub.f32        s11, s14, s10
        vmov    r2, s9  @ int
        vmov    r3, s11 @ int
        rsbs    r1, r3, #0
        orrs    r1, r1, r3
        vmov    r3, s8  @ int
        rsbs    r1, r1, #0
        eors    r3, r3, r2
        and     r3, r3, r1, lsr #31
        eors    r3, r3, r2
        vmov    s10, r3 @ int
        vsub.f32        s9, s9, s8
        vmul.f32        s14, s12, s10
        vadd.f32        s11, s15, s11
        vsub.f32        s15, s15, s9
        vmov    r2, s11
        vadd.f32        s11, s14, s13
        vmov    r3, s15
        vsub.f32        s15, s11, s13
        eors    r3, r3, r2
        vsub.f32        s2, s14, s15
        and     r3, r3, r1, lsr #31
        eors    r3, r3, r2
        vmov    s7, r3
        vfnms.f32 s14, s12, s10
        vfma.f32 s14, s13, s10
        vfma.f32 s14, s12, s7
        vadd.f32        s15, s2, s14 
        vadd.f32        s2, s11, s15
        vadd.f32        s0, s2, s12
        vsub.f32        s11, s2, s11
        vsub.f32        s12, s0, s12
        vsub.f32        s15, s15, s11
        vsub.f32        s2, s2, s12
        vadd.f32        s1, s15, s2
        vsub.f32        s2, s1, s2
        vsub.f32        s2, s15, s2
        bx      lr

	.align	2
	.global tw_div_fast_ct_asm	
	.thumb
	.thumb_func
	.type	tw_div_fast_ct_asm, %function
tw_div_fast_ct_asm:
        vldr.32 s15, =1065353217
        vldr.32 s12, =1065353214
        vdiv.f32        s13, s15, s3
        nop
        push    {lr}
        vfnms.f32 s15, s13, s3
        vneg.f32        s15, s15
        vmul.f32        s14, s13, s12
        vfms.f32 s15, s13, s4
        vmov.f32        s10, s14
        vfnms.f32 s10, s13, s12
        vfma.f32 s10, s13, s15
        vmov r0, s16
        vadd.f32        s9, s10, s14
        vsub.f32        s14, s9, s14
        vmul.f32        s15, s4, s9
        vsub.f32        s10, s10, s14
        vmul.f32        s16, s9, s3
        vmul.f32        s11, s10, s3
        vmov.f32        s8, s16
        vadd.f32        s14, s15, s11
        vfnms.f32 s8, s9, s3
        vsub.f32        s13, s14, s11
        vadd.f32        s7, s8, s14
        vsub.f32        s12, s14, s13
        vsub.f32        s6, s15, s13
        vsub.f32        s12, s11, s12
        vfnms.f32 s11, s10, s3
        vadd.f32        s12, s12, s6
        vsub.f32        s6, s7, s14
        vfma.f32 s11, s9, s5
        vsub.f32        s5, s7, s6
        vfnms.f32 s15, s9, s4
        vsub.f32        s8, s8, s6
        vadd.f32        s13, s15, s11
        vsub.f32        s14, s14, s5
        vmov.f32        s15, s12
        vfma.f32 s15, s10, s4
        vadd.f32        s14, s14, s8
        vadd.f32        s13, s13, s15
        vmul.f32        s4, s9, s1
        vneg.f32        s13, s13
        vsub.f32        s8, s13, s14
        vmul.f32        s12, s10, s0
        vsub.f32        s11, s8, s7
        vadd.f32        s14, s14, s8
        vsub.f32        s3, s11, s16
        vadd.f32        s6, s7, s11
        vadd.f32        s3, s3, s16
        vsub.f32        s8, s8, s6
        vsub.f32        s11, s11, s3
        vmul.f32        s7, s9, s0
        vadd.f32        s16, s8, s11
        vsub.f32        s13, s13, s14
        vsub.f32        s11, s16, s11
        vadd.f32        s14, s4, s12
        vmov.f32        s5, s7
        vsub.f32        s6, s8, s11
        vfnms.f32 s5, s9, s0
        vmov.f32        s15, s4
        vadd.f32        s8, s5, s14
        vmov.f32        s11, s12
        vfnms.f32 s15, s9, s1
        vfma.f32 s11, s9, s2
        vsub.f32        s2, s8, s14
        vfnms.f32 s11, s10, s0
        vsub.f32        s0, s8, s2
        vadd.f32        s15, s15, s11
        vsub.f32        s11, s14, s12
        vadd.f32        s3, s13, s16
        vsub.f32        s4, s4, s11
        vsub.f32        s11, s14, s11
        vsub.f32        s14, s14, s0
        vsub.f32        s12, s12, s11
        vsub.f32        s11, s5, s2
        vmov    r2, s6  @ int
        rsbs    r3, r2, #0
        vmov    r1, s3  @ int
        vadd.f32        s11, s14, s11
        vsub.f32        s3, s3, s16
        vadd.f32        s12, s12, s4
        orrs    r3, r3, r2
        vfma.f32 s12, s10, s1
        vmov    r2, s16 @ int
        vmov s16, r0
        vadd.f32        s15, s15, s12
        vadd.f32        s10, s13, s6
        vadd.f32        s12, s11, s15
        vsub.f32        s13, s13, s3
        rsbs    r3, r3, #0
        eors    r2, r2, r1
        and     r2, r2, r3, lsr #31
        eors    r2, r2, r1
        vmov    r1, s13
        vadd.f32        s13, s8, s12
        vmov    lr, s10
        vadd.f32        s5, s7, s13
        vsub.f32        s8, s13, s8
        vsub.f32        s7, s5, s7
        vsub.f32        s11, s12, s11
        vsub.f32        s10, s13, s7
        vsub.f32        s12, s12, s8
        vsub.f32        s15, s15, s11
        vadd.f32        s13, s12, s10
        vmov    s9, r2  @ int
        vsub.f32        s10, s13, s10
        vadd.f32        s11, s15, s13
        vsub.f32        s12, s12, s10
        vmov    r0, s11 @ int
        vsub.f32        s11, s11, s13
        vmov    r2, s12 @ int
        vadd.f32        s12, s15, s12
        vsub.f32        s15, s15, s11
        rsb     ip, r2, #0
        eor     r1, lr, r1
        orr     ip, ip, r2
        and     r1, r1, r3, lsr #31
        vmov    r2, s12
        vmov    r3, s15
        rsb     ip, ip, #0
        eors    r3, r3, r2
        and     r3, r3, ip, lsr #31
        eors    r3, r3, r2
        vmov    r2, s13 @ int
        eors    r2, r2, r0
        and     r2, r2, ip, lsr #31
        eors    r2, r2, r0
        vmul.f32        s14, s5, s9
        vmov    s13, r2 @ int
        vmov.f32        s15, s14
        vmov.f32        s12, s13
        vfnms.f32 s15, s5, s9
        vadd.f32        s8, s14, s13
        vfma.f32 s15, s13, s9
        eor     r1, r1, lr
        vsub.f32        s12, s8, s12
        vmov    s11, r1
        vsub.f32        s14, s14, s12
        vfma.f32 s15, s5, s11
        vadd.f32        s14, s14, s15
        vmov    s15, r3
        vadd.f32        s14, s14, s15
        vadd.f32        s2, s8, s14
        vadd.f32        s0, s2, s5
        vsub.f32        s8, s2, s8
        vsub.f32        s15, s0, s5
        vsub.f32        s7, s14, s8
        vsub.f32        s2, s2, s15
        vadd.f32        s1, s7, s2
        vsub.f32        s2, s1, s2
        vsub.f32        s2, s7, s2
        ldr     pc, [sp], #4
