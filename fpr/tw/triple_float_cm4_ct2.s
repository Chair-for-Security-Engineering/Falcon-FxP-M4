@ Copyright 2025 NXP
@ SPDX-License-Identifier: MIT

.syntax	unified
.cpu	cortex-m4
.file	"triple_float_cm4_ct2.s"
.text

.align	2
.global tw_add_sub_ct_asm 
.thumb
.thumb_func
.type	tw_add_sub_ct_asm, %function
tw_add_sub_ct_asm:
	vmov	r5, s3
	sub	sp, #52	@ 0x34
	nop
	vmov	r6, s4
	vmov	r0, s5
	vmov r7, s0
	eor.w	r8, r5, #2147483648	@ 0x80000000
	bic.w	r3, r7, #2147483648	@ 0x80000000
	bic.w	r1, r5, #2147483648	@ 0x80000000
	subs.w	r3, r3, r1
	sbc.w	r1, r1, r1
	uadd8	r3, r1, r1
	sel	r3, r5, r7
	sel	r1, r7, r5
	sel	r2, r8, r7
	sel	r4, r7, r8
	vmov	s7, r3
	vmov	s0, r2
	vmov r5, s2
	eor.w	r2, r0, #2147483648	@ 0x80000000
	bic.w	r3, r5, #2147483648	@ 0x80000000
	bic.w	r7, r0, #2147483648	@ 0x80000000
	subs.w	r3, r3, r7
	sbc.w	r7, r7, r7
	uadd8	r3, r7, r7
	sel	r3, r0, r5
	sel	r7, r5, r0
	sel	r8, r5, r2
	sel	r2, r2, r5
	vmov	s10, r7
	vmov	s9, r8
	eor.w	sl, r6, #2147483648	@ 0x80000000
	bic.w	r5, r3, #2147483648	@ 0x80000000
	bic.w	r7, r1, #2147483648	@ 0x80000000
	subs.w	r5, r5, r7
	sbc.w	r7, r7, r7
	uadd8	r5, r7, r7
	sel	r5, r1, r3
	sel	r1, r3, r1
	sel	r0, r4, r2
	sel	r4, r2, r4
	vmov r7, s1
	bic.w	r3, r7, #2147483648	@ 0x80000000
	bic.w	r8, r6, #2147483648	@ 0x80000000
	subs.w	r3, r3, r8
	sbc.w	r8, r8, r8
	uadd8	r3, r8, r8
	sel	r3, r6, r7
	sel	r8, r7, r6
	sel	r2, sl, r7
	sel	r6, r7, sl 
	bic.w	r7, r1, #2147483648	@ 0x80000000
	bic.w	sl, r8, #2147483648	@ 0x80000000
	subs.w	r7, r7, sl
	sbc.w	sl, sl, sl
	uadd8	r7, sl, sl
	sel	r7, r8, r1
	vmov	s12, r7
	sel	r1, r1, r8
	sel	r7, r6, r4
	sel	r4, r4, r6 
	vmov	s11, r1
	vadd.f32	s8, s11, s10
	vmov	s5, r4
	vadd.f32	s13, s8, s12
	vsub.f32	s2, s8, s10
	vmov	s3, r7
	bic.w	r7, r3, #2147483648	@ 0x80000000
	bic.w	r1, r5, #2147483648	@ 0x80000000
	subs.w	r7, r7, r1
	sbc.w	r1, r1, r1
	uadd8	r7, r1, r1
	sel	r7, r5, r3
	sel	r1, r3, r5
	sel	r6, r0, r2
	sel	r4, r2, r0
	vsub.f32	s1, s11, s2
	vmov	s15, r1
	vsub.f32	s11, s13, s8
	vsub.f32	s4, s8, s2
	vsub.f32	s12, s12, s11
	vsub.f32	s2, s13, s11
	vadd.f32	s11, s13, s15
	vsub.f32	s10, s10, s4
	vmov	s4, r7
	vsub.f32	s8, s8, s2
	vsub.f32	s2, s11, s13
	vadd.f32	s8, s8, s12
	vadd.f32	s10, s10, s1
	vadd.f32	s12, s11, s4
	vsub.f32	s1, s11, s2
	vsub.f32	s15, s15, s2
	vsub.f32	s13, s13, s1
	vsub.f32	s2, s12, s11
	vadd.f32	s1, s12, s7
	vadd.f32	s13, s13, s15
	vsub.f32	s6, s12, s2
	vsub.f32	s15, s1, s12
	vsub.f32	s11, s11, s6
	vsub.f32	s4, s4, s2
	vsub.f32	s7, s7, s15
	vsub.f32	s2, s1, s15
	vadd.f32	s15, s5, s9
	vadd.f32	s4, s11, s4
	vsub.f32	s12, s12, s2
	vadd.f32	s11, s15, s3
	vsub.f32	s2, s15, s9
	vadd.f32	s7, s12, s7
	vmov	s14, r4
	vsub.f32	s12, s15, s2
	vsub.f32	s5, s5, s2
	vsub.f32	s2, s11, s15
	vsub.f32	s9, s9, s12
	vsub.f32	s6, s11, s2
	vadd.f32	s12, s11, s14
	vsub.f32	s15, s15, s6
	vmov	s6, r6
	vadd.f32	s9, s9, s5
	vsub.f32	s5, s12, s11
	vsub.f32	s3, s3, s2
	vsub.f32	s14, s14, s5
	vadd.f32	s3, s15, s3
	vsub.f32	s5, s12, s5
	vadd.f32	s15, s12, s6
	vsub.f32	s11, s11, s5
	vsub.f32	s2, s15, s12
	vadd.f32	s5, s11, s14
	vsub.f32	s14, s15, s2
	vsub.f32	s6, s6, s2
	vadd.f32	s2, s15, s0
	vsub.f32	s12, s12, s14
	vsub.f32	s14, s2, s15
	vadd.f32	s12, s12, s6
	vsub.f32	s11, s2, s14
	vsub.f32	s0, s0, s14
	vsub.f32	s15, s15, s11
	vadd.f32	s11, s1, s7
	vadd.f32	s15, s15, s0
	vsub.f32	s14, s11, s1
	vmov	r7, s11
	vsub.f32	s14, s7, s14
	movs	r4, #0
	vmov	r2, s14
	vadd.f32	s14, s2, s15
	negs	r6, r2
	vsub.f32	s2, s14, s2
	orrs	r6, r2
	eors	r2, r7
	vsub.f32	s15, s15, s2
	and.w	r2, r2, r6, asr #31
	eors	r2, r7
	vmov	r3, s15
	vmov	s15, r2
	vmov	r5, s14
	vadd.f32	s11, s15, s4
	rsb	r1, r3, #0
	orr.w	r1, r1, r3
	eors	r3, r5
	vsub.f32	s14, s11, s15
	and.w	r3, r3, r1, asr #31
	eors	r3, r5
	vmov	s15, r3
	vsub.f32	s14, s4, s14
	vadd.f32	s7, s15, s12
	vmov	r3, s14
	and.w	r7, r7, r6, asr #31
	str	r7, [sp, #24]
	vmov	r7, s11
	vsub.f32	s15, s7, s15
	rsb	r0, r3, #0
	orr.w	r0, r0, r3
	eors	r3, r7
	vsub.f32	s15, s12, s15
	lsrs	r6, r6, #31
	and.w	r3, r3, r0, asr #31
	eors	r3, r7
	and.w	r7, r7, r0, asr #31
	add.w	r0, r6, r0, lsr #31
	add.w	r6, sp, r6, lsl #2
	vmov	r2, s15
	vmov	s15, r3
	str	r4, [sp, #32]
	str	r7, [r6, #24]
	vmov	r6, s7
	vadd.f32	s11, s15, s13
	negs	r7, r2
	orrs	r7, r2
	eor.w	r3, r2, r6
	vsub.f32	s15, s11, s15
	and.w	r3, r3, r7, asr #31
	eors	r3, r6
	vmov	s12, r3
	vsub.f32	s15, s13, s15
	and.w	r5, r5, r1, asr #31
	add	r3, sp, #48	@ 0x30
	mov.w	r1, r1, lsr #31
	vadd.f32	s14, s12, s5
	and.w	r6, r6, r7, asr #31
	add.w	r7, r1, r7, lsr #31
	add.w	r11, r3, r1, lsl #2
	vmov	r3, s15
	vmov	r8, s11
	vsub.f32	s12, s14, s12
	negs	r2, r3
	orrs	r2, r3
	eor.w	r3, r3, r8
	vsub.f32	s15, s5, s12
	and.w	r3, r3, r2, asr #31
	eor.w	r3, r3, r8
	vmov	sl, s15
	vmov	s15, r3
	and.w	r3, r8, r2, asr #31
	add.w	r2, r0, r2, lsr #31
	add.w	r0, sp, r0, lsl #2
	vmov	r1, s14 @ start r1
	str.w	r3, [r0, #24]
	vadd.f32	s13, s15, s8
	rsb	r0, sl, #0
	orr.w	r0, r0, sl
	eor.w	r3, sl, r1
	vsub.f32	s15, s13, s15
	and.w	r3, r3, r0, asr #31
	eor.w	r3, r3, r1
	vmov	s14, r3
	vsub.f32	s15, s8, s15
	vadd.f32	s12, s14, s3
	vmov	r3, s15
	vmov	sl, s13
	vsub.f32	s14, s12, s14
	rsb	r8, r3, #0
	orr.w	r8, r8, r3
	eor.w	r3, r3, sl
	vsub.f32	s15, s3, s14
	and.w	r3, r3, r8, asr #31
	eor.w	r3, r3, sl
	vmov	r9, s15
	vmov	s15, r3
	and.w	r3, sl, r8, asr #31
	add.w	r8, r2, r8, lsr #31
	add.w	r2, sp, r2, lsl #2
	vmov	sl, s12
	str	r3, [r2, #24]
	vadd.f32	s13, s15, s10
	rsb	r2, r9, #0
	orr.w	r2, r2, r9
	eor.w	r3, r9, sl
	vsub.f32	s15, s13, s15
	and.w	r3, r3, r2, asr #31
	eor.w	r3, r3, sl
	vsub.f32	s10, s10, s15
	vmov	s15, r3
	vadd.f32	s14, s15, s9
	and.w	r1, r1, r0, asr #31 @ end r1
	vsub.f32	s15, s14, s15
	add.w	r0, r7, r0, lsr #31
	and.w	r3, sl, r2, asr #31
	vsub.f32	s9, s9, s15
	add.w	sl, sp, r8, lsl #2
	add.w	r2, r0, r2, lsr #31
	add.w	r8, sp, r8, lsl #2
	add.w	r7, sp, r7, lsl #2
	add.w	r0, sp, r0, lsl #2
	vstr	s13, [sl, #24]
	vstr	s10, [r8, #28]
	add.w	r8, sp, r2, lsl #2
	vldr.32	s0, [sp, #24]
	vldr.32	s2, [sp, #32]
	vldr.32	s1, [sp, #28]
	str	r5, [sp, #24]
	str	r4, [sp, #32]
	str.w	r6, [r11, #-24]
	str.w	r1, [r7, #24] @end r1
	str.w	r3, [r0, #24]
	vstr	s14, [r8, #24]
	vstr	s9, [r8, #28]
	vldr.32	s3, [sp, #24]
	vldr.32	s4, [sp, #28]
	vldr.32	s5, [sp, #32]
	add	sp, #52	@ 0x34
	bx lr


	.align	2
	.global tw_rint_fast_ct_asm
	.thumb
	.thumb_func
	.type	tw_rint_fast_ct_asm, %function
tw_rint_fast_ct_asm:
        vmov    r2, s0  @ int
        lsls    r3, r2, #7
        ubfx    ip, r2, #23, #8
        bic     r3, r3, #-1073741824
        vmov s14, lr
        orr     r3, r3, #1073741824
        rsb     lr, ip, #157
        sub     r1, ip, #126
        rsb     ip, ip, #125
        and     ip, r3, ip, asr #16
        orr     r1, lr, r1, lsr #16
        eor     r3, r3, r2, asr #31
        sub     r3, r3, r2, asr #31
        and     r1, r1, #31
        asrs    r3, r3, r1
        vmov    s15, r3 @ int
        vcvt.f32.s32    s15, s15
        vmov    r0, s1  @ int
        vsub.f32        s15, s0, s15
        orrs    r0, r0, r2
        vmov    r3, s15
        lsrs    r0, r0, #31
        add     r1, r3, #-1056964608
        rsb     r3, r3, #1056964608
        orrs    r1, r1, r3
        lsls    r0, r0, #1
        and     lr, lr, #31
        asrs    r3, r1, #31
        rsb     r0, r0, #1
        bic     r0, r0, r3
        rsb     r3, lr, #31
        lsl     r3, ip, r3
        lsr     ip, ip, lr
        bic     lr, r3, #-1073741824
        add     lr, lr, #1073741824
        add     lr, lr, #-1
        orr     r3, r3, lr, lsr #1
        lsr     lr, r3, #29
        movs    r3, #200
        asr     r3, r3, lr
        vcvt.s32.f32    s15, s0
        and     r3, r3, #1
        add     r3, r3, ip
        eor     r3, r3, r2, asr #31
        vmov    ip, s15 @ int
        sub     r3, r3, r2, asr #31
        eor     r3, r3, ip
        and     r3, r3, r1, asr #31
        eor     r3, r3, ip
        vmov lr, s14
        add     r0, r0, r3
        bx lr

	.align	2
	.global tw_floor_fast_asm
	.thumb
	.thumb_func
	.type	tw_floor_fast_asm, %function
tw_floor_fast_asm:
        vmrs    r1, fpscr
        bic.w   r1, r1, #12582912
        orrs    r0, r1, #8388608
        vmsr    fpscr, r0
        vadd.f32        s0, s1
        vmsr    fpscr, r1
        vmov r1, s0
        lsls    r3, r1, #7
        ubfx    r0, r1, #23, #8
        bic     r3, r3, #-1073741824
        sub     r2, r0, #126
        orr     r3, r3, #1073741824
        rsb     r0, r0, #157
        eor     r3, r3, r1, asr #31
        orr     r0, r0, r2, lsr #16
        sub     r3, r3, r1, asr #31
        and     r0, r0, #31
        asr     r0, r3, r0
        bx lr

	.align	2
	.global tw_from_int_ct_asm 
	.thumb
	.thumb_func
	.type	tw_from_int_ct_asm, %function
tw_from_int_ct_asm:
        mov r1, #0
        vmov    s0, r0  @ int
        vcvt.f32.s32    s0, s0
        vmov s1, s2, r1, r1
        bx      lr

	.align	2
	.global tw_trunc_fast_ct_asm 
	.thumb
	.thumb_func
	.type	tw_trunc_fast_ct_asm, %function
tw_trunc_fast_ct_asm:
        vmrs    r1, fpscr
        bic.w   r1, r1, #12582912
        orrs    r0, r1, #12582912
        vmsr    fpscr, r0
        vadd.f32        s0, s0, s1
        mov    r0, r1
        vmsr    fpscr, r0
        vcvt.u32.f32    s0, s0
        vmov    r0, s0 @ int
        bx lr

	.align	2
	.global tw_trunc_full_ct_asm 
	.thumb
	.thumb_func
	.type	tw_trunc_full_ct_asm, %function
tw_trunc_full_ct_asm:
        vmov s4, r4
        vmov s5, r5
        vmov s6, r6
        vmov s8, r8
        vmov s9, r9
        vmov s10, r10
        vmov s11, fp
        vmov s12, lr
        vmov r3, s0
        vmov r5, s1
        vmov r4, s2
        ubfx    r3, r3, #23, #8
        ubfx    r1, r5, #23, #8
        ubfx    r0, r5, #0, #23
        asrs    r5, r5, #31
        sub     lr, r3, r1
        vmov s3, r3
        orr     r2, r5, #1
        vmov s15, r2
        rsb     r8, lr, #30
        rsb     r2, lr, #37
        asr     r9, r2, #31
        add     r0, r0, #8388608
        and     r2, r2, r8
        rsb     r3, lr, #38
        asr     r8, r8, #31
        ubfx    fp, r4, #23, #8
        and     r10, r3, r2, asr #31
        mov     r6, r3
        sub     r1, r1, fp
        bic     fp, r3, r8
        lsl     r3, r0, r3
        bic     r3, r3, r8
        bic     ip, r8, r9
        bic     r3, r3, #-2147483648
        mov     r5, r3
        and     r3, r6, ip
        lsl     r6, r0, r6
        and     r6, r6, ip
        vmov s14, r6
        orr     fp, r10, fp
        rsb     r6, lr, #100
        sub     r10, lr, #70
        sub     ip, lr, #69
        and     r6, r6, r10
        lsr     ip, r0, ip
        and     ip, ip, r6, lsr #31
        and     r10, r10, r2
        and     ip, ip, r2, asr #31
        rsb     r2, lr, #69
        lsl     r2, r0, r2
        and     r2, r2, r10, asr #31
        bic     r2, r2, #-2147483648
        orr     fp, fp, r3
        add     ip, ip, r2
        subs    r1, r1, #24
        vmov r2, s15
        sub     r1, r1, fp
        ubfx    r3, r4, #0, #23
        add     r3, r3, #8388608
        mul     ip, r2, ip
        rsb     r6, r1, #7
        sub     r2, lr, #38
        sub     lr, lr, #7
        asr     fp, r1, #31
        lsl     r6, r3, r6
        lsr     r2, r0, r2
        lsr     r0, r0, lr
        and     r6, r6, fp
        and     r2, r2, r10, asr #31
        mvn     fp, fp
        bic     r0, r0, r8
        sub     lr, r1, #31
        lsr     r8, r3, r1
        bic     r10, r8, r9
        bic     r6, r6, r9
        add     r2, r2, r5
        and     lr, fp, lr, asr #31
        vmov r5, s14
        and     lr, lr, r10
        bic     r6, r6, #-2147483648
        add     r2, r2, r5
        add     r6, r6, lr
        and     r8, r8, r9
        asrs    r4, r4, #31
        orr     r4, r4, #1
        add     r2, r2, ip, lsr #31
        add     r6, r6, r8
        bic     ip, ip, #-2147483648
        mla     r6, r4, r6, ip
        add     ip, r1, #24
        vmov r5, s15
        lsr     r3, r3, ip
        bic     r3, r3, r9
        subs    r1, r1, #8
        mul     r2, r5, r2
        and     r3, r3, r1, asr #31
        add     r0, r0, r2, lsr #31
        add     r3, r3, r6, lsr #31
        bic     r2, r2, #-2147483648
        mla     r3, r4, r3, r2
        vmov r1, s0
        ldr     r2, =1073741696
        and     r2, r2, r1, lsl #7
        add     r2, r2, #1073741824
        mla     r2, r5, r0, r2
        lsrs    r1, r3, #31
        mla     r2, r4, r1, r2
        lsrs    r0, r2, #30
        eor     r0, r0, #1
        bic     r3, r3, #-2147483648
        lsl     r4, r3, r0
        vmov r3, s3
        lsls    r2, r2, r0
        subs    r3, r3, r0
        bic     r0, r4, #-2147483648
        movs    r5, #0
        ubfx    r6, r6, #30, #1
        adds    r0, r0, r0
        adc     r1, r5, r5
        adds    r0, r0, r6
        uxth    r3, r3
        adc     r1, r1, #0
        add     r2, r2, r4, lsr #31
        add     r2, r2, r1
        rsb     r6, r3, #189
        adds    r0, r0, r5
        sub     r5, r3, #157
        lsl     r5, r2, r5
        rsb     r4, r3, #157
        lsrs    r0, r0, r6
        lsr     r4, r2, r4
        sub     r1, r3, #127
        orrs    r0, r0, r5
        asrs    r1, r1, #31
        orrs    r0, r0, r4
        lsrs    r2, r2, r6
        bic     r0, r0, r1
        bic     r1, r2, r1
        vmov r4, s4
        vmov r5, s5
        vmov r6, s6
        vmov r8, s8
        vmov r9, s9
        vmov r10, s10
        vmov lr, s12
        vmov fp, s11
        bx lr


	.align	2
	.global tw_sqrt_fast_ct_asm
	.thumb
	.thumb_func
	.type	tw_sqrt_fast_ct_asm, %function
tw_sqrt_fast_ct_asm:
        vsqrt.f32       s15, s0
        vldr.32		s14, =1065353218
        vdiv.f32        s13, s14, s15
        vmov.f32        s7, #5.0e-1
        vmul.f32        s9, s13, s7
        vmul.f32        s15, s13, s0
        vmov.f32        s10, #1.5e+0
        vmul.f32        s14, s9, s15
        vmov.f32        s12, s15
        vsub.f32        s10, s10, s14
        vfnms.f32 s12, s13, s0
        vmul.f32        s11, s13, s10
        vfma.f32 s12, s13, s1
        vfnms.f32 s14, s9, s15
        vmov.f32        s15, s11
        vfma.f32 s14, s9, s12
        vfnms.f32 s15, s13, s10
        vneg.f32        s14, s14
        vfma.f32 s15, s13, s14
        vadd.f32        s10, s15, s11
        vsub.f32        s11, s10, s11
        vmul.f32        s14, s1, s10
        vsub.f32        s15, s15, s11
        vmul.f32        s5, s10, s0
        vmul.f32        s12, s15, s0
        vmov.f32        s8, s5
        vadd.f32        s13, s14, s12
        vfnms.f32 s8, s10, s0
        vsub.f32        s4, s13, s12
        vadd.f32        s9, s8, s13
        vsub.f32        s11, s13, s4
        vsub.f32        s3, s9, s13
        vsub.f32        s4, s14, s4
        vsub.f32        s6, s9, s3
        vsub.f32        s11, s12, s11
        vsub.f32        s13, s13, s6
        vsub.f32        s8, s8, s3
        vadd.f32        s11, s11, s4
        vfnms.f32 s14, s10, s1
	vfma.f32 s12, s10, s2
	vadd.f32        s2, s13, s8
        vfnms.f32 s12, s15, s0
	vmov.f32        s13, s11
        vadd.f32        s14, s14, s12
        vfma.f32 s13, s15, s1
	vmul.f32        s6, s10, s7
        vadd.f32        s14, s14, s13
	vmul.f32        s13, s15, s7
        vadd.f32        s11, s2, s14
        vadd.f32        s10, s9, s11
        vsub.f32        s15, s11, s2
        vadd.f32        s12, s5, s10
        vsub.f32        s9, s10, s9
        vsub.f32        s5, s12, s5
        vsub.f32        s11, s11, s9
        vsub.f32        s10, s10, s5
        vsub.f32        s14, s14, s15
        vadd.f32        s4, s11, s10
        vmul.f32        s7, s6, s12
        vsub.f32        s10, s4, s10
        vadd.f32        s5, s14, s4
        vsub.f32        s8, s11, s10
        vmov    r2, s5  @ int
        vmov    r3, s8  @ int
        rsbs    r1, r3, #0
        orrs    r1, r1, r3
        vmov    r3, s4  @ int
        rsbs    r1, r1, #0
        eors    r3, r3, r2
        and     r3, r3, r1, lsr #31
        eors    r3, r3, r2
        vmov    s9, r3  @ int
        vmul.f32        s10, s13, s12
        vmul.f32        s15, s6, s9
        vsub.f32        s1, s5, s4
        vadd.f32        s11, s15, s10
        vadd.f32        s8, s14, s8
        vsub.f32        s14, s14, s1
        vsub.f32        s3, s11, s10
        vmov.f32        s4, s7
        vfnms.f32 s4, s6, s12
	vmov    r2, s8
        vadd.f32        s5, s4, s11
        vmov    r3, s14
        vsub.f32        s1, s5, s11
        vsub.f32        s8, s11, s3
        eors    r3, r3, r2
        vsub.f32        s14, s5, s1
        vsub.f32        s3, s15, s3
        vsub.f32        s8, s10, s8
        and     r3, r3, r1, lsr #31
        eors    r3, r3, r2
        vadd.f32        s8, s8, s3
        vfnms.f32 s10, s13, s12
        vsub.f32        s11, s11, s14
        vfnms.f32 s15, s6, s9
        vmov    s14, r3 @ int
        vfma.f32 s10, s6, s14
	vsub.f32        s4, s4, s1
        vadd.f32        s15, s15, s10
        vmov.f32        s10, s8
        vfnma.f32 s10, s13, s9
	vadd.f32        s11, s11, s4
        vsub.f32        s15, s10, s15
        vsub.f32        s13, s15, s11
        vsub.f32        s10, s13, s5
        vadd.f32        s11, s11, s13
        vsub.f32        s8, s10, s7
        vadd.f32        s5, s5, s10
        vadd.f32        s8, s8, s7
        vsub.f32        s13, s13, s5
        vsub.f32        s10, s10, s8
        vsub.f32        s15, s15, s11
        vadd.f32        s7, s13, s10
        vsub.f32        s10, s7, s10
        vadd.f32        s8, s15, s7
        vsub.f32        s11, s13, s10
        vmov    r2, s8  @ int
        vmov    r3, s11 @ int
        rsbs    r1, r3, #0
        orrs    r1, r1, r3
        vmov    r3, s7  @ int
        rsbs    r1, r1, #0
        eors    r3, r3, r2
        and     r3, r3, r1, lsr #31
        eors    r3, r3, r2
        vmov    s10, r3 @ int
        vsub.f32        s8, s8, s7
        vmul.f32        s13, s10, s12
        vadd.f32        s11, s15, s11
        vsub.f32        s15, s15, s8
        vmov    r2, s11
        vadd.f32        s11, s13, s9
        vmov    r3, s15
        vsub.f32        s15, s11, s9
        eors    r3, r3, r2
        vsub.f32        s2, s13, s15
        and     r3, r3, r1, lsr #31
        eors    r3, r3, r2
        vmov.f32        s15, s13
        vmov    s13, r3
        vfnms.f32 s15, s12, s10
        vfma.f32 s15, s9, s10
        vfma.f32 s15, s12, s13
        vadd.f32        s15, s2, s15
        vadd.f32        s15, s15, s14
        vadd.f32        s2, s11, s15
        vadd.f32        s0, s2, s12
        vsub.f32        s11, s2, s11
        vsub.f32        s12, s0, s12
        vsub.f32        s15, s15, s11
        vsub.f32        s2, s2, s12
        vadd.f32        s1, s15, s2
        vsub.f32        s2, s1, s2
        vsub.f32        s2, s15, s2
	bx lr
