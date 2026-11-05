	.syntax	unified
	.cpu	cortex-m4
	.file	"sign_sampler_cm4.s"
	.text

@ =======================================================================
@ void fndsa_ffsamp_fft_inner(sampler_state *ss, unsigned logn, fpr *tmp)
@ =======================================================================

	.align	2
	.global	fndsa_ffsamp_fft_inner
	.thumb
	.thumb_func
	.type	fndsa_ffsamp_fft_inner, %function
fndsa_ffsamp_fft_inner:
	@ We mimic the C implementation (sign_sampler.c); this function
	@ is in assembly so that we can optimize its stack allocation
	@ (GCC uses about 80 bytes per call level, which is a lot, since
	@ this is the main recursive function).

	@ If logn == 1, we tail-call fndsa_ffsamp_fft_deepest
	cmp	r1, #1
	bne	fndsa_ffsamp_fft_inner__L1
	mov.w	r1, r2
	b.w	fndsa_ffsamp_fft_deepest

fndsa_ffsamp_fft_inner__L1:
	push	{ r4, r5, r6, lr }

	@ r4 <- ss
	@ r5 <- logn
	@ r6 <- tmp
	movs	r4, r0
	movs	r5, r1
	movs	r6, r2

	@ Write into rd the address tmp + off*(n/4), expressed in 8-byte
	@ chunk (hence, tmp + 2*off*n, when counting in bytes).
.macro QC  rd, off
	.if ((\off) == 0)
	mov.w	\rd, r6
	.else
	movs	\rd, #(2 * (\off))
	lsls	\rd, r5
	add.w	\rd, \rd, r6
	.endif
.endm

	@ Decompose G into LDL; the decomposed matrix replaces G. The
	@ fifth argument (l10, aliasing g01) goes on the stack.
	sub	sp, sp, #8
	QC	r0, 8
	str	r0, [sp]
	mov.w	r0, r5
	QC	r1, 12
	QC	r2, 8
	QC	r3, 14
	bl	fndsa_fpoly_LDL_fft
	add	sp, sp, #8

	@ Split d11 into the right sub-tree (right_00, right_01);
	@ right_11 is a copy of right_00.
	mov.w	r0, r5
	QC	r1, 20
	QC	r2, 18
	QC	r3, 14
	bl	fndsa_fpoly_split_selfadj_fft
	QC	r0, 21
	QC	r1, 20
	movs	r2, #2
	lsls	r2, r5
	bl	memcpy

	@ Split t1 and make the first recursive call on the two
	@ halves, using the right sub-tree, then merge the result
	@ into 18..21
	mov.w	r0, r5
	QC	r1, 14
	QC	r2, 16
	QC	r3, 4
	bl	fndsa_fpoly_split_fft
	movs	r0, r4
	subs	r1, r5, #1
	QC	r2, 14
	bl	fndsa_ffsamp_fft_inner
	mov.w	r0, r5
	QC	r1, 18
	QC	r2, 14
	QC	r3, 16
	bl	fndsa_fpoly_merge_fft

	@ Compute tb0 = t0 + (t1 - z1)*l10 (into t0) and move z1 into t1.
	QC	r0, 14
	QC	r1, 4
	movs	r2, #8
	lsls	r2, r5
	bl	memcpy
	mov.w	r0, r5
	QC	r1, 14
	QC	r2, 18
	bl	fndsa_fpoly_sub
	QC	r0, 4
	QC	r1, 18
	movs	r2, #8
	lsls	r2, r5
	bl	memcpy
	mov.w	r0, r5
	QC	r1, 14
	QC	r2, 8
	bl	fndsa_fpoly_mul_fft
	mov.w	r0, r5
	QC	r1, 0
	QC	r2, 14
	bl	fndsa_fpoly_add

	@ Split d00 to obtain the left-subtree.
	mov.w	r0, r5
	QC	r1, 20
	QC	r2, 18
	QC	r3, 12
	bl	fndsa_fpoly_split_selfadj_fft
	QC	r0, 21
	QC	r1, 20
	movs	r2, #2
	lsls	r2, r5
	bl	memcpy

	@ Split tb0 and perform the second recursive call on the
	@ split output; the final merge produces z0, which we write
	@ into t0.
	mov.w	r0, r5
	QC	r1, 14
	QC	r2, 16
	QC	r3, 0
	bl	fndsa_fpoly_split_fft
	movs	r0, r4
	subs	r1, r5, #1
	QC	r2, 14
	bl	fndsa_ffsamp_fft_inner
	mov.w	r0, r5
	QC	r1, 0
	QC	r2, 14
	QC	r3, 16
	bl	fndsa_fpoly_merge_fft

	pop	{ r4, r5, r6, pc }
	.size	fndsa_ffsamp_fft_inner,.-fndsa_ffsamp_fft_inner
