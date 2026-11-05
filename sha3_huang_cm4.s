	.syntax	unified
	.cpu	cortex-m4
	.file	"sha3_huang_cm4.s"
	.text

	@ Credit: Henry S. Warren, Hacker's Delight, Addison-Wesley, 2002
	@ Inverse of the bit-interleaving of pqm4/common/keccakf1600.S.
	.macro	fromBitInterleaving x0, x1, t
	movs	\t, \x0
	bfi	\x0, \x1, #16, #16
	bfc	\x1, #0, #16
	orr	\x1, \x1, \t, LSR #16

	eor	\t, \x0, \x0, LSR #8
	and	\t, #0x0000FF00
	eors	\x0, \x0, \t
	eor	\x0, \x0, \t, LSL #8

	eor	\t, \x0, \x0, LSR #4
	and	\t, #0x00F000F0
	eors	\x0, \x0, \t
	eor	\x0, \x0, \t, LSL #4

	eor	\t, \x0, \x0, LSR #2
	and	\t, #0x0C0C0C0C
	eors	\x0, \x0, \t
	eor	\x0, \x0, \t, LSL #2

	eor	\t, \x0, \x0, LSR #1
	and	\t, #0x22222222
	eors	\x0, \x0, \t
	eor	\x0, \x0, \t, LSL #1

	eor	\t, \x1, \x1, LSR #8
	and	\t, #0x0000FF00
	eors	\x1, \x1, \t
	eor	\x1, \x1, \t, LSL #8

	eor	\t, \x1, \x1, LSR #4
	and	\t, #0x00F000F0
	eors	\x1, \x1, \t
	eor	\x1, \x1, \t, LSL #4

	eor	\t, \x1, \x1, LSR #2
	and	\t, #0x0C0C0C0C
	eors	\x1, \x1, \t
	eor	\x1, \x1, \t, LSL #2

	eor	\t, \x1, \x1, LSR #1
	and	\t, #0x22222222
	eors	\x1, \x1, \t
	eor	\x1, \x1, \t, LSL #1
	.endm

@ Glue for FNDSA_SHAKE_HUANG_BUF: the permutation itself is pqm4's
@ KeccakF1600_StatePermute, on a bit-interleaved state.

@ void fndsa_sha3_huang_process_block(uint64_t *A, unsigned r)
	.align	2
	.global	fndsa_sha3_huang_process_block
	.thumb
	.thumb_func
	.type	fndsa_sha3_huang_process_block, %function
fndsa_sha3_huang_process_block:
	b.w	KeccakF1600_StatePermute
	.size	fndsa_sha3_huang_process_block,.-fndsa_sha3_huang_process_block

@ void fndsa_sha3_extract_rate136(uint64_t *A, uint8_t *buf)
@ Decode the 136-byte SHAKE256 rate from the bit-interleaved state.
	.align	2
	.global	fndsa_sha3_extract_rate136
	.thumb
	.thumb_func
	.type	fndsa_sha3_extract_rate136, %function
fndsa_sha3_extract_rate136:
	push	{r4, r5}
	movs	r3, #17
1:
	ldrd	r4, r5, [r0], #8
	fromBitInterleaving r4, r5, r2
	str	r4, [r1], #4
	subs	r3, r3, #1
	str	r5, [r1], #4
	bne	1b
	pop	{r4, r5}
	bx	lr
	.size	fndsa_sha3_extract_rate136,.-fndsa_sha3_extract_rate136
