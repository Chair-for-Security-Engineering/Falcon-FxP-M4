@ =======================================================================
@ Cycle-count benchmarks for the fpr contract.
@
@ Each BENCH function calls one target with the caller's arguments and
@ returns the DWT cycle count of that call, e.g.
@     extern uint32_t bench_fpoly_FFT(unsigned logn, fpr *f);
@ Only contract symbols are named, so the same file benchmarks every
@ backend. Inline ops go through the wrappers in fpr_bench_wrap.c.
@ Adapted from the BENCH macro in c-fn-dsa's bench_cm4/bench_cm4.s.
@ =======================================================================
#include "knobs.h"

	.syntax	unified
	.cpu	cortex-m4
	.file	"fpr_bench_cm4.S"
	.text

	@ DWT_CYCCNT; the driver must enable it first.
	.equ	SYSCNT_ADDR, 0xE0001004

@ =======================================================================
@ BENCH name, target
@
@ Times one call to 'target'. Four registers are pushed to keep the stack
@ 8-byte aligned, since some targets are C functions.
@ =======================================================================
.macro	BENCH	name, target
	.align	2
	.global	\name
	.thumb
	.thumb_func
	.type	\name, %function
\name:
	push.w	{ r4, r10, r11, lr }
	movw	r11, #(SYSCNT_ADDR & 0xFFFF)
	movt	r11, #(SYSCNT_ADDR >> 16)
	ldr.w	r10, [r11]
	bl	\target
	ldr.w	r0, [r11]
	sub.w	r0, r0, r10
	pop	{ r4, r10, r11, pc }
	.size	\name,.-\name
.endm

@ Empty function for calibration. The fndsa_bench_nop* aliases have
@ different C prototypes; the baseline wrappers call them where a real
@ wrapper calls the operation, so subtracting leaves only the op body.
	.align	2
	.thumb
	.thumb_func
	.type	do_nothing, %function
	.global	fndsa_bench_nop
	.type	fndsa_bench_nop, %function
	.global	fndsa_bench_nop_f1
	.type	fndsa_bench_nop_f1, %function
	.global	fndsa_bench_nop_f2
	.type	fndsa_bench_nop_f2, %function
	.global	fndsa_bench_nop_p3
	.type	fndsa_bench_nop_p3, %function
	.global	fndsa_bench_nop_p4
	.type	fndsa_bench_nop_p4, %function
	.global	fndsa_bench_nop_sqrt_sigma
	.type	fndsa_bench_nop_sqrt_sigma, %function
do_nothing:
fndsa_bench_nop:
fndsa_bench_nop_f1:
fndsa_bench_nop_f2:
fndsa_bench_nop_p3:
fndsa_bench_nop_p4:
fndsa_bench_nop_sqrt_sigma:
	bx	lr
	.size	do_nothing,.-do_nothing

@ --- calibration --------------------------------------------------------
BENCH	bench_none,          do_nothing
BENCH	bench_base0,         fndsa_bw_none
BENCH	bench_base1,         fndsa_bw_base1
BENCH	bench_base1nc,       fndsa_bw_base1nc
BENCH	bench_base2,         fndsa_bw_base2
BENCH	bench_base2nc,       fndsa_bw_base2nc
BENCH	bench_base2s2,       fndsa_bw_base2s2
BENCH	bench_base2s2nc,     fndsa_bw_base2s2nc
BENCH	bench_base_addsub,   fndsa_bw_base_addsub
BENCH	bench_base_recip_x,  fndsa_bw_base_recip_x
BENCH	bench_base_mul_recip, fndsa_bw_base_mul_recip
BENCH	bench_base_sqrt_sigma, fndsa_bw_base_sqrt_sigma

@ --- scalar fpr operations ----------------------------------------------
@ Through wrappers, since an fpr does not always fit in registers.
BENCH	bench_fpr_add,       fndsa_bw_add
BENCH	bench_fpr_mul,       fndsa_bw_mul
BENCH	bench_fpr_div,       fndsa_bw_div
BENCH	bench_fpr_sqrt,      fndsa_bw_sqrt
BENCH	bench_fpr_sqrt_sigma, fndsa_bw_sqrt_sigma
BENCH	bench_fpr_add_sub,   fndsa_bw_add_sub
BENCH	bench_fpr_sub,       fndsa_bw_sub
BENCH	bench_fpr_neg,       fndsa_bw_neg
BENCH	bench_fpr_half,      fndsa_bw_half
BENCH	bench_fpr_double,    fndsa_bw_double
BENCH	bench_fpr_sqr,       fndsa_bw_sqr
BENCH	bench_fpr_inv,       fndsa_bw_inv
BENCH	bench_fpr_recip_x,   fndsa_bw_recip_x
BENCH	bench_fpr_mul_recip, fndsa_bw_mul_recip
BENCH	bench_fpr_div2e,     fndsa_bw_div2e
BENCH	bench_fpr_mul2e,     fndsa_bw_mul2e
BENCH	bench_fpr_of,        fndsa_bw_of
BENCH	bench_fpr_of32,      fndsa_bw_of32
BENCH	bench_fpr_rint,      fndsa_bw_rint
BENCH	bench_fpr_floor,     fndsa_bw_floor
BENCH	bench_fpr_trunc,     fndsa_bw_trunc
BENCH	bench_fpc_mul,       fndsa_bw_fpc_mul

@ fxp32.64 only; weak so other backends still link.
	.weak	fndsa_bw_fpc_mul_gm
BENCH	bench_fpc_mul_gm,    fndsa_bw_fpc_mul_gm
	.weak	fndsa_bw_fpc_mul_gmj
BENCH	bench_fpc_mul_gmj,   fndsa_bw_fpc_mul_gmj
	.weak	fndsa_bw_mul_gm
BENCH	bench_fpr_mul_gm,    fndsa_bw_mul_gm

@ --- direct upstream-method section -------------------------------------
@ The upstream bench_cm4 method: the target is called directly with
@ operands in registers. Targets are weak; absent ones are skipped.
.macro	BENCHD	name, target
	.align	2
	.global	\name
	.thumb
	.thumb_func
	.type	\name, %function
\name:
	push.w	{ r10, r11, lr }
	movw	r11, #(SYSCNT_ADDR & 0xFFFF)
	movt	r11, #(SYSCNT_ADDR >> 16)
	ldr.w	r10, [r11]
	bl	\target
	ldr.w	r0, [r11]
	sub.w	r0, r0, r10
	pop	{ r10, r11, pc }
	.size	\name,.-\name
.endm

	.weak	fndsa_fpr_mul
	.weak	fndsa_fpr_sqr
	.weak	fndsa_fpr_div
	.weak	fndsa_fpr_sqrt
	.weak	fndsa_fpr_of32
	.weak	fndsa_fpr_scaled

BENCHD	bench_direct_sqr,    fndsa_fpr_sqr
BENCHD	bench_direct_div,    fndsa_fpr_div
BENCHD	bench_direct_sqrt,   fndsa_fpr_sqrt
BENCHD	bench_direct_of32,   fndsa_fpr_of32
BENCHD	bench_direct_scaled, fndsa_fpr_scaled

@ For fxp32.64, fndsa_fpr_mul takes x = r0:r1:r2, y = r3:r4:r5 and
@ returns r0:r1:r2. From C, y1/y2 arrive on the stack and are loaded
@ into r4/r5 first.
	.align	2
	.global	bench_direct_mul
	.thumb
	.thumb_func
	.type	bench_direct_mul, %function
bench_direct_mul:
	push.w	{ r4-r11, lr }
	ldr.w	r4, [sp, #36]
	ldr.w	r5, [sp, #40]
	movw	r11, #(SYSCNT_ADDR & 0xFFFF)
	movt	r11, #(SYSCNT_ADDR >> 16)
	ldr.w	r10, [r11]
	bl	fndsa_fpr_mul
	ldr.w	r0, [r11]
	sub.w	r0, r0, r10
	pop.w	{ r4-r11, pc }
	.size	bench_direct_mul,.-bench_direct_mul

@ Address of each weak target (0 if absent).
	.section .rodata
	.align	2
	.global	fndsa_bench_dir_present
fndsa_bench_dir_present:
	.word	fndsa_fpr_mul
	.word	fndsa_fpr_sqr
	.word	fndsa_fpr_div
	.word	fndsa_fpr_sqrt
	.word	fndsa_fpr_of32
	.word	fndsa_fpr_scaled
	.text

@ --- polynomial / FFT level ---------------------------------------------
@ Called directly; LDL, gram and mulconst need a wrapper because they
@ would pass an argument on the stack.
BENCH	bench_fpoly_FFT,             fndsa_fpoly_FFT
/* The compact transforms exist only in FFT_NEW builds. */
#if FNDSA_FFT_NEW
BENCH	bench_fpoly_FFT_int32,       fndsa_fpoly_FFT_int32
BENCH	bench_fpoly_FFT_int32_selfadj, fndsa_fpoly_FFT_int32_selfadj
BENCH	bench_fpoly_FFT_int8,        fndsa_fpoly_FFT_int8
BENCH	bench_fpoly_FFT_uint16,      fndsa_fpoly_FFT_uint16
BENCH	bench_expand_int8,            fndsa_bench_expand_int8
BENCH	bench_expand_uint16,          fndsa_bench_expand_uint16
BENCH	bench_expand_uint16_fpr,      fndsa_bench_expand_uint16_fpr
#endif
BENCH	bench_fpoly_iFFT,            fndsa_fpoly_iFFT
BENCH	bench_fpoly_set_small,       fndsa_fpoly_set_small
BENCH	bench_fpoly_add,             fndsa_fpoly_add
BENCH	bench_fpoly_sub,             fndsa_fpoly_sub
BENCH	bench_fpoly_neg,             fndsa_fpoly_neg
BENCH	bench_fpoly_mul_fft,         fndsa_fpoly_mul_fft
BENCH	bench_fpoly_LDL_fft,         fndsa_bw_LDL_fft
BENCH	bench_fpoly_split_fft,       fndsa_fpoly_split_fft
BENCH	bench_fpoly_split_selfadj,   fndsa_fpoly_split_selfadj_fft
BENCH	bench_fpoly_merge_fft,       fndsa_fpoly_merge_fft
BENCH	bench_fpoly_gram_fft,        fndsa_bw_gram_fft
BENCH	bench_fpoly_mulconst,        fndsa_bw_mulconst
