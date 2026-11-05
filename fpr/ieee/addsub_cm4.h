#ifndef FNDSA_FPR_IEEE_ADDSUB_CM4_H__
#define FNDSA_FPR_IEEE_ADDSUB_CM4_H__

/* Cortex-M4 shim for FPR_ADD_SUB, from c-fn-dsa's sign_inner.h. Included
   by fpr_ieee.h. */

/* On the ARM Cortex-M4, we have a dedicated function, but since it does
   not follow the AAPCS requirements (it returns two 64-bit values), it
   uses a custom convention which requires inline assembly for invocation.
   In particular, about all registers are clobbered. */
#define FPR_ADD_SUB(a, b, x, y)   do { \
		fpr t_add_sub_x = (x); \
		fpr t_add_sub_y = (y); \
		register uint32_t t_add_sub_x0 __asm__("r0") = \
			(uint32_t)t_add_sub_x; \
		register uint32_t t_add_sub_x1 __asm__("r1") = \
			(uint32_t)(t_add_sub_x >> 32); \
		register uint32_t t_add_sub_y0 __asm__("r2") = \
			(uint32_t)t_add_sub_y; \
		register uint32_t t_add_sub_y1 __asm__("r3") = \
			(uint32_t)(t_add_sub_y >> 32); \
		__asm__( \
			"bl	fndsa_fpr_add_sub" \
			: "+r" (t_add_sub_x0), "+r" (t_add_sub_x1), \
			  "+r" (t_add_sub_y0), "+r" (t_add_sub_y1) \
			: \
			: "r4", "r5", "r6", "r7", "r8", \
			  "r10", "r11", "r12", "r14", "cc"); \
		(a) = (uint64_t)t_add_sub_x0 | ((uint64_t)t_add_sub_x1 << 32); \
		(b) = (uint64_t)t_add_sub_y0 | ((uint64_t)t_add_sub_y1 << 32); \
	} while (0)

#endif
