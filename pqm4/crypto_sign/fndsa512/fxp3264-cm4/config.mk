# Fixed-point 32.64, Cortex-M4 assembly.
elf/crypto_sign_fndsa512_fxp3264-cm4_%.elf: CPPFLAGS += -DFNDSA_ASM_CORTEXM4=1 -DFNDSA_FPR_REP=FNDSA_FPR_REP_FXP3264
obj/libcrypto_sign_fndsa512_fxp3264-cm4.a: CPPFLAGS += -DFNDSA_ASM_CORTEXM4=1 -DFNDSA_FPR_REP=FNDSA_FPR_REP_FXP3264
