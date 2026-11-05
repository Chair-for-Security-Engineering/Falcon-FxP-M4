# Triple-word (TWFalcon), Cortex-M4 assembly; ffSampling stays C.
elf/crypto_sign_fndsa512_tw-m4f_%.elf: CPPFLAGS += -DFNDSA_ASM_CORTEXM4=1 -DFNDSA_FPR_REP=FNDSA_FPR_REP_TW
obj/libcrypto_sign_fndsa512_tw-m4f.a: CPPFLAGS += -DFNDSA_ASM_CORTEXM4=1 -DFNDSA_FPR_REP=FNDSA_FPR_REP_TW
