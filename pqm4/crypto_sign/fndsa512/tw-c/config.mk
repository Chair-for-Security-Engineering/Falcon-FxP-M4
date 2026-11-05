# Triple-word (TWFalcon), portable C.
elf/crypto_sign_fndsa512_tw-c_%.elf: CPPFLAGS += -DFNDSA_ASM_CORTEXM4=0 -DFNDSA_FPR_REP=FNDSA_FPR_REP_TW
obj/libcrypto_sign_fndsa512_tw-c.a: CPPFLAGS += -DFNDSA_ASM_CORTEXM4=0 -DFNDSA_FPR_REP=FNDSA_FPR_REP_TW
