# ieee, portable C.
elf/crypto_sign_fndsa512_ieee-c_%.elf: CPPFLAGS += -DFNDSA_ASM_CORTEXM4=0
obj/libcrypto_sign_fndsa512_ieee-c.a: CPPFLAGS += -DFNDSA_ASM_CORTEXM4=0
