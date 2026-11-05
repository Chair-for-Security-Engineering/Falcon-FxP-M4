# ieee, portable C.
elf/crypto_sign_fndsa1024_ieee-c_%.elf: CPPFLAGS += -DFNDSA_ASM_CORTEXM4=0
obj/libcrypto_sign_fndsa1024_ieee-c.a: CPPFLAGS += -DFNDSA_ASM_CORTEXM4=0
