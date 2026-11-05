# ieee, Cortex-M4 assembly.
elf/crypto_sign_fndsa512_ieee-m4f_%.elf: CPPFLAGS += -DFNDSA_ASM_CORTEXM4=1
obj/libcrypto_sign_fndsa512_ieee-m4f.a: CPPFLAGS += -DFNDSA_ASM_CORTEXM4=1
