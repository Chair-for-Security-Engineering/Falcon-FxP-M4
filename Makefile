# FN-DSA with a choice of fpr backends, natively and on the Cortex-M4 (pqm4).
#
#   make [BACKEND=fxp32.64|ieee|tw|prec] [CONFIG=<name>] [knobs] [target]
#
# Native targets: all, test, test_fpr, test_roundtrip, eval_precision,
# precision, ctgrind, gen-const, clean. Board targets: pqm4, pqm4-test,
# pqm4-speed, pqm4-signbench[-run], pqm4-bench[-run], pqm4-fprtest[-run],
# pqm4-ctbench[-run], pqm4-run, pqm4-clean. CONFIG picks a named paper
# configuration from configs.mk. Knob defaults live in knobs.h; a knob
# set on the command line (or by CONFIG) is passed as -D, anything else
# is left to the header.

# Knobs: Makefile variable and macro.
KNOBS = \
	SPLIT:FNDSA_SIGN_SPLIT \
	SYMP_ROOT:FNDSA_SYMP_ROOT \
	SYMP_TREE:FNDSA_SYMP_TREE \
	INT_GRAM:FNDSA_SIGN_INT_GRAM \
	FFT_NEW:FNDSA_FFT_NEW \
	RECIP_FUSED:FNDSA_SIGN_RECIP_FUSED \
	SQRT_FUSED:FNDSA_SIGN_SQRT_FUSED \
	INVQ_WIDE:FNDSA_SIGN_INVQ_WIDE \
	L10_64:FNDSA_SIGN_L10_64 \
	TIES_UP:FNDSA_FXP_TIES_UP \
	CMUL_FUSED:FNDSA_FXP_CMUL_FUSED \
	TWIDDLE_U64:FNDSA_TWIDDLE_U64 \
	FFT_ASM:FNDSA_FFT_ASM \
	SHAKE_HUANG:FNDSA_SHAKE_HUANG_BUF \
	GAMMA:FNDSA_KGEN_GAMMA
# SPLIT        split target, signature computed mod q'
# SYMP_ROOT    Symplecroot (d11 = q^2/g00 at the root)
# SYMP_TREE    Symplectree (excludes SYMP_ROOT)
# INT_GRAM     exact integer Gram matrix (needs SPLIT)
# FFT_NEW      compact FFT of c-fn-dsa-alt
# RECIP_FUSED  fused reciprocal at the LDL divisions (fxp32.64)
# SQRT_FUSED   fused sqrt*INV_SIGMA at the sampler leaves (fxp32.64)
# INVQ_WIDE    1/q as a normalized 64-bit word (fxp32.64)
# L10_64       non-root l10 in 0.63 format (fxp32.64)
# TIES_UP      products rounded half-up (fxp32.64)
# CMUL_FUSED   complex products rounded once (fxp32.64)
# TWIDDLE_U64  twiddles as unsigned 64-bit magnitudes (fxp32.64)
# FFT_ASM      FFT butterflies in assembly (fxp3264-cm4, needs TWIDDLE_U64)
# SHAKE_HUANG  Keccak of Huang et al. (Cortex-M4), 0 = c-fn-dsa's own
# GAMMA        keygen gamma checks

include configs.mk
ifdef CONFIG
cfg := $(CONFIG_$(CONFIG))
ifeq ($(cfg),)
$(error unknown CONFIG '$(CONFIG)'; see configs.mk)
endif
BACKEND ?= $(word 1,$(cfg))
BENCH_IMPL ?= $(word 2,$(cfg))
PREC_TARGET ?= $(word 1,$(cfg))
$(foreach kv,$(wordlist 3,$(words $(cfg)),$(cfg)),$(eval $(kv)))
endif

knob_var = $(word 1,$(subst :, ,$(1)))
knob_def = $(word 2,$(subst :, ,$(1)))
KNOB_FLAGS = $(strip $(foreach k,$(KNOBS),\
	$(if $($(call knob_var,$(k))),-D$(call knob_def,$(k))=$($(call knob_var,$(k))))))

BACKEND ?= fxp32.64
EXTRA_CFLAGS ?=
CC = clang
WARNFLAGS = -W -Wextra -Wundef -Wshadow
LIBS =

ifeq ($(BACKEND),ieee)
REPFLAGS =
OBJ_FPR = fpr/ieee/fpr_const_ieee.o fpr/ieee/sign_fpr.o
BINS = test_fndsa speed_fndsa
TESTBINS = test_fndsa
else ifeq ($(BACKEND),tw)
REPFLAGS = -DFNDSA_FPR_REP=FNDSA_FPR_REP_TW
OBJ_FPR = fpr/tw/fpr_const_tw.o fpr/tw/ops_c.o fpr/tw/triple_float.o \
	fpr/tw/triple_float_round.o
LIBS = -lm
BINS = test_fndsa test_roundtrip
TESTBINS = test_fndsa
else ifeq ($(BACKEND),fxp32.64)
REPFLAGS = -DFNDSA_FPR_REP=FNDSA_FPR_REP_FXP3264
OBJ_FPR = fpr/fxp32.64/fpr_const_fxp.o fpr/fxp32.64/ops_c.o
OBJ_FPRTEST = fpr/test_fpr.o fpr/fxp32.64/fxp_test_table.o fpr/fxp_mpfr.o
BINS = test_fpr test_fndsa test_roundtrip
TESTBINS = test_fpr test_fndsa
NEED_MPFR = 1
else ifeq ($(BACKEND),prec)
# MPFR reference plus an embedded target, whose fpr layer is compiled
# again into prec-owned objects with the target representation.
PREC_TARGET ?= none
ifeq ($(PREC_TARGET),none)
PREC_TGT_REP = 0
PREC_TGT_SRC =
else ifeq ($(PREC_TARGET),ieee)
PREC_TGT_REP = FNDSA_FPR_REP_IEEE
PREC_TGT_SRC = fpr/ieee/sign_fpr.c fpr/ieee/fpr_const_ieee.c
else ifeq ($(PREC_TARGET),tw)
PREC_TGT_REP = FNDSA_FPR_REP_TW
PREC_TGT_SRC = fpr/tw/ops_c.c fpr/tw/triple_float.c \
	fpr/tw/triple_float_round.c fpr/tw/fpr_const_tw.c
else ifeq ($(PREC_TARGET),fxp32.64)
PREC_TGT_REP = FNDSA_FPR_REP_FXP3264
PREC_TGT_SRC = fpr/fxp32.64/ops_c.c fpr/fxp32.64/fpr_const_fxp.c
else
$(error unknown PREC_TARGET '$(PREC_TARGET)' (none, ieee, tw or fxp32.64))
endif
PREC_TGT_SRC += fpr/prec/target_glue.c
PREC_TGT_OBJS = $(patsubst %.c,fpr/prec/tgt-$(PREC_TARGET)/%.o,$(notdir $(PREC_TGT_SRC)))
REPFLAGS = -DFNDSA_FPR_REP=FNDSA_FPR_REP_PREC -DFNDSA_PREC_TARGET=$(PREC_TGT_REP) \
	-DFNDSA_PREC_CHECKPOINTS=1
OBJ_FPR = fpr/prec/fpr_const_prec.o fpr/prec/ops_mpfr.o $(PREC_TGT_OBJS)
LIBS = -lmpfr -lgmp -lgc -lm
BINS = test_roundtrip eval_precision
TESTBINS = test_roundtrip
NEED_MPFR = 1
else
$(error unknown BACKEND '$(BACKEND)' (ieee, tw, fxp32.64 or prec))
endif

CFLAGS = $(WARNFLAGS) -O2 -I. $(REPFLAGS) $(KNOB_FLAGS) $(EXTRA_CFLAGS)
PREC_TGT_CFLAGS = $(WARNFLAGS) -O2 -I. $(if $(filter none,$(PREC_TARGET)),\
	-DPREC_TGT_NONE=1,-DFNDSA_FPR_REP=$(PREC_TGT_REP)) \
	-DFNDSA_PREC_TARGET=$(PREC_TGT_REP) $(KNOB_FLAGS) $(EXTRA_CFLAGS)

OBJ = codec.o mq.o sha3.o sysrng.o util.o \
	kgen.o kgen_fxp.o kgen_gauss.o kgen_mp31.o kgen_ntru.o kgen_poly.o \
	kgen_zint31.o \
	sign.o sign_core.o sign_fpoly.o sign_mq2.o sign_sampler.o vrfy.o \
	$(OBJ_FPR)
HDRS = $(wildcard *.h fpr/*.h fpr/*/*.h)

# Pinned nixpkgs, also for nix-shell's own bash (no channel needed).
NIX_SHELL = nix-shell -I nixpkgs=$(CURDIR)/nixpkgs.nix

# mpfr/gmp (test_fpr, prec, gen-const) and gc (prec) come from shell.nix
# when missing.
ifneq ($(filter gen-const,$(MAKECMDGOALS)),)
NEED_MPFR = 1
endif
ifdef NEED_MPFR
ifneq ($(shell $(CC) -E -x c -include mpfr.h $(if $(filter prec,$(BACKEND)),-include gc.h) \
	/dev/null >/dev/null 2>&1 && echo 1),1)
ifneq ($(shell command -v nix-shell 2>/dev/null),)
NATIVE_NIX = 1
endif
endif
endif

# Objects of another backend or knob set are removed before a build.
BSTAMP = .backend-stamp
BSTAMP_SIG = $(BACKEND) $(PREC_TARGET) $(KNOB_FLAGS) $(EXTRA_CFLAGS)
.DEFAULT_GOAL := all
.PHONY: FORCE all test clean gen-const ctgrind precision
FORCE:
$(BSTAMP): FORCE
	@if [ "`cat $@ 2>/dev/null`" != "$(BSTAMP_SIG)" ]; then \
		find . -name '*.o' -not -path './pqm4/*' -delete; \
		echo "$(BSTAMP_SIG)" > $@; \
	fi

NATIVE_GOALS = all test test_fndsa speed_fndsa test_roundtrip test_fpr \
	test_ct eval_precision gen-const
ifdef NATIVE_NIX
.PHONY: $(NATIVE_GOALS)
$(NATIVE_GOALS):
	+$(NIX_SHELL) shell.nix --run '$(MAKE) $@ $(MAKEOVERRIDES) NATIVE_NIX='
else
all: $(BINS)

test: $(TESTBINS)
	@set -e; for t in $(TESTBINS); do echo "== ./$$t"; ./$$t; done

test_fndsa: $(OBJ) test_fndsa.o test_sampler.o test_sign.o
speed_fndsa: $(OBJ) speed_fndsa.o
test_roundtrip: $(OBJ) test_roundtrip.o
test_ct: $(OBJ) test_ct.o
eval_precision: $(OBJ) eval_precision.o
test_fndsa speed_fndsa test_roundtrip test_ct eval_precision:
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)
test_fpr: $(OBJ) $(OBJ_FPRTEST)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS) -lmpfr -lgmp

# Host tool, built in one step so it never picks up the backend flags.
GEN_CC ?= cc
gen-const:
	$(GEN_CC) -O2 -W -Wextra -o fpr/gen/gen_const fpr/gen/gen_const.c \
		-lmpfr -lgmp
	./fpr/gen/gen_const --twiddle-u64 fpr/fxp32.64 32 64
endif

%.o: %.c $(HDRS) $(BSTAMP)
	$(CC) $(CFLAGS) -c -o $@ $<
define prec_tgt_rule
fpr/prec/tgt-$(PREC_TARGET)/$(basename $(notdir $(1))).o: $(1) $(HDRS) $(BSTAMP)
	@mkdir -p $$(@D)
	$$(CC) $$(PREC_TGT_CFLAGS) -c -o $$@ $$<
endef
$(foreach s,$(PREC_TGT_SRC),$(eval $(call prec_tgt_rule,$(s))))
test_sampler.o: sign_sampler.c
test_sign.o: sign_sampler.c sign_core.c

clean:
	find . -name '*.o' -not -path './pqm4/*' -delete
	rm -rf $(BSTAMP) fpr/prec/tgt-* fpr/gen/gen_const test_fndsa speed_fndsa \
		test_roundtrip test_ct test_fpr eval_precision

# Constant-time check: test_ct under Valgrind memcheck (see test_ct.c).
CT_KEYS ?= 4
CT_SIGS ?= 2
CT_VALGRIND ?= valgrind --error-exitcode=1 --track-origins=yes --num-callers=40
ctgrind:
	$(MAKE) test_ct EXTRA_CFLAGS='$(EXTRA_CFLAGS) -DFNDSA_CTGRIND=1 -g'
	$(CT_VALGRIND) ./test_ct $(CT_KEYS) $(CT_SIGS)

# ----------------------------------------------------------------------
# Precision evaluation (prec backend): RUNS keys x SIGS signatures,
# sharded over JOBS processes (the data does not depend on JOBS), merged
# into $(PREC_RESULTS)/<name>.{stats,worst.csv,perkey.csv,summary.txt}
# (and .csv with CSV=1). The name carries the compiled-in knobs.
#   make precision [CONFIG=..|PREC_TARGET=..] [knobs] [LOGN=10 RUNS=25 SIGS=1
#        SEED=1 JOBS=N KEY_FIRST= TAG= CSV=1 KEEP_SHARDS=1]
PREC_TARGETS ?= ieee tw
LOGN ?= 10
RUNS ?= 25
SIGS ?= 1
SEED ?= 1
JOBS ?= $(shell nproc 2>/dev/null || echo 4)
PREC_RESULTS ?= results/precision
TREE_ID ?= $(shell cat TREE_ID 2>/dev/null || git rev-parse --short HEAD 2>/dev/null || echo unknown)

precision:
	@set -e; \
	if [ -n "$(EXTRA_CFLAGS)" ] && [ "$(ALLOW_EXTRA_CFLAGS)" != 1 ]; then \
		echo "precision: EXTRA_CFLAGS is not in the result name;" \
			"set ALLOW_EXTRA_CFLAGS=1 and a TAG" >&2; exit 1; \
	fi; \
	mkdir -p $(PREC_RESULTS)/shards; \
	for t in $(or $(PREC_TARGET),$(PREC_TARGETS)); do \
		$(MAKE) BACKEND=prec PREC_TARGET=$$t eval_precision; \
		name="$$t$$(./eval_precision --variant-tokens)-logn$(LOGN)-r$(RUNS)"; \
		[ "$(SIGS)" = 1 ] || name="$$name-m$(SIGS)"; \
		[ -z "$(KEY_FIRST)" ] || name="$$name-k$(KEY_FIRST)"; \
		name="$$name-s$(SEED)$(TAG)"; \
		base="$(PREC_RESULTS)/$$name"; shp="$(PREC_RESULTS)/shards/$$name"; \
		echo "== precision: $$name jobs=$(JOBS) tree=$(TREE_ID) =="; \
		rm -f "$$shp".*; \
		seq 0 $$(( $(JOBS) - 1 )) | xargs -P $(JOBS) -n 1 sh -c \
			'./eval_precision --logn $(LOGN) --runs $(RUNS) --sigs $(SIGS) \
				$(if $(KEY_FIRST),--key-first $(KEY_FIRST)) --seed $(SEED) \
				--tree $(TREE_ID) --shard "$$1" --jobs $(JOBS) \
				$(if $(filter 1,$(CSV)),--csv) \
				--out "$$(printf "%s.%04d" "$$0" "$$1")"' "$$shp"; \
		./eval_precision --worst-header > "$$base.worst.csv"; \
		cat $$(ls "$$shp".*.worst | sort) >> "$$base.worst.csv"; \
		if [ "$(CSV)" = 1 ]; then \
			./eval_precision --csv-header > "$$base.csv"; \
			cat $$(ls "$$shp".*.csv | sort) >> "$$base.csv"; \
		fi; \
		./eval_precision --out "$$base" --tree $(TREE_ID) \
			--merge "$$shp".*.stats | tee "$$base.summary.txt"; \
		[ "$(KEEP_SHARDS)" = 1 ] || rm -f "$$shp".*; \
	done

# ----------------------------------------------------------------------
# pqm4 (pqm4/, vendored). Every implementation dir under
# pqm4/crypto_sign/<scheme>/ is built; BENCH_IMPL (or CONFIG) picks one,
# SCHEME=512|1024 one scheme. The knobs reach the compiler through
# MUPQ_EXTRA_CPPFLAGS, a retained pqm4 variable (a change rebuilds).
#   make pqm4-signbench-run CONFIG=fxp-rec KEYS=2 ITER=25
#   make pqm4-run BENCH_IMPL=tw-m4f BIN=speed
PLATFORM ?= nucleo-l4r5zi
ifeq ($(shell uname),Darwin)
PQM4_UART ?= $(firstword $(wildcard /dev/tty.usbmodem*))
else
PQM4_UART ?= /dev/ttyACM0
endif
SCHEME ?=
ITER ?= 1
KEYS ?= 2
BIN ?= test
PQM4_SCHEMES = $(if $(SCHEME),fndsa$(SCHEME),fndsa512 fndsa1024)
PQM4_ALL = $(patsubst pqm4/%/config.mk,%,$(sort $(wildcard \
	$(foreach s,$(PQM4_SCHEMES),pqm4/crypto_sign/$(s)/*/config.mk))))
PQM4_IMPLS = $(if $(BENCH_IMPL),$(filter %/$(BENCH_IMPL),$(PQM4_ALL)),$(PQM4_ALL))
# The fpr benchmark and test are the same code in both schemes: 512 only.
PQM4_IMPLS_512 = $(sort $(patsubst crypto_sign/fndsa1024/%,crypto_sign/fndsa512/%,\
	$(PQM4_IMPLS)))
PQM4_EXT = $(if $(filter nucleo-l4r5zi,$(PLATFORM)),hex,bin)
# CONFIG is also a pqm4 variable (its retained-settings file): reset it.
# Recipes that run it start with '+', so make shares its job slots.
PQM4_MAKE = $(MAKE) -C pqm4 CONFIG=obj/.config.mk PLATFORM=$(PLATFORM) MUPQ_ITERATIONS=$(ITER) \
	MUPQ_EXTRA_CPPFLAGS='$(strip $(KNOB_FLAGS) -DSIGNBENCH_KEYS=$(KEYS) $(EXTRA_CFLAGS))'
# $(call pqm4_bin,<impl>,<binary>): the image path below pqm4/.
pqm4_bin = bin/$$(echo $(1) | tr / _)_$(2).$(PQM4_EXT)
# $(call pqm4_flash,<image>,<iterations>): flash, print the board output.
pqm4_flash = cd pqm4 && python3 -c "from mupq import platforms; \
p = platforms.OpenOCD('st_nucleo_l4r5.cfg', '$(PQM4_UART)'); p.__enter__(); \
print(p.run('$(1)', $(2))); p.__exit__(None, None, None)"
PQM4_UART_CHECK = @[ -e "$(PQM4_UART)" ] || { echo "no board at" \
	"'$(PQM4_UART)'; set PQM4_UART=/dev/..." >&2; exit 1; }

PQM4_GOALS = pqm4 pqm4-test pqm4-speed pqm4-signbench pqm4-signbench-run \
	pqm4-bench pqm4-bench-run pqm4-fprtest pqm4-fprtest-run pqm4-ctbench \
	pqm4-ctbench-run pqm4-run
.PHONY: $(PQM4_GOALS) pqm4-clean
ifneq ($(filter pqm4%,$(MAKECMDGOALS)),)
ifeq ($(PQM4_IMPLS),)
$(error no pqm4 implementation matches BENCH_IMPL='$(BENCH_IMPL)')
endif
ifeq ($(shell command -v arm-none-eabi-gcc >/dev/null && python3 -c 'import tqdm, serial' 2>/dev/null && echo ok),)
ifneq ($(shell command -v nix-shell 2>/dev/null),)
PQM4_NIX = 1
endif
endif
endif

ifdef PQM4_NIX
$(PQM4_GOALS):
	+$(NIX_SHELL) shell-cm4.nix --run '$(MAKE) $@ $(MAKEOVERRIDES) PQM4_NIX='
else
pqm4:
	+@set -e; for i in $(PQM4_IMPLS); do \
		$(PQM4_MAKE) IMPLEMENTATION_PATH=$$i tests-bin; done

# On-device self-test: 30 OK and no ERROR per implementation.
pqm4-test:
	$(PQM4_UART_CHECK)
	+@fail=0; for i in $(PQM4_IMPLS); do \
		img=$(call pqm4_bin,$$i,test); \
		if ! $(PQM4_MAKE) IMPLEMENTATION_PATH=$$i $$img >/dev/null 2>&1; then \
			echo "$$i: BUILD-FAILED"; fail=1; continue; fi; \
		out=`($(call pqm4_flash,$$img,1)) 2>&1`; \
		ok=`echo "$$out" | grep -c OK`; err=`echo "$$out" | grep -c ERROR`; \
		if [ $$ok = 30 ] && [ $$err = 0 ]; then echo "$$i: PASS"; \
		else echo "$$i: FAIL (OK=$$ok ERROR=$$err)"; fail=1; fi; \
	done; [ $$fail = 0 ]

# One binary per implementation. pqm4-signbench, pqm4-bench and
# pqm4-fprtest only build it; the other targets also flash it and print
# the board output.
pqm4-speed: PQM4_BIN = speed
pqm4-speed: PQM4_COUNT = $(ITER)
pqm4-signbench pqm4-signbench-run: PQM4_BIN = signbench
pqm4-signbench-run: PQM4_COUNT = $(KEYS) * $(ITER)
pqm4-bench pqm4-bench-run: PQM4_BIN = fprbench
pqm4-fprtest pqm4-fprtest-run: PQM4_BIN = fprtest
pqm4-run: PQM4_BIN = $(BIN)
pqm4-run: PQM4_COUNT = $(ITER)
pqm4-bench pqm4-bench-run: PQM4_SET = $(PQM4_IMPLS_512)
pqm4-fprtest pqm4-fprtest-run: PQM4_SET = $(filter %-m4f %-cm4,$(PQM4_IMPLS_512))
pqm4-run: PQM4_SET = $(filter crypto_sign/fndsa$(or $(SCHEME),512)/$(BENCH_IMPL),\
	$(PQM4_IMPLS))
PQM4_SET ?= $(PQM4_IMPLS)
PQM4_COUNT ?= 1
pqm4-signbench pqm4-bench pqm4-fprtest:
	+@set -e; for i in $(PQM4_SET); do \
		$(PQM4_MAKE) IMPLEMENTATION_PATH=$$i $(call pqm4_bin,$$i,$(PQM4_BIN)); done
pqm4-speed pqm4-signbench-run pqm4-bench-run pqm4-fprtest-run pqm4-run:
	$(PQM4_UART_CHECK)
	@[ "$@" != pqm4-run ] || [ -n "$(BENCH_IMPL)" ] || { \
		echo "pqm4-run needs BENCH_IMPL or CONFIG" >&2; exit 1; }
	+@set -e; for i in $(PQM4_SET); do \
		img=$(call pqm4_bin,$$i,$(PQM4_BIN)); \
		$(PQM4_MAKE) IMPLEMENTATION_PATH=$$i $$img; \
		echo; echo "===== $$i [$(PQM4_BIN)] ====="; \
		($(call pqm4_flash,$$img,$(PQM4_COUNT))); done

# Constant-time rows of the fpr benchmark (FPRBENCH_CT in fprbench.c).
pqm4-ctbench pqm4-ctbench-run:
	$(MAKE) $(@:pqm4-ctbench%=pqm4-bench%) EXTRA_CFLAGS='$(EXTRA_CFLAGS) -DFPRBENCH_CT=1'
endif

pqm4-clean:
	rm -rf pqm4/obj pqm4/elf pqm4/bin
