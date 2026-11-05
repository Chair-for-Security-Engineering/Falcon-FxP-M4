#!/usr/bin/env bash
# KAT gate over the backend x knob matrix and the paper configurations:
# 'make test' per configuration below, which must reach the test_fndsa
# KAT suite (fxp32.64 runs test_fpr first). Runs in this tree, one
# configuration after the other (the build stamp rebuilds on every knob
# change); logs under $OUT, exit status 1 if any fails.
#
#   [OUT=results/kat-matrix] [J=16] scripts/kat-matrix.sh [filter]
#
# filter is a grep -E pattern on the configuration line.
set -u
cd "$(dirname "$0")/.."
OUT=${OUT:-results/kat-matrix}
J=${J:-16}
FILTER=${1:-}
mkdir -p "$OUT"

configs() {
	cat <<'EOF'
BACKEND=ieee
BACKEND=ieee GAMMA=0
BACKEND=ieee EXTRA_CFLAGS=-DFNDSA_SHAKE256X4=1
BACKEND=ieee GAMMA=0 EXTRA_CFLAGS=-DFNDSA_SHAKE256X4=1
BACKEND=ieee "EXTRA_CFLAGS=-DFNDSA_AVX2=0 -DFNDSA_SSE2=0"
BACKEND=ieee SPLIT=1
BACKEND=ieee SPLIT=1 INT_GRAM=0
BACKEND=ieee SPLIT=1 SYMP_ROOT=1
BACKEND=ieee SPLIT=1 SYMP_ROOT=1 FFT_NEW=1
BACKEND=ieee SYMP_ROOT=1
BACKEND=ieee SYMP_TREE=1
BACKEND=ieee FFT_NEW=1
BACKEND=tw
BACKEND=tw GAMMA=0
BACKEND=tw SPLIT=1
BACKEND=tw SPLIT=1 INT_GRAM=0
BACKEND=tw SYMP_ROOT=1
BACKEND=tw SYMP_TREE=1
BACKEND=fxp32.64
BACKEND=fxp32.64 GAMMA=0
BACKEND=fxp32.64 EXTRA_CFLAGS=-DFNDSA_SHAKE256X4=1
BACKEND=fxp32.64 SPLIT=0
BACKEND=fxp32.64 SPLIT=0 SYMP_ROOT=0
BACKEND=fxp32.64 SPLIT=0 INVQ_WIDE=0
BACKEND=fxp32.64 INT_GRAM=0
BACKEND=fxp32.64 SYMP_ROOT=0
BACKEND=fxp32.64 SYMP_TREE=1
BACKEND=fxp32.64 L10_64=0
BACKEND=fxp32.64 TIES_UP=0
BACKEND=fxp32.64 FFT_NEW=1
BACKEND=fxp32.64 FFT_NEW=1 TIES_UP=0
BACKEND=fxp32.64 TWIDDLE_U64=0
BACKEND=fxp32.64 CMUL_FUSED=0
BACKEND=fxp32.64 RECIP_FUSED=0
BACKEND=fxp32.64 SQRT_FUSED=0
BACKEND=fxp32.64 INVQ_WIDE=0
BACKEND=fxp32.64 INVQ_WIDE=0 TIES_UP=0
EOF
	sed -n 's/^CONFIG_\([^ ]*\) .*/CONFIG=\1/p' configs.mk
}

pass=0 fail=0 failed=
while IFS= read -r cfg; do
	[ -n "$FILTER" ] && ! grep -Eq -- "$FILTER" <<<"$cfg" && continue
	name=$(tr -c 'A-Za-z0-9.=_-' '_' <<<"$cfg" | sed 's/_*$//')
	log=$OUT/$name.log
	cmd="make -s -j$J $cfg test"
	want='^Test KAT: .* done\.$'
	if eval "$cmd" > "$log" 2>&1 && grep -Eq "$want" "$log"; then
		echo "PASS  $cfg"
		pass=$((pass + 1))
	else
		echo "FAIL  $cfg   ($log)"
		fail=$((fail + 1))
		failed="$failed
  $cfg"
	fi
done < <(configs)
echo "$pass passed, $fail failed$failed"
[ "$fail" = 0 ]
