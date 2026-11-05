#!/usr/bin/env bash
# Constant-time check on the board: 'make pqm4-ctbench-run' for each
# signing configuration of the paper. Every fpr operation and FFT-domain
# routine is timed on many random inputs; a row passes iff its cycle count
# is the same on all of them. Summary in $OUT/summary.txt.
#
#   [PQM4_UART=/dev/ttyACM0] [STLINK_SERIAL=..] [OUT=results/ct-board] \
#       scripts/ct-board-matrix.sh [grep -E filter]
# STLINK_SERIAL picks the probe when several boards are connected.
set -u
cd "$(dirname "$0")/.."
export PQM4_UART=${PQM4_UART:-/dev/ttyACM0}
[ -z "${STLINK_SERIAL:-}" ] || export STLINK_SERIAL
OUT=${OUT:-results/ct-board}
FILTER=${1:-}
mkdir -p "$OUT"

# The paper's board configurations (configs.mk) except ieee-base-origshake,
# whose fpr code is that of ieee-base.
configs() {
	cat <<EOF
ieee-m4f|ieee-base|CONFIG=ieee-base
ieee-m4f|ieee-base-fftnew|CONFIG=ieee-base-fftnew
ieee-m4f|ieee-split|CONFIG=ieee-split
ieee-m4f|ieee-rec|CONFIG=ieee-rec
ieee-m4f|ieee-rec-fftnew|CONFIG=ieee-rec-fftnew
fxp3264-cm4|fxp-split|CONFIG=fxp-split
fxp3264-cm4|fxp-unfused|CONFIG=fxp-unfused
fxp3264-cm4|fxp-rec|CONFIG=fxp-rec
fxp3264-cm4|fxp-rec-fftnew|CONFIG=fxp-rec-fftnew
fxp3264-cm4|fxp-nonsplit|CONFIG=fxp-nonsplit
fxp3264-cm4|fxp-nonsplit-tiesup|CONFIG=fxp-nonsplit-tiesup
fxp3264-cm4|fxp-nonsplit-tiesup-fftnew|CONFIG=fxp-nonsplit-tiesup-fftnew
EOF
}

fail=0
summary=$OUT/summary.txt
echo "ct-board-matrix $(date -Is) tree=$(cat TREE_ID 2>/dev/null || git rev-parse --short HEAD 2>/dev/null)" > "$summary"
while IFS='|' read -r impl name knobs; do
	[ -n "$FILTER" ] && ! grep -Eq -- "$FILTER" <<<"$impl|$name|$knobs" && continue
	log=$OUT/$impl-$name.log
	# shellcheck disable=SC2086
	timeout 2700 make pqm4-ctbench-run BENCH_IMPL="$impl" \
		PQM4_UART="$PQM4_UART" $knobs > "$log" 2>&1
	res=$(grep -E '^ct:' "$log")
	case "$res" in
	*" 0 FAIL"*) st=PASS ;;
	*) st=FAIL; fail=$((fail + 1)) ;;
	esac
	echo "$st $impl-$name ${res:-no result}" | tee -a "$summary"
	grep -E 'jit [0-9]+ +FAIL$' "$log" | sed 's/^/    /' >> "$summary"
done < <(configs)
echo "$fail failing" | tee -a "$summary"
[ "$fail" = 0 ]
