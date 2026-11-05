#!/usr/bin/env bash
# Sign-cycle table of the paper: pqm4-signbench-run for every row of
# CONFIG (default scripts/board-configs.txt, the presets of configs.mk),
# both schemes, PASSES repeats. The defaults are the paper's size, 100
# keys x 100 signatures; finished runs are skipped, so a re-run resumes.
# Tables: scripts/board_tables.py OUT.
#
#   [PQM4_UART=/dev/ttyACM0] [STLINK_SERIAL=..] [OUT=..] [KEYS=100 ITER=100]
#       [PASSES=1] [CONFIG_FILE=..] scripts/board.sh
set -u

cd "$(dirname "$0")/.."

OUT=${OUT:-results/board/$(date +%F)}
PASSES=${PASSES:-1}
KEYS=${KEYS:-100}
ITER=${ITER:-100}
SCHEMES=${SCHEMES:-"512 1024"}
# Generous: the slowest row (ieee-c, Falcon-1024) needs about 5 s per
# signature.
RUN_TIMEOUT=${RUN_TIMEOUT:-$((KEYS * ITER * 10 + 600))}
export PQM4_UART=${PQM4_UART:-/dev/ttyACM0}
# STLINK_SERIAL picks the probe when several boards are connected.
[ -z "${STLINK_SERIAL:-}" ] || export STLINK_SERIAL
DRY_RUN=${DRY_RUN:-0}
CONFIG_FILE=${CONFIG_FILE:-scripts/board-configs.txt}

# impl|row name|make knobs. Row names double as log names.
CONFIGS=$(cat "$CONFIG_FILE")

mkdir -p "$OUT"
printf '%s\n' "$CONFIGS" | grep -v '^$' > "$OUT/configs.txt"
{
	echo "tree: $(cat TREE_ID 2>/dev/null || git rev-parse --short HEAD 2>/dev/null)"
	echo "start: $(date)"
	echo "host: $(hostname)  keys: $KEYS  iter: $ITER  passes: $PASSES"
} >> "$OUT/README.txt"

if [ "$DRY_RUN" = 1 ]; then
	printf '%s\n' "$CONFIGS" | grep -v '^$'
	exit 0
fi

want=$((KEYS * ITER))
status="$OUT/status.tsv"

for pass in $(seq 1 "$PASSES"); do
	for scheme in $SCHEMES; do
		printf '%s\n' "$CONFIGS" | grep -v '^$' | while IFS='|' read -r impl name knobs; do
			log="$OUT/p$pass/$impl-$scheme-$name.log"
			mkdir -p "$OUT/p$pass"
			have=$(grep -c '^sign cycles:' "$log" 2>/dev/null || true)
			if [ "${have:-0}" -ge "$want" ]; then
				echo "skip $log"
				continue
			fi
			echo "== $(date +%T) pass $pass $scheme $impl $name [$knobs]"
			# shellcheck disable=SC2086
			timeout -k 30 "$RUN_TIMEOUT" nice -n 5 make pqm4-signbench-run \
				BENCH_IMPL="$impl" SCHEME="$scheme" KEYS="$KEYS" ITER="$ITER" \
				PQM4_UART="$PQM4_UART" $knobs > "$log" 2>&1
			rc=$?
			have=$(grep -c '^sign cycles:' "$log" 2>/dev/null || true)
			printf '%s\t%s\t%s\t%s\trc=%s\tsamples=%s\n' \
				"$pass" "$scheme" "$impl" "$name" "$rc" "${have:-0}" >> "$status"
			echo "   rc=$rc samples=${have:-0}"
		done
	done
done

echo "end: $(date)" >> "$OUT/README.txt"
python3 scripts/board_tables.py "$OUT" > "$OUT/tables.md" 2> "$OUT/tables.err"
touch "$OUT/DONE"
