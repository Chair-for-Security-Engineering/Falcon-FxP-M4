#!/usr/bin/env sh
# The precision table of the paper: 'make precision' for each of its rows,
# ieee and fxp32.64, both degrees, then scripts/precision_tables.py. The
# default size is a short run of 10 keys x 10 signatures; the paper used
# KEYS=100000 SIGS=100 (about 50 / 90 core-hours per row at Falcon-512 /
# 1024). Finished rows are skipped, so a re-run resumes.
#
#   [KEYS=10] [SIGS=10] [LOGNS="9 10"] [JOBS=N] [RES=results/precision-paper]
#       [DRY_RUN=1] scripts/precision.sh
set -eu
cd "$(dirname "$0")/.."

KEYS=${KEYS:-10}
SIGS=${SIGS:-10}
LOGNS=${LOGNS:-"9 10"}
JOBS=${JOBS:-$(nproc 2>/dev/null || echo 4)}
SEED=${SEED:-1}
RES=${RES:-results/precision-paper}
DRY_RUN=${DRY_RUN:-0}

# One line per table row: number|label|ieee knobs|fxp32.64 knobs, with
# '-' where the paper has no cell. Every knob not listed is 0, except
# CMUL_FUSED and GAMMA; L10_64 and ties-up do not change the precision.
rows() {
	cat <<'EOF'
1|baseline (non-split)||
2|Symplecroot|SYMP_ROOT=1|SYMP_ROOT=1
3|split|SPLIT=1|SPLIT=1
4|split, int. Gram|SPLIT=1 INT_GRAM=1|SPLIT=1 INT_GRAM=1
5|split, Symplecroot|SPLIT=1 SYMP_ROOT=1|SPLIT=1 SYMP_ROOT=1
6|+ int. Gram|SPLIT=1 SYMP_ROOT=1 INT_GRAM=1|SPLIT=1 SYMP_ROOT=1 INT_GRAM=1
7|+ FFT-new|SPLIT=1 SYMP_ROOT=1 INT_GRAM=1 FFT_NEW=1|SPLIT=1 SYMP_ROOT=1 INT_GRAM=1 FFT_NEW=1
8|non-split, fused reciprocal|-|RECIP_FUSED=1
9|split, fused reciprocal|-|SPLIT=1 RECIP_FUSED=1
10|+ fused sqrt|-|SPLIT=1 RECIP_FUSED=1 SQRT_FUSED=1
11|+ wide 1/q|-|SPLIT=1 RECIP_FUSED=1 SQRT_FUSED=1 INVQ_WIDE=1
12|+ Symplecroot|-|SPLIT=1 RECIP_FUSED=1 SQRT_FUSED=1 INVQ_WIDE=1 SYMP_ROOT=1
13|+ int. Gram (recommended)|-|SPLIT=1 RECIP_FUSED=1 SQRT_FUSED=1 INVQ_WIDE=1 SYMP_ROOT=1 INT_GRAM=1
14|non-split, fused recip./sqrt, wide 1/q, Symplecroot|-|RECIP_FUSED=1 SQRT_FUSED=1 INVQ_WIDE=1 SYMP_ROOT=1
15|recommended, FFT-new|-|SPLIT=1 RECIP_FUSED=1 SQRT_FUSED=1 INVQ_WIDE=1 SYMP_ROOT=1 INT_GRAM=1 FFT_NEW=1
EOF
}

# All knobs explicit: the header defaults differ between the targets.
ZERO="SPLIT=0 SYMP_ROOT=0 SYMP_TREE=0 INT_GRAM=0 FFT_NEW=0 RECIP_FUSED=0
	SQRT_FUSED=0 INVQ_WIDE=0 L10_64=0 TIES_UP=0 CMUL_FUSED=1 GAMMA=1"

mkdir -p "$RES"
fail=0
# A file, not a pipe: the loop must not run in a subshell (fail).
ROWS=$(mktemp); trap 'rm -f "$ROWS"' EXIT
rows > "$ROWS"
for logn in $LOGNS; do
	while IFS='|' read -r num label ieee fxp; do
		for t in ieee fxp32.64; do
			knobs=$ieee
			[ "$t" = ieee ] || knobs=$fxp
			[ "$knobs" = - ] && continue
			tag=-row$num
			# The result name make precision gives this run.
			m=; [ "$SIGS" = 1 ] || m=-m$SIGS
			done_sum=$(ls "$RES"/"$t"*-logn"$logn"-r"$KEYS""$m"-s"$SEED""$tag".summary.txt \
				2>/dev/null || true)
			if [ -n "$done_sum" ]; then
				echo "skip  row $num $t logn $logn"
				continue
			fi
			flags="PREC_TARGET=$t $ZERO $knobs LOGN=$logn RUNS=$KEYS SIGS=$SIGS
				SEED=$SEED JOBS=$JOBS PREC_RESULTS=$RES TAG=$tag"
			if [ "$DRY_RUN" = 1 ]; then
				echo "make precision" $flags
				continue
			fi
			echo "$(date '+%T') row $num $t logn $logn"
			# shellcheck disable=SC2086
			if ! make precision $flags > "$RES/row$num-$t-logn$logn.log" 2>&1; then
				echo "  FAILED, see $RES/row$num-$t-logn$logn.log"
				fail=1
			fi
		done
	done < "$ROWS"
done
[ "$DRY_RUN" = 1 ] || python3 scripts/precision_tables.py --dir "$RES"
exit $fail
