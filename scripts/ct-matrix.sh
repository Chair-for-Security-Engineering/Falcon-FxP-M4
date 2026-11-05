#!/usr/bin/env bash
# Constant-time check of signing over the paper's configurations: 'make
# ctgrind' (test_ct under Valgrind memcheck) per configuration, each in its
# own copy of the tree, P at a time. Defaults are the paper's run: gcc -O3
# without vectorization (like the M4 build), 100 keys x 100 signatures per
# degree. The paper's second run is CT_OPT=-O0, where every C condition
# becomes a jump that memcheck checks. ieee runs its integer emulation.
#
#   [OUT=results/ct] [P=14] [CT_KEYS=100] [CT_SIGS=100] [CT_CC=gcc] \
#       [CT_OPT=-O3] scripts/ct-matrix.sh [grep -E filter]
set -u
cd "$(dirname "$0")/.."
if ! command -v valgrind >/dev/null; then
	exec nix-shell -I nixpkgs="$PWD/nixpkgs.nix" shell.nix --run "$(printf '%q ' "$0" "$@")"
fi
OUT=${OUT:-results/ct}
P=${P:-14}
CT_KEYS=${CT_KEYS:-100}
CT_SIGS=${CT_SIGS:-100}
CT_CC=${CT_CC:-gcc}
CT_OPT=${CT_OPT:--O3}
case $CT_CC in
*clang*) CT_NOVEC='-fno-vectorize -fno-slp-vectorize' ;;
*) CT_NOVEC='-fno-tree-vectorize' ;;
esac
FILTER=${1:-}
ROOT=$PWD
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$OUT"
rm -f "$OUT/configs.txt"

# The paper's board configurations (configs.mk) except ieee-base-origshake
# (natively the same build as ieee-base), plus the fxp split row without
# the 64-bit l10 and ties-up.
configs() {
	cat <<EOF
ieee-base|CONFIG=ieee-base
ieee-base-fftnew|CONFIG=ieee-base-fftnew
ieee-split|CONFIG=ieee-split
ieee-rec|CONFIG=ieee-rec
ieee-rec-fftnew|CONFIG=ieee-rec-fftnew
fxp-split|CONFIG=fxp-split
fxp-unfused|CONFIG=fxp-unfused
fxp-rec|CONFIG=fxp-rec
fxp-rec-fftnew|CONFIG=fxp-rec-fftnew
fxp-nonsplit|CONFIG=fxp-nonsplit
fxp-nonsplit-tiesup|CONFIG=fxp-nonsplit-tiesup
fxp-nonsplit-tiesup-fftnew|CONFIG=fxp-nonsplit-tiesup-fftnew
fxp-split-nol10|BACKEND=fxp32.64 L10_64=0 TIES_UP=0
EOF
}

run() { # run <name> <knobs>
	local name=$1 knobs=$2 dir=$WORK/$1 extra=
	rsync -a --exclude .git --exclude pqm4 --exclude results \
		--exclude '*.o' "$ROOT"/ "$dir"/
	case $knobs in
	*CONFIG=ieee*) extra='-DFNDSA_AVX2=0 -DFNDSA_SSE2=0' ;;
	esac
	extra="$extra $CT_OPT $CT_NOVEC"
	# shellcheck disable=SC2086
	(cd "$dir" && make -s ctgrind CC="$CT_CC" $knobs EXTRA_CFLAGS="$extra" \
		CT_KEYS="$CT_KEYS" CT_SIGS="$CT_SIGS") > "$OUT/$name.log" 2>&1
	echo $? > "$OUT/$name.rc"
}

n=0
while IFS='|' read -r name knobs; do
	[ -n "$FILTER" ] && ! grep -Eq -- "$FILTER" <<<"$name|$knobs" && continue
	echo "$name|$knobs" >> "$OUT/configs.txt"
	run "$name" "$knobs" &
	n=$((n + 1))
	if [ "$n" -ge "$P" ]; then
		wait -n
		n=$((n - 1))
	fi
done < <(configs)
wait

# Summary: per configuration, pass/fail and the distinct source sites of
# the reported errors (first frame inside the library).
fail=0
{
	echo "ct-matrix $(date -Is) keys=$CT_KEYS sigs=$CT_SIGS cc=$CT_CC" \
		"flags=$CT_OPT $CT_NOVEC" \
		"tree=$(cat TREE_ID 2>/dev/null || git rev-parse --short HEAD)"
	while IFS='|' read -r name knobs; do
		log=$OUT/$name.log
		rc=$(cat "$OUT/$name.rc" 2>/dev/null || echo 99)
		errs=$(sed -n 's/.*ERROR SUMMARY: \([0-9]*\) errors.*/\1/p' \
			"$log" | tail -1)
		if [ "$rc" = 0 ] && [ "${errs:-x}" = 0 ]; then
			echo "PASS  $name  ($knobs)"
		else
			echo "FAIL  $name  ($knobs) rc=$rc errors=${errs:-?}"
			fail=$((fail + 1))
			awk '/^==[0-9]+== (Conditional|Use of uninit|Syscall|Invalid)/ {
				kind = $0; sub(/^==[0-9]+== /, "", kind); getline;
				site = $0; sub(/^==[0-9]+== +(at|by) 0x[0-9A-F]+: /,
					"", site);
				print "      " kind " @ " site }' "$log" |
				sort | uniq -c | sort -rn
		fi
	done < "$OUT/configs.txt"
} | tee "$OUT/summary.txt"
[ "$fail" = 0 ]
