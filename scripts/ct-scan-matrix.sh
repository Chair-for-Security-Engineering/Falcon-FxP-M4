#!/usr/bin/env bash
# Constant-time scan of the Cortex-M4 signing images: builds the pqm4 image
# of each paper configuration and runs scripts/ct_scan.py on its signing
# call graph (hard checks, plus triage against scripts/ct_scan.allow).
# Fails on any hard finding or unreviewed site; summary in $OUT/summary.txt.
#
#   [OUT=results/ct-scan] [SCHEMES="512 1024"] scripts/ct-scan-matrix.sh [grep -E filter]
set -u
cd "$(dirname "$0")/.."
OUT=${OUT:-results/ct-scan}
SCHEMES=${SCHEMES:-512 1024}
FILTER=${1:-}
mkdir -p "$OUT"
if ! command -v arm-none-eabi-gcc >/dev/null; then
	exec nix-shell -I nixpkgs="$PWD/nixpkgs.nix" shell-cm4.nix --run "$(printf '%q ' "$0" "$@")"
fi

# The paper's board configurations (configs.mk), plus the fxp split row without
# the 64-bit l10 and ties-up.
configs() {
	cat <<EOF
ieee-m4f|ieee-base|CONFIG=ieee-base
ieee-m4f|ieee-base-fftnew|CONFIG=ieee-base-fftnew
ieee-m4f|ieee-base-origshake|CONFIG=ieee-base-origshake
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
fxp3264-cm4|fxp-split-nol10|L10_64=0 TIES_UP=0
EOF
}

scan() { # scan <elf> <tag> <root...>
	local elf=$1 tag=$2 r
	shift 2
	for r in "$@"; do
		python3 scripts/ct_scan.py "$elf" --root "$r" --src . \
			--src fpr/fxp32.64 --src fpr/ieee --src pqm4/common \
			--allow scripts/ct_scan.allow > "$OUT/$tag-$r.txt" 2>&1
		echo "$? $(grep -E '^roots' "$OUT/$tag-$r.txt")" \
			"$(grep -E '^unreviewed' "$OUT/$tag-$r.txt")"
	done
}

fail=0
summary=$OUT/summary.txt
echo "ct-scan-matrix $(date -Is) tree=$(cat TREE_ID 2>/dev/null || git rev-parse --short HEAD 2>/dev/null)" > "$summary"
while IFS='|' read -r impl name knobs; do
	[ -n "$FILTER" ] && ! grep -Eq -- "$FILTER" <<<"$impl|$name|$knobs" && continue
	for s in $SCHEMES; do
		tag=$impl-$s-$name
		# shellcheck disable=SC2086
		if ! make -s pqm4 BENCH_IMPL="$impl" SCHEME="$s" $knobs \
			> "$OUT/$tag.build.log" 2>&1; then
			echo "BUILD-FAIL $tag" | tee -a "$summary"
			fail=$((fail + 1))
			continue
		fi
		elf=pqm4/elf/crypto_sign_fndsa${s}_${impl}_speed.elf
		while read -r rc rest; do
			st=PASS
			case "$rest" in *"unreviewed sites: 0"*) ;; *) st=FAIL ;; esac
			[ "$rc" = 0 ] || st=FAIL
			[ "$st" = PASS ] || fail=$((fail + 1))
			echo "$st $tag $rest" | tee -a "$summary"
		done < <(scan "$elf" "$tag" crypto_sign_signature)
	done
done < <(configs)

echo "$fail failing" | tee -a "$summary"
[ "$fail" = 0 ]
