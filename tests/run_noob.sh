#!/bin/sh
# Level-3 (B): noob. Unit-tests normExpSignal + MASS::huber against the C
# primitives (must agree to a few ULP), then compares full betas after --prep B
# to getBetas(noob(sdf)) (lineage-bounded, like Q/P/D). See compare_noob.R.
#
# Needs a fetched store with the platform's .cm, the R oracle, and the
# normexp_test unit harness (built by `make test-noob`). Skips cleanly if the
# store or IDATs are missing.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
unit="$root/normexp_test"
store=${SESAME_INDEX_DIR:-$root/data}
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ]  || { echo "FAIL: $bin not built"; exit 1; }
[ -x "$unit" ] || { echo "FAIL: $unit not built (make test-noob)"; exit 1; }

PASS=0; FAIL=0
run_one() {
    plat=$1; rel=$2
    pfx="$idats/$rel"
    [ -d "$store/$plat" ] || { echo "SKIP $plat B: no store $store/$plat"; return; }
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "SKIP $plat B: no IDAT $pfx"; return; fi

    # raw and noob betas, so the oracle can isolate B on raw-identical probes
    SESAME_INDEX_DIR="$store" "$bin" betas          "$pfx" 2>/dev/null > "$work/c_raw.txt"
    SESAME_INDEX_DIR="$store" "$bin" betas --prep B "$pfx" 2>/dev/null > "$work/c_noob.txt"

    if Rscript --vanilla "$here/compare_noob.R" "$plat" "$pfx" \
         "$work/c_noob.txt" "$unit" "$work/c_raw.txt" 2>"$work/r.err"; then
        PASS=$((PASS+1))
    else
        sed 's/^/    /' "$work/r.err" | head -8
        FAIL=$((FAIL+1))
    fi
}

run_one MSA MSA/207760740030_R01C03

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
