#!/bin/sh
# Level-3 (P): pOOBAH. Every C-vs-R disagreement must be a boundary flip or a D2
# NA case (see compare_poobah.R); an unexplained one is an algorithm bug.
#
# Needs a fetched store with the platform's .cm (for the background mask) and the
# R oracle. Skips cleanly if missing.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
store=${SESAME_INDEX_DIR:-$root/data}
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }

PASS=0; FAIL=0
run_one() {
    plat=$1; rel=$2
    pfx="$idats/$rel"
    [ -d "$store/$plat" ] || { echo "SKIP $plat P: no store $store/$plat"; return; }
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "SKIP $plat P: no IDAT $pfx"; return; fi

    SESAME_INDEX_DIR="$store" "$bin" betas --prep P "$pfx" 2>/dev/null \
      | python3 -c "import sys; print('\n'.join(l.split('\t')[0] for l in sys.stdin if l.rstrip('\n').split('\t')[1]=='NA'))" \
      > "$work/c_masked.txt"

    if Rscript --vanilla "$here/compare_poobah.R" "$plat" "$pfx" "$work/c_masked.txt" \
         2>"$work/r.err"; then
        PASS=$((PASS+1))
    else
        sed 's/^/    /' "$work/r.err" | head -4
        FAIL=$((FAIL+1))
    fi
}

run_one MSA MSA/207760740030_R01C03

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
