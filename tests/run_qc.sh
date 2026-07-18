#!/bin/sh
# Level-3 (QC): the sesameQC panel vs R's sesameQC_calcStats. Every metric must
# agree to within data-lineage scale; a formula bug would be far off. See
# compare_qc.R. Needs a fetched store with the platform's .cm and the R oracle;
# skips cleanly if missing.
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
    [ -d "$store/$plat" ] || { echo "SKIP $plat QC: no store $store/$plat"; return; }
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "SKIP $plat QC: no IDAT $pfx"; return; fi

    SESAME_INDEX_DIR="$store" "$bin" qc "$pfx" 2>/dev/null > "$work/c_qc.tsv"

    if Rscript --vanilla "$here/compare_qc.R" "$plat" "$pfx" "$work/c_qc.tsv" \
         2>"$work/r.err"; then
        PASS=$((PASS+1))
    else
        sed 's/^/    /' "$work/r.err" | head -30
        FAIL=$((FAIL+1))
    fi
}

run_one MSA MSA/207760740030_R01C03

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
