#!/bin/sh
# GCT (bisConversionControl): the qc.tsv GCT column must match R's
# bisConversionControl(sdf). Uses the store's <plat>.typeI_ext.tsv.gz extension
# lists (EPIC/HM450 only; NA elsewhere, as in R). Needs a store with the ordering
# + mask + typeI_ext asset, the R oracle, and IDATs; skips clean otherwise.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
store=${SESAME_INDEX_DIR:-$root/data}
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }
command -v Rscript >/dev/null 2>&1 || { echo "SKIP gct: no Rscript"; exit 0; }

PASS=0; FAIL=0
run_one() {
    plat=$1; rel=$2
    ord="$store/$plat/$plat.ordering.tsv.gz"
    ext="$store/$plat/$plat.typeI_ext.tsv.gz"
    pfx="$idats/$rel"
    [ -f "$ord" ] || { echo "SKIP $plat GCT: no $ord"; return; }
    [ -f "$ext" ] || { echo "SKIP $plat GCT: no $ext (run export_typeI_ext.R)"; return; }
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "SKIP $plat GCT: no IDAT $pfx"; return; fi

    SESAME_INDEX_DIR="$store" "$bin" preprocess --platform "$plat" --output qc \
        --out "$work/pp" "$pfx" 2>/dev/null
    cgct=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="GCT")c=i} NR==2{print $c}' "$work/pp/qc.tsv")

    if Rscript --vanilla - "$plat" "$pfx" "$cgct" <<'PY' 2>"$work/r.err"
suppressMessages(library(sesame))
a <- commandArgs(TRUE); plat<-a[1]; pfx<-a[2]; cg<-as.numeric(a[3])
sdf <- readIDATpair(pfx, platform=plat)
rg  <- bisConversionControl(sdf)
d   <- abs(rg - cg)
cat(sprintf("ok   %-6s GCT: R=%.7f C=%.7f |diff|=%.2e\n", plat, rg, cg, d))
if (d > 1e-4) { cat("FAIL: GCT diverges\n"); quit(status=1) }
PY
    then PASS=$((PASS+1)); else sed 's/^/    /' "$work/r.err" | head -3; FAIL=$((FAIL+1)); fi
}

run_one EPIC  EPIC/GSM2995280_201868590258_R01C01
run_one HM450 HM450/3999492009_R01C01

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
