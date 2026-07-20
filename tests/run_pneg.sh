#!/bin/sh
# detectionPnegEcdf: the negative-control ECDF detection p-value. Compared on the
# RAW signal (prep ""), where it is channel-agnostic (pmax over each colour), so
# the handful of pre-C raw Infinium-I channel disagreements do not enter. Gate:
# bit-identical (<=1 ULP) to R's detectionPnegEcdf(sdf, return.pval=TRUE) over
# every probe R reports. C additionally emits p=1 for all-NA control probes that
# R's SigDF omits (positional ordering carries them) -- those are allowed.
#
# Needs a store with the platform ordering, the R oracle, and IDATs; skips clean.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
dump="$root/pipeline_dump"
store=${SESAME_INDEX_DIR:-$root/data}
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$dump" ] || { echo "FAIL: $dump not built (make pipeline_dump)"; exit 1; }
command -v Rscript >/dev/null 2>&1 || { echo "SKIP pneg: no Rscript"; exit 0; }

PASS=0; FAIL=0
run_one() {
    plat=$1; rel=$2
    ord="$store/$plat/$plat.ordering.tsv.gz"
    pfx="$idats/$rel"
    [ -f "$ord" ] || { echo "SKIP $plat pneg: no $ord"; return; }
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "SKIP $plat pneg: no IDAT $pfx"; return; fi

    SESAME_INDEX_DIR="$store" "$dump" --prep "" --what pval --detection pneg "$pfx" \
        2>/dev/null > "$work/c_pneg.txt"

    if Rscript --vanilla - "$plat" "$pfx" "$work/c_pneg.txt" <<'PY' 2>"$work/r.err"
suppressMessages(library(sesame))
a <- commandArgs(TRUE); plat<-a[1]; pfx<-a[2]; cf<-a[3]
sdf <- readIDATpair(pfx, platform=plat)
rp  <- detectionPnegEcdf(sdf, return.pval=TRUE)
cb  <- read.table(cf, colClasses=c("character","numeric")); cP<-setNames(cb$V2,cb$V1)
ids <- intersect(names(rp), names(cP))
both<- ids[!is.na(rp[ids]) & !is.na(cP[ids])]
mx  <- if (length(both)) max(abs(rp[both]-cP[both])) else 0
# any probe R reports (non-NA) that C left NA is a real mismatch:
conly <- ids[!is.na(rp[ids]) & is.na(cP[ids])]
cat(sprintf("ok   %-8s pneg: max|diff|=%.2e over %d probes (C-extra NA=%d)\n",
            plat, mx, length(both), length(conly)))
if (mx > 1e-9 || length(conly) != 0) { cat("FAIL: pneg diverges\n"); quit(status=1) }
PY
    then PASS=$((PASS+1)); else sed 's/^/    /' "$work/r.err" | head -4; FAIL=$((FAIL+1)); fi
}

run_one EPICv2 EPICv2/206909630040_R03C01
run_one MSA    MSA/207760740030_R01C03
run_one EPIC   EPIC/GSM2995280_201868590258_R01C01
run_one HM450  HM450/3999492009_R01C01

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
