#!/bin/sh
# Level-3 (E): dyeBiasL (linear dye bias). E is applied in its natural pipeline
# position -- AFTER channel inference (prep "CE"), where C and R agree on every
# Infinium-I channel, so the fR/fG median pools match and betas are bit-identical
# (up to float rounding). Applying E to the RAW SigDF instead would expose the
# handful of pre-C raw channel disagreements, which C resolves in the C step and
# which are unrelated to E. Compares full betas to R's
# getBetas(dyeBiasL(inferInfiniumIChannel(sdf))).
#
# Needs a store with the platform ordering, the R oracle, and IDATs; skips clean.
set -eu

## The R oracle. Overridable because the binary is not called the same
## thing everywhere: on the lab HPC plain `Rscript` is a different R
## without a usable sesame, and the one to use is Rscript-4.6.0.
RSCRIPT=${RSCRIPT:-Rscript}

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
dump="$root/pipeline_dump"
## The shared store yame fetch fills. It is keyed on the browser
## hierarchy -- species/platform -- so a platform's assets sit directly
## at $store/$plat/ and a genome build's at $store/<build>/.
yhome=${YAME_DATA_HOME:-${XDG_DATA_HOME:-$HOME/.local/share}/yame}
store=$yhome
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$dump" ] || { echo "FAIL: $dump not built (make pipeline_dump)"; exit 1; }
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "SKIP dyebiasL: no Rscript"; exit 0; }

PASS=0; FAIL=0
run_one() {
    plat=$1; rel=$2
    ord="$store/$plat/$plat.ordering.tsv.gz"
    pfx="$idats/$rel"
    [ -f "$ord" ] || { echo "SKIP $plat E: no $ord"; return; }
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "SKIP $plat E: no IDAT $pfx"; return; fi

    YAME_DATA_HOME="$yhome" "$dump" --prep CE --what beta "$pfx" \
        2>/dev/null > "$work/c_CE.txt"

    if "$RSCRIPT" --vanilla - "$plat" "$pfx" "$work/c_CE.txt" <<'PY' 2>"$work/r.err"
suppressMessages(library(sesame))
a <- commandArgs(TRUE); plat<-a[1]; pfx<-a[2]; cf<-a[3]
sdf <- readIDATpair(pfx, platform=plat)
rB  <- getBetas(dyeBiasL(inferInfiniumIChannel(sdf)))
cb  <- read.table(cf, colClasses=c("character","numeric")); cB<-setNames(cb$V2,cb$V1)
ids <- intersect(names(rB), names(cB))
d   <- abs(rB[ids]-cB[ids]); d <- d[!is.na(d)]
nam <- sum(xor(is.na(rB[ids]), is.na(cB[ids])))
mx  <- if (length(d)) max(d) else 0
cat(sprintf("ok   %-8s E: max|diff|=%.2e over %d betas, NA-mismatch=%d\n",
            plat, mx, length(d), nam))
if (mx > 1e-9 || nam != 0) { cat("FAIL: E diverges from dyeBiasL\n"); quit(status=1) }
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
