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
    ## Raw betas too, to find probes that already disagree before E runs --
    ## the ordering-lineage set, excluded the same way compare_noob.R does.
    YAME_DATA_HOME="$yhome" "$dump" --prep "" --what beta "$pfx" \
        2>/dev/null > "$work/c_raw.txt"

    if "$RSCRIPT" --vanilla - "$plat" "$pfx" "$work/c_CE.txt" "$work/c_raw.txt" \
         <<'PY' 2>"$work/r.err"
suppressMessages(library(sesame))
a <- commandArgs(TRUE); plat<-a[1]; pfx<-a[2]; cf<-a[3]; rawf<-a[4]
rd  <- function(f) { x <- read.table(f, colClasses=c("character","numeric"))
                     setNames(x$V2, x$V1) }
sdf <- readIDATpair(pfx, platform=plat)
rB  <- getBetas(dyeBiasL(inferInfiniumIChannel(sdf)))
cB  <- rd(cf)
ids <- intersect(names(rB), names(cB))

## The store ships InfiniumAnnotation's ordering; R reads sesameData's older
## manifest. Where the two disagree on a probe's DESIGN TYPE (Infinium-I vs II)
## M and U come from different addresses, so the betas differ from the raw stage
## on and no downstream step can reconcile them -- e.g. EPIC/HM450 cg07162498
## and cg09334382, which v8.1 calls Infinium-I and sesameData Infinium-II.
## That is an annotation-lineage difference, not a dye-bias one, so exclude it
## and test E on the probes whose raw betas already agree (cf. compare_noob.R).
## A G<->R channel flip is NOT in this set: inferInfiniumIChannel re-derives the
## channel from the data, so both sides converge on it.
rRaw <- getBetas(sdf); cRaw <- rd(rawf)
rawd <- abs(rRaw[ids] - cRaw[ids])
lin  <- (!is.na(rawd) & rawd >= 1e-9) | xor(is.na(rRaw[ids]), is.na(cRaw[ids]))
ids  <- ids[!lin]

d   <- abs(rB[ids]-cB[ids]); d <- d[!is.na(d)]
nam <- sum(xor(is.na(rB[ids]), is.na(cB[ids])))
mx  <- if (length(d)) max(d) else 0
cat(sprintf("ok   %-8s E: max|diff|=%.2e over %d betas, NA-mismatch=%d (%d lineage-divergent excluded)\n",
            plat, mx, length(d), nam, sum(lin)))
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
