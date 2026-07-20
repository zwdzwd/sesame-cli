#!/bin/sh
# imputeBetasByGenomicNeighbors: `impute --method neighbors` fills each missing
# MAPPED probe with the mean of its nearest same-strand non-missing genomic
# neighbours (manifest-mapped coords). Compared to R's imputeBetasByGenomicNeighbors
# on C's own betas: every MAPPED missing probe must match (imputed-or-not, value
# within float32), C must never impute one R leaves NA, and the only probes R
# imputes that C does not must be UNMAPPED ('*') -- which C deliberately leaves NA
# (R's shared-0-position artifact). Needs a store with the ordering + mapcoord
# asset, yame, pipeline_dump, R, and an IDAT.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
dump="$root/pipeline_dump"
yame="$root/YAME/yame"
store=${SESAME_INDEX_DIR:-$root/data}
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ]  || { echo "FAIL: $bin not built"; exit 1; }
[ -x "$dump" ] || { echo "FAIL: $dump not built (make pipeline_dump)"; exit 1; }
[ -x "$yame" ] || { echo "SKIP neighbors: no $yame"; exit 0; }
command -v Rscript >/dev/null 2>&1 || { echo "SKIP neighbors: no Rscript"; exit 0; }

plat=EPICv2
ord="$store/$plat/$plat.ordering.tsv.gz"
mc="$store/$plat/$plat.hg38.mapcoord.tsv.gz"
pfx="$idats/$plat/206909630040_R03C01"
[ -f "$ord" ] || { echo "SKIP neighbors: no $ord"; exit 0; }
[ -f "$mc" ]  || { echo "SKIP neighbors: no $mc (run export_manifest_coords.R)"; exit 0; }
if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
    echo "SKIP neighbors: no IDAT $pfx"; exit 0; fi

SESAME_INDEX_DIR="$store" "$bin" preprocess --platform $plat --output beta \
    --out "$work/pp" "$pfx" 2>/dev/null
SESAME_INDEX_DIR="$store" "$dump" --prep QCDPB --what beta "$pfx" 2>/dev/null > "$work/src.txt"
"$bin" impute --method neighbors --platform $plat --coords "$mc" \
    "$work/pp/beta.cg" "$work/out.cg" 2>/dev/null
"$yame" unpack "$work/out.cg" 2>/dev/null > "$work/vals.txt"
zcat < "$ord" | tail -n +2 | cut -f1 > "$work/ids.txt"
paste "$work/ids.txt" "$work/vals.txt" > "$work/c.txt"

if Rscript --vanilla - "$plat" "$work/src.txt" "$work/c.txt" <<'PY'
suppressMessages(library(sesame)); suppressMessages(library(GenomicRanges))
a <- commandArgs(TRUE); plat<-a[1]
sv <- read.table(a[2], colClasses=c("character","numeric")); b<-setNames(sv$V2, sv$V1)
ri <- imputeBetasByGenomicNeighbors(b, plat)
cc <- read.table(a[3], colClasses=c("character","numeric")); ci<-setNames(cc$V2, cc$V1)
mft <- sesameData_getManifestGRanges(plat)
mapped <- setNames(as.character(seqnames(mft)) != "*", names(mft))
miss <- names(b)[is.na(b)]; miss <- miss[miss %in% names(ci) & miss %in% names(ri)]
mm <- miss[mapped[miss]]                          # mapped missing probes
conly <- sum(!is.na(ci[miss]) & is.na(ri[miss]))  # C imputes, R doesn't (must be 0)
ron   <- miss[is.na(ci[miss]) & !is.na(ri[miss])] # R imputes, C doesn't
ron_unmapped_only <- all(!mapped[ron])
dm <- abs(ri[mm]-ci[mm]); dm <- dm[!is.na(dm)]
mapped_namis <- sum(xor(is.na(ri[mm]), is.na(ci[mm])))
cat(sprintf("neighbors: mapped-miss=%d max|diff|=%.2e mapped-NAmis=%d Conly=%d Ronly=%d(all-unmapped=%s)\n",
    length(mm), if(length(dm))max(dm)else 0, mapped_namis, conly, length(ron), ron_unmapped_only))
if ((length(dm) && max(dm) > 2e-3) || mapped_namis != 0 || conly != 0 || !ron_unmapped_only) {
    cat("FAIL: neighbors diverges on mapped probes\n"); quit(status=1) }
PY
then echo; echo "passed 1, failed 0"; else echo; echo "passed 0, failed 1"; exit 1; fi
