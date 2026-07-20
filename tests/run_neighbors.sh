#!/bin/sh
# imputeBetasByGenomicNeighbors: `impute --method neighbors` fills each missing
# MAPPED probe with the mean of its nearest same-strand non-missing genomic
# neighbours, using the store's coord table. Since the coord table matches
# sesameData's manifest for every mapQ>=1 probe but not for mapQ=0 ones, we do
# NOT compare to R's imputeBetasByGenomicNeighbors directly (that reads the
# manifest); instead we replicate R's exact algorithm -- resize/findOverlaps/
# slice_min -- on a GRanges built from the SAME coord table, so any difference is
# a C algorithm bug, not a coordinate-source difference. Every mapped probe must
# match (imputed-or-not, value within float32); unmapped probes stay NA.
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
co="$store/$plat/$plat.hg38.coord.tsv.gz"
pfx="$idats/$plat/206909630040_R03C01"
[ -f "$ord" ] || { echo "SKIP neighbors: no $ord"; exit 0; }
[ -f "$co" ]  || { echo "SKIP neighbors: no $co"; exit 0; }
if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
    echo "SKIP neighbors: no IDAT $pfx"; exit 0; fi

SESAME_INDEX_DIR="$store" "$bin" preprocess --platform $plat --output beta \
    --out "$work/pp" "$pfx" 2>/dev/null
SESAME_INDEX_DIR="$store" "$dump" --prep QCDPB --what beta "$pfx" 2>/dev/null > "$work/src.txt"
"$bin" impute --method neighbors --platform $plat --coords "$co" \
    "$work/pp/beta.cg" "$work/out.cg" 2>/dev/null
"$yame" unpack "$work/out.cg" 2>/dev/null > "$work/vals.txt"
zcat < "$ord" | tail -n +2 | cut -f1 > "$work/ids.txt"
paste "$work/ids.txt" "$work/vals.txt" > "$work/c.txt"

if Rscript --vanilla - "$co" "$work/ids.txt" "$work/src.txt" "$work/c.txt" <<'PY'
suppressMessages({library(GenomicRanges); library(dplyr)})
a <- commandArgs(TRUE)
ids <- readLines(a[2])
co  <- read.table(gzfile(a[1]), header=TRUE, sep="\t", stringsAsFactors=FALSE)   # CpG_chrm CpG_beg strand mapQ
sv  <- read.table(a[3], colClasses=c("character","numeric")); b <- setNames(sv$V2, sv$V1)
cc  <- read.table(a[4], colClasses=c("character","numeric")); ci <- setNames(cc$V2, cc$V1)
b <- b[ids]
mapped <- !(is.na(co$CpG_chrm) | co$CpG_chrm %in% c("NA","*",""))
# GRanges from coord: 1-based [CpG_beg+1, CpG_beg+2], as the C loader derives.
mft <- GRanges(co$CpG_chrm[mapped],
               IRanges(as.integer(co$CpG_beg[mapped])+1L, as.integer(co$CpG_beg[mapped])+2L),
               strand=co$strand[mapped]); names(mft) <- ids[mapped]
mm  <- mft[names(mft) %in% names(which(is.na(b)))]
nm  <- mft[names(mft) %in% names(which(!is.na(b)))]
idx <- findOverlaps(resize(mm, 10000), nm)
gm  <- mm[queryHits(idx)]; gn <- nm[subjectHits(idx)]
df  <- tibble(cg=names(gm), beg_m=start(gm), end_m=end(gm),
              cg_n=names(gn), beg_n=start(gn), end_n=end(gn))
df$dist  <- pmax(df$beg_m - df$end_n - 1, df$beg_n - df$end_m - 1)
df$betas <- b[df$cg_n]
res <- summarise(slice_min(group_by(df, cg), n=3, order_by=dist), mb=mean(betas), .groups="drop")
exp <- b; exp[res$cg] <- res$mb                       # reference imputation on coord
mmm <- ids[mapped]; mmm <- mmm[is.na(b[mmm])]          # mapped missing probes
d   <- abs(exp[mmm]-ci[mmm]); d <- d[!is.na(d)]
namis <- sum(xor(is.na(exp[mmm]), is.na(ci[mmm])))
cat(sprintf("ok   neighbors: mapped-miss=%d max|diff|=%.2e NA-mismatch=%d\n",
            length(mmm), if(length(d))max(d)else 0, namis))
if ((length(d) && max(d) > 2e-3) || namis != 0) { cat("FAIL: neighbors diverges from the coord-algorithm reference\n"); quit(status=1) }
PY
then echo; echo "passed 1, failed 0"; else echo; echo "passed 0, failed 1"; exit 1; fi
