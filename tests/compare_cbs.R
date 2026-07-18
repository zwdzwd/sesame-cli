#!/usr/bin/env Rscript
# Oracle for sesame__cbs (deterministic CBS): run DNAcopy::segment -- the exact
# routine sesame's cnSegmentation calls -- on real K562 copy-number bin signals,
# with sesame's parameters and a FIXED seed (cnSegmentation itself never seeds, so
# its output is otherwise irreproducible). For each chromosome dumps
#   <out>.sig   : chrom <TAB> bin_signal          (the C-side input, per chrom)
#   <out>.seg   : chrom <TAB> end_bin             (R's segment ends, 1-based, per chrom)
# so the C harness can segment each chromosome's signal and compare breakpoints.
#
#   Rscript compare_cbs.R <out_prefix>

suppressMessages({ library(sesame); library(sesameData); library(DNAcopy) })

args <- commandArgs(trailingOnly = TRUE)
outp <- if (length(args) >= 1) args[1] else "cbs"

sdfs <- sesameDataGet("EPICv2.8.SigDF")
seg  <- cnSegmentation(sdfs[["K562_206909630040_R01C01"]])   # deterministic binning
bs   <- seg$bin.signals
bc   <- seg$bin.coords
chrom <- as.character(GenomicRanges::seqnames(bc[names(bs)]))

sigf <- file(paste0(outp, ".sig"), "w")
segf <- file(paste0(outp, ".seg"), "w")
for (ch in unique(chrom)) {
    v <- as.numeric(bs[chrom == ch])
    if (length(v) < 2) next
    for (val in v) cat(ch, "\t", val, "\n", sep = "", file = sigf)

    cna <- CNA(genomdat = v, chrom = rep(ch, length(v)),
               maploc = seq_along(v), data.type = "logratio")
    set.seed(42)
    s <- segment(cna, min.width = 5, nperm = 10000, alpha = 0.001,
                 undo.splits = "sdundo", undo.SD = 2.2, verbose = 0)
    ends <- cumsum(s$output$num.mark)                 # 1-based bin ends within chrom
    for (e in ends) cat(ch, "\t", e, "\n", sep = "", file = segf)
}
close(sigf); close(segf)
cat(sprintf("wrote %s.sig and %s.seg\n", outp, outp))
