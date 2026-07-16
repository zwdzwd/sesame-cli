#!/usr/bin/env Rscript
# Level-3 differential test: one prep step, applied in isolation, C vs R.
#
#   Rscript tests/compare_prep.R <platform> <prefix> <code> <c_col.tsv>
#
# Compares the resulting Infinium colour channel per probe. The channel is a
# discrete label, so the gate is exact equality -- there is no tolerance to hide
# behind. A single differing call means the port is wrong.

suppressMessages(library(sesame))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: compare_prep.R <platform> <prefix> <code> <c_col.tsv>")
platform <- args[1]; prefix <- args[2]; code <- args[3]; cfile <- args[4]

sdf <- suppressWarnings(readIDATpair(prefix, platform = platform, min_beads = NULL))

for (ch in strsplit(code, "")[[1]]) {
    sdf <- switch(ch,
        "C" = inferInfiniumIChannel(sdf),
        stop(sprintf("compare_prep.R has no oracle for code '%s'", ch)))
}

rc <- as.character(sdf$col)
cc <- read.table(cfile, sep = "\t", header = FALSE, colClasses = "character")
tag <- sprintf("%s/%s prep=%s", platform, basename(prefix), code)

if (!identical(sdf$Probe_ID, cc[[1]])) {
    cat(sprintf("FAIL %s: probe IDs/order differ\n", tag)); quit(status = 1)
}

d <- which(rc != cc[[2]])
if (length(d) == 0) {
    n <- table(rc)
    cat(sprintf("ok   %-42s %d probes, channels identical (G=%s R=%s II=%s)\n",
                tag, length(rc), n["G"], n["R"], n["2"]))
    quit(status = 0)
}

cat(sprintf("FAIL %s: %d of %d channel calls differ\n", tag, length(d), length(rc)))
for (i in head(d, 5))
    cat(sprintf("     %s: R=%s C=%s  MG=%s UG=%s MR=%s UR=%s\n",
                sdf$Probe_ID[i], rc[i], cc[[2]][i],
                sdf$MG[i], sdf$UG[i], sdf$MR[i], sdf$UR[i]))
quit(status = 1)
