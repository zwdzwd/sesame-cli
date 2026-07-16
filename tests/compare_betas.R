#!/usr/bin/env Rscript
# Differential test: sesame betas vs R, with no preprocessing (prep="").
#
#   Rscript tests/compare_betas.R <platform> <idat_prefix> <sesame_out.f64>
#
# The R side is readIDATpair(min_beads=NULL) + getBetas(), which is what
# openSesame(prefix, prep="", func=getBetas) reduces to with preprocessing off.
# Both sides compute
#     Inf-I G : max(MG,1)/max(MG+UG,2)
#     Inf-I R : max(MR,1)/max(MR+UR,2)
#     Inf-II  : max(UG,1)/max(UG+UR,2)
# then NA out masked probes. Identical uint16 inputs and a correctly-rounded
# division mean the results must be *bit-identical*, not merely close. We
# assert exactly that.
#
# The C side must be written with --f64 (raw little-endian float64).
# Do NOT compare via text: R's parser does not correctly round 17-significant-
# digit decimals -- as.numeric("0.96236179722418258") returns a double one ULP
# away from the value C printed, even though C's own strtod round-trips it
# exactly. Comparing text therefore manufactures a ~1e-16 "divergence" across
# ~19% of probes that does not exist in the computation.

suppressMessages(library(sesame))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: compare_betas.R <platform> <prefix> <c_out.f64>")
platform <- args[1]; prefix <- args[2]; cfile <- args[3]

sdf <- suppressWarnings(readIDATpair(prefix, platform = platform, min_beads = NULL))
rb  <- unname(getBetas(sdf))   # mask=TRUE by default
cb  <- readBin(cfile, "double", n = length(rb), size = 8, endian = "little")

ok <- TRUE
tag <- sprintf("%s/%s", platform, basename(prefix))

if (length(rb) != length(cb)) {
    cat(sprintf("FAIL %s: length R=%d C=%d\n", tag, length(rb), length(cb)))
    quit(status = 1)
}

# NA pattern must agree exactly.
if (!identical(is.na(rb), is.na(cb))) {
    n <- sum(is.na(rb) != is.na(cb))
    cat(sprintf("FAIL %s: NA pattern differs at %d probes (R %d NA, C %d NA)\n",
                tag, n, sum(is.na(rb)), sum(is.na(cb))))
    ok <- FALSE
}

# Values must be bit-identical.
both <- !is.na(rb) & !is.na(cb)
if (ok && !identical(rb[both], cb[both])) {
    d <- abs(rb[both] - cb[both])
    cat(sprintf("FAIL %s: %d/%d betas differ, max|diff|=%.3g\n",
                tag, sum(d > 0), sum(both), max(d)))
    ok <- FALSE
}

if (ok)
    cat(sprintf("ok   %-38s %7d betas bit-identical (%d NA)\n",
                tag, length(rb), sum(is.na(rb))))

quit(status = if (ok) 0 else 1)
