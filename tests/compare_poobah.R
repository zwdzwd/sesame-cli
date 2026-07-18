#!/usr/bin/env Rscript
# Level-3 (P): sesame's pOOBAH vs R's, on a fixed raw SigDF.
#
#   Rscript tests/compare_poobah.R <platform> <prefix> <c_masked.txt>
#
# Not bit-identical to R, and cannot be: the published .cm background mask and
# ordering are a newer lineage than sesameData (see NUMERICS.md). What this gate
# pins is that pOOBAH has NO ALGORITHM error -- every probe where C and R
# disagree must be *explainable*:
#
#   (a) a boundary flip: R's detection p is within `tol` of the 0.05 cutoff, so
#       a tiny lineage-induced shift in the background pool moves it across; or
#   (b) a D2 fix: R's p is NA (a channel is all-NA), which R silently leaves
#       unmasked and sesame deliberately masks (p := 1).
#
# A disagreement with R p far from 0.05 and not NA would be a real bug and fails.

suppressMessages(library(sesame))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: compare_poobah.R <platform> <prefix> <c_masked.txt>")
platform <- args[1]; prefix <- args[2]; cfile <- args[3]
tol <- 1e-3

sdf <- suppressWarnings(readIDATpair(prefix, platform = platform, min_beads = NULL))
pv  <- pOOBAH(sdf, return.pval = TRUE)
rmask <- names(pv)[!is.na(pv) & pv > 0.05]
cmask <- readLines(cfile)
tag <- sprintf("%s/%s", platform, basename(prefix))

dis <- union(setdiff(rmask, cmask), setdiff(cmask, rmask))
# a disagreeing probe must be at the boundary or a D2 NA case
p_dis <- pv[dis]
explained <- is.na(p_dis) | abs(p_dis - 0.05) < tol
unexplained <- dis[!explained]

nb <- sum(!is.na(p_dis) & abs(p_dis - 0.05) < tol)
nd2 <- sum(is.na(p_dis))

if (length(unexplained) == 0) {
    cat(sprintf("ok   %-34s R=%d C=%d masked; %d disagree (all explained: %d boundary, %d D2)\n",
                tag, length(rmask), length(cmask), length(dis), nb, nd2))
    quit(status = 0)
}

cat(sprintf("FAIL %s: %d disagreements NOT explained by boundary/D2\n",
            tag, length(unexplained)))
for (p in head(unexplained, 5))
    cat(sprintf("     %s: R p = %s\n", p, format(pv[p])))
quit(status = 1)
