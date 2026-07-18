#!/usr/bin/env Rscript
# QC panel: sesame-cli `qc` vs R sesameQC_calcStats on a raw SigDF.
#
#   Rscript tests/compare_qc.R <platform> <prefix> <c_qc.tsv>
#
# Not bit-identical, and cannot be: the published ordering/mask is a newer
# lineage than sesameData (fewer probes; different background mask), so probe
# counts, detection, and beta stats all shift a little -- the same lineage the P
# and B tests carry. What this gate pins is that every metric agrees to within
# lineage scale; a real formula bug would be far off.
#
# Two metrics are reported but NOT gated, both by design:
#   num_dtna / frac_dtna -- the D2 fix makes R's pOOBAH never return NA, so R
#   reports 0 here while sesame-cli counts probes with no signal in either
#   channel (its intended meaning). See NUMERICS.md.

suppressMessages(library(sesame))
options(warn = -1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: compare_qc.R <platform> <prefix> <c_qc.tsv>")
platform <- args[1]; prefix <- args[2]; ctsv <- args[3]
tag <- sprintf("%s/%s", platform, basename(prefix))

sdf <- suppressWarnings(readIDATpair(prefix, platform = platform, min_beads = NULL))
rq  <- unlist(sesameQC_calcStats(sdf)@stat)          # named numeric vector

ctab <- read.table(ctsv, sep = "\t", header = TRUE, check.names = FALSE,
                   stringsAsFactors = FALSE)
cv <- as.list(ctab[1, , drop = TRUE])

# tolerances by metric kind. Counts get an absolute slack: the ordering lineage
# plus channel/mask boundary moves only ~a dozen probes, which is >1% for the
# small categories (e.g. rs NAs, G->R switches) yet nowhere near a real bug,
# which would miscount by hundreds.
count_re <- "^(num_|na_intensity|InfI_switch)"        # integer-valued
exempt   <- c("num_dtna", "frac_dtna")                # D2 divergence (reported)
tol_frac <- 2e-2                                       # fractions / means / ratios
tol_cnt  <- 1e-2                                       # counts, relative ...
abs_cnt  <- 25                                         # ... or this absolute slack

metrics <- intersect(names(rq), names(cv))
fail <- 0L; worst <- ""; worst_re <- 0

cat(sprintf("%-20s %14s %14s %10s\n", "metric", "R", "C", "rel"))
for (m in metrics) {
    r <- as.numeric(rq[[m]]); c <- suppressWarnings(as.numeric(cv[[m]]))
    if (is.na(r) && is.na(c)) next
    denom <- max(abs(r), 1e-9)
    re <- abs(c - r) / denom
    is_cnt <- grepl(count_re, m)
    ok_metric <- if (is_cnt) (re <= tol_cnt || abs(c - r) <= abs_cnt)
                 else        (is.finite(re) && re <= tol_frac)
    status <- "ok"
    if (m %in% exempt) {
        status <- "D2"                                # reported, not gated
    } else if (!ok_metric) {
        status <- "FAIL"; fail <- fail + 1L
    }
    if (status == "FAIL" || (status == "ok" && re > worst_re)) {
        if (status != "D2") { worst_re <- re; worst <- m }
    }
    # only print the notable rows to keep output short
    if (status != "ok" || re > 1e-3)
        cat(sprintf("%-20s %14.6g %14.6g %10.2e  %s\n", m, r, c, re, status))
}

cat(sprintf("\n%-4s %s  %d metrics, worst gated rel=%.2e (%s)\n",
            ifelse(fail == 0, "ok", "FAIL"), tag, length(metrics), worst_re, worst))
quit(status = ifelse(fail == 0, 0, 1))
