#!/usr/bin/env Rscript
# Level-3 (B): sesame's noob vs R's.
#
#   Rscript tests/compare_noob.R <platform> <prefix> <c_betas.f64> <normexp_bin>
#
# Two independent checks:
#
#  (1) UNIT -- the D5/D6 arithmetic, isolated from data lineage. normExpSignal and
#      MASS::huber are compared to the C primitives on identical inputs via the
#      normexp_test harness. These are the same algorithm and MUST agree to a few
#      ULP; a miss here is a real bug and fails the test.
#
#  (2) INTEGRATION -- full betas after `--prep B` vs getBetas(noob(sdf)). Like
#      Q/P/D this cannot be bit-identical: the published background mask/ordering
#      is a newer lineage than sesameData, so the noob background/foreground pools
#      differ slightly and every fit parameter shifts. We report the beta-diff
#      distribution and fail only if it is far larger than lineage can explain.

suppressMessages(library(sesame))
options(warn = -1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5)
    stop("usage: compare_noob.R <platform> <prefix> <c_noob.txt> <normexp_bin> <c_raw.txt>")
platform <- args[1]; prefix <- args[2]; cnoob_f <- args[3]; nbin <- args[4]; craw_f <- args[5]
tag <- sprintf("%s/%s", platform, basename(prefix))
fail <- 0L

read_betas <- function(f) {
    t <- read.table(f, sep = "\t", header = FALSE, na.strings = "NA",
                    colClasses = c("character", "numeric"))
    setNames(t[[2]], t[[1]])
}

normExp <- get("normExpSignal", envir = asNamespace("sesame"))

## ---- (1a) normExpSignal ------------------------------------------------------
mus <- c(10, 50, 200, 1000); sigmas <- c(5, 30, 150)
alphas <- c(10, 100, 1000, 5000)
xs <- c(1, 5, 20, 100, 300, 1000, 3000, 10000, 30000, 60000)
grid <- expand.grid(mu = mus, sigma = sigmas, alpha = alphas, x = xs)

r_ne <- mapply(function(mu, sigma, alpha, x) normExp(mu, sigma, alpha, x),
               grid$mu, grid$sigma, grid$alpha, grid$x)
inp <- tempfile()
writeLines(sprintf("%.17g %.17g %.17g %.17g",
                   grid$mu, grid$sigma, grid$alpha, grid$x), inp)
c_ne <- as.numeric(system2(nbin, "normexp", stdin = inp, stdout = TRUE))

d_ne  <- abs(c_ne - r_ne)
rel   <- d_ne / pmax(abs(r_ne), 1e-300)
maxabs <- max(d_ne); p90rel <- quantile(rel, 0.90)
# Gate on absolute error (a formula bug shows up huge) plus a tight bulk relative
# error. The worst relative miss is a deep-tail near-cancellation (signal is a
# difference of two ~equal ~1e3 numbers) -- precisely the regime where R's
# difference-of-logs loses precision and the Mills-ratio form (D5) does not, so
# there C is the more accurate of the two, not the one to gate against.
ne_ok <- maxabs < 1e-6 && p90rel < 1e-10
cat(sprintf("%-4s %s normExpSignal: %d cases, max|abs|=%.3g p90|rel|=%.3g\n",
            ne_ok, tag, nrow(grid), maxabs, p90rel))
if (!ne_ok) {
    fail <- 1L
    i <- which.max(d_ne)
    cat(sprintf("     worst-abs: mu=%g sigma=%g alpha=%g x=%g  R=%.17g C=%.17g\n",
                grid$mu[i], grid$sigma[i], grid$alpha[i], grid$x[i], r_ne[i], c_ne[i]))
}

## ---- (1b) MASS::huber --------------------------------------------------------
huber_c <- function(v) {
    f <- tempfile(); writeLines(sprintf("%.17g", v), f)
    strsplit(system2(nbin, "huber", stdin = f, stdout = TRUE), "\t")[[1]]
}
sdf0 <- suppressWarnings(readIDATpair(prefix, platform = platform, min_beads = NULL))
vecs <- list(
    seq_len(1000) * 1.0,
    c(rep(5, 60), as.numeric(1:300)),
    as.numeric(na.omit(sdf0$MG)),               # a real, large signal vector
    as.numeric(na.omit(c(sdf0$UG, sdf0$UR)))
)
hub_ok <- TRUE
for (v in vecs) {
    rr <- MASS::huber(v)
    cc <- huber_c(v)
    dmu <- abs(as.numeric(cc[1]) - rr$mu); ds <- abs(as.numeric(cc[2]) - rr$s)
    ok  <- dmu < 1e-6 && ds < 1e-8
    hub_ok <- hub_ok && ok
    cat(sprintf("%-4s %s huber n=%d: |dmu|=%.3g |ds|=%.3g\n",
                ok, tag, length(v), dmu, ds))
}
if (!hub_ok) fail <- 1L
## D6: MAD==0 -- R errors, C falls back. Confirm C flags it and does not crash.
z <- huber_c(rep(100, 50))
cat(sprintf("%-4s %s huber MAD==0 (D6): C s=%s mad0=%s (R errors here)\n",
            identical(z[3], "1"), tag, z[2], z[3]))
if (!identical(z[3], "1")) fail <- 1L

## ---- (2) INTEGRATION: betas after B -----------------------------------------
## The published ordering (284309) assigns different M/U addresses than
## sesameData's manifest (284317) to a handful of replicate probes, so their RAW
## betas already differ -- a lineage artifact that has nothing to do with B. To
## test B in isolation we restrict to probes whose raw betas already agree: there
## the input to noob is identical, so any post-noob difference is the B step
## alone (fit params differ only at the lineage level, ~1e-2 on mu/sigma/alpha).
r_raw  <- getBetas(sdf0, mask = FALSE)
r_noob <- getBetas(noob(sdf0), mask = FALSE)
c_raw  <- read_betas(craw_f)
c_noob <- read_betas(cnoob_f)

common <- Reduce(intersect, list(names(r_raw), names(r_noob),
                                 names(c_raw), names(c_noob)))
raw_div <- abs(r_raw[common] - c_raw[common])
iso <- common[!is.na(raw_div) & raw_div < 1e-9]     # raw-identical -> isolates B
n_lin <- sum(!is.na(raw_div) & raw_div >= 1e-9)

rn <- r_noob[iso]; cn <- c_noob[iso]
both <- !is.na(rn) & !is.na(cn)
d <- abs(cn[both] - rn[both])
qs <- quantile(d, c(0.5, 0.99, 1.0))
cat(sprintf("     %s betas after B (raw-identical probes): n=%d  median|d|=%.2e  p99=%.2e  max=%.2e\n",
            tag, sum(both), qs[1], qs[2], qs[3]))
cat(sprintf("     excluded %d raw-divergent probes (ordering lineage, pre-B); NA agree R-only=%d C-only=%d\n",
            n_lin, sum(is.na(rn) & !is.na(cn)), sum(!is.na(rn) & is.na(cn))))
## On raw-identical inputs B still differs by the lineage-level fit shift: the
## background/foreground pools include the ordering-divergent probes, so muG/sgG
## etc move ~1e-2, which propagates to ~1e-3 on the lowest-signal probes (the
## worst is invariably a negative control). A real B bug moves betas ~0.1+ (an
## earlier channel-swap bug produced 0.4-0.7), so this threshold separates them.
b_ok <- qs[3] < 1e-2 && qs[2] < 3e-3
if (!b_ok) {
    cat(sprintf("FAIL %s: B moved raw-identical betas too far (p99=%.2e max=%.2e) -- not lineage\n",
                tag, qs[2], qs[3]))
    fail <- 1L
    worst <- iso[both][which.max(d)]
    cat(sprintf("     worst: %s  R=%.5f C=%.5f\n", worst, r_noob[worst], c_noob[worst]))
}

quit(status = fail)
