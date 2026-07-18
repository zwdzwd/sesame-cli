#!/usr/bin/env Rscript
# DML: sesame-cli `dml` vs R's DML/summaryExtractTest.
#
#   Rscript tests/compare_dml.R <sesame-binary>
#
# Unlike the preprocessing steps this has NO data-lineage caveat: DML consumes a
# betas matrix, so the same matrix + metadata in gives the same OLS out. We
# generate a small dataset, fit it both ways with ~ group + age, and require the
# estimates, t p-values, holdout F p-values, effect sizes, and BH adjustment to
# agree to ~1e-8 (a real bug would be far off). A couple of NA betas exercise the
# per-probe dropping.

suppressMessages(library(sesame))
options(warn = -1)
args <- commandArgs(trailingOnly = TRUE)
bin <- if (length(args) >= 1) args[1] else "./sesame"

set.seed(42)
nsamp <- 24L; nprobe <- 400L
group <- rep(c("A", "B"), each = nsamp %/% 2L)
age   <- round(runif(nsamp, 20, 80), 1)
meta  <- data.frame(sample = paste0("S", seq_len(nsamp)),
                    group = group, age = age, stringsAsFactors = FALSE)

B <- matrix(rnorm(nprobe * nsamp, 0.5, 0.08), nprobe, nsamp)
B[1:120, ]   <- B[1:120, ]   + outer(rep(1, 120), ifelse(group == "B", 0.12, 0)) # group signal
B[100:180, ] <- B[100:180, ] + outer(rep(0.004, 81), age - 50)                    # age signal
B <- pmin(pmax(B, 1e-3), 1 - 1e-3)
rownames(B) <- sprintf("cg%08d", seq_len(nprobe))
colnames(B) <- meta$sample
B[5, 3]  <- NA; B[5, 4] <- NA; B[200, 10] <- NA        # exercise NA dropping

dir <- tempfile("dml"); dir.create(dir)
bfile <- file.path(dir, "betas.tsv"); mfile <- file.path(dir, "meta.tsv")
write.table(cbind(Probe_ID = rownames(B), B), bfile, sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(meta, mfile, sep = "\t", quote = FALSE, row.names = FALSE)

## ---- sesame-cli ---- (system2 pastes through a shell, so quote every arg)
cout <- system2(bin, c("dml", "--betas", shQuote(bfile), "--meta", shQuote(mfile),
                       "--formula", shQuote("~ group + age"), "--threads", "2"),
                stdout = TRUE)
C <- read.table(text = paste(cout, collapse = "\n"), header = TRUE,
                sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
rownames(C) <- C$Probe_ID

## ---- R ----
smry <- DML(B, ~ group + age, meta = meta)
R <- summaryExtractTest(smry)
rownames(R) <- R$Probe_ID

probes <- intersect(rownames(R), rownames(C))
# align coefficient columns by position within each block (names differ only in
# how R make.names-mangles "(Intercept)")
Rest  <- R[probes, grep("^Est_",   colnames(R)), drop = FALSE]
Cest  <- C[probes, grep("^Est_",   colnames(C)), drop = FALSE]
Rp    <- R[probes, grep("^Pval_",  colnames(R)), drop = FALSE]
Cp    <- C[probes, grep("^Pval_",  colnames(C)), drop = FALSE]
Rf    <- R[probes, grep("^FPval_", colnames(R)), drop = FALSE]
Cf    <- C[probes, grep("^FPval_", colnames(C)), drop = FALSE]
Re    <- R[probes, grep("^Eff_",   colnames(R)), drop = FALSE]
Ce    <- C[probes, grep("^Eff_",   colnames(C)), drop = FALSE]

relmax <- function(a, b) {
    a <- as.matrix(a); b <- as.matrix(b)
    ok <- !is.na(a) & !is.na(b)
    if (!any(ok)) return(0)
    max(abs(a[ok] - b[ok]) / pmax(abs(a[ok]), 1e-8))
}
d_est <- relmax(Rest, Cest); d_p <- relmax(Rp, Cp)
d_f   <- relmax(Rf, Cf);     d_e <- relmax(Re, Ce)

## BH adjustment: C's FPadj must match p.adjust of the F p-values
Rfadj_ref <- apply(as.matrix(Cf), 2, p.adjust, method = "BH")
Cfadj     <- as.matrix(C[probes, grep("^FPadj_", colnames(C)), drop = FALSE])
d_adj <- relmax(Rfadj_ref, Cfadj)

tol <- 1e-8
worst <- max(d_est, d_p, d_f, d_e, d_adj)
cat(sprintf("DML  probes=%d coefs=%d  relmax: Est=%.2e Pval=%.2e FPval=%.2e Eff=%.2e BH=%.2e\n",
            length(probes), ncol(Cest), d_est, d_p, d_f, d_e, d_adj))
cat(sprintf("%-4s worst=%.2e (tol=%.0e); R kept %d probes, C %d\n",
            ifelse(worst < tol, "ok", "FAIL"), worst, tol, nrow(R), nrow(C)))
quit(status = ifelse(worst < tol, 0, 1))
