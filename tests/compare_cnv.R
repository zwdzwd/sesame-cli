#!/usr/bin/env Rscript
# Tight oracle for `sesame cnv`'s numeric core (the OLS fit + log2 ratio), in the
# spirit of compare_dml.R: feed R the SAME per-probe totals C uses -- extracted
# from C's own .cg files with `attach-probe` -- so there is no data-lineage or
# channel-assignment ambiguity, and R's lm must agree with C's Householder QR to
# ~machine precision. (The separate question of whether those totals match R's
# readIDATpair is a documented col-lineage divergence, not a numeric one.)
#
# Inputs (all Probe_ID-keyed TSVs from `sesame attach-probe`):
#   target : Probe_ID <TAB> total                    (total_intensity.cg)
#   normals: Probe_ID <TAB> A_M <TAB> A_U <TAB> ...   (cnvnormals.cg, --all)
#   coords : Probe_ID <TAB> chrom <TAB> pos <TAB> ... (coord.tsv.gz)
#
#   Rscript compare_cnv.R <target.tsv> <normals.tsv> <coords.tsv> <out.tsv>

args <- commandArgs(trailingOnly = TRUE)
tgtf <- args[1]; nrmf <- args[2]; crdf <- args[3]
outfile <- if (length(args) >= 4) args[4] else "cnv_oracle.tsv"

rd <- function(f) read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                             quote = "", comment.char = "", check.names = FALSE)

tg <- rd(tgtf); nm <- rd(nrmf); cd <- rd(crdf)
names(tg)[1:2] <- c("Probe_ID", "total")
names(cd)[1:2] <- c("Probe_ID", "chrom")

## normals: sum each sample's _M/_U pair into a total per normal
ncols <- names(nm)[-1]
mcol <- ncols[grepl("_M$", ncols)]
Ntot <- sapply(mcol, function(mc) {
    uc <- sub("_M$", "_U", mc)
    as.numeric(nm[[mc]]) + as.numeric(nm[[uc]])
})
rownames(Ntot) <- nm$Probe_ID

## align everything on Probe_ID; keep probes with target, all normals, and a coord
target <- setNames(as.numeric(tg$total), tg$Probe_ID)
coord  <- setNames(cd$chrom, cd$Probe_ID)
pb <- Reduce(intersect, list(names(target)[!is.na(target)],
                             rownames(Ntot)[complete.cases(Ntot)],
                             names(coord)[!is.na(coord) & coord != "" & coord != "*"]))

y <- target[pb]
X <- Ntot[pb, , drop = FALSE]
fit <- lm(y ~ ., data = data.frame(y = y, X = X))
signal <- log2(y / pmax(predict(fit), 1))

df <- data.frame(Probe_ID = pb, cnv = as.numeric(signal), stringsAsFactors = FALSE)
write.table(df, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %s: %d probe signals (%d normals)\n", outfile, nrow(df), ncol(Ntot)))
