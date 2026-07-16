#!/usr/bin/env Rscript
# Export a platform's ordering table from sesameData to a flat TSV.
#
#   Rscript tools/export_ordering.R EPICv2 out/EPICv2.ordering.tsv.gz
#
# This is the bootstrap path: it takes the ordering table *straight out of the
# object R itself uses* (sesameDataGet("<platform>.address")$ordering), so any
# beta mismatch against R is arithmetic in the C code, never a data difference.
# That property is what makes differential testing meaningful.
#
# (The eventual production path builds the same table from the public manifest
# TSVs via sesameAnno_buildAddressFile(); see the plan. Same columns either way.)
#
# Output columns (header, tab-separated):
#   Probe_ID  M  U  col  mask
# where
#   M, U  = bead addresses; NA (Infinium-II has no M address)
#   col   = G | R | 2      (2 = Infinium-II, matching sesame's col factor)
#   mask  = 0 | 1          (the default design mask, create_default_mask()$ref_issue)
#
# ROW ORDER IS LOAD-BEARING (R/sesame.R:504 reorders the SigDF to manifest
# order). It is preserved exactly and must not be sorted downstream.

suppressMessages(library(sesame))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("usage: export_ordering.R <platform> <out.tsv[.gz]>")
}
platform <- args[1]
outfile  <- args[2]

obj <- sprintf("%s.address", platform)
if (!sesameDataHas(obj)) {
    stop(sprintf("sesameData has no '%s'", obj))
}

ord <- sesameDataGet(obj)$ordering
stopifnot(all(c("Probe_ID", "M", "U", "col") %in% colnames(ord)))

# col: factor(levels=c("G","R")) with NA marking Infinium-II. Encode NA as "2"
# to match sesame's internal SigDF col factor levels c("G","R","2").
col <- as.character(ord$col)
col[is.na(col)] <- "2"
stopifnot(all(col %in% c("G", "R", "2")))

# mask may be absent on some objects; default FALSE.
msk <- if ("mask" %in% colnames(ord)) as.logical(ord$mask) else rep(FALSE, nrow(ord))
msk[is.na(msk)] <- FALSE

out <- data.frame(
    Probe_ID = ord$Probe_ID,
    M        = ord$M,
    U        = ord$U,
    col      = col,
    mask     = as.integer(msk),
    stringsAsFactors = FALSE
)

dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
con <- if (grepl("\\.gz$", outfile)) gzfile(outfile, "wt") else file(outfile, "wt")
write.table(out, con, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE, na = "NA")
close(con)

message(sprintf("%s: %d probes -> %s", platform, nrow(out), outfile))
message(sprintf("  Infinium-I  G: %d", sum(col == "G")))
message(sprintf("  Infinium-I  R: %d", sum(col == "R")))
message(sprintf("  Infinium-II  : %d", sum(col == "2")))
message(sprintf("  masked       : %d", sum(msk)))
message(sprintf("  controls     : %d", sum(grepl("^ctl", out$Probe_ID))))
