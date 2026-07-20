#!/usr/bin/env Rscript
# Export the Type-I extension-base class per probe, positional to a platform's
# ordering table, for the GCT (bisulfite-conversion) QC metric.
#
#   Rscript tools/export_typeI_ext.R <platform> <ordering.tsv.gz> <out.tsv.gz>
#
# Output: a gzipped single-column file, header "ext", then one token per probe in
# ordering order -- "C" for a Type-I extension-C probe (typeI.extC), "T" for an
# extension-A/T probe (typeI.extT), "." otherwise. Only EPIC/HM450 carry the
# probeInfo extension lists in sesameData; other platforms emit all "." (GCT
# unsupported, exactly as R's bisConversionControl is out of the box).
suppressMessages(library(sesame))

a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 3) stop("usage: export_typeI_ext.R <platform> <ordering.tsv.gz> <out.tsv.gz>")
platform <- a[1]; ordp <- a[2]; outp <- a[3]

ord <- read.table(gzfile(ordp), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ids <- ord$Probe_ID
ext <- rep(".", length(ids))

pinfo <- tryCatch(sesameDataGet(paste0(platform, ".probeInfo")), error = function(e) NULL)
if (!is.null(pinfo) && !is.null(pinfo$typeI.extC)) {
    ext[ids %in% pinfo$typeI.extC] <- "C"
    ext[ids %in% pinfo$typeI.extT] <- "T"
}

con <- gzfile(outp, "w")
writeLines(c("ext", ext), con)
close(con)
cat(sprintf("wrote %s: %d C, %d T, %d none (%d probes)\n",
            outp, sum(ext == "C"), sum(ext == "T"), sum(ext == "."), length(ids)))
