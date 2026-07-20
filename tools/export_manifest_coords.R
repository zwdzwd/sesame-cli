#!/usr/bin/env Rscript
# Export the empirically-mapped manifest coordinates per probe, positional to a
# platform's ordering table, for imputeBetasByGenomicNeighbors (which in R uses
# sesameData_getManifestGRanges -- a mapping that can differ from the design
# coordinates, e.g. remapped nv- variant probes).
#
#   Rscript tools/export_manifest_coords.R <platform> <ordering.tsv.gz> <out.tsv.gz>
#
# Output: gzipped TSV, header "chrm beg end strand", one row per probe in ordering
# order. beg/end are 1-based inclusive (GRanges start/end); unmapped probes get
# chrm "NA", beg/end "NA", strand "*".
suppressMessages(library(sesame))
suppressMessages(library(GenomicRanges))

a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 3) stop("usage: export_manifest_coords.R <platform> <ordering.tsv.gz> <out.tsv.gz>")
platform <- a[1]; ordp <- a[2]; outp <- a[3]

ord <- read.table(gzfile(ordp), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ids <- ord$Probe_ID
mft <- sesameData_getManifestGRanges(platform)
idx <- match(ids, names(mft))

df <- data.frame(
    chrm   = ifelse(is.na(idx), "NA", as.character(seqnames(mft))[idx]),
    beg    = ifelse(is.na(idx), "NA", as.character(start(mft)[idx])),
    end    = ifelse(is.na(idx), "NA", as.character(end(mft)[idx])),
    strand = ifelse(is.na(idx), "*",  as.character(strand(mft))[idx]),
    stringsAsFactors = FALSE)

con <- gzfile(outp, "w")
write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
cat(sprintf("wrote %s: %d mapped, %d unmapped (%d probes)\n",
            outp, sum(!is.na(idx)), sum(is.na(idx)), length(ids)))
