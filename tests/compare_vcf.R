#!/usr/bin/env Rscript
# Oracle for `sesame vcf`: R's formatVCF genotypes the SNP probes into a VCF. Run
# it on the same raw IDAT + SNP annotation the C side uses, and dump the parsed
# per-probe genotype for comparison (Probe_ID, GT, GS, PVF).
#
#   Rscript compare_vcf.R <idat_prefix> <platform> <snp.tsv.gz> <out.tsv>

suppressMessages({ library(sesame) })

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]; platform <- args[2]; snpf <- args[3]
outfile <- if (length(args) >= 4) args[4] else "vcf_oracle.tsv"

sdf  <- readIDATpair(prefix, platform = platform)          # raw
anno <- read.table(gzfile(snpf), header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
v <- formatVCF(sdf, anno)                                  # data.frame, INFO col carries GT/GS/PVF

info <- as.character(v$INFO)
pull <- function(tag) sub(paste0(".*", tag, "=([^;]*).*"), "\\1", info)
df <- data.frame(Probe_ID = pull("Probe_ID"),
                 GT = pull("GT"),
                 GS = as.integer(v$QUAL),
                 PVF = as.numeric(pull("PVF")),
                 stringsAsFactors = FALSE)
write.table(df, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %s: %d genotyped probes\n", outfile, nrow(df)))
