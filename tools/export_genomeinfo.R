#!/usr/bin/env Rscript
# Export a genome's genome-level annotation from sesameData, into one folder:
#   <outdir>/<genome>/seqinfo.tsv.gz    chrom<TAB>length            (binning: NEEDED)
#   <outdir>/<genome>/gaps.tsv.gz       chrom<TAB>start<TAB>end     (binning: NEEDED)
#   <outdir>/<genome>/cytoband.tsv.gz   chrom start end band stain (plot-only)
#
# seqinfo + gaps drive getBinCoordinates (tile the genome, carve out gaps);
# cytoband is only for the ideogram in visualizeSegments and is bundled for
# completeness. Genome-level: shared across all platforms of the build.
#
# This output is NOT bundled in sesame-cli. It populates the standalone
# zhou-lab/genomes repo (one <genome>/ folder per build, so plotting tools can
# reuse it), and `sesame fetch genome <build>` pulls it into the local store.
#
#   Rscript tools/export_genomeinfo.R [genome=hg38] [outdir=data/genome]

suppressMessages(library(sesameData))

args   <- commandArgs(trailingOnly = TRUE)
genome <- if (length(args) >= 1) args[1] else "hg38"
root   <- if (length(args) >= 2) args[2] else "data/genome"
outdir <- file.path(root, genome)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

gi <- sesameData_getGenomeInfo(genome)
writegz <- function(df, name) {
    gz <- gzfile(file.path(outdir, name), "w")
    write.table(df, gz, sep = "\t", quote = FALSE, row.names = FALSE)
    close(gz)
}

sl <- gi$seqLength
if (is.null(names(sl))) {          # some builds (e.g. mm39) leave seqLength unnamed;
    si <- GenomicRanges::seqinfo(gi$gapInfo)   # recover chrom names from the gap GRanges
    sl <- stats::setNames(si@seqlengths, si@seqnames)
}
writegz(data.frame(chrom = names(sl), length = as.integer(sl)), "seqinfo.tsv.gz")

gap <- gi$gapInfo
gapdf <- data.frame(
    chrom = as.character(GenomicRanges::seqnames(gap)),
    start = GenomicRanges::start(gap),
    end   = GenomicRanges::end(gap))
gapdf <- gapdf[order(match(gapdf$chrom, names(sl)), gapdf$start), ]
writegz(gapdf, "gaps.tsv.gz")

nb <- 0L
if (!is.null(gi$cytoBand)) {
    cb <- gi$cytoBand
    cbdf <- data.frame(chrom = as.character(cb$chrom), start = cb$chromStart,
                       end = cb$chromEnd, band = cb$name, stain = cb$gieStain)
    writegz(cbdf, "cytoband.tsv.gz")
    nb <- nrow(cbdf)
}

cat(sprintf("wrote %s/: seqinfo (%d chr), gaps (%d), cytoband (%d)\n",
            outdir, length(sl), nrow(gapdf), nb))
