#!/usr/bin/env Rscript
# A2 step 0: per-peak GC content + width for the matched-null background (GSE268609 ATAC peaks).
# GC is the field-standard covariate the prior matched-null already used; A2 adds mean-accessibility
# (computed in Python from the full matrix) + distance-to-TSS as the further covariates.
# Out: processed/per_dataset/GSE268609_peak_gc.csv  (peak, gc, width)  -- gc=NA for non-genome contigs.
suppressMessages({library(BSgenome.Hsapiens.UCSC.hg38); library(GenomicRanges)})
P <- Sys.getenv("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
G <- BSgenome.Hsapiens.UCSC.hg38
peaks <- readLines(file.path(P, "processed/per_dataset/GSE268609_peaks_names.txt"))
# parse "chr-start-end": start/end are the last two '-'-delimited fields (chr may itself contain '-')
parts <- strsplit(peaks, "-")
chr <- vapply(parts, function(z) paste(z[seq_len(length(z) - 2)], collapse = "-"), character(1))
st  <- as.integer(vapply(parts, function(z) z[length(z) - 1], character(1)))
en  <- as.integer(vapply(parts, function(z) z[length(z)], character(1)))
gc <- rep(NA_real_, length(peaks)); width <- en - st
ok <- chr %in% seqnames(G)
gr <- GRanges(chr[ok], IRanges(st[ok] + 1L, en[ok]))            # BED 0-based -> 1-based
seqs <- getSeq(G, gr)
gc[ok] <- as.numeric(letterFrequency(seqs, "GC", as.prob = TRUE)[, 1])
out <- data.frame(peak = peaks, gc = gc, width = width)
write.csv(out, file.path(P, "processed/per_dataset/GSE268609_peak_gc.csv"), row.names = FALSE)
cat(sprintf("wrote GC+width for %d peaks (%d on-genome, %d NA contigs)\n",
            length(peaks), sum(ok), sum(!ok)))
