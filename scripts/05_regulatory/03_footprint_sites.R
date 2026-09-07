#!/usr/bin/env Rscript
# FP step 2: exact genomic positions of the short-list TF motifs within accessible peaks, for
# donor-level differential footprinting (Stage 2). Short list locked in PHASE_REPORT.md.
suppressMessages({
  library(TFBSTools); library(JASPAR2020); library(motifmatchr)
  library(BSgenome.Hsapiens.UCSC.hg38); library(GenomicRanges); library(Signac)
})
set.seed(42)
P <- Sys.getenv("PROJ", unset = "/home/neurofuture/Bioanalysis/01.SGZ_aging_project"); PD <- file.path(P, "processed/per_dataset")
genome <- BSgenome.Hsapiens.UCSC.hg38
TFs <- c("RELA", "NFKB1", "CEBPB", "STAT3", "IRF1", "TEAD1", "TEAD4")

pfm <- getMatrixSet(JASPAR2020, list(species = 9606, collection = "CORE"))
nm <- vapply(as.list(pfm), TFBSTools::name, character(1))
sel <- pfm[nm %in% TFs]
cat("selected motifs:", paste(sprintf("%s(%s)", vapply(as.list(sel), TFBSTools::name, character(1)), names(sel)), collapse = ", "), "\n")

peaks <- StringToGRanges(readLines(file.path(PD, "GSE268609_peaks_names.txt")), sep = c("-", "-"))
peaks <- peaks[seqnames(peaks) %in% paste0("chr", c(1:22, "X", "Y"))]
seqlevels(peaks) <- seqlevelsInUse(peaks)
cat("accessible peaks:", length(peaks), "\n")

pos <- matchMotifs(sel, peaks, genome = genome, out = "positions")   # GRangesList by motif id
# collapse multiple motif-ids per TF into one GRanges of centers per TF
sites <- list()
for (tf in TFs) {
  ids <- names(sel)[vapply(as.list(sel), TFBSTools::name, character(1)) == tf]
  gr <- unlist(GRangesList(pos[ids]))
  if (length(gr) == 0) { cat("  ", tf, ": 0 sites\n"); next }
  ctr <- resize(gr, width = 1, fix = "center")
  ctr <- unique(ctr)
  sites[[tf]] <- ctr
  cat(sprintf("  %-6s : %d motif sites within peaks\n", tf, length(ctr)))
}
saveRDS(sites, file.path(P, "processed/footprint_tmp/motif_sites.rds"))
cat("wrote processed/footprint_tmp/motif_sites.rds\n")
