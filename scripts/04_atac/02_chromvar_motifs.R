#!/usr/bin/env Rscript
# chromVAR: TF-motif accessibility deviations YA vs HA in niche senders (268609 ATAC).
# Signac 1.17.1 lacks RunChromVAR -> call chromVAR::computeDeviations directly.
# Caches the post-AddMotifs object so reruns skip the slow 8.9G readMM + motif scan.
suppressMessages({
  library(Signac); library(Seurat); library(Matrix)
  library(JASPAR2020); library(TFBSTools); library(motifmatchr)
  library(BSgenome.Hsapiens.UCSC.hg38); library(chromVAR); library(SummarizedExperiment)
})
set.seed(42)
PROJ <- Sys.getenv("PROJ", unset = "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
PD <- file.path(PROJ, "processed/per_dataset")
cache <- file.path(PD, "GSE268609_peaks_motifs.rds")

if (file.exists(cache)) {
  cat("loading cached post-AddMotifs object...\n"); obj <- readRDS(cache)
} else {
  pk <- as(readMM(file.path(PD, "GSE268609_peaks_counts.mtx")), "CsparseMatrix")
  rownames(pk) <- readLines(file.path(PD, "GSE268609_peaks_names.txt"))
  colnames(pk) <- readLines(file.path(PD, "GSE268609_peaks_barcodes.txt"))
  meta <- read.csv(file.path(PD, "GSE268609_metadata.csv"), row.names = 1)
  meta <- meta[colnames(pk), ]
  CL <- c(Astrocytes = "Astro", Microglia = "Micro", mOli = "Oligo", OPCs = "OPC", Endothelial = "Endo")
  meta$ct <- CL[as.character(meta$Cluster)]
  keep <- meta$Group %in% c("YA", "HA") & !is.na(meta$ct)
  pk <- pk[, keep]; meta <- meta[keep, ]
  pk <- pk[grepl("^chr[0-9XY]+-", rownames(pk)), ]
  cat("cells:", ncol(pk), " peaks:", nrow(pk), "\n")
  obj <- CreateSeuratObject(CreateChromatinAssay(counts = pk, sep = c("-", "-")),
                            assay = "peaks", meta.data = meta)
  pfm <- getMatrixSet(JASPAR2020, opts = list(species = 9606, collection = "CORE"))
  obj <- AddMotifs(obj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
  saveRDS(obj, cache)
}

# --- chromVAR directly ---
counts <- GetAssayData(obj, assay = "peaks", layer = "counts")         # peaks x cells
ranges <- StringToGRanges(rownames(counts), sep = c("-", "-"))
se <- SummarizedExperiment(assays = list(counts = counts), rowRanges = ranges)
se <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)
mm <- obj[["peaks"]]@motifs@data[rownames(counts), ]                   # peaks x motifs
bg <- getBackgroundPeaks(se)
dev_obj <- computeDeviations(object = se, annotations = mm, background_peaks = bg)
dev <- assays(dev_obj)$z                                              # motif x cell
if (is.null(rownames(dev))) rownames(dev) <- colnames(mm)
mname <- unlist(obj[["peaks"]]@motifs@motif.names)[rownames(dev)]

res <- list()
for (ctype in c("Astro", "Micro", "Oligo", "OPC", "Endo")) {
  ci <- which(obj$ct == ctype)
  dn <- obj$SampleNumber[ci]; gr <- obj$Group[ci]; donors <- unique(dn)
  M <- sapply(donors, function(d) Matrix::rowMeans(dev[, ci[dn == d], drop = FALSE]))
  dg <- vapply(donors, function(d) gr[dn == d][1], character(1))
  ya <- which(dg == "YA"); ha <- which(dg == "HA")
  if (length(ya) < 3 || length(ha) < 3) next
  dz <- rowMeans(M[, ha, drop = FALSE]) - rowMeans(M[, ya, drop = FALSE])
  pv <- apply(M, 1, function(x) tryCatch(wilcox.test(x[ha], x[ya])$p.value, error = function(e) NA))
  res[[ctype]] <- data.frame(motif = rownames(dev), tf = mname, celltype = ctype,
                             dz = dz, p = pv, row.names = NULL)
}
out <- do.call(rbind, res); out$padj <- p.adjust(out$p, "BH")
write.csv(out, file.path(PROJ, "results/de/chromvar_motifs_YAvsHA.csv"), row.names = FALSE)
cat("\n=== top aging-changed TF motifs per niche celltype (YA->HA; dz>0 = up with age) ===\n")
for (ctype in c("Astro", "Micro", "Oligo", "OPC", "Endo")) {
  s <- out[out$celltype == ctype, ]; s <- s[order(s$p), ]
  cat("\n--", ctype, "--\n"); print(head(s[, c("tf", "dz", "p", "padj")], 10), row.names = FALSE)
}
cat("\nwrote results/de/chromvar_motifs_YAvsHA.csv\n")
