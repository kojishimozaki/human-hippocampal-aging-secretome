#!/usr/bin/env Rscript
# Strengthen #3 motif layer: compute chromVAR deviations on ALL niche groups
# (YA/HA/MCI/AD/SA), then donor-pseudobulk t-tests for BOTH the aging (HA vs YA)
# and the higher-powered disease (AD vs HA) contrasts. Reuses cached motif matches.
suppressMessages({
  library(Signac); library(Seurat); library(Matrix); library(chromVAR)
  library(BSgenome.Hsapiens.UCSC.hg38); library(SummarizedExperiment)
})
set.seed(42)
PROJ <- Sys.getenv("PROJ", unset = "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
PD <- file.path(PROJ, "processed/per_dataset")
devcache <- file.path(PD, "GSE268609_chromvar_allniche.rds")

if (file.exists(devcache)) {
  cat("loading cached all-niche deviations...\n")
  S <- readRDS(devcache); dev <- S$dev; meta <- S$meta; mnames <- S$mnames
} else {
  cobj <- readRDS(file.path(PD, "GSE268609_peaks_motifs.rds"))
  mm <- cobj[["peaks"]]@motifs@data
  mnames <- unlist(cobj[["peaks"]]@motifs@motif.names)
  peakset <- rownames(mm); rm(cobj); gc()
  pk <- as(readMM(file.path(PD, "GSE268609_peaks_counts.mtx")), "CsparseMatrix")
  rownames(pk) <- readLines(file.path(PD, "GSE268609_peaks_names.txt"))
  colnames(pk) <- readLines(file.path(PD, "GSE268609_peaks_barcodes.txt"))
  meta <- read.csv(file.path(PD, "GSE268609_metadata.csv"), row.names = 1)
  meta <- meta[colnames(pk), ]
  CL <- c(Astrocytes = "Astro", Microglia = "Micro", mOli = "Oligo", OPCs = "OPC", Endothelial = "Endo")
  meta$ct <- CL[as.character(meta$Cluster)]
  keep <- !is.na(meta$ct) & meta$Group %in% c("YA", "HA", "MCI", "AD", "SA")
  pk <- pk[peakset, keep]; meta <- meta[keep, ]
  cat("all-niche cells:", ncol(pk), " peaks:", nrow(pk), "\n")
  se <- SummarizedExperiment(assays = list(counts = pk), rowRanges = StringToGRanges(rownames(pk), sep = c("-", "-")))
  se <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)
  bg <- getBackgroundPeaks(se)
  dev_obj <- computeDeviations(object = se, annotations = mm, background_peaks = bg)
  dev <- assays(dev_obj)$z
  if (is.null(rownames(dev))) rownames(dev) <- colnames(mm)
  saveRDS(list(dev = dev, meta = meta, mnames = mnames), devcache)
  cat("saved all-niche deviations\n")
}

test_contrast <- function(ctype, g1, g2) {
  ci <- which(meta$ct == ctype)
  dn <- meta$SampleNumber[ci]; gr <- meta$Group[ci]; donors <- unique(dn)
  M <- sapply(donors, function(d) Matrix::rowMeans(dev[, ci[dn == d], drop = FALSE]))
  dg <- vapply(donors, function(d) gr[dn == d][1], character(1))
  i1 <- which(dg == g1); i2 <- which(dg == g2)
  if (length(i1) < 3 || length(i2) < 3) return(NULL)
  dz <- rowMeans(M[, i1, drop = FALSE]) - rowMeans(M[, i2, drop = FALSE])
  p <- apply(M, 1, function(x) tryCatch(t.test(x[i1], x[i2])$p.value, error = function(e) NA))
  data.frame(tf = mnames[rownames(dev)], celltype = ctype, contrast = paste0(g1, "v", g2),
             dz = dz, p = p, padj = p.adjust(p, "BH"), row.names = NULL)
}

all <- list()
for (ct in c("Astro", "Micro", "Oligo", "OPC", "Endo")) {
  for (cc in list(c("HA", "YA"), c("AD", "HA"))) {
    r <- test_contrast(ct, cc[1], cc[2]); if (!is.null(r)) all[[length(all) + 1]] <- r
  }
}
out <- do.call(rbind, all)
write.csv(out, file.path(PROJ, "results/de/chromvar_motifs_strengthened.csv"), row.names = FALSE)

for (cc in c("HAvYA", "ADvHA")) {
  cat("\n######### contrast:", cc, "(dz>0 = up in first group) #########\n")
  for (ct in c("Astro", "Micro", "Oligo", "OPC", "Endo")) {
    s <- out[out$celltype == ct & out$contrast == cc, ]
    s <- s[order(s$padj), ]
    nsig <- sum(s$padj < 0.1, na.rm = TRUE)
    cat(sprintf("\n-- %s %s : FDR<0.1 = %d --\n", ct, cc, nsig))
    print(head(s[, c("tf", "dz", "p", "padj")], 8), row.names = FALSE)
  }
}
cat("\nwrote results/de/chromvar_motifs_strengthened.csv\n")
