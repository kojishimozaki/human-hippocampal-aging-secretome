#!/usr/bin/env Rscript
# Phase 2 (regulatory) — matched-null TF-motif enrichment in the aging niche secretome's linked
# peaks, with 3-layer convergence. Plan locked in results/regulatory/PHASE1_5_PRE_REGISTRATION.md.
# Tests the pre-specified hypothesis: UP arm ~ SASP TFs (NF-kB/C-EBP/AP-1/STAT/IRF/GATA/TEAD/SMAD),
# DOWN arm ~ homeostatic/vascular TFs (ETS/Wnt-TCF/SOX/KLF/FOXO). Reports actual arm + null honestly.
suppressMessages({
  library(Signac); library(Matrix); library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
})
set.seed(42)
P <- Sys.getenv("PROJ", unset = "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
PD <- file.path(P, "processed/per_dataset")
dir.create(file.path(P, "results/regulatory"), showWarnings = FALSE, recursive = TRUE)
genome <- BSgenome.Hsapiens.UCSC.hg38

# ---- pre-specified candidate panel ----
UP_TF <- c("NFKB1","NFKB2","RELA","RELB","REL","CEBPA","CEBPB","CEBPD","CEBPG",
           "JUN","JUNB","JUND","FOS","FOSB","FOSL1","FOSL2","BATF","BATF3",
           "STAT1","STAT2","STAT3","STAT1::STAT2","IRF1","IRF2","IRF3","IRF8","SPI1",
           "GATA2","GATA4","GATA6","TEAD1","TEAD3","TEAD4",
           "SMAD2","SMAD3","SMAD4","SMAD2::SMAD3::SMAD4","RUNX2","RUNX3")
DOWN_TF <- c("ETS1","ETV1","ETV2","ETV4","ETV6","ERG","FLI1","ELK1","ELK3","ELK4","ELF1","EHF","FEV",
             "TCF7","TCF7L1","TCF7L2","LEF1","SOX18","KLF2","KLF4","KLF5","SP1","FOXO3","FOXO4")
pred_dir <- c(setNames(rep("UP", length(UP_TF)), UP_TF), setNames(rep("DOWN", length(DOWN_TF)), DOWN_TF))
PANEL_TF <- c(UP_TF, DOWN_TF)

# ---- motif x peak matrix ----
cat("loading cached motif object...\n")
obj <- readRDS(file.path(PD, "GSE268609_peaks_motifs.rds"))
mm <- obj[["peaks"]]@motifs@data                       # peaks x motifs (binary)
mnames <- unlist(obj[["peaks"]]@motifs@motif.names)    # motif-id -> TF name
panel_ids <- names(mnames)[mnames %in% PANEL_TF]
cat("panel motifs available:", length(panel_ids), "for", length(unique(mnames[panel_ids])), "TFs\n")
peaks <- StringToGRanges(rownames(mm), sep = c("-", "-"))
peaks <- trim(peaks)

# ---- per-peak GC + width (for matched null) ----
cat("computing per-peak GC...\n")
sq <- getSeq(genome, peaks)
gc <- as.numeric(letterFrequency(sq, "GC", as.prob = TRUE))
wd <- width(peaks)
gc_bin <- cut(gc, quantile(gc, seq(0, 1, length.out = 21), na.rm = TRUE), include.lowest = TRUE, labels = FALSE)
w_bin  <- cut(log10(wd), quantile(log10(wd), seq(0, 1, length.out = 6), na.rm = TRUE), include.lowest = TRUE, labels = FALSE)
stratum <- paste(gc_bin, w_bin, sep = "_")
strata_idx <- split(seq_along(peaks), stratum)

# ---- secretome gene bodies +/-2kb ----
sec <- read.csv(file.path(P, "refs/secretome_union.csv"))$gene
fe <- read.delim(gzfile(file.path(P, "raw/GSE268609_lazarov2024/GSE268609_features.tsv.gz")), header = FALSE)
fe <- fe[fe$V3 == "Gene Expression" & fe$V2 %in% sec & grepl("^chr[0-9XY]+$", fe$V4), ]
g <- GRanges(fe$V4, IRanges(pmax(1, fe$V5 - 2000), fe$V6 + 2000), symbol = fe$V2)
ov <- findOverlaps(peaks, g)
peak2gene <- data.frame(peak = queryHits(ov), gene = g$symbol[subjectHits(ov)])
gene2peaks <- split(peak2gene$peak, peak2gene$gene)

# ---- RNA DE + chromVAR ----
de <- read.csv(file.path(P, "results/de/de_GSE268609_per_celltype.csv"))
cv <- read.csv(file.path(P, "results/de/chromvar_motifs_YAvsHA.csv"))   # motif, tf, celltype, dz, p, padj
NICHE <- c("Astro", "Micro", "Oligo", "OPC", "Endo")

matched_bg <- function(fg_idx, k = 10) {
  st <- table(stratum[fg_idx])
  bg <- unlist(lapply(names(st), function(s) {
    pool <- setdiff(strata_idx[[s]], fg_idx)
    if (length(pool) == 0) return(integer(0))
    sample(pool, min(length(pool), k * st[[s]]), replace = FALSE)
  }))
  unique(bg)
}

rows <- list()
for (ct in NICHE) {
  d <- de[de$celltype == ct & !is.na(de$padj) & de$gene %in% sec, ]
  for (arm in c("UP", "DOWN")) {
    genes <- if (arm == "UP") d$gene[d$padj < 0.1 & d$log2FoldChange > 0] else d$gene[d$padj < 0.1 & d$log2FoldChange < 0]
    genes <- intersect(genes, names(gene2peaks))
    evaluable <- length(genes) >= 4
    fg <- unique(unlist(gene2peaks[genes]))
    if (!evaluable || length(fg) < 5) {
      rows[[length(rows) + 1]] <- data.frame(celltype = ct, arm = arm, n_genes = length(genes),
        n_fg_peaks = length(fg), motif = NA, tf = NA, OR = NA, p = NA, padj = NA,
        dz = NA, dz_sign_ok = NA, n_targets = NA, verdict = "NOT_EVALUABLE", pred_match = NA)
      next
    }
    bg <- matched_bg(fg, 10)
    fg_has <- mm[fg, panel_ids, drop = FALSE]
    bg_has <- mm[bg, panel_ids, drop = FALSE]
    a <- Matrix::colSums(fg_has); cc <- Matrix::colSums(bg_has)
    nfg <- length(fg); nbg <- length(bg)
    res <- data.frame(motif = panel_ids, tf = mnames[panel_ids], a = a, c = cc)
    res$OR <- (a / (nfg - a)) / pmax(cc, 0.5) * (nbg - cc)
    res$p <- mapply(function(ai, ci) fisher.test(matrix(c(ai, nfg - ai, ci, nbg - ci), 2),
                                                 alternative = "greater")$p.value, a, cc)
    res$padj <- p.adjust(res$p, "BH")
    # L2 chromVAR dz (this celltype)
    cvc <- cv[cv$celltype == ct, ]; res$dz <- cvc$dz[match(res$motif, cvc$motif)]
    # `arm` is a length-1 loop variable, so ifelse() returned a length-1 logical that was
    # recycled across all 59 motifs: every UP arm read all-TRUE and every DOWN arm all-FALSE,
    # which is why the pre-registered down-arm hypothesis "found nothing". Audit BLOCKER B-7 /
    # finding F-2-001; fixed 2026-08-26. Vectorised comparison, arm chosen once.
    res$dz_sign_ok <- if (arm == "UP") res$dz > 0 else res$dz < 0
    # L3 targets: # arm genes whose linked peaks contain the motif
    res$n_targets <- vapply(res$motif, function(mid) {
      sum(vapply(genes, function(gg) {
        pks <- gene2peaks[[gg]]; any(mm[pks, mid] > 0)
      }, logical(1)))
    }, integer(1))
    L1 <- res$padj < 0.10 & !is.na(res$padj)
    L2 <- res$dz_sign_ok %in% TRUE
    L3 <- res$n_targets >= 3
    npass <- L1 + L2 + L3
    res$verdict <- ifelse(npass == 3, "CONVERGENT", ifelse(npass == 2, "SUGGESTIVE", "ns"))
    res$pred_match <- pred_dir[res$tf] == arm
    res$celltype <- ct; res$arm <- arm; res$n_genes <- length(genes); res$n_fg_peaks <- nfg
    rows[[length(rows) + 1]] <- res[, c("celltype","arm","n_genes","n_fg_peaks","motif","tf",
                                        "OR","p","padj","dz","dz_sign_ok","n_targets","verdict","pred_match")]
  }
}
out <- do.call(rbind, rows)
out <- out[order(out$celltype, out$arm, out$padj), ]
write.csv(out, file.path(P, "results/regulatory/tf_motif_enrichment.csv"), row.names = FALSE)

conv <- out[out$verdict %in% c("CONVERGENT", "SUGGESTIVE") & !is.na(out$motif), ]
conv <- conv[order(factor(conv$verdict, c("CONVERGENT","SUGGESTIVE")), conv$padj), ]
write.csv(conv, file.path(P, "results/regulatory/tf_convergence_summary.csv"), row.names = FALSE)

cat("\n=== evaluable arms (genes>=4) ===\n")
ev <- unique(out[out$verdict != "NOT_EVALUABLE", c("celltype","arm","n_genes","n_fg_peaks")])
print(ev, row.names = FALSE)
cat("\n=== CONVERGENT / SUGGESTIVE (L1 enrichFDR<0.1, L2 chromVAR dz-sign, L3 >=3 targets) ===\n")
if (nrow(conv) == 0) cat("  NONE -> H0: no pre-specified TF-program enrichment detected at this power.\n") else
  print(conv[, c("celltype","arm","tf","OR","padj","dz","n_targets","verdict","pred_match")], row.names = FALSE)
cat("\nwrote results/regulatory/tf_motif_enrichment.csv + tf_convergence_summary.csv\n")
