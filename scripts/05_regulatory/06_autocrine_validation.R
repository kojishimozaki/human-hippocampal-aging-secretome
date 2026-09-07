#!/usr/bin/env Rscript
# Strengthen the positive autocrine lead (CCN2/CTGF, SPON1) BEYOND the prior-based NicheNet activity:
#  (1) receptor side — does the receiver niche cell EXPRESS (and age-regulate) the ligand's receptors?
#  (2) targets — are the ligand's top NicheNet targets actually aging-UP in the receiver?
#  (3) robustness — is the ligand itself aging-UP in an independent cohort (GSE199243)?
# Expression proxy = DESeq2 baseMean (gene tested => expressed in >=50% donors).
suppressMessages({ library(Matrix) })
P <- Sys.getenv("PROJ", unset = "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
lr  <- readRDS(file.path(P, "refs/nichenet/lr_network_human_21122021.rds"))
ltm <- readRDS(file.path(P, "refs/nichenet/ligand_target_matrix_nsga2r_final.rds"))   # targets x ligands
de  <- read.csv(file.path(P, "results/de/de_GSE268609_per_celltype.csv"))
de2 <- read.csv(file.path(P, "results/de/de_GSE199243_per_celltype.csv"))             # cross-cohort glia
NICHE <- c("Astro", "Micro", "Oligo", "OPC", "Endo")
LIG <- list(CCN2 = c("Astro", "Micro"), SPON1 = c("Endo"), APOD = c("Astro", "Oligo"), IL15 = c("Micro"))

sig <- function(d, g) { r <- d[d$gene == g, ]; if (!nrow(r) || is.na(r$padj)) return(NA); r }
fmt <- function(r) if (length(r) == 1 && is.na(r)) "n.t." else sprintf("lfc=%+.2f padj=%.3f base=%.0f", r$log2FoldChange, r$padj, r$baseMean)

for (L in names(LIG)) {
  cat("\n=================== LIGAND:", L, "===================\n")
  # (3) ligand aging direction, primary + cross-cohort
  for (ct in LIG[[L]]) {
    cat(sprintf("  [ligand@%s 268609] %s\n", ct, fmt(sig(de[de$celltype==ct,], L))))
    r2 <- sig(de2[de2$celltype==ct,], L); cat(sprintf("  [ligand@%s 199243 x-cohort] %s\n", ct, fmt(r2)))
  }
  # (1) receptors expressed/regulated in the receiver(s)
  recs <- sort(unique(lr$to[lr$from == L]))
  cat("  receptors (lr_network):", paste(recs, collapse=", "), "\n")
  for (ct in LIG[[L]]) {
    expr <- character(0)
    for (R in recs) { r <- sig(de[de$celltype==ct,], R); if (!(length(r)==1 && is.na(r))) expr <- c(expr, sprintf("%s(%s)", R, fmt(r))) }
    cat(sprintf("  [receptors expressed@%s] %s\n", ct, if (length(expr)) paste(expr, collapse="; ") else "none tested/expressed"))
  }
  # (2) top NicheNet targets that are aging-UP-sig in the receiver
  if (L %in% colnames(ltm)) {
    top <- names(sort(ltm[, L], decreasing = TRUE))[1:250]
    for (ct in LIG[[L]]) {
      d <- de[de$celltype == ct & !is.na(de$padj), ]
      upsig <- d$gene[d$padj < 0.1 & d$log2FoldChange > 0]
      hit <- intersect(top, upsig)
      cat(sprintf("  [top-250 targets that are aging-UP-sig @%s] %d: %s\n", ct, length(hit),
                  paste(utils::head(hit, 15), collapse=", ")))
    }
  }
}
cat("\n(autocrine support = ligand aging-UP [+x-cohort] AND receptor expressed in receiver AND targets enriched for aging-UP genes)\n")
