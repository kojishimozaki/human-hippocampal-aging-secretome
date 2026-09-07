#!/usr/bin/env Rscript
# N3 L3-B (SUPPLEMENTAL, honest exploration) — NicheNet niche->NEURON receiver.
# Adapts scripts/05_regulatory/02_nichenet.R (ligand activity + 2000x geneset-perm null) and
# 06_autocrine_validation.R (receptor expression/DE), pointing the RECEIVER at neurons. Prior-based,
# correlative NOT causal; computes NicheNet ligand activity WITHOUT nichenetr (readRDS + base R + Matrix).
#
# Question: do the niche aging-UP secretome LIGANDS signal into NEURONS? (i) are their receptors expressed
# in neurons; (ii) are those receptors age-regulated in neurons; (iii) do the ligands predict the neuron
# aging response better than random (perm null), binary (padj<0.1 response) + continuous (vs log2FC, the
# powered variant for the sparse neuron DE). Inherits the niche autocrine NicheNet's weak calibration
# (null passed ~82-91% of ligands). The receptor readout (i,ii) is the most interpretable for
# ambient-vs-response: a genuine paracrine response needs expressed receptors; ambient contamination does not.
suppressMessages({ library(Matrix) })
set.seed(42)
P <- Sys.getenv("PROJ", unset = "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
NICHE  <- c("Astro", "Micro", "Oligo", "OPC", "Endo")
NEURON <- c("DG_GC", "CA_ExN", "InN")
NPERM  <- 2000

lr  <- readRDS(file.path(P, "refs/nichenet/lr_network_human_21122021.rds"))
ltm <- readRDS(file.path(P, "refs/nichenet/ligand_target_matrix_nsga2r_final.rds"))   # targets x ligands
sec <- read.csv(file.path(P, "refs/secretome_union.csv"))$gene
frozen <- read.csv(file.path(P, "results/validation/frozen_primary_signature.csv"))
de  <- read.csv(file.path(P, "results/de/de_GSE268609_neuron_aging_per_celltype.csv"))  # niche + neuron, ~arm+grp
ligands_all <- intersect(colnames(ltm), unique(lr$from))

# ---- candidate niche aging-UP secretome ligands (frozen 60 UP + broader N0 niche UP) ----
frozen_up <- unique(frozen$gene[frozen$direction == "UP"])
niche_up  <- unique(de$gene[de$celltype %in% NICHE & !is.na(de$padj) & de$padj < 0.1 &
                            de$log2FoldChange > 0 & de$gene %in% sec])
cand        <- intersect(union(frozen_up, niche_up), ligands_all)
cand_frozen <- intersect(frozen_up, ligands_all)
cat(sprintf("niche aging-UP secretome ligands that are NicheNet ligands: %d (of which frozen-60-UP: %d)\n",
            length(cand), length(cand_frozen)))
cat("  candidates:", paste(sort(cand), collapse = ", "), "\n\n")

# ================= (iii) ligand activity in NEURON receivers =================
la_rows <- list()
for (rc in NEURON) {
  d <- de[de$celltype == rc & !is.na(de$padj), ]
  background <- intersect(d$gene, rownames(ltm))
  if (length(background) < 50) { cat("skip", rc, "(<50 background)\n"); next }
  M <- ltm[background, , drop = FALSE]
  lfc <- d$log2FoldChange[match(background, d$gene)]
  run_perm <- function(resp, label) {
    pcc <- as.numeric(cor(M, resp)); names(pcc) <- colnames(M)
    perm <- matrix(0, NPERM, ncol(M))
    for (i in seq_len(NPERM)) perm[i, ] <- cor(M, sample(resp))
    pp <- (colSums(sweep(perm, 2, pcc, ">=")) + 1) / (NPERM + 1)
    data.frame(receiver = rc, response = label, ligand = names(pcc), pearson = pcc,
               z = (pcc - mean(pcc)) / sd(pcc), perm_p = pp, perm_padj = p.adjust(pp, "BH"),
               is_candidate = names(pcc) %in% cand, is_frozen_up = names(pcc) %in% cand_frozen)
  }
  # binary response (padj<0.1, both directions) -- only if >=5 mapped
  gs <- intersect(d$gene[d$padj < 0.1], background)
  if (length(gs) >= 5) la_rows[[paste0(rc, "_bin")]] <- run_perm(as.integer(background %in% gs), "binary_padj0.1")
  else cat(sprintf("  %s: binary response NOT_EVALUABLE (only %d padj<0.1 genes mapped)\n", rc, length(gs)))
  # continuous response (log2FC) -- powered variant for sparse neuron DE
  okv <- !is.na(lfc)
  Mc <- M[okv, , drop = FALSE]; lfcv <- lfc[okv]
  pcc <- as.numeric(cor(Mc, lfcv)); names(pcc) <- colnames(Mc)
  perm <- matrix(0, NPERM, ncol(Mc)); for (i in seq_len(NPERM)) perm[i, ] <- cor(Mc, sample(lfcv))
  pp <- (colSums(sweep(perm, 2, pcc, ">=")) + 1) / (NPERM + 1)
  la_rows[[paste0(rc, "_cont")]] <- data.frame(receiver = rc, response = "continuous_lfc", ligand = names(pcc),
      pearson = pcc, z = (pcc - mean(pcc)) / sd(pcc), perm_p = pp, perm_padj = p.adjust(pp, "BH"),
      is_candidate = names(pcc) %in% cand, is_frozen_up = names(pcc) %in% cand_frozen)
}
la <- do.call(rbind, la_rows)
dir.create(file.path(P, "results/regulatory"), showWarnings = FALSE, recursive = TRUE)
write.csv(la[order(la$receiver, la$response, -la$pearson), ],
          file.path(P, "results/regulatory/neuron_nichenet_ligand_activity.csv"), row.names = FALSE)

cat("\n=== (iii) candidate niche ligands predicting the NEURON aging response (perm null) ===\n")
for (rc in NEURON) for (resp in c("binary_padj0.1", "continuous_lfc")) {
  s <- la[la$receiver == rc & la$response == resp & la$is_candidate, ]
  if (!nrow(s)) next
  nsig <- sum(s$perm_padj < 0.10, na.rm = TRUE)
  cat(sprintf("-- %s / %s: %d candidate ligands | perm_padj<0.10: %d | best: %s (perm_p=%.3f)\n",
              rc, resp, nrow(s), nsig, s$ligand[which.min(s$perm_p)], min(s$perm_p)))
}

# ================= (i)+(ii) receptors of candidate ligands, on NEURONS =================
rec_rows <- list()
for (L in cand) {
  recs <- sort(unique(lr$to[lr$from == L]))
  for (rc in NEURON) {
    for (R in recs) {
      r <- de[de$celltype == rc & de$gene == R, ]
      expressed <- nrow(r) > 0 && !is.na(r$baseMean[1])
      rec_rows[[length(rec_rows) + 1]] <- data.frame(ligand = L, receiver = rc, receptor = R,
          expressed = expressed,
          lfc = if (expressed) r$log2FoldChange[1] else NA,
          padj = if (expressed && !is.na(r$padj[1])) r$padj[1] else NA,
          age_regulated = expressed && !is.na(r$padj[1]) && r$padj[1] < 0.1)
    }
  }
}
rec <- do.call(rbind, rec_rows)
write.csv(rec, file.path(P, "results/regulatory/neuron_nichenet_receptors.csv"), row.names = FALSE)

cat("\n=== (i)+(ii) niche-ligand RECEPTORS on neurons (expression + age-regulation) ===\n")
for (rc in NEURON) {
  s <- rec[rec$receiver == rc, ]
  uniqR <- unique(s$receptor)
  exprR <- unique(s$receptor[s$expressed])
  aregR <- unique(s$receptor[s$age_regulated])
  cat(sprintf("-- %s: %d unique receptors of candidate ligands | expressed(tested): %d (%.0f%%) | age-regulated padj<0.1: %d %s\n",
              rc, length(uniqR), length(exprR), 100 * length(exprR) / max(1, length(uniqR)),
              length(aregR), if (length(aregR)) paste0("[", paste(aregR, collapse = ", "), "]") else ""))
}
cat("\nNOTE: prior-based, CORRELATIVE not causal; inherits weak NicheNet calibration; snRNA receptors sparse.\n")
cat("wrote results/regulatory/neuron_nichenet_ligand_activity.csv + neuron_nichenet_receptors.csv\n")
cat("DONE (N3 L3-B).\n")
