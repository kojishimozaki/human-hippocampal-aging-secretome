#!/usr/bin/env Rscript
# Phase 2 (regulatory, orthogonal) — NicheNet-style ligand activity: do aging-UP niche secretome
# LIGANDS best predict the aging transcriptional response in niche cells (autocrine/paracrine)?
# Prior-based, HYPOTHESIS-GENERATING (not causal). Computes NicheNet's standard ligand-activity
# (Pearson of ligand->target regulatory potential vs the binary aging response) WITHOUT nichenetr.
#
# DEVIATION from the PHASE1_5 pre-specified plan (logged): receiver (ii) NSC/Neuroblast is DROPPED because
# the niche DE output has no NSC/Neuroblast contrast (rare cells, underpowered, not load-bearing by
# design). Receivers restricted to the 5 niche cell types. Reason = data availability, not result-chasing.
suppressMessages({ library(Matrix) })
set.seed(42)
P <- Sys.getenv("PROJ", unset = "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
lr  <- readRDS(file.path(P, "refs/nichenet/lr_network_human_21122021.rds"))
ltm <- readRDS(file.path(P, "refs/nichenet/ligand_target_matrix_nsga2r_final.rds"))   # targets x ligands
de  <- read.csv(file.path(P, "results/de/de_GSE268609_per_celltype.csv"))
sec <- read.csv(file.path(P, "refs/secretome_union.csv"))$gene
NICHE <- c("Astro", "Micro", "Oligo", "OPC", "Endo")
ligands_all <- intersect(colnames(ltm), unique(lr$from))

# candidate ligands = aging-UP secretome ligands (sig up in >=1 niche celltype) that are NicheNet ligands
up <- de[!is.na(de$padj) & de$padj < 0.1 & de$log2FoldChange > 0 & de$gene %in% sec, ]
cand <- intersect(unique(up$gene), ligands_all)
cat("aging-UP secretome ligands (NicheNet senders):", length(cand), "\n  ", paste(sort(cand), collapse = ", "), "\n\n")

rows <- list()
for (rc in NICHE) {
  d <- de[de$celltype == rc & !is.na(de$padj), ]
  background <- intersect(d$gene, rownames(ltm))
  geneset <- intersect(d$gene[d$padj < 0.1], background)    # aging response (both directions)
  if (length(geneset) < 5) { cat("skip receiver", rc, "(<5 response genes mapped)\n"); next }
  resp <- as.integer(background %in% geneset)
  M <- ltm[background, , drop = FALSE]
  pcc <- as.numeric(cor(M, resp, method = "pearson")); names(pcc) <- colnames(ltm)
  la <- data.frame(ligand = names(pcc), pearson = pcc)
  la$z <- (la$pearson - mean(la$pearson)) / sd(la$pearson)  # vs all-ligand distribution (hub-biased)
  la$pct <- rank(la$pearson) / nrow(la)
  # geneset-PERMUTATION null (controls hub-ligand bias): does the ligand predict the REAL aging
  # response better than random same-size genesets? empirical p over 2000 label permutations.
  nperm <- 2000; ge <- colSums(matrix(0, 0, ncol(M)))       # placeholder
  perm_ge <- matrix(0, nrow = nperm, ncol = ncol(M))
  for (i in seq_len(nperm)) perm_ge[i, ] <- cor(M, sample(resp), method = "pearson")
  la$perm_p <- (colSums(sweep(perm_ge, 2, pcc, FUN = ">=")) + 1) / (nperm + 1)
  la$perm_padj <- p.adjust(la$perm_p, "BH")                 # over all 1226 ligands
  la$is_aging_up_secretome_ligand <- la$ligand %in% cand
  la$receiver <- rc
  rows[[rc]] <- la
}
out <- do.call(rbind, rows)
write.csv(out[order(out$receiver, -out$pearson), ], file.path(P, "results/regulatory/nichenet_ligand_activity.csv"), row.names = FALSE)

cat("=== aging-UP secretome ligands: activity + geneset-permutation null, per receiver ===\n")
cat("   (z = vs all 1226 ligands [hub-biased]; perm_p/perm_padj = vs random same-size genesets)\n")
for (rc in NICHE) {
  s <- out[out$receiver == rc & out$is_aging_up_secretome_ligand, ]
  if (!nrow(s)) next
  s <- s[order(s$perm_p), ]
  cat(sprintf("\n-- receiver %s --\n", rc))
  print(utils::head(s[, c("ligand", "pearson", "z", "perm_p", "perm_padj")], 8), row.names = FALSE)
}
sig <- out[out$is_aging_up_secretome_ligand & out$perm_padj < 0.10, ]
cat("\n=== aging-UP secretome ligands with PERMUTATION FDR<0.10 (specific, not hub-bias) ===\n")
if (nrow(sig)) print(sig[order(sig$perm_padj), c("receiver","ligand","pearson","z","perm_p","perm_padj")], row.names = FALSE) else cat("  NONE survive permutation FDR.\n")
cat("\nwrote results/regulatory/nichenet_ligand_activity.csv\n")
