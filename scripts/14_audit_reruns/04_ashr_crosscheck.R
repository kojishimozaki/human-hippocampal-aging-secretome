#!/usr/bin/env Rscript
# [F-5-006 cross-check, step 2/2] Independent shrinkage of the frozen 60 with DESeq2 + ashr.
#
# WHY. pyDESeq2 0.5.4's lfc_shrink collapsed 17 of the frozen 60 to ~0 while reporting
# "converged" for all of them; fitting the mathematically equivalent parameterisation showed
# exactly one branch failing per gene, always the branch where that gene's coefficient is
# negative. scripts/14_audit_reruns/01_covariates_and_shrinkage.py repairs this by taking the
# non-failing branch. That repair needs an implementation that does not share pyDESeq2's
# optimiser. apeglm is not installed in any env on this host, so this uses ashr -- which the
# project already adopted for the AD leg (scripts/11_trajectory/02_deseq2_ashr.R, and the
# manuscript justifies the choice at Methods L82).
#
# WHAT THIS CAN AND CANNOT SETTLE. ashr uses a different prior family, so it cannot validate
# apeGLM's numbers. It can settle whether the collapse to ~0 is implementation-dependent.
# The criteria below were fixed before this script was run
# (audit_log/2026-08-25_p0_cheap_reruns/RESOLUTION.md, section 7).
#
# env: trajR (DESeq2 1.50.2, ashr 2.2.63). apeglm is deliberately NOT installed.
# Inputs : results/audit_reruns/pseudobulk/{counts,metadata}_<ct>  (written by 03_export_pseudobulk.py)
#          results/audit_reruns/lfc_shrinkage_frozen60.csv          (pyDESeq2 MLE + repaired apeGLM)
#          results/validation/frozen_primary_signature.csv
# Out    : results/audit_reruns/ashr_crosscheck_frozen60.csv
suppressMessages({library(DESeq2); library(ashr)})

PROJ  <- Sys.getenv("PROJ", getwd())
PBDIR <- file.path(PROJ, "results/audit_reruns/pseudobulk")
OUT   <- file.path(PROJ, "results/audit_reruns")
NICHE <- c("Astro", "Micro", "Oligo", "OPC", "Endo")

py <- read.csv(file.path(OUT, "lfc_shrinkage_frozen60.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(py) == 60)
py$key <- paste(py$celltype, py$gene)

rows <- list()
for (ct in NICHE) {
  cts <- as.matrix(read.delim(file.path(PBDIR, sprintf("counts_%s.tsv.gz", ct)),
                              row.names = 1, check.names = FALSE))
  md  <- read.delim(file.path(PBDIR, sprintf("metadata_%s.tsv", ct)),
                    row.names = 1, check.names = FALSE)
  md$grp <- factor(md$grp, levels = c("YA", "HA"))     # YA is the reference -> grp_HA_vs_YA
  cts <- cts[, rownames(md), drop = FALSE]
  dds <- DESeq(DESeqDataSetFromMatrix(cts, md, ~grp), quiet = TRUE)
  res <- results(dds, contrast = c("grp", "HA", "YA"))
  sh  <- tryCatch(lfcShrink(dds, res = res, type = "ashr", quiet = TRUE),
                  error = function(e) stop(sprintf("ashr failed for %s: %s", ct,
                                                   conditionMessage(e))))
  conv <- as.logical(mcols(dds)$betaConv); conv[is.na(conv)] <- FALSE
  names(conv) <- rownames(dds)

  fz <- py[py$celltype == ct, ]
  for (i in seq_len(nrow(fz))) {
    g <- fz$gene[i]
    present <- g %in% rownames(res)
    rows[[paste(ct, g)]] <- data.frame(
      celltype = ct, gene = g, present_in_R = present,
      baseMean_R   = if (present) res[g, "baseMean"]       else NA_real_,
      lfc_MLE_R    = if (present) res[g, "log2FoldChange"] else NA_real_,
      lfcSE_MLE_R  = if (present) res[g, "lfcSE"]          else NA_real_,
      padj_R       = if (present) res[g, "padj"]           else NA_real_,
      betaConv_R   = if (present) conv[[g]]                else NA,
      lfc_ashr     = if (present) sh[g, "log2FoldChange"]  else NA_real_,
      lfcSE_ashr   = if (present) sh[g, "lfcSE"]           else NA_real_,
      stringsAsFactors = FALSE)
  }
  cat(sprintf("  %s: %d donors (%s), %d genes, ashr finite %d/%d\n", ct, ncol(cts),
              paste(table(md$grp), collapse = "/"), nrow(dds),
              sum(is.finite(sh$log2FoldChange)), nrow(sh)))
}
R <- do.call(rbind, rows)
M <- merge(py, R, by = c("celltype", "gene"), all.x = TRUE)
stopifnot(nrow(M) == 60)
M$ratio_ashr    <- M$lfc_ashr / M$lfc_MLE
M$ratio_apeglm  <- M$lfc_shrunk_repaired / M$lfc_MLE
M$suspect       <- M$altparam_absdiff > 0.05           # the 17 pyDESeq2 disagreed on
write.csv(M[order(M$celltype, M$gene), ], file.path(OUT, "ashr_crosscheck_frozen60.csv"),
          row.names = FALSE)

cat("\n=================== pre-registered checks (section 7) ===================\n")
ok <- is.finite(M$lfc_MLE_R) & is.finite(M$lfc_MLE)
g9_r    <- cor(M$lfc_MLE_R[ok], M$lfc_MLE[ok])
g9_max  <- max(abs(M$lfc_MLE_R[ok] - M$lfc_MLE[ok]))
g9 <- (g9_r > 0.999) && (g9_max < 0.05) && all(ok)
cat(sprintf("G-9 engine gate      : n=%d/60 present, Pearson r = %.6f, max |diff| = %.5f -> %s\n",
            sum(ok), g9_r, g9_max, ifelse(g9, "PASS", "FAIL")))
cat(sprintf("                       R betaConv %d/60\n", sum(M$betaConv_R %in% TRUE)))

sus <- M[M$suspect, ]
c1_bad <- sus[abs(sus$lfc_ashr) <= 0.5 * abs(sus$lfc_MLE), ]
cat(sprintf("C1 suspect genes     : %d suspect; ashr collapses (|ashr| <= 0.5|MLE|) for %d -> %s\n",
            nrow(sus), nrow(c1_bad), ifelse(nrow(c1_bad) == 0, "PASS", "FAIL")))
if (nrow(c1_bad)) print(c1_bad[, c("celltype", "gene", "lfc_MLE", "lfc_shrunk_repaired", "lfc_ashr")])

c2_r  <- cor(M$lfc_shrunk_repaired, M$lfc_ashr, use = "complete.obs")
disc  <- M[(abs(M$ratio_apeglm) > 0.8 & abs(M$ratio_ashr) < 0.4) |
           (abs(M$ratio_ashr) > 0.8 & abs(M$ratio_apeglm) < 0.4), ]
c2 <- (c2_r > 0.95) && (nrow(disc) == 0)
cat(sprintf("C2 consistency       : Pearson r(apeGLM_repaired, ashr) = %.4f, discordant genes = %d -> %s\n",
            c2_r, nrow(disc), ifelse(c2, "PASS", "FAIL")))
if (nrow(disc)) print(disc[, c("celltype", "gene", "lfc_MLE", "ratio_apeglm", "ratio_ashr")])

med <- median(M$ratio_ashr, na.rm = TRUE)
c3 <- (med > 0.6) && (med < 1.0)
cat(sprintf("C3 shrinkage level   : median ashr ratio = %.4f (range %.3f-%.3f) -> %s\n",
            med, min(M$ratio_ashr, na.rm = TRUE), max(M$ratio_ashr, na.rm = TRUE),
            ifelse(c3, "PASS", "FAIL")))
nsign <- sum(sign(M$lfc_ashr) == sign(M$lfc_MLE), na.rm = TRUE)
cat(sprintf("sign preserved       : %d/60 -> %s\n", nsign, ifelse(nsign == 60, "PASS", "FAIL")))

cat("\n--- suspect genes: does ashr agree with the repair or with the collapse? ---\n")
print(sus[order(sus$celltype, sus$gene),
          c("celltype", "gene", "lfc_MLE", "lfc_shrunk", "lfc_shrunk_altparam",
            "lfc_shrunk_repaired", "lfc_ashr", "ratio_ashr")], row.names = FALSE, digits = 4)

named <- c("SERPINE1", "IL15", "CCN2", "ANXA1", "TNC", "COL21A1", "APOD", "WNT5B", "KDR",
           "PLAT", "LAMA4", "APOC1", "LGALS9", "NPC2", "APOE", "ECM2", "IGFBP5",
           "COL12A1", "KITLG")
cat("\n--- genes named in the Results prose ---\n")
print(M[M$gene %in% named, c("celltype", "gene", "baseMean", "lfc_MLE",
                             "lfc_shrunk_repaired", "lfc_ashr", "ratio_apeglm", "ratio_ashr")],
      row.names = FALSE, digits = 4)

cat(sprintf("\nVERDICT: G-9 %s | C1 %s | C2 %s | C3 %s | sign %s\n",
            ifelse(g9, "PASS", "FAIL"), ifelse(nrow(c1_bad) == 0, "PASS", "FAIL"),
            ifelse(c2, "PASS", "FAIL"), ifelse(c3, "PASS", "FAIL"),
            ifelse(nsign == 60, "PASS", "FAIL")))
cat(sprintf("wrote %s/ashr_crosscheck_frozen60.csv\n", OUT))
