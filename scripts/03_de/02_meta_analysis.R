#!/usr/bin/env Rscript
# Meta-analyze per-cohort pseudobulk age-DE across cohorts (REML random-effects),
# restricted to secretome genes. This is the statistical BACKBONE: no single
# cohort has the donor count for a credible age effect, so power comes from
# pooling effect sizes here.
#
# Input : results/de/de_<COHORT>_per_celltype.csv
#         cols: gene, celltype, log2FoldChange, lfcSE, pvalue, padj, cohort, n_donors
# Refs  : refs/secretome_union.csv  (gene, is_core, category)
# Output: results/de/secretome_meta.csv
#
# Run: PROJ=<dir> MIN_COHORT=2 conda run -n sgz_r Rscript scripts/03_de/02_meta_analysis.R
suppressMessages({library(tidyverse); library(metafor)})

PROJ   <- Sys.getenv("PROJ", unset = "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
de_dir <- file.path(PROJ, "results", "de")
files  <- list.files(de_dir, pattern = "^de_.*_per_celltype\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No per-cohort DE files in ", de_dir)
message("Cohort DE files: ", paste(basename(files), collapse = ", "))

de  <- map_dfr(files, read_csv, show_col_types = FALSE)
sec <- read_csv(file.path(PROJ, "refs", "secretome_union.csv"), show_col_types = FALSE)
de  <- de %>% filter(gene %in% sec$gene, !is.na(log2FoldChange), !is.na(lfcSE), lfcSE > 0)

MIN_COHORT <- as.integer(Sys.getenv("MIN_COHORT", unset = "2"))  # triage: 2 ; final: 4

meta <- de %>%
  group_by(gene, celltype) %>%
  filter(n_distinct(cohort) >= MIN_COHORT) %>%
  group_modify(~{
    fit <- tryCatch(rma(yi = .x$log2FoldChange, sei = .x$lfcSE, method = "REML"),
             error = function(e) tryCatch(rma(yi = .x$log2FoldChange, sei = .x$lfcSE, method = "DL"),
                                           error = function(e) NULL))
    if (is.null(fit)) return(tibble())
    n_pos <- sum(.x$log2FoldChange > 0); n_neg <- sum(.x$log2FoldChange < 0)
    tibble(meta_b = fit$b[1, 1], meta_se = fit$se, meta_p = fit$pval,
           meta_ci_lb = fit$ci.lb, meta_ci_ub = fit$ci.ub, meta_I2 = fit$I2,
           n_cohort = n_distinct(.x$cohort),
           n_sign_concordant = max(n_pos, n_neg),
           direction = ifelse(fit$b[1, 1] > 0, "up", "down"))
  }) %>% ungroup() %>%
  mutate(meta_padj = p.adjust(meta_p, method = "BH")) %>%
  left_join(sec %>% select(gene, is_core, category), by = "gene") %>%
  arrange(meta_padj)

write_csv(meta, file.path(de_dir, "secretome_meta.csv"))

robust <- meta %>% filter(meta_padj < 0.05, meta_I2 < 50, n_sign_concordant == n_cohort)
cat(sprintf("meta rows: %d | padj<0.05: %d | robust(I2<50 & fully concordant): %d\n",
            nrow(meta), sum(meta$meta_padj < 0.05, na.rm = TRUE), nrow(robust)))
cat("\nTop 15 by meta_padj:\n")
print(head(meta %>% select(gene, celltype, meta_b, meta_p, meta_padj, meta_I2, n_cohort, category), 15))
