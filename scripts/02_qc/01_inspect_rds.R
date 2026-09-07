#!/usr/bin/env Rscript
# Generic inspector for a Seurat/SCE RDS: prints structure + metadata schema,
# writes meta.data to CSV so Python can join annotations onto the MTX counts.
# Usage: Rscript 01_inspect_rds.R <in.rds> [out_prefix]
suppressMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 01_inspect_rds.R <in.rds> [out_prefix]")
in_rds <- args[1]
out_prefix <- if (length(args) >= 2) args[2] else sub("\\.[Rr][Dd][Ss]$", "", in_rds)

cat("Reading:", in_rds, "\n")
obj <- readRDS(in_rds)
cat("Top-level class:", paste(class(obj), collapse = ", "), "\n")

if (inherits(obj, "Seurat")) {
  cat("Cells:", ncol(obj), " Features:", nrow(obj), "\n")
  cat("Assays:", paste(Assays(obj), collapse = ", "), "\n")
  cat("Default assay:", DefaultAssay(obj), "\n")
  md <- obj@meta.data
  cat("\n--- meta.data colnames ---\n"); print(colnames(md))
  cat("\n--- head(meta.data, 3) ---\n"); print(utils::head(md, 3))
  for (cn in colnames(md)) {
    v <- md[[cn]]
    if (is.factor(v) || is.character(v) || (is.numeric(v) && length(unique(v)) <= 30)) {
      u <- length(unique(v))
      if (u <= 60) {
        cat(sprintf("\n--- %s (%d levels) ---\n", cn, u))
        print(sort(table(v), decreasing = TRUE))
      }
    } else if (is.numeric(v)) {
      cat(sprintf("\n--- %s (numeric) ---\n", cn)); print(summary(v))
    }
  }
  write.csv(md, paste0(out_prefix, "_metadata.csv"), row.names = TRUE)
  cat("\nWrote:", paste0(out_prefix, "_metadata.csv"), "\n")
} else {
  cat("Not a Seurat object; structure (depth 2):\n")
  utils::str(obj, max.level = 2)
}
