#!/usr/bin/env Rscript
# receiver-map STEP 00b — extract ligand->receptor edges from the SHA-pinned NicheNet rds and
# CROSS-CHECK the rds-derived receptor set against the manifest-authoritative 76 (step 00).
# Pre-registration: docs/RECEIVER_MAP_PLAN.md, tag prereg/receiver-map-v1.  env: sgz_r.
# The manifest remains authoritative for the 76; this step provides edges (n_ligands + LOO readiness)
# and an independent re-derivation as a provenance cross-check. Mismatch => quit(status=3) (STOP).
suppressMessages({})
P <- Sys.getenv("PROJ", unset = "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
setwd(P)
OUT <- "results/receiver"

lig <- read.csv(file.path(OUT, "ligands_22.csv"), stringsAsFactors = FALSE)$ligand
man76 <- sort(read.csv(file.path(OUT, "_universe76_from_manifest.csv"), stringsAsFactors = FALSE)$receptor)

lr <- readRDS("refs/nichenet/lr_network_human_21122021.rds")
edges <- unique(lr[lr$from %in% lig, c("from", "to")])
colnames(edges) <- c("ligand", "receptor")
rds_recs <- sort(unique(edges$receptor))

cat("STEP 00b receptor-edges\n")
cat("  ligands in:", length(lig), " edges:", nrow(edges), " distinct rds receptors:", length(rds_recs), "\n")

if (!setequal(rds_recs, man76)) {
  cat("  !! MISMATCH rds-derived vs manifest-76\n")
  cat("    in rds not manifest:", paste(setdiff(rds_recs, man76), collapse = ", "), "\n")
  cat("    in manifest not rds:", paste(setdiff(man76, rds_recs), collapse = ", "), "\n")
  writeLines(c(paste("rds_only:", paste(setdiff(rds_recs, man76), collapse = ",")),
               paste("manifest_only:", paste(setdiff(man76, rds_recs), collapse = ","))),
             file.path(OUT, "_step00b_DISCREPANCY.txt"))
  quit(status = 3)
}
cat("  cross-check PASS: rds-derived receptor set == manifest 76 (setequal)\n")

nlig <- as.data.frame(table(edges$receptor), stringsAsFactors = FALSE)
colnames(nlig) <- c("receptor", "n_ligands_connecting")
nlig <- nlig[order(-nlig$n_ligands_connecting, nlig$receptor), ]

write.csv(edges[order(edges$ligand, edges$receptor), ], file.path(OUT, "ligand_receptor_edges.csv"), row.names = FALSE)
write.csv(nlig, file.path(OUT, "frozen_receptor_universe_76.csv"), row.names = FALSE)
cat("  wrote:", file.path(OUT, "ligand_receptor_edges.csv"), "and frozen_receptor_universe_76.csv\n")
cat("  receptors mapping >=3 ligands:", paste(head(nlig$receptor[nlig$n_ligands_connecting >= 3], 12), collapse = ", "), "\n")
