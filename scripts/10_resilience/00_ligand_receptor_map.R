#!/usr/bin/env Rscript
# CONSEQUENCE (Phase A) receiver-competence resilience pipeline — STEP 00: ligand->receptor map.
# Provenance: recovered VERBATIM from origin session 027c5fdb (transcript line 773); this is the code that
#   produced the committed results/step0/ + downstream results/resilience/* (agent-verified byte-for-byte).
# env: sgz_r.  DATA-DEPENDENT (reads refs/nichenet/*) -> run from $PROJ; NOT part of `make figures`.
# Out: results/step0/consequence_ligand_receptor_edge_audit.tsv + _up_ligands/_nichenet_ligands/_receptors_prefilter.txt
# Frozen UP secreted (34) -> NicheNet ligands (22) -> distinct receptors pre-filter (76). See scripts/10_resilience/README.md.
suppressMessages({library(data.table)})
sig <- read.csv("results/validation/frozen_primary_signature.csv", stringsAsFactors=FALSE)
up <- sort(unique(sig$gene[sig$direction=="UP"]))
cat("frozen UP secreted unique symbols:", length(up), "\n")
lr <- readRDS("refs/nichenet/lr_network_human_21122021.rds")
lig_in <- intersect(up, unique(lr$from))
cat("UP recognized as NicheNet ligands:", length(lig_in), "\n")
edges <- lr[lr$from %in% lig_in, c("from","to","database","source")]
edges <- edges[!duplicated(edges), ]
recs <- sort(unique(edges$to))
cat("distinct receptors (pre feature/expr filter):", length(recs), "\n")
# save edge audit + lists
colnames(edges) <- c("ligand","receptor","database","source")
write.table(edges, "results/step0/consequence_ligand_receptor_edge_audit.tsv", sep="\t", quote=FALSE, row.names=FALSE)
writeLines(up,      "results/step0/_up_ligands.txt")
writeLines(lig_in, "results/step0/_nichenet_ligands.txt")
writeLines(recs,   "results/step0/_receptors_prefilter.txt")
cat("n edges:", nrow(edges), "\n")
cat("UP genes NOT NicheNet ligands:", paste(setdiff(up, lig_in), collapse=", "), "\n")
cat("saved: consequence_ligand_receptor_edge_audit.tsv (+ _up_ligands/_nichenet_ligands/_receptors_prefilter)\n")
