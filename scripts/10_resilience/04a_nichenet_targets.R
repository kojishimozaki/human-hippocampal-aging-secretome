#!/usr/bin/env Rscript
# CONSEQUENCE (Phase A) — STEP 04a: NicheNet ranked target lists (top-100 per ligand) for the downstream test.
# Provenance: VERBATIM from origin session 027c5fdb (transcript line 892). Deviation D2 (see deviation_log.md):
#   NicheNet enrichment statistic not fully pre-fixed -> both competitive Wilcoxon (04b/04c) and fgsea (04b) reported.
# env: sgz_r.  DATA-DEPENDENT (refs/nichenet/ligand_target_matrix) -> run from $PROJ.
# Out: results/resilience/_nichenet_targets_top100_ranked.tsv (input to 04b).  (stderr -> model_warnings.log.)
lig <- readLines("results/step0/_nichenet_ligands.txt")
M <- readRDS("refs/nichenet/ligand_target_matrix_nsga2r_final.rds")
ligp <- intersect(lig, colnames(M))
rows <- list(); zero<-0; tie<-0
for(l in ligp){
  v <- M[,l]; v <- v[v>0]                       # exclude 0-weight (not a target)
  o <- order(v, decreasing=TRUE); top <- names(v)[o][1:min(100,length(v))]
  pot <- v[o][1:min(100,length(top))]
  if(length(v)>=100){ if(v[o][100]==v[o][min(101,length(v))]) tie<-tie+1 }
  rows[[l]] <- data.frame(ligand=l, target=top, potential=as.numeric(pot), rank=seq_along(top))
}
zero <- sum(sapply(ligp,function(l) sum(M[,l]==0)))
long <- do.call(rbind, rows)
write.table(long,"results/resilience/_nichenet_targets_top100_ranked.tsv",sep="\t",quote=FALSE,row.names=FALSE)
cat("top-100 ranked targets per ligand saved; ligands:",length(ligp),
    "| ties at rank-100 boundary (logged):",tie,"| 0-weight entries excluded from targets (context):",zero,"\n")
