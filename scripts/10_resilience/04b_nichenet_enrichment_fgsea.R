#!/usr/bin/env Rscript
# CONSEQUENCE (Phase A) — STEP 04b (part 2/2): fgsea sensitivity (ranked-list) for the NicheNet union sets.
# Provenance: VERBATIM from origin session 027c5fdb (transcript line 896, R part). Run AFTER 04b_nichenet_enrichment_build.py.
# env: sgz_r.  Run from $PROJ. seed 42 (rule 10).  Out: results/resilience/nichenet_method_sensitivity.tsv.  (stderr -> model_warnings.log.)
suppressMessages(library(fgsea))
rk <- read.delim("results/resilience/_sadminusres_ranking.tsv", row.names=1); stats<-setNames(rk[[1]],rownames(rk))
det <- readLines("results/resilience/_ctrl_detectable_genes.txt")
T <- read.delim("results/resilience/_nichenet_targets_top100_ranked.tsv")
out<-list()
for(N in c(20,50,100)){
  tN<-T[T$rank<=N & T$target %in% det,]
  union<-unique(tN$target[tN$target %in% names(stats)])
  set.seed(42)
  fg<-fgsea(pathways=list(UNION=union), stats=stats, minSize=5, maxSize=2000, nproc=1)
  out[[as.character(N)]]<-data.frame(topN=N,n=length(union),NES=fg$NES,pval=fg$pval,padj=fg$padj)
}
res<-do.call(rbind,out)
write.table(res,"results/resilience/nichenet_method_sensitivity.tsv",sep="\t",quote=FALSE,row.names=FALSE)
cat("=== fgsea UNION (ranked-list) sensitivity ===\n"); print(res)
