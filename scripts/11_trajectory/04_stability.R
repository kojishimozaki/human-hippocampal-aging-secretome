#!/usr/bin/env Rscript
# Phase B (TRAJECTORY) B7 stability — leave-one-donor-out + arm-adjusted-vs-not category stability.
# NOTE: per-refit ashr is computationally prohibitive (~130 refits), so the stability checks classify on the
# RAW results() LFC (fast; a documented proxy for the ashr classification — raw-vs-ashr category concordance
# is 0.83 on the full data). The PRIMARY classification (file 03) is unchanged and ashr-based. env: trajR.
suppressMessages({library(DESeq2)})
PROJ<-Sys.getenv("PROJ",getwd()); OUT<-file.path(PROJ,"results/trajectory")
NICHE<-c("Astro","Micro","Endo","OPC","Oligo")
sig<-read.csv(file.path(PROJ,"results/validation/frozen_primary_signature.csv"),stringsAsFactors=FALSE)
prim<-read.delim(file.path(OUT,"frozen60_GSE268609_hit_trajectory.tsv"))
classify<-function(A,D,T,rel=0.5,abso=0.25){
  if(!is.finite(A)) return("not_testable"); if(A<0) return("aging_direction_not_recapitulated")
  if(A<abso) return("weak_indeterminate_aging_projection"); d<-max(abso,rel*A)
  if(D>=d) return("amplified"); if(D<=-d){ if(T>0) return("AD_attenuated") else return("AD_divergent") }
  return("maintained")}
rawcat<-function(dds,genes,s){
  HY<-results(dds,contrast=c("grp","HA","YA"))[genes,"log2FoldChange"]
  AH<-results(dds,contrast=c("grp","AD","HA"))[genes,"log2FoldChange"]
  A<-s*HY; D<-s*AH; mapply(classify,A,D,A+D)}
loo<-list(); arm<-list()
for(ct in NICHE){
  cts<-as.matrix(read.delim(file.path(OUT,sprintf("pseudobulk_counts_%s.tsv.gz",ct)),row.names=1,check.names=FALSE))
  md<-read.delim(file.path(OUT,sprintf("pseudobulk_metadata_%s.tsv",ct)),row.names=1,check.names=FALSE)
  md$grp<-factor(md$grp,levels=c("YA","HA","AD")); md$arm<-factor(md$arm); cts<-cts[,rownames(md),drop=FALSE]
  fz<-sig[sig$celltype==ct,]; genes<-fz$gene[fz$gene %in% rownames(cts)]; s<-sign(fz$log2FoldChange[match(genes,fz$gene)])
  pc<-setNames(prim$trajectory_category[prim$celltype==ct],prim$gene[prim$celltype==ct])[genes]
  # ashr-PRIMARY (arm-adjusted) category vs RAW no-arm category. NOTE: this conflates TWO differences --
  # the arm covariate AND ashr-vs-raw shrinkage -- so it is NOT a clean arm-vs-no-arm comparison (report as such).
  d0<-DESeq(DESeqDataSetFromMatrix(cts,md,~grp),quiet=TRUE)
  ac<-rawcat(d0,genes,s); arm[[ct]]<-data.frame(celltype=ct,gene=genes,primary_cat=pc,noarm_rawcat=ac,agree=(pc==ac))
  # leave-one-donor-out (raw)
  ns<-setNames(rep(0L,length(genes)),genes); sa<-setNames(rep(0L,length(genes)),genes)
  for(dn in rownames(md)){
    sub<-md[setdiff(rownames(md),dn),,drop=FALSE]; if(any(table(sub$grp)<4)) next
    dds<-tryCatch(DESeq(DESeqDataSetFromMatrix(cts[,rownames(sub),drop=FALSE],sub,~arm+grp),quiet=TRUE),error=function(e)NULL)
    if(is.null(dds)) next
    cc<-rawcat(dds,genes,s); sa<-sa+as.integer(cc==pc); ns<-ns+1L
  }
  loo[[ct]]<-data.frame(celltype=ct,gene=genes,primary_cat=pc,n_loo=ns,n_same=sa,frac_stable=sa/pmax(ns,1))
  cat(sprintf("  %s: LOO %d donors; ashr-primary-vs-raw-no-arm agree %d/%d (conflates arm + shrinkage)\n",ct,as.integer(ns[1]),sum(arm[[ct]]$agree),length(genes)))
}
LOO<-do.call(rbind,loo); ARM<-do.call(rbind,arm)
write.table(LOO,file.path(OUT,"frozen60_leave_one_donor_stability.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(ARM,file.path(OUT,"frozen60_arm_adjusted_vs_not.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat(sprintf("\nLOO raw-category stability: median frac_stable=%.2f ; hits 100%%-stable=%d/%d\n",
            median(LOO$frac_stable), sum(LOO$frac_stable==1), nrow(LOO)))
cat(sprintf("ashr-primary (arm-adj) vs RAW no-arm: category agreement=%d/%d (conflates arm covariate + ashr-vs-raw; NOT clean arm-vs-no-arm)\n", sum(ARM$agree), nrow(ARM)))
