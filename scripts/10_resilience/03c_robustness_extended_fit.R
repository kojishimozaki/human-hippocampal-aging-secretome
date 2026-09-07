#!/usr/bin/env Rscript
# CONSEQUENCE (Phase A) — STEP 03c (part 2/2): extended robustness (Run-fixed / no-Run / leave-one-Run /
#   leave-one-ligand-family / subtype-specific).
# Provenance: VERBATIM from origin session 027c5fdb (transcript line 883, R part). Run AFTER 03b_robustness_build.py
#   (reads _robustness_scores_2stage/3stage.tsv) and 03c_robustness_extended_build.py (reads _famdrop_2stage.tsv).
# env: sgz_r.  Run from $PROJ.  Out: results/resilience/{robustness_run_fixed, robustness_no_run,
#   robustness_leave_one_run, robustness_leave_one_ligand_family, robustness_subtype_specific}.tsv.  (stderr -> model_warnings.log.)
suppressMessages(library(nlme))
prim2 <- read.delim("results/resilience/_robustness_scores_2stage.tsv"); prim <- prim2[prim2$variant=="primary",]
fit2 <- function(d, random=~1|Run/donor){
  d$group<-factor(d$group,levels=c("RES","SAD")); d$maturation_stage<-factor(d$maturation_stage,levels=c("immature","mature"))
  contrasts(d$maturation_stage)<-contr.sum(2); d$Run<-factor(d$Run); d$donor<-factor(d$donor)
  m<-tryCatch(lme(receiver_score~group*maturation_stage,random=random,data=d,method="REML"),error=function(e)NULL)
  if(is.null(m)) return(c(est=NA,se=NA,df=NA,p=NA))
  tt<-summary(m)$tTable; g<-grep("^groupSAD$",rownames(tt)); c(est=tt[g,"Value"],se=tt[g,"Std.Error"],df=tt[g,"DF"],p=tt[g,"p-value"])
}
# no-Run
nr<-fit2(prim, random=~1|donor); write.table(data.frame(model="no_Run_donoronly",t(nr)),"results/resilience/robustness_no_run.tsv",sep="\t",quote=FALSE,row.names=FALSE)
# Run fixed (check estimability)
d<-prim; d$group<-factor(d$group,levels=c("RES","SAD")); d$maturation_stage<-factor(d$maturation_stage,levels=c("immature","mature")); contrasts(d$maturation_stage)<-contr.sum(2); d$Run<-factor(d$Run); d$donor<-factor(d$donor)
mf<-tryCatch(lme(receiver_score~group*maturation_stage+Run,random=~1|donor,data=d,method="REML"),error=function(e)NULL)
if(is.null(mf)){ rf<-data.frame(model="Run_fixed",est=NA,se=NA,df=NA,p=NA,note="not estimable / rank-deficient")
} else { tt<-summary(mf)$tTable; g<-grep("^groupSAD$",rownames(tt)); rf<-data.frame(model="Run_fixed",est=tt[g,"Value"],se=tt[g,"Std.Error"],df=tt[g,"DF"],p=tt[g,"p-value"],note="group estimable from within-Run variation (Runs 2/4/5 mixed)") }
write.table(rf,"results/resilience/robustness_run_fixed.tsv",sep="\t",quote=FALSE,row.names=FALSE)
# leave-one-Run-out
lor<-do.call(rbind,lapply(sort(unique(prim$Run)),function(rn){
  d<-prim[prim$Run!=rn,]; ng<-length(unique(d$group))
  if(ng<2) return(data.frame(dropped_Run=rn,est=NA,se=NA,df=NA,p=NA,note="only one group remains -> not estimable"))
  r<-fit2(d); data.frame(dropped_Run=rn,t(r),note=paste0("RES+SAD both remain; n_donor=",length(unique(d$donor))))}))
write.table(lor,"results/resilience/robustness_leave_one_run.tsv",sep="\t",quote=FALSE,row.names=FALSE)
# leave-one-ligand-family-out
fd<-read.delim("results/resilience/_famdrop_2stage.tsv")
lf<-do.call(rbind,lapply(split(fd,fd$variant),function(x){r<-fit2(x[,setdiff(names(x),"variant")]); data.frame(variant=x$variant[1],t(r))}))
write.table(lf,"results/resilience/robustness_leave_one_ligand_family.tsv",sep="\t",quote=FALSE,row.names=FALSE)
# subtype-specific (score ~ group + Run, donor independent; one obs/donor/subtype)
d3<-read.delim("results/resilience/_robustness_scores_3stage.tsv")
ss<-do.call(rbind,lapply(sort(unique(d3$maturation_stage)),function(s){
  x<-d3[d3$maturation_stage==s,]; x$group<-factor(x$group,levels=c("RES","SAD")); x$Run<-factor(x$Run)
  m<-tryCatch(lm(receiver_score~group+Run,data=x),error=function(e)NULL)
  if(is.null(m)) return(data.frame(subtype=s,est=NA,se=NA,df=NA,p=NA))
  co<-summary(m)$coefficients; g<-grep("^groupSAD$",rownames(co))
  if(length(g)==0) return(data.frame(subtype=s,est=NA,se=NA,df=NA,p=NA,note="groupSAD not estimable"))
  data.frame(subtype=s,est=co[g,1],se=co[g,2],df=m$df.residual,p=co[g,4])}))
write.table(ss,"results/resilience/robustness_subtype_specific.tsv",sep="\t",quote=FALSE,row.names=FALSE)
cat("=== A5 completion ===\n")
cat("no-Run:",sprintf("est=%.3f p=%.3f",nr["est"],nr["p"]),"\n")
cat("Run-fixed:",if(is.na(rf$est))"NOT ESTIMABLE" else sprintf("est=%.3f p=%.3f",rf$est,rf$p),"\n")
cat("leave-one-Run-out p-range:",sprintf("[%.3f, %.3f]",min(lor$p,na.rm=TRUE),max(lor$p,na.rm=TRUE)),"all>0.05:",all(lor$p>0.05,na.rm=TRUE),"\n")
print(lf); cat("\nsubtype-specific:\n"); print(ss)
