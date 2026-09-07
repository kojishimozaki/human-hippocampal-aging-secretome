#!/usr/bin/env Rscript
# CONSEQUENCE (Phase A) — STEP 03b (part 2/2): fit nlme on the robustness-variant score tables.
# Provenance: VERBATIM from origin session 027c5fdb (transcript line 849, R part). Run AFTER 03b_robustness_build.py.
# env: sgz_r.  Run from $PROJ.  Out: results/resilience/robustness_summary.tsv.  (stderr -> model_warnings.log.)
suppressMessages(library(nlme))
fit2 <- function(d){
  d$group<-factor(d$group,levels=c("RES","SAD")); d$maturation_stage<-factor(d$maturation_stage,levels=c("immature","mature"))
  contrasts(d$maturation_stage)<-contr.sum(2); d$Run<-factor(d$Run); d$donor<-factor(d$donor)
  m<-tryCatch(lme(receiver_score~group*maturation_stage,random=~1|Run/donor,data=d,method="REML"),error=function(e)NULL)
  if(is.null(m)) m<-lme(receiver_score~group*maturation_stage,random=~1|donor,data=d,method="REML")
  tt<-summary(m)$tTable; g<-grep("groupSAD$",rownames(tt)); c(est=tt[g,"Value"],se=tt[g,"Std.Error"],df=tt[g,"DF"],p=tt[g,"p-value"])
}
L<-read.delim("results/resilience/_robustness_scores_2stage.tsv")
res<-do.call(rbind,lapply(split(L,L$variant),fit2)); res<-data.frame(variant=rownames(res),res,row.names=NULL)
# 3-stage
d3<-read.delim("results/resilience/_robustness_scores_3stage.tsv")
d3$group<-factor(d3$group,levels=c("RES","SAD")); d3$maturation_stage<-factor(d3$maturation_stage); d3$Run<-factor(d3$Run); d3$donor<-factor(d3$donor)
m3<-lme(receiver_score~group*maturation_stage,random=~1|Run/donor,data=d3,method="REML"); tt3<-summary(m3)$tTable; g3<-grep("groupSAD$",rownames(tt3))
# raw pooled
dr<-read.delim("results/resilience/_robustness_rawpooled.tsv"); rp<-fit2(dr)
main<-res[res$variant %in% c("primary","ligand_balanced","drop_ECM_GF_specific"),]
loo<-res[grepl("^LOO_",res$variant),]
cat("=== robustness: stage-averaged SAD-RES ===\n")
print(format(main,digits=3))
cat(sprintf("3-stage model groupSAD: est=%.4f p=%.3f\n", tt3[g3,"Value"], tt3[g3,"p-value"]))
cat(sprintf("raw-pooled mature: est=%.4f p=%.3f\n", rp["est"], rp["p"]))
cat(sprintf("leave-one-receptor-out (n=%d): estimate range [%.3f, %.3f], all p>0.05: %s\n",
    nrow(loo), min(loo$est), max(loo$est), all(loo$p>0.05)))
summ<-rbind(main[,c("variant","est","se","df","p")],
  data.frame(variant="three_stage",est=tt3[g3,"Value"],se=tt3[g3,"Std.Error"],df=tt3[g3,"DF"],p=tt3[g3,"p-value"]),
  data.frame(variant="raw_pooled_mature",est=rp["est"],se=rp["se"],df=rp["df"],p=rp["p"]),
  data.frame(variant="LOO_receptor_min",est=min(loo$est),se=NA,df=NA,p=max(loo$p)),
  data.frame(variant="LOO_receptor_max",est=max(loo$est),se=NA,df=NA,p=min(loo$p)))
write.table(summ,"results/resilience/robustness_summary.tsv",sep="\t",quote=FALSE,row.names=FALSE)
cat("saved robustness_summary.tsv\n")
