#!/usr/bin/env Rscript
# CONSEQUENCE (Phase A) — STEP 03: PRIMARY stage-averaged SAD-RES mixed model (+ secondary contrasts).
# Provenance: recovered VERBATIM from origin session 027c5fdb (transcript line 840, the FINAL nlme version;
#   line 826 was an lme4 attempt that failed because lme4 is absent -> deviation D1). Produced the committed
#   primary_model_results.tsv (SAD-RES -0.069609, p 0.544259) + secondary_contrasts.tsv byte-for-byte.
# Deviation D1 (see results/resilience/deviation_log.md): nlme::lme(random=~1|Run/donor) is used instead of
#   lme4 (1|Run)+(1|donor) because each donor is nested in exactly one Run -> identical model; lme4/lmerTest
#   not installed. contr.sum on stage makes the groupSAD coefficient the pre-specified stage-averaged estimand.
# env: sgz_r.  Run from $PROJ (reads results/resilience/donor_stage_receptor_scores.tsv).
# Out: results/resilience/{primary_model_results.tsv, secondary_contrasts.tsv}.  (Original redirected stderr -> model_warnings.log.)
suppressMessages(library(nlme))
d <- read.delim("results/resilience/donor_stage_receptor_scores.tsv", stringsAsFactors=FALSE)
d$group <- factor(d$group, levels=c("RES","SAD"))
d$maturation_stage <- factor(d$maturation_stage, levels=c("immature","mature"))
contrasts(d$maturation_stage) <- contr.sum(2)   # so groupSAD main effect = stage-averaged SAD-RES
d$Run <- factor(d$Run); d$donor <- factor(d$donor)
cat("n obs:", nrow(d), "| donors:", nlevels(d$donor), "| Runs:", nlevels(d$Run), "\n")

m <- lme(receiver_score ~ group*maturation_stage, random = ~1|Run/donor, data=d, method="REML")
tt <- summary(m)$tTable
cat("\n=== fixed effects (nlme tTable) ===\n"); print(round(tt,4))
cat("\n=== variance components ===\n"); print(VarCorr(m))
b <- fixef(m); V <- vcov(m); nm <- names(b)
gi <- grep("groupSAD$", nm); ii <- grep("groupSAD:", nm); si <- grep("maturation_stage1$", nm)
df_grp <- tt[gi,"DF"]; df_within <- tt[si,"DF"]
ci <- function(est,se,df) c(est-qt(.975,df)*se, est+qt(.975,df)*se)
lincon <- function(cv,df){est<-sum(cv*b); se<-sqrt(t(cv)%*%V%*%cv)[1]; t<-est/se; data.frame(estimate=est,SE=se,df=df,t=t,p=2*pt(-abs(t),df),ci_lo=ci(est,se,df)[1],ci_hi=ci(est,se,df)[2])}
cv0<-setNames(numeric(length(b)),nm)
# primary = stage-averaged SAD-RES = groupSAD main (sum-coded stage)
cvP<-cv0; cvP[gi]<-1; prim<-lincon(cvP,df_grp); prim$test<-"primary_stage_averaged_SAD_minus_RES"
# stage-specific: immature(stage1=+1): groupSAD + (+1)*interaction ; mature(stage1=-1): groupSAD + (-1)*interaction
cvI<-cv0; cvI[gi]<-1; cvI[ii]<- 1; sI<-lincon(cvI,df_grp); sI$test<-"SAD_minus_RES_immature"
cvM<-cv0; cvM[gi]<-1; cvM[ii]<--1; sM<-lincon(cvM,df_grp); sM$test<-"SAD_minus_RES_mature"
# interaction term
intx<-data.frame(estimate=b[ii],SE=tt[ii,"Std.Error"],df=tt[ii,"DF"],t=tt[ii,"t-value"],p=tt[ii,"p-value"],ci_lo=NA,ci_hi=NA,test="interaction_group_x_stage")
out<-rbind(prim,sI,sM,intx)[,c("test","estimate","SE","df","t","p","ci_lo","ci_hi")]
cat("\n=== PRIMARY + key secondary (two-sided) ===\n"); print(format(out,digits=3))

# fallback donor-only (singular policy)
mf <- lme(receiver_score ~ group*maturation_stage, random = ~1|donor, data=d, method="REML")
ttf<-summary(mf)$tTable; gif<-grep("groupSAD$",rownames(ttf))
cat("\n=== fallback donor-only: stage-averaged SAD-RES ===\n")
cat(sprintf("estimate=%.4f SE=%.4f df=%.0f t=%.3f p=%.4f\n", ttf[gif,"Value"],ttf[gif,"Std.Error"],ttf[gif,"DF"],ttf[gif,"t-value"],ttf[gif,"p-value"]))

write.table(data.frame(test="primary_stage_averaged_SAD_minus_RES",estimate=prim$estimate,SE=prim$SE,df=prim$df,
  ci_lo=prim$ci_lo,ci_hi=prim$ci_hi,t=prim$t,p_two_sided=prim$p,model="group*stage, random=~1|Run/donor (REML, nlme)",
  units="z-composite SD"), "results/resilience/primary_model_results.tsv",sep="\t",quote=FALSE,row.names=FALSE)
write.table(out,"results/resilience/secondary_contrasts.tsv",sep="\t",quote=FALSE,row.names=FALSE)
# append fallback row
write.table(data.frame(test="fallback_donoronly_SAD_minus_RES",estimate=ttf[gif,"Value"],SE=ttf[gif,"Std.Error"],
  df=ttf[gif,"DF"],t=ttf[gif,"t-value"],p=ttf[gif,"p-value"],ci_lo=NA,ci_hi=NA),
  "results/resilience/secondary_contrasts.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
cat("\nsaved primary_model_results.tsv + secondary_contrasts.tsv\n")
