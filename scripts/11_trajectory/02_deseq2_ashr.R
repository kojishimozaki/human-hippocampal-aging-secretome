#!/usr/bin/env Rscript
# Phase B (TRAJECTORY) B2-B5 — single R DESeq2 fit/celltype + ashr; v1.4 two-layer identity gate.
# Gate 1 = global coefficient additive identity on coef(dds) (<=1e-10). Gate 2 = frozen-60 results-level
# identity (<=1e-6) + convergence + ashr validity. The global results()-LFC non-additivity is reported
# descriptively (not a gate, not excluded). env: trajR (DESeq2 1.50.2, ashr 2.2.63). FAIL-CLOSED.
suppressMessages({library(DESeq2); library(ashr)})
PROJ <- Sys.getenv("PROJ", getwd()); OUT <- file.path(PROJ,"results/trajectory")
NICHE <- c("Astro","Micro","Endo","OPC","Oligo")
sig <- read.csv(file.path(PROJ,"results/validation/frozen_primary_signature.csv"), stringsAsFactors=FALSE)
frozenkey <- paste(sig$gene, sig$celltype)
warnlog <- file(file.path(OUT,"model_warnings.log"), open="wt")
hit_rows <- list(); coef_rows <- list(); res_rows <- list(); shr_rows <- list()

for (ct in NICHE) {
  cts <- as.matrix(read.delim(file.path(OUT,sprintf("pseudobulk_counts_%s.tsv.gz",ct)), row.names=1, check.names=FALSE))
  md  <- read.delim(file.path(OUT,sprintf("pseudobulk_metadata_%s.tsv",ct)), row.names=1, check.names=FALSE)
  md$grp <- factor(md$grp, levels=c("YA","HA","AD")); md$arm <- factor(md$arm); cts <- cts[, rownames(md), drop=FALSE]
  design <- if (nlevels(md$arm) >= 2) ~arm + grp else ~grp
  mm <- model.matrix(design, md); full_rank <- (qr(mm)$rank == ncol(mm))
  dds <- DESeq(DESeqDataSetFromMatrix(cts, md, design), quiet=TRUE)
  g <- rownames(dds); bm <- mcols(dds)$baseMean; names(bm) <- g
  conv <- as.logical(mcols(dds)$betaConv); conv[is.na(conv)] <- FALSE; names(conv) <- g
  cf <- coef(dds); bHA <- cf[,"grp_HA_vs_YA"]; bAD <- cf[,"grp_AD_vs_YA"]
  coef_HY <- bHA; coef_AY <- bAD; coef_AH <- bAD - bHA
  coef_iderr <- coef_AY - (coef_HY + coef_AH)
  coef_rows[[ct]] <- data.frame(celltype=ct, gene=g, coef_HA_YA=coef_HY, coef_AD_HA=coef_AH, coef_AD_YA=coef_AY,
                                coefficient_identity_error=coef_iderr,
                                coef_finite=is.finite(bHA)&is.finite(bAD), betaConv=conv, full_rank=full_rank)
  cons <- list(HA_YA=c("grp","HA","YA"), AD_HA=c("grp","AD","HA"), AD_YA=c("grp","AD","YA"))
  raw <- list(); ashrL <- list()
  for (nm in names(cons)) {
    r <- results(dds, contrast=cons[[nm]]); raw[[nm]] <- r
    s <- tryCatch(lfcShrink(dds, res=r, type="ashr", quiet=TRUE),
                  error=function(e){ writeLines(sprintf("[%s %s] ashr ERROR: %s", ct, nm, conditionMessage(e)), warnlog); stop(sprintf("ashr failed %s %s", ct, nm)) })
    ashrL[[nm]] <- s
    shr_rows[[paste(ct,nm)]] <- data.frame(celltype=ct, contrast=nm, n_genes=nrow(s), n_ashr_finite=sum(is.finite(s$log2FoldChange)), shrinkage_ok=TRUE)
  }
  rHY<-raw$HA_YA$log2FoldChange; rAH<-raw$AD_HA$log2FoldChange; rAY<-raw$AD_YA$log2FoldChange
  res_iderr <- rAY - (rHY + rAH)
  res_rows[[ct]] <- data.frame(celltype=ct, gene=g, results_HA_YA=rHY, results_AD_HA=rAH, results_AD_YA=rAY,
                               results_identity_error=res_iderr, baseMean=bm, betaConv=conv,
                               results_finite=is.finite(rHY)&is.finite(rAH)&is.finite(rAY), is_frozen_hit=(paste(g,ct)%in%frozenkey))
  fz <- sig[sig$celltype==ct,]
  for (i in seq_len(nrow(fz))) {
    gn <- fz$gene[i]; present <- gn %in% g
    G <- function(r,col) if (present) r[gn,col] else NA_real_
    hit_rows[[paste(ct,gn)]] <- data.frame(gene=gn, celltype=ct, frozen_direction=fz$direction[i], frozen_log2FC=fz$log2FoldChange[i],
      present_in_de=present, baseMean=if(present) bm[gn] else NA, betaConv=if(present) conv[gn] else FALSE,
      HA_YA_LFC_coef=if(present) coef_HY[gn] else NA, AD_HA_LFC_coef=if(present) coef_AH[gn] else NA, AD_YA_LFC_coef=if(present) coef_AY[gn] else NA,
      coefficient_identity_error=if(present) coef_iderr[gn] else NA,
      HA_YA_LFC_raw=G(raw$HA_YA,"log2FoldChange"), HA_YA_SE=G(raw$HA_YA,"lfcSE"), HA_YA_LFC_ashr=if(present) ashrL$HA_YA[gn,"log2FoldChange"] else NA,
      AD_HA_LFC_raw=G(raw$AD_HA,"log2FoldChange"), AD_HA_SE=G(raw$AD_HA,"lfcSE"), AD_HA_LFC_ashr=if(present) ashrL$AD_HA[gn,"log2FoldChange"] else NA,
      AD_YA_LFC_raw=G(raw$AD_YA,"log2FoldChange"), AD_YA_SE=G(raw$AD_YA,"lfcSE"), AD_YA_LFC_ashr=if(present) ashrL$AD_YA[gn,"log2FoldChange"] else NA,
      results_identity_error=if(present) (G(raw$AD_YA,"log2FoldChange")-(G(raw$HA_YA,"log2FoldChange")+G(raw$AD_HA,"log2FoldChange"))) else NA,
      pvalue_HA_YA=G(raw$HA_YA,"pvalue"), padj_HA_YA=G(raw$HA_YA,"padj"),
      pvalue_AD_HA=G(raw$AD_HA,"pvalue"), padj_AD_HA=G(raw$AD_HA,"padj"),
      pvalue_AD_YA=G(raw$AD_YA,"pvalue"), padj_AD_YA=G(raw$AD_YA,"padj"))
  }
  cat(sprintf("  %s: donors %d, genes %d, full_rank %s, coef-identity max|err| %.2e, results-identity max|err| %.2e\n",
              ct, ncol(cts), nrow(dds), full_rank, max(abs(coef_iderr[is.finite(coef_iderr)])), max(abs(res_iderr[is.finite(res_iderr)]))))
}
hits<-do.call(rbind,hit_rows); coefA<-do.call(rbind,coef_rows); resA<-do.call(rbind,res_rows); shr<-do.call(rbind,shr_rows)
write.table(hits, file.path(OUT,"_frozen60_contrasts_R.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(coefA, file.path(OUT,"singlefit_coefficient_identity_audit.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(resA,  file.path(OUT,"results_lfc_additivity_audit.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(shr,   file.path(OUT,"shrinkage_audit.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# GATE 1: coefficient identity over finite+converged+full-rank
c1u <- coefA[coefA$coef_finite & coefA$betaConv & coefA$full_rank,]
g1max <- max(abs(c1u$coefficient_identity_error)); gate1 <- g1max <= 1e-10
# GATE 2: frozen 60 mandatory
need <- hits[hits$present_in_de,]
need$raw_finite  <- is.finite(need$HA_YA_LFC_raw)&is.finite(need$AD_HA_LFC_raw)&is.finite(need$AD_YA_LFC_raw)
need$se_finite   <- is.finite(need$HA_YA_SE)&is.finite(need$AD_HA_SE)&is.finite(need$AD_YA_SE)
need$ashr_finite <- is.finite(need$HA_YA_LFC_ashr)&is.finite(need$AD_HA_LFC_ashr)&is.finite(need$AD_YA_LFC_ashr)
need$res_id_ok   <- abs(need$results_identity_error) <= 1e-6
gate2 <- nrow(need)==60 && all(need$betaConv & is.finite(need$baseMean) & need$raw_finite & need$se_finite & need$ashr_finite & need$res_id_ok)

# coefficient identity summary
sc<-file(file.path(OUT,"singlefit_coefficient_identity_summary.txt"),open="wt"); wl<-function(s){cat(s,"\n");writeLines(s,sc)}
wl("=== Gate 1: single-fit coefficient additive identity (amendment v1.4) ===")
wl(sprintf("universe (finite coef + converged + full-rank): %d genes", nrow(c1u)))
wl(sprintf("max |coefficient_identity_error| = %.3e  (threshold 1e-10)", g1max))
wl(sprintf("GATE 1: %s", ifelse(gate1,"PASS","FAIL"))); close(sc)
# results additivity descriptive summary
sr<-file(file.path(OUT,"results_lfc_additivity_summary.txt"),open="wt"); wr<-function(s){cat(s,"\n");writeLines(s,sr)}
rf<-resA[resA$results_finite,]; v<-abs(rf$results_identity_error); viol<-rf[v>1e-6,]
wr("=== descriptive QC: results()-LFC additivity (NOT a gate; nothing excluded) ===")
wr(sprintf("total finite genes: %d", nrow(rf)))
wr(sprintf("identity <=1e-6: %d ; >1e-6: %d", sum(v<=1e-6), sum(v>1e-6)))
wr(sprintf("|error| median=%.2e 95th=%.2e 99th=%.2e max=%.2e", median(v), quantile(v,.95), quantile(v,.99), max(v)))
wr(sprintf("violator baseMean: min=%.3f median=%.3f max=%.3f", min(viol$baseMean), median(viol$baseMean), max(viol$baseMean)))
wr(sprintf("frozen-hit violators: %d ; coefficient-level identity max|err|: %.2e (exact)", sum(viol$is_frozen_hit), g1max)); close(sr)
# execution audit
ae<-file(file.path(OUT,"amendment_v1.4_execution_audit.txt"),open="wt"); wa<-function(s)writeLines(s,ae)
wa(sprintf("v1.4 tag: prereg/secretome-resilience-trajectory-v1.4"))
wa(sprintf("DESeq2 %s ; ashr %s ; R %s", packageVersion("DESeq2"), packageVersion("ashr"), getRversion()))
wa(sprintf("GATE 1 coefficient identity: max|err|=%.3e -> %s", g1max, ifelse(gate1,"PASS","FAIL")))
wa(sprintf("GATE 2 frozen-60: present %d/60, betaConv %d, results-id<=1e-6 %d, ashr-finite %d -> %s",
           nrow(need), sum(need$betaConv), sum(need$res_id_ok), sum(need$ashr_finite), ifelse(gate2,"PASS","FAIL")))
wa(sprintf("global results()-LFC non-additive genes (>1e-6): %d (descriptive; not excluded)", sum(v>1e-6)))
wa(sprintf("B6 decision: %s", ifelse(gate1&&gate2,"GO","STOP"))); close(ae); close(warnlog)
cat(sprintf("\nGATE 1 (coef <=1e-10): %s | GATE 2 (frozen-60): %s\n", ifelse(gate1,"PASS","FAIL"), ifelse(gate2,"PASS","FAIL")))
if(!gate1) stop(sprintf("GATE 1 FAIL: coef identity %.3e > 1e-10", g1max))
if(!gate2) stop("GATE 2 FAIL: a frozen hit failed convergence/results-identity/ashr/finite")
cat("[B2-B5 v1.4] BOTH GATES PASS -> proceed to B6.\n")
