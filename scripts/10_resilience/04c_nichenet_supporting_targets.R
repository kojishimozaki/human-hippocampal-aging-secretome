#!/usr/bin/env Rscript
# CONSEQUENCE (Phase A) — STEP 04c (part 1/2): NicheNet supporting top-50 target lists.
# Provenance: VERBATIM from origin session 027c5fdb (transcript line 853, R part); feeds 04c_nichenet_supporting_wilcox.py.
# CAVEAT (agent-flagged, honest): the committed nichenet_supporting_results.tsv UNION row is numerically identical to
#   STEP 04b's nichenet_top50_results.tsv UNION (n=228, median_diff -0.0207, p 0.4378); line 853 is the only writer of
#   that exact filename in the transcript and uses the same competitive-Wilcoxon method, but whether the committed file
#   was written by 853 or later overwritten by a top-50 copy could not be 100% confirmed. Verify per-ligand rows on re-run.
# env: sgz_r.  DATA-DEPENDENT (refs/nichenet + results/de/de_GSE325391_RESvSAD.csv) -> run from $PROJ.
# Out: results/resilience/{_nichenet_targets_top50.tsv, _nichenet_targets_union.txt}.  (stderr -> model_warnings.log.)
suppressMessages({})
cat("fgsea present:", requireNamespace("fgsea",quietly=TRUE), "\n")
de <- read.csv("results/de/de_GSE325391_RESvSAD.csv", stringsAsFactors=FALSE)
cat("RESvSAD cols:", paste(colnames(de),collapse=","), "\n")
# orient: contrast 'RESvSAD' = test RES vs ref SAD -> log2FC>0 = up in RES. Check a row.
cat("example rows (gene,log2FC,stat):\n"); print(head(de[order(-abs(de$stat)),c("gene","log2FoldChange","stat")],3))
lig <- readLines("results/step0/_nichenet_ligands.txt")
M <- readRDS("refs/nichenet/ligand_target_matrix_nsga2r_final.rds")   # rows=targets, cols=ligands
cat("ligand_target_matrix dim:", paste(dim(M),collapse=" x "), "| ligands present:", length(intersect(lig,colnames(M))),"/",length(lig),"\n")
ligp <- intersect(lig, colnames(M))
top <- lapply(ligp, function(l){ v<-sort(M[,l],decreasing=TRUE); names(v)[1:50] })
names(top) <- ligp
long <- do.call(rbind, lapply(ligp, function(l) data.frame(ligand=l, target=top[[l]])))
write.table(long, "results/resilience/_nichenet_targets_top50.tsv", sep="\t", quote=FALSE, row.names=FALSE)
writeLines(sort(unique(long$target)), "results/resilience/_nichenet_targets_union.txt")
cat("union top-50 targets:", length(unique(long$target)), "across", length(ligp), "ligands\n")
