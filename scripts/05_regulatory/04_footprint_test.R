#!/usr/bin/env Rscript
# FP step 3: donor-level differential TF footprinting (YA vs HA), Astro/Micro/Endo.
# For each motif, aggregate the Tn5-insertion meta-profile around its sites PER DONOR, derive a
# footprint-depth score (flank/center protection) + an accessibility score, then test YA vs HA at
# the DONOR level (Wilcoxon, n=8/9). The Tn5 sequence bias is identical at the same sites in both
# groups, so it CANCELS in the YA-vs-HA differential (no bias model needed). Avoids pooling
# pseudoreplication (donor is the unit). Short list locked in PHASE_REPORT.md.
suppressMessages({ library(data.table); library(GenomicRanges) })
P <- Sys.getenv("PROJ", unset = "/home/neurofuture/Bioanalysis/01.SGZ_aging_project"); TMP <- file.path(P, "processed/footprint_tmp")
W <- 200L; CENTER <- 15L; FL1 <- 50L; FL2 <- 150L      # window / center / flank (bp)
sites <- readRDS(file.path(TMP, "motif_sites.rds"))
TFs <- names(sites)

# windows DT (all TFs): chr, start, end, center, tf  (+/- W around motif centers)
win <- rbindlist(lapply(TFs, function(tf) {
  s <- sites[[tf]]
  data.table(chr = as.character(seqnames(s)), center = start(s),
             start = start(s) - W, end = start(s) + W, tf = tf)
}))
win[, start := pmax(1L, start)]
nsites <- win[, .(n = uniqueN(center)), by = tf]
setkey(win, chr, start, end)

prof_all <- list(); score_all <- list(); test_all <- list()
for (ct in c("Astro", "Micro", "Endo")) {
  f <- file.path(TMP, sprintf("ins_%s.tsv.gz", ct))
  if (!file.exists(f)) { cat("missing", f, "\n"); next }
  ins <- fread(cmd = sprintf("zcat %s", f), header = FALSE, col.names = c("chr", "pos", "group", "donor"),
               colClasses = list(character = 1, integer = 2, character = 3, integer = 4))
  dtot <- ins[, .(dtot = .N), by = .(donor, group)]            # per-donor total insertions (depth)
  ins[, `:=`(start = pos, end = pos)]; setkey(ins, chr, start, end)
  hits <- foverlaps(ins, win, type = "within", nomatch = 0L)   # insertions falling in any motif window
  hits[, rel := pos - center]
  cat(sprintf("%s: %d insertions, %d in motif windows\n", ct, nrow(ins), nrow(hits)))
  # meta-profile per (tf, group, rel) averaged across donors -> for figure
  prof <- hits[, .N, by = .(tf, group, rel, donor)]
  prof <- merge(prof, dtot, by = c("donor", "group"))
  prof[, cpm := N / dtot * 1e6]
  prof_all[[ct]] <- prof[, .(cpm = mean(cpm)), by = .(tf, group, rel)][, celltype := ct]
  # per (tf, donor) footprint-depth + accessibility
  sc <- hits[, .(center_n = sum(abs(rel) <= CENTER),
                 flank_n = sum(abs(rel) >= FL1 & abs(rel) <= FL2),
                 win_n = .N), by = .(tf, donor, group)]
  sc <- merge(sc, dtot, by = c("donor", "group"))
  sc <- merge(sc, nsites, by = "tf")
  cw <- 2 * CENTER + 1; fw <- 2 * (FL2 - FL1 + 1)
  sc[, fp_depth := log2(((flank_n / fw) + 1) / ((center_n / cw) + 1))]      # higher = deeper footprint
  sc[, access := win_n / n / dtot * 1e6]                                    # CPM-like accessibility/site
  sc[, celltype := ct]
  score_all[[ct]] <- sc
  # donor-level YA vs HA test per TF
  tt <- sc[, {
    ya <- fp_depth[group == "YA"]; ha <- fp_depth[group == "HA"]
    aya <- access[group == "YA"]; aha <- access[group == "HA"]
    .(celltype = ct, n_YA = length(ya), n_HA = length(ha),
      fp_YA = mean(ya), fp_HA = mean(ha), fp_dz = mean(ha) - mean(ya),
      fp_p = tryCatch(wilcox.test(ha, ya)$p.value, error = function(e) NA_real_),
      acc_dz = mean(aha) - mean(aya),
      acc_p = tryCatch(wilcox.test(aha, aya)$p.value, error = function(e) NA_real_))
  }, by = tf]
  test_all[[ct]] <- tt
  rm(ins, hits); gc()
}
test <- rbindlist(test_all)
test[, fp_padj := p.adjust(fp_p, "BH")][, acc_padj := p.adjust(acc_p, "BH")]
fwrite(test, file.path(P, "results/regulatory/footprint_test.csv"))
fwrite(rbindlist(score_all), file.path(P, "results/regulatory/footprint_scores_per_donor.csv"))
fwrite(rbindlist(prof_all), file.path(P, "results/regulatory/footprint_metaprofiles.csv"))

cat("\n=== donor-level differential footprint (YA vs HA), short-list TFs ===\n")
cat("    fp_dz>0 = deeper footprint (more TF binding) with age; fp_padj = BH over TFxcelltype\n")
print(test[order(fp_p), .(celltype, tf, fp_YA = round(fp_YA,3), fp_HA = round(fp_HA,3),
                          fp_dz = round(fp_dz,3), fp_p = round(fp_p,4), fp_padj = round(fp_padj,3),
                          acc_dz = round(acc_dz,3), acc_p = round(acc_p,4))], row.names = FALSE)
cat("\nFDR<0.1 footprint changes:", test[fp_padj < 0.1, .N], "| nominal p<0.05:", test[fp_p < 0.05, .N], "\n")
cat("wrote results/regulatory/footprint_test.csv (+ per-donor scores, meta-profiles)\n")
