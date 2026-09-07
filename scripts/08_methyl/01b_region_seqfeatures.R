#!/usr/bin/env Rscript
# M1 methylome — per-region sequence features for the matched-null (GC, CpG density, CpG O/E).
# CpG density is the dominant methylation-specific confounder (added vs A2's GC+width+acc+TSS match).
# hg38 via BSgenome.Hsapiens.UCSC.hg38. Vectorised over all regions. seed-free.
# In: results/methyl/regions_meta.csv   Out: results/methyl/region_seqfeatures.csv
suppressMessages({library(BSgenome.Hsapiens.UCSC.hg38); library(GenomicRanges); library(Biostrings)})
PROJ <- Sys.getenv("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
g <- BSgenome.Hsapiens.UCSC.hg38

reg <- read.csv(file.path(PROJ, "results/methyl/regions_meta.csv"), stringsAsFactors = FALSE)
gr <- GRanges(reg$chrom, IRanges(pmax(1L, reg$start), reg$end))
seqlevels(gr) <- seqlevels(g)[match(seqlevels(gr), seqlevels(g))]
seqlengths(gr) <- seqlengths(g)[seqlevels(gr)]
gr <- trim(gr)                                  # clamp ends to chrom length

seqs <- getSeq(g, gr)
lf <- letterFrequency(seqs, c("A", "C", "G", "T", "N"))
di <- dinucleotideFrequency(seqs)
w  <- width(gr)
len_eff <- pmax(1, w - lf[, "N"])               # effective (non-N) length
gc   <- (lf[, "C"] + lf[, "G"]) / pmax(1, rowSums(lf))
ncpg <- di[, "CG"]
cpg_density <- ncpg / pmax(1, len_eff)
# CpG observed/expected = (CpG * L) / (C * G)   (Gardiner-Garden)
cpg_oe <- (as.numeric(ncpg) * len_eff) / pmax(1, as.numeric(lf[, "C"]) * as.numeric(lf[, "G"]))

out <- data.frame(region_id = reg$region_id, gc = round(gc, 5), n_cpg = ncpg,
                  cpg_density = round(cpg_density, 6), cpg_oe = round(cpg_oe, 5),
                  n_frac = round(lf[, "N"] / pmax(1, w), 4))
write.csv(out, file.path(PROJ, "results/methyl/region_seqfeatures.csv"), row.names = FALSE)
cat(sprintf("wrote region_seqfeatures.csv for %d regions | GC %.3f | CpG/bp %.4f | CpGo/e %.3f (medians)\n",
            nrow(out), median(gc), median(cpg_density), median(cpg_oe, na.rm = TRUE)))
