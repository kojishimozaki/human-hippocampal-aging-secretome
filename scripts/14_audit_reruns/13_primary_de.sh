#!/usr/bin/env bash
# [B-12 / F-0-006] The launch command for the PRIMARY differential expression table.
#
# results/de/de_GSE268609_per_celltype.csv is the origin of all 60 frozen hits and therefore of
# every downstream leg. The 2026-08-26 audit found no Makefile target invoked the engine that
# produces it, and its git history is a single initial commit -- the same gap as the
# de_GSE325391 tables (script 12).
#
# VERIFIED 2026-08-27, and the verification is not a clean pass -- read this before using the
# output. Running this command reproduces 99,543 of the committed table's 99,544 rows. The one
# missing row is LINC01238 in oligodendrocytes, which has a total pseudobulk count of 2 across
# 2 of 17 donors and so fails the current gene filter (sum >= 10 AND non-zero in >= max(3, 0.5n)
# = 8). The committed table was made with a looser filter than the code now applies. Dropping
# that one gene also moves GOLGA8M (oligodendrocytes) from -0.305 to +0.077 through the
# recomputed size factors and independent filtering; both genes have padj > 0.9.
#
# What is unaffected, and it is what the paper rests on: all 60 frozen hits are present, all 60
# keep padj < 0.1, and their largest log2FC difference is 4.3e-05.
#
# The committed table is deliberately NOT overwritten by this script -- it writes beside it. See
# audit_log/2026-08-27_provenance_b12_b13_b14/.
set -euo pipefail
: "${PROJ:=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
export PROJ
PY="${PY:-$HOME/miniforge3/envs/bio/bin/python}"
OUT="${OUT:-$PROJ/results/de/de_GSE268609_per_celltype.regenerated.csv}"

"$PY" "$PROJ/scripts/03_de/01_pseudobulk_de.py" \
  --h5ad "$PROJ/processed/per_dataset/GSE268609_anchor.h5ad" \
  --cohort GSE268609 --mode group --group-col Group \
  --test-level HA --ref-level YA \
  --celltypes Astro Micro Oligo OPC Endo \
  --out "$OUT"
echo "wrote $OUT (committed table left untouched at results/de/de_GSE268609_per_celltype.csv)"
