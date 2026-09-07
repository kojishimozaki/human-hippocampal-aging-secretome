#!/usr/bin/env bash
# [B-8 / F-0-001] Recovered launch commands for results/de/de_GSE325391_*.csv
#
# The audit found these three files have ZERO write sites anywhere in the repository and only
# an initial-commit git history, so "RES vs CTRL is 2 genes" (Fig 6) had no reproducible origin.
# The commands below were reconstructed from the committed CSVs' own `cohort` column
# (GSE325391_RESvCTRL etc.), the h5ad's `group` levels (RES / SAD / CTRL / MAD) and
# `celltype_l1` (GC_lineage), then VERIFIED: all three outputs reproduce the committed files
# with a maximum absolute difference of 0 across baseMean, log2FoldChange, lfcSE, stat,
# pvalue and padj (29,688 / 29,521 / 29,563 rows respectively, identical gene sets).
#
# ROW ORDER DIFFERS: this emits baseMean-descending order (01_pseudobulk_de.py:L139-L144, the
# 2026-06-03 duplicate-collapse fix) while the committed CSVs predate that block and keep the
# h5ad var_names order. No duplicate (gene,celltype) rows exist in this cohort, so the block
# only reorders -- content is identical. Running this REWRITES the three committed files.
#
# Verified 2026-08-26. See audit_log/2026-08-26_remaining_still_open/RESOLUTION.md.
#
# NOTE this closes only the FIRST half of B-8. The Seurat RDS -> MTX export that produced
# GSE268609_{rna,peaks}_counts.mtx is still absent (`writeMM` appears 0 times in scripts/)
# and is NOT recovered by this script.
set -euo pipefail
: "${PROJ:=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
export PROJ
PY="${PY:-$HOME/miniforge3/envs/bio/bin/python}"
OUTDIR="${OUTDIR:-$PROJ/results/de}"
H5="$PROJ/processed/per_dataset/GSE325391_resilience.h5ad"

for pair in "RES CTRL" "RES SAD" "SAD CTRL"; do
  set -- $pair; T="$1"; R="$2"
  echo "==> GSE325391 ${T} vs ${R}"
  "$PY" "$PROJ/scripts/03_de/01_pseudobulk_de.py" \
    --h5ad "$H5" \
    --cohort "GSE325391_${T}v${R}" \
    --mode group --group-col group \
    --test-level "$T" --ref-level "$R" \
    --out "$OUTDIR/de_GSE325391_${T}v${R}.csv"
done
echo "done: $OUTDIR/de_GSE325391_{RESvCTRL,RESvSAD,SADvCTRL}.csv"
