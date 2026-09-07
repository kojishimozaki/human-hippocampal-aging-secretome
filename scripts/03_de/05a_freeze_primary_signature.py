#!/usr/bin/env python
"""Freeze the PRIMARY niche-secretome aging signature as an immutable artifact, so the
A1 external-replication test (scripts/03_de/05_external_replication.py) *tests* a fixed
hypothesis rather than re-discovering one. Phase-1.5 pre-specification (locked before any
external-cohort DE is run).

Inputs (committed):
  results/de/de_GSE268609_per_celltype.csv  -- primary donor-pseudobulk DE (Lazarov 2024,
      GSE268609, PMID 41741649; YA n=8 vs HA n=9, design ~grp). cols incl. gene, celltype,
      log2FoldChange, lfcSE, padj, baseMean.
  refs/secretome_union.csv                  -- analysis universe (HPA predicted-secreted ∪
      UniProt reviewed-secreted; 2,224 genes), cols gene,n_sources,sources,is_core,category.

Frozen definition (matches the manuscript's Fig-2 "60 hits"):
  niche celltypes = {Astro, Micro, Endo, OPC, Oligo}
  frozen hits      = secretome_union ∩ {padj < 0.1} within each niche celltype
                   = Astro 26 / Micro 12 / Endo 11 / OPC 6 / Oligo 5 = 60 gene×celltype pairs
                     (57 unique genes; 36 UP / 24 DOWN).

Outputs (immutable; do NOT regenerate after the external test is run):
  results/validation/frozen_primary_signature.csv      -- the 60 pre-specified hits (candidate set)
  results/validation/frozen_primary_secretome_logfc.csv-- ALL secretome gene×niche-celltype primary
        log2FC/padj (the fixed vector used by metric (ii) Spearman correlation and metric (iii)
        background-concordance; frozen so the external test cannot drift the primary side).

Deterministic (no RNG). seed not required.
"""
import os
import pandas as pd

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
NICHE = ["Astro", "Micro", "Endo", "OPC", "Oligo"]
PADJ = 0.1

de = pd.read_csv(f"{PROJ}/results/de/de_GSE268609_per_celltype.csv")
sec = pd.read_csv(f"{PROJ}/refs/secretome_union.csv")[["gene", "is_core", "category"]]

# restrict to niche celltypes and the secretome universe
de_niche = de[de.celltype.isin(NICHE)].merge(sec, on="gene", how="inner")

cols = ["gene", "celltype", "log2FoldChange", "lfcSE", "pvalue", "padj", "baseMean",
        "is_core", "category"]

# (1) full per-celltype secretome log2FC vector (fixed primary side for metrics ii/iii)
full = de_niche[cols].copy()
full["direction"] = full.log2FoldChange.apply(lambda x: "UP" if x > 0 else "DOWN")
full = full.sort_values(["celltype", "padj"]).reset_index(drop=True)
out_full = f"{PROJ}/results/validation/frozen_primary_secretome_logfc.csv"
os.makedirs(os.path.dirname(out_full), exist_ok=True)
full.to_csv(out_full, index=False)

# (2) the 60 pre-specified hits (candidate set)
hits = full[full.padj < PADJ].copy().reset_index(drop=True)
out_hits = f"{PROJ}/results/validation/frozen_primary_signature.csv"
hits.to_csv(out_hits, index=False)

# report (for the lock document)
print(f"secretome universe genes: {sec.gene.nunique()}")
print(f"frozen full secretome×niche rows: {len(full)} -> {out_full}")
print(f"\nFROZEN SIGNATURE (padj<{PADJ}) per niche celltype:")
print(hits.groupby("celltype").size().reindex(NICHE).to_string())
print(f"TOTAL frozen hits: {len(hits)} | unique genes: {hits.gene.nunique()} | "
      f"UP {(hits.direction=='UP').sum()} / DOWN {(hits.direction=='DOWN').sum()}")
print(f"-> {out_hits}")
