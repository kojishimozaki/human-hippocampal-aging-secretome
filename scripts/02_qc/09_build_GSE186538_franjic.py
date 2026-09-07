#!/usr/bin/env python
"""Build Franjic GSE186538 HUMAN hippocampal-entorhinal snRNA AnnData for A1-secondary
external replication (continuous mid-old age slope) + B2 subregion substrate.

Dataset: GSE186538 (PMID 34798047, Franjic et al. 2022, Neuron; verified via !Series_pubmed_id).
Human only (pig/macaque matrices ignored). 6 donors; ages from per-GSM titles (GEO SOFT provenance):
  HSB179=48, HSB181=44, HSB231=79, HSB237=51, HSB282=48, HSB628=50  -> mid-to-old, NO young adult.
  CAVEAT (locked, RESOLUTION.md §A1): the age slope is leverage-dominated by the single 79 y donor
  (five cluster at 44-51); this cohort is SUPPORTING, not a YA-vs-HA replication.

Author cell labels live in meta column 'cluster' (e.g. "Astro AQP4 GFAP", "Oligo OPALIN ...").
Mapped to the 5 niche celltypes by label prefix; everything else -> 'other' (excluded from DE).
'region' (DG/CA1/CA24/EC/SUB) retained for B2.

Out: processed/per_dataset/GSE186538_franjic_human.h5ad  (X=counts, layers['counts']=counts,
     obs: donor_id, age_years, region, cluster, celltype_l1, sex='U')
"""
import os
import numpy as np, pandas as pd, scanpy as sc, scipy.io, scipy.sparse as sp

P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
D = f"{P}/raw/GSE186538_franjic2022"
AGE = {"HSB179": 48, "HSB181": 44, "HSB231": 79, "HSB237": 51, "HSB282": 48, "HSB628": 50}

genes = pd.read_csv(f"{D}/GSE186538_Human_genes.txt.gz", header=None)[0].astype(str).values
meta = pd.read_csv(f"{D}/GSE186538_Human_cell_meta.txt.gz", sep="\t")
print(f"genes={len(genes)}  cells={len(meta)}  cols={list(meta.columns)}")

print("reading MTX (genes x cells, 482M nnz) ...", flush=True)
X = scipy.io.mmread(f"{D}/GSE186538_Human_counts.mtx.gz").tocsr()   # 33939 genes x 219058 cells
assert X.shape == (len(genes), len(meta)), f"shape {X.shape}"
X = X.T.tocsr()                                                     # -> cells x genes

A = sc.AnnData(X=X)
A.var_names = genes
A.obs_names = meta["cell_name"].astype(str).values
A.var_names_make_unique()
A.layers["counts"] = A.X.copy()
A.obs["donor_id"] = meta["samplename"].astype(str).values
A.obs["region"] = meta["region"].astype(str).values
A.obs["cluster"] = meta["cluster"].astype(str).values
A.obs["age_years"] = A.obs["donor_id"].map(AGE).astype(float)
A.obs["sex"] = "U"


def mapct(c):
    c = str(c)
    if c.startswith("Astro"): return "Astro"
    if c.startswith("Micro"): return "Micro"
    if c.startswith("Oligo"): return "Oligo"
    if c.startswith("OPC"):  return "OPC"
    if c.startswith("Endo") or c.startswith("aEndo"): return "Endo"
    return "other"


A.obs["celltype_l1"] = A.obs["cluster"].map(mapct).astype("category")

os.makedirs(f"{P}/processed/per_dataset", exist_ok=True)
A.write(f"{P}/processed/per_dataset/GSE186538_franjic_human.h5ad")

print("\ndonor ages:", {d: AGE[d] for d in sorted(AGE)})
print("\nper-donor x niche-celltype cell counts:")
niche = A.obs.celltype_l1 != "other"
ct = pd.crosstab(A.obs.loc[niche, "donor_id"], A.obs.loc[niche, "celltype_l1"])
print(ct.to_string())
print("\nfrozen-signature gene presence:")
fz = pd.read_csv(f"{P}/results/validation/frozen_primary_signature.csv")
present = fz.gene.isin(set(A.var_names))
print(f"  {present.sum()}/{len(fz)} frozen hit rows have their gene in Franjic var_names "
      f"({fz.loc[present,'gene'].nunique()}/{fz.gene.nunique()} unique genes)")
print(f"wrote processed/per_dataset/GSE186538_franjic_human.h5ad  ({A.n_obs} cells x {A.n_vars} genes)")
