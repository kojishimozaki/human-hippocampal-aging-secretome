#!/usr/bin/env python
"""Re-annotate GSE199243 glia with Leiden clustering + CLUSTER-LEVEL marker scoring
(replaces the crude per-cell argmax) for a fair cross-cohort re-test vs 268609.
Adds obs['celltype_l1_v2']; keeps raw counts in layers['counts'].
"""
import os
import scanpy as sc, pandas as pd
P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
A = sc.read_h5ad(f"{P}/processed/per_dataset/GSE199243_glia.h5ad")
A.X = A.layers["counts"].copy()
sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
A.raw = A
sc.pp.highly_variable_genes(A, n_top_genes=2000, batch_key="donor_id")
A2 = A[:, A.var.highly_variable].copy()
sc.pp.scale(A2, max_value=10)
sc.tl.pca(A2, n_comps=30)
sc.pp.neighbors(A2, n_neighbors=15, n_pcs=30)
sc.tl.leiden(A2, resolution=1.0, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
A.obs["leiden"] = A2.obs["leiden"].values

markers = {
    "Astro": ["AQP4", "GFAP", "SLC1A2", "SLC1A3", "GJA1", "ALDH1L1"],
    "Micro": ["CSF1R", "C1QB", "C1QA", "P2RY12", "CX3CR1", "AIF1"],
    "Oligo": ["PLP1", "MOG", "MOBP", "MBP"],
    "OPC": ["PDGFRA", "CSPG4", "OLIG1", "OLIG2"],
    "Endo": ["CLDN5", "FLT1", "PECAM1", "VWF"],
    "Vascular": ["PDGFRB", "RGS5", "DCN", "COL1A2"],
    "ExN": ["SLC17A7", "SATB2", "RBFOX3", "PROX1"],
    "InN": ["GAD1", "GAD2"],
}
score_cols = []
for ct, gs in markers.items():
    gs = [g for g in gs if g in A.var_names]
    sc.tl.score_genes(A, gs, score_name=f"sc_{ct}", use_raw=True)
    score_cols.append(f"sc_{ct}")
cl_scores = A.obs.groupby("leiden", observed=True)[score_cols].mean()
cl_assign = cl_scores.idxmax(axis=1).str.replace("sc_", "")
A.obs["celltype_l1_v2"] = A.obs["leiden"].map(cl_assign).astype("category")
A.X = A.layers["counts"].copy()
A.write(f"{P}/processed/per_dataset/GSE199243_glia.h5ad")

print("leiden clusters:", A.obs.leiden.nunique())
print("celltype_l1_v2:\n", A.obs.celltype_l1_v2.value_counts().to_string())
print("\ncluster -> assignment:\n", cl_assign.to_string())
print("\nv1(per-cell argmax) vs v2(cluster) crosstab:\n",
      pd.crosstab(A.obs.celltype_l1, A.obs.celltype_l1_v2).to_string())
