#!/usr/bin/env python
"""Build GSE325391 resilience h5ad (granule-cell lineage) from R-exported RNA MTX + RDS meta.
donor_id = sample; group = RES/CTRL/MAD/SAD; cell_type = GC-lineage subtype.
All cells are ExN_GC (granule lineage) -> celltype_l1 = 'GC_lineage' (one bucket for DE).
Raw counts in layers['counts']. (Resilience DE = RES vs CTRL / RES vs SAD at DE time.)
"""
import os
import pandas as pd, anndata as ad, scipy.io

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
pre = f"{PROJ}/processed/per_dataset/GSE325391_rna"

X = scipy.io.mmread(f"{pre}_counts.mtx").tocsr()
genes = [l.strip() for l in open(f"{pre}_genes.txt")]
barcodes = [l.strip() for l in open(f"{pre}_barcodes.txt")]
if X.shape == (len(genes), len(barcodes)):
    X = X.T.tocsr()

A = ad.AnnData(X=X)
A.var_names = pd.Index(genes, dtype=str)
A.obs_names = pd.Index(barcodes, dtype=str)
A.var_names_make_unique()

m = pd.read_csv(f"{PROJ}/processed/per_dataset/GSE325391_adultgc_metadata.csv", index_col=0)
m.index = m.index.astype(str)
A = A[A.obs_names.isin(m.index)].copy()
m = m.loc[A.obs_names]
A.obs["donor_id"] = m["sample"].astype(str).values
A.obs["group"] = m["group"].astype(str).values
A.obs["cell_type"] = m["cell_type"].astype(str).values
A.obs["celltype_l1"] = "GC_lineage"
A.obs["sex"] = "U"
A.layers["counts"] = A.X.copy()

# remap to symbols if genes look like Ensembl (fallback via 268609 features map)
if A.var_names[0].startswith("ENSG"):
    import gzip
    e2s = {}
    with gzip.open(f"{PROJ}/raw/GSE268609_lazarov2024/GSE268609_features.tsv.gz", "rt") as f:
        for line in f:
            p = line.rstrip().split("\t")
            if len(p) >= 3 and p[2] == "Gene Expression":
                e2s[p[0]] = p[1]
    A.var["symbol"] = [e2s.get(g, g) for g in A.var_names]
    A.var_names = A.var["symbol"].values
    A.var_names_make_unique()
    print("remapped Ensembl -> symbol")

out = f"{PROJ}/processed/per_dataset/GSE325391_resilience.h5ad"
A.write(out)
print("wrote", out, A.shape)
print("group -> n_donors:\n", A.obs.groupby("group", observed=True)["donor_id"].nunique().to_string())
print("cell_type:\n", A.obs["cell_type"].value_counts().to_string())
