#!/usr/bin/env python
"""Build GSE268609 anchor h5ad from R-exported RNA MTX (counts) + RDS meta.data.
donor_id = SampleNumber; Group = YA/HA/MCI/AD/SA; celltype_l1 from Lazarov Cluster.
Raw counts in layers['counts']. (Normal-aging subset = YA+HA done at DE time.)
"""
import os
import numpy as np, pandas as pd, anndata as ad, scipy.io, scipy.sparse as sp

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
pre = f"{PROJ}/processed/per_dataset/GSE268609_rna"

X = scipy.io.mmread(f"{pre}_counts.mtx").tocsr()
genes = [l.strip() for l in open(f"{pre}_genes.txt")]
barcodes = [l.strip() for l in open(f"{pre}_barcodes.txt")]
if X.shape == (len(genes), len(barcodes)):
    X = X.T.tocsr()
print("counts cells x genes:", X.shape)

A = ad.AnnData(X=X)
A.var_names = pd.Index(genes, dtype=str)
A.obs_names = pd.Index(barcodes, dtype=str)
A.var_names_make_unique()

m = pd.read_csv(f"{PROJ}/processed/per_dataset/GSE268609_metadata.csv", index_col=0)
m.index = m.index.astype(str)
A = A[A.obs_names.isin(m.index)].copy()
m = m.loc[A.obs_names]
A.obs["Cluster"] = m["Cluster"].astype(str).values
A.obs["donor_id"] = m["SampleNumber"].astype(str).values
A.obs["Group"] = m["Group"].astype(str).values
A.obs["Neurogenic"] = m["Neurogenic"].astype(str).values
A.obs["sex"] = "U"

CL = {"Astrocytes": "Astro", "Microglia": "Micro", "mOli": "Oligo", "OPCs": "OPC",
      "Endothelial": "Endo", "Ependymal": "Ependymal", "mGC": "DG_GC", "NSC": "NSC",
      "Neuroblast": "Neuroblast", "Immature": "Immature", "CA_neurons": "CA_ExN",
      "CA2-4_neurons": "CA_ExN", "GABA_neurons": "InN"}
A.obs["celltype_l1"] = pd.Categorical([CL.get(x, "Other") for x in A.obs["Cluster"]])

# ENSG -> HGNC symbol (Gene Expression features), so DE/figures key on symbols. Symbol
# assignment belongs to the build step; it was previously done inside the Fig 1 script,
# which also OVERWROTE this file (reproducibility hazard — moved here 2026-06-03 audit fix).
# 10 symbols carry >1 Ensembl ID (read-through genes, e.g. GOLGA8M); collapse each to its
# highest-total-count copy rather than var_names_make_unique(), so one symbol == one row
# downstream and no duplicate (gene,celltype) DE rows are produced.
if str(A.var_names[0]).startswith("ENSG"):
    import gzip
    e2s = {}
    with gzip.open(f"{PROJ}/raw/GSE268609_lazarov2024/GSE268609_features.tsv.gz", "rt") as f:
        for line in f:
            p = line.rstrip().split("\t")
            if len(p) >= 3 and p[2] == "Gene Expression":
                e2s[p[0]] = p[1]
    A.var_names = pd.Index([e2s.get(g, g) for g in A.var_names], dtype=str)
    n_before = A.n_vars
    tot = np.asarray(A.X.sum(0)).ravel()                       # total counts per column
    order = np.argsort(-tot, kind="stable")                    # highest-count first
    keep = order[~pd.Index(A.var_names[order]).duplicated()]   # first (highest) per symbol
    A = A[:, np.sort(keep)].copy()                             # restore original column order
    print(f"remapped ENSG->symbol; collapsed {n_before} -> {A.n_vars} unique symbols "
          f"(dropped {n_before - A.n_vars} duplicate-symbol copies)")

A.layers["counts"] = A.X.copy()

out = f"{PROJ}/processed/per_dataset/GSE268609_anchor.h5ad"
A.write(out)
print("wrote", out, A.shape)
print("celltype_l1:\n", A.obs["celltype_l1"].value_counts().to_string())
print("\nGroup -> n_donors:\n",
      A.obs.groupby("Group", observed=True)["donor_id"].nunique().to_string())
print("\nniche senders x Group (n donors with >=10 cells):")
for ct in ["Astro", "Micro", "Oligo", "OPC", "Endo"]:
    s = A.obs[A.obs.celltype_l1 == ct]
    nd = s.groupby("Group", observed=True)["donor_id"].apply(
        lambda d: (d.value_counts() >= 10).sum())
    print(f"  {ct}: {nd.to_dict()}")
