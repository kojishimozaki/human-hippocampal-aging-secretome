#!/usr/bin/env python
"""GSE199243 (glia lifespan): per-sample DENSE genes x cells txt.gz -> concatenated h5ad.

- maps age/sex via refs/metadata_master.csv (by GSM parsed from filename)
- merges sequencing-depth variants (Sample53_1K/200/500) into one donor
- marker-score argmax -> celltype_l1 (niche senders + neurons; no clustering needed)
- raw counts in layers['counts']; restricts to adults at DE time
Output: processed/per_dataset/GSE199243_glia.h5ad
"""
import os, glob, re
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import scipy.sparse as sp

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
EXT = f"{PROJ}/processed/per_dataset/GSE199243_extracted"
meta = pd.read_csv(f"{PROJ}/refs/metadata_master.csv")
meta = meta[meta.dataset == "GSE199243"].set_index("gsm")

files = sorted(f for f in glob.glob(f"{EXT}/*.txt.gz") if "barcode" not in os.path.basename(f))
ads = []
for f in files:
    base = os.path.basename(f)                       # GSM5967896_Sample43.txt.gz
    gsm = base.split("_")[0]
    sample = re.sub(r"\.txt\.gz$", "", base.split("_", 1)[1])   # Sample43
    donor = re.sub(r"_(1K|1k|200|500)$", "", sample)            # merge depth variants
    df = pd.read_csv(f, delim_whitespace=True, index_col=0)     # genes x cells (cols=barcodes)
    a = ad.AnnData(X=sp.csr_matrix(df.values.T.astype(np.float32)))   # cells x genes
    a.var_names = df.index.astype(str)
    a.obs_names = [f"{gsm}_{bc}" for bc in df.columns.astype(str)]
    a.obs["gsm"] = gsm; a.obs["sample"] = sample; a.obs["donor_id"] = donor
    age = meta.loc[gsm, "age_years"] if gsm in meta.index else np.nan
    a.obs["age_years"] = float(age) if pd.notna(age) else np.nan
    a.obs["sex"] = meta.loc[gsm, "sex"] if gsm in meta.index else "U"
    ads.append(a)
    print(f"{base}: {a.shape} donor={donor} age={a.obs['age_years'][0]}", flush=True)

A = ad.concat(ads, join="outer")
A.var_names_make_unique()
A.layers["counts"] = A.X.copy()
A.var["mt"] = A.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True, percent_top=None)
A = A[(A.obs.n_genes_by_counts >= 300) & (A.obs.pct_counts_mt < 10)].copy()

A.X = A.layers["counts"].copy()
sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
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
    sc.tl.score_genes(A, gs, score_name=f"sc_{ct}")
    score_cols.append(f"sc_{ct}")
S = A.obs[score_cols].values
A.obs["celltype_l1"] = pd.Categorical(
    [c.replace("sc_", "") for c in np.array(score_cols)[S.argmax(1)]])
A.X = A.layers["counts"].copy()

out = f"{PROJ}/processed/per_dataset/GSE199243_glia.h5ad"
A.write(out)
print("wrote", out, A.shape, flush=True)
print("celltype_l1:\n", A.obs["celltype_l1"].value_counts().to_string())
print("donor x age:\n", A.obs.groupby("donor_id")["age_years"].first().to_string())
