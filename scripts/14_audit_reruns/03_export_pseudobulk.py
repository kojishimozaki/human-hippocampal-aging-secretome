#!/usr/bin/env python
"""[F-5-006 cross-check, step 1/2] Export the YA-vs-HA donor pseudobulk for the R side.

`scripts/14_audit_reruns/04_ashr_crosscheck.R` needs the SAME matrices that
`01_covariates_and_shrinkage.py` fitted, otherwise the apeGLM-vs-ashr comparison is not
apples-to-apples. Rather than have the R script reach into an 8.8 GB h5ad, this writes the
already-filtered pseudobulk out once.

Everything here is copied from 01_covariates_and_shrinkage.py:
  * YA/HA nuclei only, donors with >= 10 nuclei of that cell type (MIN_CELLS)
  * counts layer summed per donor
  * gene filter  (total >= 10) AND (non-zero in >= max(3, n_donors/2) donors)
    -- applied here so the R side fits exactly the gene set pyDESeq2 fitted
  * group levels written as plain YA/HA (the 1_YA / 2_HA re-labelling in script 01 was a
    pyDESeq2-specific work-around; R's DESeq2 takes an explicit factor level order instead)

RUNTIME ~10-15 min, dominated by loading the anchor. seed 42 (nothing here is random).

Out: results/audit_reruns/pseudobulk/counts_<celltype>.tsv.gz    genes x donors
     results/audit_reruns/pseudobulk/metadata_<celltype>.tsv     donors x covariates
     results/audit_reruns/pseudobulk/_manifest.tsv               shapes + checksums
"""
import os, hashlib
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_cache")
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
import warnings; warnings.filterwarnings("ignore")

P = os.environ.get("PROJ", os.getcwd())
OUT = f"{P}/results/audit_reruns/pseudobulk"
os.makedirs(OUT, exist_ok=True)
NICHE = ["Astro", "Micro", "Oligo", "OPC", "Endo"]
MIN_CELLS = 10

cov = pd.read_csv(f"{P}/results/de/gse268609_donor_covariates.csv")
cov["SampleNumber"] = cov.SampleNumber.astype(str)
pmi_map = dict(zip(cov.SampleNumber, cov.pmi))
arm_map = dict(zip(cov.SampleNumber, cov.arm))
sex_map = dict(zip(cov.SampleNumber, cov.sex_inferred))

A = sc.read_h5ad(f"{P}/processed/per_dataset/GSE268609_anchor.h5ad")
A = A[A.obs.Group.isin(["YA", "HA"])].copy()
print(f"YA+HA nuclei {A.n_obs}, donors {A.obs.donor_id.nunique()}", flush=True)


def pseudobulk(sub):
    Xs = sub.layers["counts"] if "counts" in sub.layers else sub.X
    Xs = Xs.tocsr() if sp.issparse(Xs) else sp.csr_matrix(Xs)
    d = sub.obs["donor_id"].astype(str).values
    u = pd.unique(d)
    M = np.vstack([np.asarray(Xs[d == x].sum(0)).ravel() for x in u])
    return pd.DataFrame(M, index=u, columns=sub.var_names.astype(str))


rows = []
for ct in NICHE:
    sub = A[A.obs.celltype_l1.astype(str) == ct].copy()
    nc = sub.obs.groupby("donor_id", observed=True).size()
    sub = sub[sub.obs.donor_id.isin(nc[nc >= MIN_CELLS].index)].copy()
    pb = pseudobulk(sub)
    ncell = sub.obs.groupby("donor_id", observed=True).size().reindex(pb.index)
    # identical filter to 01_covariates_and_shrinkage.py::fit
    keep = (pb.sum(0) >= 10) & (pb.astype(bool).sum(0) >= max(3, int(0.5 * len(pb))))
    pbf = pb.loc[:, keep].astype(int)
    md = pd.DataFrame({
        "grp": sub.obs.groupby("donor_id", observed=True).Group.first().astype(str).reindex(pb.index).values,
        "pmi": [pmi_map.get(d, np.nan) for d in pb.index],
        "arm": [arm_map.get(d, None) for d in pb.index],
        "sex": [sex_map.get(d, None) for d in pb.index],
        "n_cells": ncell.values.astype(int),
        "logncell": np.log10(ncell.values.astype(float))}, index=pb.index)
    md.index.name = "donor"

    cpath = f"{OUT}/counts_{ct}.tsv.gz"
    mpath = f"{OUT}/metadata_{ct}.tsv"
    pbf.T.to_csv(cpath, sep="\t")          # genes x donors, as 11_trajectory/01 writes them
    md.to_csv(mpath, sep="\t")
    h = hashlib.sha256(pd.util.hash_pandas_object(pbf, index=True).values.tobytes()).hexdigest()[:16]
    rows.append(dict(celltype=ct, n_donors=len(md), n_YA=int((md.grp == "YA").sum()),
                     n_HA=int((md.grp == "HA").sum()), n_nuclei=int(sub.n_obs),
                     n_genes_before_filter=pb.shape[1], n_genes=pbf.shape[1], sha256_16=h))
    print(f"  {ct}: {len(md)} donors (YA {rows[-1]['n_YA']} / HA {rows[-1]['n_HA']}), "
          f"{sub.n_obs} nuclei, {pbf.shape[1]} genes kept of {pb.shape[1]}", flush=True)

M = pd.DataFrame(rows)
M.to_csv(f"{OUT}/_manifest.tsv", sep="\t", index=False)
print("\n" + M.to_string(index=False))
print(f"\nwrote {OUT}/counts_*.tsv.gz, metadata_*.tsv, _manifest.tsv")
