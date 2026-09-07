#!/usr/bin/env python
"""Arm-batch robustness for the published niche YA-vs-HA secretome signature.

N0 gate (docs/NEURON_NICHE_PLAN.md) found GSE268609 pools two preparation arms:
whole-hippocampus multiome (orig.ident 1-14) + DG-microdissected (orig.ident 15-39),
imbalanced across groups (DG-arm fraction YA 75% / HA 56%). The published niche DE used
~grp (+~pmi+grp sensitivity) but did NOT model arm. This script tests whether the frozen
60-hit signature survives ~arm+grp adjustment, exactly mirroring 04_aging_robustness.py's
~pmi+grp sensitivity. Donor-level pyDESeq2. Read-only inputs; writes one sensitivity CSV.

Output: results/de/aging_arm_sensitivity.csv
"""
import os, re, gzip
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_cache")
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import warnings; warnings.filterwarnings("ignore")

P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
NICHE = ["Astro", "Micro", "Oligo", "OPC", "Endo"]

# ---- (A) orig.ident -> arm, keyed to donor (SampleNumber). arm = WH (oid<=14) | DG (oid>=15) ----
meta = pd.read_csv(f"{P}/processed/per_dataset/GSE268609_metadata.csv", index_col=0)
donors = meta[["orig.ident", "SampleNumber", "Group"]].drop_duplicates().rename(columns={"orig.ident": "oid"})
donors["SampleNumber"] = donors.SampleNumber.astype(str)
donors["arm"] = np.where(donors.oid.astype(int) <= 14, "WH", "DG")
arm_map = dict(zip(donors.SampleNumber, donors.arm))
print("=== arm x Group (donor-level) ===")
print(pd.crosstab(donors.Group, donors.arm).to_string())

# ---- (B) frozen 60-hit signature (the published claim) -------------------------------------
frozen = pd.read_csv(f"{P}/results/validation/frozen_primary_signature.csv")
sec = set(pd.read_csv(f"{P}/refs/secretome_union.csv").gene)
print(f"\nfrozen hits: {len(frozen)} (per celltype: {frozen.celltype.value_counts().to_dict()})")

# ---- (C) anchor, YA/HA ---------------------------------------------------------------------
A = sc.read_h5ad(f"{P}/processed/per_dataset/GSE268609_anchor.h5ad")
A = A[A.obs.Group.isin(["YA", "HA"])].copy()
A.obs["arm"] = A.obs["donor_id"].astype(str).map(arm_map)
print(f"YA+HA cells: {A.n_obs}; donors: {A.obs.donor_id.nunique()}; arm NA: {A.obs.arm.isna().sum()}")

def pseudobulk(sub):
    X = sub.layers["counts"]; X = X.tocsr() if sp.issparse(X) else X
    don = sub.obs["donor_id"].astype(str).values; uniq = pd.unique(don)
    mat = np.vstack([np.asarray(X[don == d].sum(0)).ravel() for d in uniq])
    return pd.DataFrame(mat, index=uniq, columns=sub.var_names.astype(str))

def fit(pb, md, design):
    keep = (pb.sum(0) >= 10) & (pb.astype(bool).sum(0) >= max(3, int(0.5 * len(pb))))
    sub = pb.loc[:, keep]
    dds = DeseqDataSet(counts=sub.astype(int), metadata=md, design=design, quiet=True)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=["grp", "HA", "YA"], quiet=True); ds.summary()
    return ds.results_df

# ---- (D) ~grp (baseline refit) vs ~arm+grp, anchored on the frozen hits ---------------------
rows = []; hitrows = []
for ct in NICHE:
    sub = A[A.obs.celltype_l1.astype(str) == ct].copy()
    nc = sub.obs.groupby("donor_id", observed=True).size()
    sub = sub[sub.obs.donor_id.isin(nc[nc >= 10].index)].copy()
    pb = pseudobulk(sub)
    md = pd.DataFrame(index=pb.index)
    g = sub.obs.groupby("donor_id", observed=True).Group.first().astype(str)
    a = sub.obs.groupby("donor_id", observed=True).arm.first().astype(str)
    md["grp"] = g.reindex(pb.index).values
    md["arm"] = a.reindex(pb.index).values
    # arm x group balance within this celltype (collinearity watch)
    bal = pd.crosstab(md.grp, md.arm)
    estimable = (md.arm.nunique() >= 2) and (bal.min().min() >= 1)
    r_base = fit(pb, md, "~grp")
    r_arm = fit(pb, md, "~arm + grp") if estimable else None
    fz = frozen[frozen.celltype == ct]
    common = [gname for gname in fz.gene if gname in r_base.index and (r_arm is not None and gname in r_arm.index)]
    if r_arm is None:
        rows.append(dict(celltype=ct, n_donor=len(md), wh=int((md.arm == "WH").sum()),
                         dg=int((md.arm == "DG").sum()), n_hits=len(fz), estimable=False,
                         note="arm not estimable (collinear with group)"))
        continue
    cs_base = r_base.loc[r_base.index.isin(sec), "log2FoldChange"]
    cs_arm = r_arm.loc[r_arm.index.isin(sec), "log2FoldChange"]
    sh = cs_base.index.intersection(cs_arm.index)
    lfc_corr = float(cs_base.loc[sh].corr(cs_arm.loc[sh]))
    sign_ret = sum(np.sign(r_base.loc[gname, "log2FoldChange"]) == np.sign(r_arm.loc[gname, "log2FoldChange"]) for gname in common)
    padj_ret = sum((r_arm.loc[gname, "padj"] < 0.1) for gname in common if pd.notna(r_arm.loc[gname, "padj"]))
    rows.append(dict(celltype=ct, n_donor=len(md), wh=int((md.arm == "WH").sum()), dg=int((md.arm == "DG").sum()),
                     n_hits=len(fz), n_eval=len(common), estimable=True, lfc_corr=round(lfc_corr, 3),
                     sign_retained=sign_ret, padj01_retained=padj_ret))
    for gname in common:
        hitrows.append(dict(celltype=ct, gene=gname,
            base_lfc=round(float(r_base.loc[gname, "log2FoldChange"]), 3), base_padj=float(r_base.loc[gname, "padj"]),
            arm_lfc=round(float(r_arm.loc[gname, "log2FoldChange"]), 3), arm_padj=float(r_arm.loc[gname, "padj"]),
            sign_kept=bool(np.sign(r_base.loc[gname, "log2FoldChange"]) == np.sign(r_arm.loc[gname, "log2FoldChange"]))))

sens = pd.DataFrame(rows)
hits = pd.DataFrame(hitrows)
sens.to_csv(f"{P}/results/de/aging_arm_sensitivity.csv", index=False)
hits.to_csv(f"{P}/results/de/aging_arm_sensitivity_perhit.csv", index=False)
print("\n=== arm-adjusted sensitivity (~grp vs ~arm+grp), niche celltypes ===")
print(sens.to_string(index=False))
if len(hits):
    ne = len(hits); sr = int(hits.sign_kept.sum()); pr = int((hits.arm_padj < 0.1).sum())
    print(f"\n=== SUMMARY === evaluable frozen hits: {ne}/{len(frozen)} | sign retained {sr}/{ne} | "
          f"padj<0.1 retained {pr}/{ne} | lfc_corr range {sens.lfc_corr.min():.2f}-{sens.lfc_corr.max():.2f}")
print("\nwrote results/de/aging_arm_sensitivity.csv (+ _perhit.csv)")
