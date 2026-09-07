#!/usr/bin/env python
"""[B-1 / F-4-001, F-4b-001, F-4b-002] What are the NSC / Neuroblast / Immature nuclei really?

WHY
---
GSE268609's depositor labels three rare clusters NSC (877 nuclei), Neuroblast (1,309) and
Immature (1,366). Pass 4 showed they carry none of their lineage's canonical markers at
meaningful level — EOMES is 0 UMI in all 1,309 Neuroblast nuclei, MKI67/TOP2A are 0 in all 877
NSC nuclei, NSC's SOX2 (0.95 cp10k) sits below astrocytes' (1.53) — while retaining PLP1 at 52%
and SLC1A2 at 51% of the mature populations they should not resemble. The labels come from a
scANVI transfer with no abstain option whose reference includes fetal forebrain.

That question is already settled by an output nobody read: rare/marker_floor_reassignment.csv
(pass4b_06) re-scores every nucleus on ten lineage panels with a positive-evidence floor and
sends 99% of NSC to astrocyte+oligodendrocyte, 85% of Neuroblast to oligodendrocyte and 99.8%
of Immature to neuron, while Astro/mOli/mGC/Microglia keep their own labels at 83-99.8%.

So this script is not run to decide whether the labels are valid. It is run to decide the
MECHANISM — doublet, ambient, or mis-transfer — which determines how the manuscript should
describe them. Criteria are pre-registered in
audit_log/2026-08-26_p0_2_neurogenic_labels/RESOLUTION.md section 1.

RELATION TO THE AUDIT KIT
-------------------------
Derived from compute/pass4_07_neurogenic_lineage_test.py, which was written and never run.
Three declared changes (RESOLUTION section 1-1):

  1. BUG FIX. pass4_07:L99 does `T.X @ Xs.varm["PCs"]` with T subset to the 2,000 HVGs, but
     scanpy 1.11.5 returns varm["PCs"] at full width (n_vars, n_comps) with zeros outside the
     highly-variable mask, so that raises a dimension mismatch and the script cannot complete.
     Taking PCs[hv] fixes the shape; the projection is unchanged because the masked-out rows
     are zero.
  2. SEEDS. random_state=42 on HVG selection and PCA, and the module rng seeded 42 rather than
     0, per the project's single-seed rule. At NSIM=2000 the Monte-Carlo standard error on a
     reported fraction is at most sqrt(0.25/2000) = 1.1%.
  3. ADDED CONTROL. Real Astro and mGC nuclei are classified alongside the targets. Without a
     control showing the classifier does not call genuine single cell types heterotypic
     doublets, a positive call on NSC would not be interpretable. This is gate G-12.

Panels, NSIM, kNN(15), 30 PCs and the E-3/E-4/E-5 logic are otherwise verbatim.

RUNTIME ~40-60 min, dominated by loading three h5ads (9.4 + 8.8 + 11.6 GB). Reads processed/
read-only.

Out: results/audit_reruns/neurogenic/doublet_sim.csv        E-1 + G-12
     results/audit_reruns/neurogenic/qc_by_celltype.csv     E-3
     results/audit_reruns/neurogenic/fig1D_composition_loo.csv  E-4
     results/audit_reruns/neurogenic/external_search.csv    E-5
"""
import os
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_cache")
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, scipy.stats as sst
from sklearn.neighbors import KNeighborsClassifier
import warnings; warnings.filterwarnings("ignore")

SEED = 42
rng = np.random.default_rng(SEED)
np.random.seed(SEED)
NSIM = 2000
P = os.environ.get("PROJ", os.getcwd())
OUT = f"{P}/results/audit_reruns/neurogenic"
os.makedirs(OUT, exist_ok=True)

A = sc.read_h5ad(f"{P}/processed/per_dataset/GSE268609_anchor.h5ad")
C = (A.layers["counts"] if "counts" in A.layers else A.X).tocsr().astype(np.float32)
ct = A.obs.celltype_l1.astype(str).values
grp = A.obs.Group.astype(str).values
don = A.obs.donor_id.astype(str).values
print(f"anchor {A.shape}; cell types {sorted(set(ct))}", flush=True)

# ---------------------------------------------------------------- E-3 QC by cell type
umi = np.asarray(C.sum(1)).ravel()
gen = np.asarray((C > 0).sum(1)).ravel()
mtx = [i for i, g in enumerate(A.var_names) if str(g).startswith("MT-")]
pmt = np.asarray(C[:, mtx].sum(1)).ravel() / np.maximum(umi, 1) * 100
Q = (pd.DataFrame(dict(celltype=ct, umi=umi, genes=gen, pct_mt=pmt))
     .groupby("celltype").agg(n=("umi", "size"), umi_median=("umi", "median"),
                              umi_mean=("umi", "mean"), genes_median=("genes", "median"),
                              pct_mt_median=("pct_mt", "median")).reset_index())
Q["umi_median_rank"] = Q.umi_median.rank(ascending=False).astype(int)
Q["genes_median_rank"] = Q.genes_median.rank(ascending=False).astype(int)
Q.to_csv(f"{OUT}/qc_by_celltype.csv", index=False)
print("\n=== E-3 QC by cell type ===")
print(Q.sort_values("umi_median", ascending=False).round(2).to_string(index=False), flush=True)

# ---------------------------------------------------------------- E-1 synthetic doublets
PAIRS = {"Astro_x_Oligo": ("Astro", "Oligo"), "Neuron_x_Oligo": ("DG_GC", "Oligo"),
         "Astro_x_Astro": ("Astro", "Astro"), "Oligo_x_Oligo": ("Oligo", "Oligo")}
sim = {}
for name, (a, b) in PAIRS.items():
    ia = np.nonzero(ct == a)[0]; ib = np.nonzero(ct == b)[0]
    M = np.empty((NSIM, C.shape[1]), dtype=np.float32)
    for i in range(NSIM):
        M[i] = (np.asarray(C[rng.choice(ia)].todense()).ravel()
                + np.asarray(C[rng.choice(ib)].todense()).ravel())
    sim[name] = sp.csr_matrix(M)
    del M
    print(f"  simulated {name}: {sim[name].shape}", flush=True)

ref_names = list(sim) + ["Astro", "Oligo", "DG_GC", "Micro", "OPC"]
Xs, ys = [], []
for n in ref_names:
    if n in sim:
        Xs.append(sim[n]); ys += [n] * sim[n].shape[0]
    else:
        idx = rng.choice(np.nonzero(ct == n)[0], min(2000, int((ct == n).sum())), replace=False)
        Xs.append(C[idx]); ys += [n] * len(idx)
Xs = sc.AnnData(sp.vstack(Xs).astype(np.float32))
Xs.var_names = A.var_names
Xs.obs["y"] = ys
sc.pp.normalize_total(Xs, target_sum=1e4); sc.pp.log1p(Xs)
sc.pp.highly_variable_genes(Xs, n_top_genes=2000)
hv = Xs.var.highly_variable.values
sc.pp.scale(Xs, max_value=10)
sc.tl.pca(Xs, n_comps=30, random_state=SEED)
PCs = Xs.varm["PCs"]
if PCs.shape[0] == Xs.n_vars:          # scanpy returns full width with zeros outside the mask
    PCs = PCs[hv]
assert PCs.shape[0] == int(hv.sum()), f"PC row count {PCs.shape[0]} != n_hv {int(hv.sum())}"
knn = KNeighborsClassifier(15).fit(Xs.obsm["X_pca"], ys)
mu = Xs.var["mean"].values[hv]; sd = Xs.var["std"].values[hv]

TARGETS = ["NSC", "Neuroblast", "Immature", "Endo", "Ependymal", "Astro", "mGC"]
rows = []
for target in TARGETS:
    idx = np.nonzero(ct == target)[0]
    if not len(idx):
        print(f"  [E-1] {target}: absent from celltype_l1, skipped"); continue
    if len(idx) > 4000:                                  # controls only; keeps memory bounded
        idx = rng.choice(idx, 4000, replace=False)
    T = sc.AnnData(C[idx].astype(np.float32)); T.var_names = A.var_names
    sc.pp.normalize_total(T, target_sum=1e4); sc.pp.log1p(T)
    Tx = np.asarray(T[:, hv].X.todense())
    Tx = np.clip((Tx - mu) / sd, -10, 10)
    pred = knn.predict(Tx @ PCs)
    vc = pd.Series(pred).value_counts(normalize=True)
    het = float(vc.get("Astro_x_Oligo", 0) + vc.get("Neuron_x_Oligo", 0))
    rows.append(dict(target=target, n=len(idx), frac_heterotypic_doublet=het,
                     **{f"pred_{k}": float(vc.get(k, 0)) for k in ref_names}))
    print(f"  [E-1] {target} (n={len(idx)}): heterotypic-doublet {het:.3f} | "
          f"{vc.head(4).round(3).to_dict()}", flush=True)
D = pd.DataFrame(rows); D.to_csv(f"{OUT}/doublet_sim.csv", index=False)

# ---------------------------------------------------------------- E-4 leave-one-donor-out Fig 1D
comp = (pd.DataFrame(dict(donor=don, grp=grp, ct=ct)).query("grp in ['YA','HA']")
        .groupby(["donor", "grp", "ct"]).size().unstack(fill_value=0))
frac = comp.div(comp.sum(1), axis=0)
rows = []
donors = frac.index.get_level_values(0)
for target in frac.columns:
    for drop in ["<none>"] + sorted(set(donors)):
        f = frac[donors != drop] if drop != "<none>" else frac
        g = f.index.get_level_values(1)
        ya = np.log2(f.loc[g == "YA", target] + 1e-6); ha = np.log2(f.loc[g == "HA", target] + 1e-6)
        if len(ya) < 3 or len(ha) < 3:
            continue
        rows.append(dict(celltype=target, dropped_donor=drop, n_YA=len(ya), n_HA=len(ha),
                         log2fc=float(ha.mean() - ya.mean()),
                         p=float(sst.mannwhitneyu(ya, ha).pvalue)))
L = pd.DataFrame(rows); L.to_csv(f"{OUT}/fig1D_composition_loo.csv", index=False)
print("\n=== E-4 leave-one-donor-out on the Fig 1D composition test (raw MWU p) ===")
print(L[L.celltype.isin(["NSC", "Micro", "Neuroblast", "Immature"])]
      .groupby("celltype").p.describe()[["count", "min", "50%", "max"]].round(5).to_string(),
      flush=True)

# ---------------------------------------------------------------- E-5 external search
del Xs, sim
res = []
for coh, fn, acol in [("GSE278576", "GSE278576_aging.h5ad", "subclass"),
                      ("GSE186538", "GSE186538_franjic_human.h5ad", "cluster")]:
    B = sc.read_h5ad(f"{P}/processed/per_dataset/{fn}")
    if acol not in B.obs.columns:
        alt = [c for c in ("celltype_l1", "cluster", "subclass") if c in B.obs.columns]
        print(f"  [E-5] {coh}: '{acol}' absent, using '{alt[0] if alt else None}'")
        acol = alt[0] if alt else None
        if acol is None:
            del B; continue
    common = [g for g in np.asarray(A.var_names)[hv] if g in set(B.var_names)]
    Bc = B[:, common].copy(); del B
    Bc.X = (Bc.layers["counts"] if "counts" in Bc.layers else Bc.X).copy()
    sc.pp.normalize_total(Bc, target_sum=1e4); sc.pp.log1p(Bc)
    cidx = [list(A.var_names).index(g) for g in common]
    for target in ["NSC", "Neuroblast", "Immature"]:
        ti = np.nonzero(ct == target)[0]
        T = sc.AnnData(C[ti][:, cidx].astype(np.float32)); T.var_names = common
        sc.pp.normalize_total(T, target_sum=1e4); sc.pp.log1p(T)
        Td = np.asarray(T.X.todense())
        cen = Td.mean(0)
        own = Td @ cen / (np.linalg.norm(Td, axis=1) * np.linalg.norm(cen) + 1e-9)
        thr = float(np.percentile(own, 5))
        Xb = Bc.X
        num = np.asarray(Xb @ cen).ravel()
        den = np.sqrt(np.asarray(Xb.multiply(Xb).sum(1)).ravel()) * np.linalg.norm(cen)
        cos = num / np.maximum(den, 1e-9)
        lab = Bc.obs[acol].astype(str).values
        for lb in pd.unique(lab):
            m = lab == lb
            res.append(dict(cohort=coh, target=target, label=lb, n=int(m.sum()),
                            frac_above_primary_p5=float((cos[m] > thr).mean()),
                            median_cos=float(np.median(cos[m])), primary_p5=thr,
                            n_common_genes=len(common)))
    del Bc
    print(f"  [E-5] {coh} done ({len(common)} shared HVGs)", flush=True)
E = pd.DataFrame(res); E.to_csv(f"{OUT}/external_search.csv", index=False)
if len(E):
    print("\n=== E-5 external labels with frac_above_primary_p5 > 0.5 ===")
    hits = E[E.frac_above_primary_p5 > 0.5].sort_values("frac_above_primary_p5", ascending=False)
    print(hits.round(3).to_string(index=False) if len(hits) else "  (none)")

print("\n=== pre-registered gates (RESOLUTION section 1-2, 1-3) ===")
g12 = D[D.target.isin(["Astro", "mGC"])].frac_heterotypic_doublet.max() < 0.20 and \
      D[D.target.isin(["Endo", "Ependymal"])].frac_heterotypic_doublet.max() < 0.50
print(f"  G-12 classifier sanity -> {'PASS' if g12 else 'FAIL'}")
for t in ["Astro", "mGC", "Endo", "Ependymal"]:
    r = D[D.target == t]
    if len(r):
        print(f"       {t:11s} heterotypic-doublet {float(r.frac_heterotypic_doublet.iloc[0]):.3f}")
for t, lab in [("NSC", "N1"), ("Neuroblast", "N2")]:
    r = D[D.target == t]
    if len(r):
        h = float(r.frac_heterotypic_doublet.iloc[0])
        band = ">=50% doublet-supported" if h >= .5 else ("20-50% partial" if h >= .2 else "<20% not doublet-explained")
        print(f"  {lab} {t}: heterotypic-doublet {h:.3f} -> {band}")
print(f"\nwrote {OUT}/")
