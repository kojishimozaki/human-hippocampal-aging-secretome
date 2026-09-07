#!/usr/bin/env python
"""[B-2 / F-4-003] Decompose the 'Endothelial' bucket and re-run the YA-vs-HA DE on pure endothelium.

WHY
---
GSE268609's deposit has no mural / fibroblast / VLMC class among its 13 clusters, so pericytes,
smooth muscle and vascular leptomeningeal cells had nowhere to go but 'Endothelial'. Pass 4
(F-4-003) measured the consequence: the primary's Endo bucket carries about a third of the CLDN5
and 6.7x the PDGFRB of GSE186538's endothelium, and its two largest UP hits are markers of
subtypes the primary cannot resolve (APOD -> VLMC; DKK2 -> the aEndo DKK2 FBLN5 subtype). The
audit's own prescribed fix was to run `compute/pass4_05_endo_bucket_decomposition.py` and report
how many of the 11 frozen endothelial hits keep padj<0.1 once the bucket is purified.

That script was written and never executed: its three outputs are absent from audit/cluster/, it
has no row in audit/pass6/run_all_status.tsv, and compute/run_all.sh never references it. This
script carries it out.

RELATION TO pass4_05
--------------------
The panels, the argmax assignment rule, the frozen hit list, the pseudobulk, the gene filter, the
>=10 nuclei/donor and >=4 donors/group gates and the three DE models are verbatim from pass4_05.
Three declared changes (audit_log/2026-08-25_p0_3_endo_bucket/RESOLUTION.md, section 1-1):

  1. random_state=42 on pca / neighbors / leiden / score_genes. pass4_05 sets no seed at all,
     which is the hidden randomness the audit itself flagged as F-0-007; without it the
     partition does not reproduce.
  2. top_score is reported per subcluster, plus how many nuclei sit in a subcluster whose top
     panel score is <= 0 -- i.e. whose label has no positive evidence behind it. The primary
     rule stays argmax; this is additional reporting only.
  3. All three resolutions are reported; resolution 0.6 governs, which is pass4_05's own
     pre-existing choice of working partition.

WHAT LOSING HITS DOES AND DOES NOT MEAN
---------------------------------------
Purification removes nuclei and can drop donors below the >=10 gate, so a fall in the number of
padj<0.1 hits is not by itself evidence of a composition artefact -- the same trap F-5-007 set.
Sign retention, the secretome log2FC correlation against the full bucket, and the surviving donor
counts are therefore reported next to every hit count.

seed 42. RUNTIME ~12-15 min, dominated by loading the 8.8 GB anchor; the clustering itself is
873 nuclei and is trivial. Reads processed/ read-only; writes nothing outside results/.

Out: results/audit_reruns/endo/subclusters.csv          per resolution x subcluster
     results/audit_reruns/endo/composition_by_donor.csv per donor, resolution 0.6
     results/audit_reruns/endo/hits_under_purification.csv  the 11 hits x 3 models
     results/audit_reruns/endo/model_summary.csv       donor/nuclei counts + sign + lfc corr
"""
import os
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_cache")
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, scipy.stats as sst
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import warnings; warnings.filterwarnings("ignore")

SEED = 42
np.random.seed(SEED)
P = os.environ.get("PROJ", os.getcwd())
OUT = f"{P}/results/audit_reruns/endo"
os.makedirs(OUT, exist_ok=True)
GOV_RES = 0.6                      # pass4_05's own working partition
MIN_CELLS, MIN_DONORS = 10, 4      # verbatim from pass4_05 / the RNA leg

PANELS = {
    "Endo":      ["CLDN5", "FLT1", "PECAM1", "VWF", "CDH5", "ABCB1"],
    "Pericyte":  ["PDGFRB", "RGS5", "KCNJ8", "ABCC9", "NOTCH3", "HIGD1B"],
    "SMC":       ["ACTA2", "TAGLN", "MYH11", "MYL9"],
    "VLMC":      ["DCN", "COL1A1", "COL1A2", "LUM", "SLC6A13"],
    "Ependymal": ["FOXJ1", "PIFO", "TTR", "HDC"],
    "aEndo":     ["DKK2", "FBLN5", "SEMA3G", "GJA5"],
}
CORE_ENDO = ["CLDN5", "FLT1", "PECAM1"]        # for the G-11 marker-sanity gate
HITS = ["PRSS3", "APOD", "PLAT", "KDR", "LAMA4", "KITLG", "EFEMP1", "DKK2", "MMP28", "APP", "EDN3"]

frozen = pd.read_csv(f"{P}/results/validation/frozen_primary_signature.csv")
fz_endo = frozen[frozen.celltype == "Endo"].set_index("gene")
assert sorted(fz_endo.index) == sorted(HITS), "frozen Endo hits no longer match the hard-coded list"
sec = set(pd.read_csv(f"{P}/refs/secretome_union.csv").gene)

# ------------------------------------------------------------------ load and embed
A = sc.read_h5ad(f"{P}/processed/per_dataset/GSE268609_anchor.h5ad")
E = A[A.obs.celltype_l1.astype(str) == "Endo"].copy()
del A
E = E[E.obs.Group.astype(str).isin(["YA", "HA"])].copy()
print(f"Endo YA+HA nuclei: {E.shape}, {E.obs.Group.value_counts().to_dict()}, "
      f"{E.obs.donor_id.nunique()} donors", flush=True)

E.X = E.layers["counts"].copy()
sc.pp.normalize_total(E, target_sum=1e4); sc.pp.log1p(E); E.raw = E
sc.pp.highly_variable_genes(E, n_top_genes=2000, batch_key="donor_id")
E2 = E[:, E.var.highly_variable].copy()
sc.pp.scale(E2, max_value=10)
sc.tl.pca(E2, n_comps=20, random_state=SEED)
sc.pp.neighbors(E2, n_neighbors=15, n_pcs=20, random_state=SEED)

for pn, gs in PANELS.items():
    gs = [g for g in gs if g in E.var_names]
    sc.tl.score_genes(E, gs, score_name=f"sc_{pn}", use_raw=True, random_state=SEED)

def marker_mean(mask, genes):
    gs = [g for g in genes if g in E.var_names]
    X = E[:, gs].X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    return float(X[mask].mean()) if mask.sum() else np.nan

rows = []
for r in [0.3, GOV_RES, 1.0]:
    k = f"sub_{r}"
    sc.tl.leiden(E2, resolution=r, key_added=k, flavor="igraph", n_iterations=2,
                 directed=False, random_state=SEED)
    E.obs[k] = E2.obs[k].values
    m = E.obs.groupby(k, observed=True)[[f"sc_{p}" for p in PANELS]].mean()
    for cl in m.index:
        sel = (E.obs[k] == cl).values
        nd = int(E.obs.loc[sel, "donor_id"].nunique())
        rows.append(dict(
            resolution=r, subcluster=cl, n=int(sel.sum()), n_donors=nd,
            assigned=m.loc[cl].idxmax().replace("sc_", ""),
            top_score=float(m.loc[cl].max()),
            positive_evidence=bool(m.loc[cl].max() > 0),
            unstable=bool(sel.sum() < 30 or nd < 4),
            frac_YA=float((E.obs.loc[sel, "Group"] == "YA").mean()),
            **{p: float(m.loc[cl, f"sc_{p}"]) for p in PANELS},
            **{f"mean_{g}": marker_mean(sel, [g]) for g in CORE_ENDO}))
SUB = pd.DataFrame(rows)
SUB.to_csv(f"{OUT}/subclusters.csv", index=False)
print("\n=== subclusters (all three resolutions) ===")
print(SUB.round(3).to_string(index=False), flush=True)

# ------------------------------------------------------------------ governing partition
K = f"sub_{GOV_RES}"
gov = SUB[SUB.resolution == GOV_RES].set_index("subcluster")
E.obs["subtype"] = E.obs[K].map(gov["assigned"]).values
n_noev = int(E.obs[K].map(~gov["positive_evidence"]).sum())
print(f"\nresolution {GOV_RES}: {len(gov)} subclusters; "
      f"{n_noev}/{E.n_obs} nuclei ({100*n_noev/E.n_obs:.1f}%) sit in a subcluster whose top "
      f"panel score is <= 0 (label has no positive evidence)")

# G-11 marker sanity
endo_cl = gov.index[gov.assigned == "Endo"]
oth_cl = gov.index[gov.assigned != "Endo"]
g11 = {}
for g in CORE_ENDO:
    e = gov.loc[endo_cl, f"mean_{g}"].max() if len(endo_cl) else np.nan
    o = gov.loc[oth_cl, f"mean_{g}"].max() if len(oth_cl) else -np.inf
    g11[g] = (e, o, bool(e >= o))
print("\nG-11 marker sanity (Endo-assigned subclusters must top the core markers):")
for g, (e, o, ok) in g11.items():
    print(f"  {g:8s} max in Endo-assigned {e:.4f} vs max elsewhere {o:.4f} -> "
          f"{'PASS' if ok else 'FAIL'}")
G11 = all(v[2] for v in g11.values())

# B-2's specific mechanism: which subtype carries APOD and DKK2?
print("\nB-2 mechanism check (subcluster with the highest mean expression):")
for g in ["APOD", "DKK2"]:
    if g in E.var_names:
        mm = {cl: marker_mean((E.obs[K] == cl).values, [g]) for cl in gov.index}
        top = max(mm, key=mm.get)
        print(f"  {g:5s} highest in subcluster {top} (assigned '{gov.loc[top,'assigned']}', "
              f"n={gov.loc[top,'n']}), mean {mm[top]:.4f}; "
              f"Endo-assigned max {max([mm[c] for c in endo_cl], default=np.nan):.4f}")

# ------------------------------------------------------------------ composition x group
comp = (E.obs.groupby(["donor_id", "Group"], observed=True)["subtype"]
        .value_counts(normalize=True).unstack(fill_value=0).reset_index())
comp["nonEndo_frac"] = 1 - (comp["Endo"] if "Endo" in comp.columns else 0)
n_tot = E.obs.groupby("donor_id", observed=True).size()
comp["n_cells"] = comp.donor_id.map(n_tot).values
comp.to_csv(f"{OUT}/composition_by_donor.csv", index=False)
ya = comp[comp.Group == "YA"].nonEndo_frac
ha = comp[comp.Group == "HA"].nonEndo_frac
mwu_p = sst.mannwhitneyu(ha, ya).pvalue if len(ya) > 1 and len(ha) > 1 else np.nan
overall_nonendo = float((E.obs.subtype != "Endo").mean())
print(f"\nE1 bucket purity : non-endothelial nuclei = {overall_nonendo*100:.1f}% "
      f"({int((E.obs.subtype != 'Endo').sum())}/{E.n_obs})  -> "
      f"{'mixed bucket CONFIRMED' if overall_nonendo >= 0.20 else 'below the 20% bar'}")
print(f"E2 composition   : nonEndo_frac YA median {ya.median():.3f} (n={len(ya)}) vs "
      f"HA median {ha.median():.3f} (n={len(ha)}), donor-level MWU p = {mwu_p:.4g}")

# ------------------------------------------------------------------ DE
def pseudobulk_de(ad, extra_cov=None, tag=""):
    X = ad.layers["counts"]
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    d = ad.obs["donor_id"].astype(str).values
    vc = pd.Series(d).value_counts()
    keep = vc[vc >= MIN_CELLS].index
    dropped = sorted(set(vc.index) - set(keep), key=str)
    m = np.isin(d, keep)
    d, X = d[m], X[m]
    if not len(d):
        print(f"  [{tag}] no donor reaches {MIN_CELLS} nuclei"); return None, {}
    uq = pd.unique(d)
    pb = pd.DataFrame(np.vstack([np.asarray(X[d == u].sum(0)).ravel() for u in uq]),
                      index=uq, columns=ad.var_names.astype(str))
    meta = ad.obs[m].groupby("donor_id", observed=True).agg(grp=("Group", "first")).loc[uq]
    if extra_cov is not None:
        meta = meta.join(extra_cov).dropna()
        pb = pb.loc[meta.index]
    per = meta.grp.value_counts()
    info = dict(n_donor=len(meta), n_YA=int(per.get("YA", 0)), n_HA=int(per.get("HA", 0)),
                n_nuclei=int(m.sum()), donors_dropped=";".join(dropped))
    if min(per.get("YA", 0), per.get("HA", 0)) < MIN_DONORS:
        print(f"  [{tag}] NOT EVALUABLE: per-group donors {per.to_dict()} (<{MIN_DONORS}/group)")
        info["evaluable"] = False
        return None, info
    kg = (pb.sum(0) >= 10) & (pb.astype(bool).sum(0) >= max(3, int(0.5 * len(pb))))
    pb = pb.loc[:, kg]
    design = "~grp" if extra_cov is None else f"~{extra_cov.columns[0]} + grp"
    dds = DeseqDataSet(counts=pb.astype(int), metadata=meta, design=design, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["grp", "HA", "YA"], quiet=True); st.summary()
    info.update(evaluable=True, design=design, n_genes=pb.shape[1])
    print(f"  [{tag}] {design}: {len(meta)} donors (YA {info['n_YA']} / HA {info['n_HA']}), "
          f"{info['n_nuclei']} nuclei, {pb.shape[1]} genes", flush=True)
    return st.results_df.copy(), info


print("\n=== DE models ===")
res, meta_info = {}, {}
res["full_bucket"], meta_info["full_bucket"] = pseudobulk_de(E, tag="full_bucket")
pure = E[E.obs.subtype == "Endo"].copy()
print(f"  endothelial-only nuclei: {pure.n_obs} ({pure.obs.donor_id.nunique()} donors)")
res["endo_only"], meta_info["endo_only"] = pseudobulk_de(pure, tag="endo_only")
cov = comp.set_index("donor_id")[["nonEndo_frac"]]
res["nonEndo_adjusted"], meta_info["nonEndo_adjusted"] = pseudobulk_de(E, extra_cov=cov,
                                                                      tag="nonEndo_adjusted")

# G-10 reproduction gate on full_bucket
base = res["full_bucket"]
if base is not None:
    common = [g for g in HITS if g in base.index]
    kept = int(sum(base.loc[g, "padj"] < 0.1 for g in common if pd.notna(base.loc[g, "padj"])))
    lf_b = np.array([base.loc[g, "log2FoldChange"] for g in common])
    lf_f = np.array([fz_endo.loc[g, "log2FoldChange"] for g in common])
    r10 = float(np.corrcoef(lf_b, lf_f)[0, 1]); d10 = float(np.abs(lf_b - lf_f).max())
    G10 = (kept == 11) and (r10 > 0.999) and (d10 < 0.05)
    print(f"\nG-10 reproduction : full_bucket keeps {kept}/11 at padj<0.1, "
          f"log2FC vs frozen r = {r10:.6f}, max |diff| = {d10:.5f} -> {'PASS' if G10 else 'FAIL'}")
else:
    G10 = False
    print("\nG-10 reproduction : FAIL (full_bucket not evaluable)")

hit_rows, summ_rows = [], []
for name, r in res.items():
    info = meta_info.get(name, {})
    if r is None:
        summ_rows.append(dict(model=name, **info, n_hits_padj01=np.nan,
                              sign_retained=np.nan, lfc_corr_secretome=np.nan))
        continue
    common = [g for g in HITS if g in r.index]
    kept = int(sum(r.loc[g, "padj"] < 0.1 for g in common if pd.notna(r.loc[g, "padj"])))
    sign_ok = int(sum(np.sign(r.loc[g, "log2FoldChange"]) == np.sign(fz_endo.loc[g, "log2FoldChange"])
                      for g in common))
    if base is not None and name != "full_bucket":
        sb = base.loc[base.index.isin(sec), "log2FoldChange"]
        sr = r.loc[r.index.isin(sec), "log2FoldChange"]
        sh = sb.index.intersection(sr.index)
        corr = float(sb.loc[sh].corr(sr.loc[sh]))
    else:
        corr = 1.0
    summ_rows.append(dict(model=name, **info, n_hits_evaluated=len(common),
                          n_hits_padj01=kept, sign_retained=sign_ok, lfc_corr_secretome=round(corr, 4)))
    for g in common:
        hit_rows.append(dict(model=name, gene=g,
                             frozen_lfc=float(fz_endo.loc[g, "log2FoldChange"]),
                             log2FoldChange=float(r.loc[g, "log2FoldChange"]),
                             lfcSE=float(r.loc[g, "lfcSE"]), baseMean=float(r.loc[g, "baseMean"]),
                             pvalue=float(r.loc[g, "pvalue"]), padj=float(r.loc[g, "padj"])))
    print(f"  {name}: {kept}/{len(common)} frozen Endo hits keep padj<0.1; "
          f"sign {sign_ok}/{len(common)}; secretome lfc corr {corr:.4f}")

H = pd.DataFrame(hit_rows); H.to_csv(f"{OUT}/hits_under_purification.csv", index=False)
S = pd.DataFrame(summ_rows); S.to_csv(f"{OUT}/model_summary.csv", index=False)
print("\n" + S.to_string(index=False))
if len(H):
    print("\n=== the 11 frozen endothelial hits, per model ===")
    print(H.pivot_table(index="gene", columns="model", values="padj").round(4).to_string())

print("\n=== pre-registered verdicts (section 1-3) ===")
print(f"  G-10 {'PASS' if G10 else 'FAIL'} | G-11 {'PASS' if G11 else 'FAIL'}")
print(f"  E1 non-endothelial fraction = {overall_nonendo*100:.1f}% "
      f"({'>=20% -> mixed bucket confirmed' if overall_nonendo >= 0.20 else '<20%'})")
print(f"  E2 composition x group MWU p = {mwu_p:.4g} "
      f"({'p<0.05 -> composition contrast' if mwu_p < 0.05 else 'no group shift detected'})")
eo = S[S.model == "endo_only"].iloc[0]
if pd.isna(eo.n_hits_padj01):
    print("  E3 endo_only = NOT EVALUABLE -> section 1-5 applies "
          "(manuscript CLOSED-CHANGE, science STILL-OPEN)")
else:
    n = int(eo.n_hits_padj01)
    band = ">=8 CLOSED-OK" if n >= 8 else ("4-7 CLOSED-CHANGE" if n >= 4 else "<=3 CLOSED-CHANGE (major)")
    print(f"  E3 endo_only keeps {n}/11 -> {band}")
print(f"\nwrote {OUT}/")
