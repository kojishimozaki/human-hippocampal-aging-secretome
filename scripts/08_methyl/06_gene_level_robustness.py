#!/usr/bin/env python
"""M1 methylome — GENE-LEVEL robustness (codex audit findings 2,3,4).

The primary 03 test unit is the REGION; secretome genes contribute many correlated regions (Astro
frozen-hit = 92 regions / 25 genes; broad sets = thousands of regions over fewer genes), so region-level
p-values are not gene-independent (finding 2). Here we (3) use the pre-specified covariate-adjusted age
slope (mCG~age+sex+log cov; the locked model) and (2) collapse to ONE value per gene, then run the
matched-null + a GENE-level sign-permutation. The sign-perm null is a **Rademacher sign-flip, two-sided**
(each gene's oriented value × random ±1; H0: gene-level oriented values are symmetric about 0) — a
one-sample test that the mean is non-zero, distinct from the region-level RNA-label shuffle in 03. (4) Micro is re-run using Micro2
only (present in all 40 donors; Micro1 is missing in 6 donors, all old 82-95 -> state-composition x age
confound). Out: results/methyl/gene_level_robustness.csv
"""
import os, sys, numpy as np, pandas as pd
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(__file__))
from importlib import import_module
S = import_module("00_setup")
PROJ, RES, CNT = S.PROJ, S.RES, f"{S.RES}/counts"
MIN_COV, NPERM, KNN = S.MIN_CG_COV, S.NPERM, 30
SEED = 42; rng = np.random.default_rng(SEED)

reg = pd.read_csv(f"{RES}/regions_meta.csv")
sf = pd.read_csv(f"{RES}/region_seqfeatures.csv").set_index("region_id").reindex(reg.region_id)
man = pd.read_csv(f"{RES}/allc_manifest.csv")
don = pd.read_csv(f"{RES}/donor_table.csv").set_index("donor")
rna = pd.read_csv(f"{PROJ}/results/de/de_GSE268609_per_celltype.csv")
froz = pd.read_csv(f"{PROJ}/results/validation/frozen_primary_signature.csv")
n = len(reg)
gt = pd.read_csv(f"{RES}/gene_tss.csv"); mid = ((reg.start + reg.end)//2).values
tss = np.full(n, np.nan)
for ch, sub in gt.groupby("chrom"):
    idx = np.where(reg.chrom.values == ch)[0]
    if not len(idx): continue
    ts = np.sort(sub.tss.values); pos = mid[idx]; j = np.clip(np.searchsorted(ts, pos), 0, len(ts)-1)
    tss[idx] = np.minimum(np.abs(ts[j]-pos), np.abs(ts[np.clip(j-1,0,len(ts)-1)]-pos))


def covadj_slope(srcs):
    """covariate-adjusted (age+sex+logcov) masked-WLS age slope per region, summing the given src suffixes."""
    sub = man[man.src_celltype.isin(srcs)]
    donors = sorted(sub.donor.unique())
    mc = np.zeros((n, len(donors))); cov = np.zeros((n, len(donors))); gcov = np.zeros(len(donors))
    for j, d in enumerate(donors):
        for _, r in sub[sub.donor == d].iterrows():
            f = f"{CNT}/{r.gsm}_{r.donor}_{r.src_celltype}.npz"
            if os.path.exists(f):
                z = np.load(f, allow_pickle=True); mc[:, j] += z["mc_cg"]; cov[:, j] += z["cov_cg"]; gcov[j] += int(z["global_cov"])
    donors = np.array(donors); age = don.loc[donors, "age"].values.astype(float)
    sex = (don.loc[donors, "sex"].values == "M").astype(float); logc = np.log(gcov)
    X = np.column_stack([np.ones(len(donors)), age, sex, logc])
    mCG = np.where(cov >= MIN_COV, mc/np.maximum(cov, 1), np.nan); w = (cov >= MIN_COV).astype(float)
    ev = (w[:, age < 40].sum(1) >= 4) & (w[:, age >= 60].sum(1) >= 8) & (w.sum(1) >= 20)
    XtWX = np.einsum("rd,dij->rij", w, np.einsum("di,dj->dij", X, X))
    XtWy = np.where(cov >= MIN_COV, mCG, 0.0) @ X
    slope = np.full(n, np.nan); s = np.where(ev)[0]
    slope[s] = np.linalg.solve(XtWX[s] + 1e-6*np.eye(4), XtWy[s][..., None])[..., 0][:, 1]
    base = np.nanmean(np.where(age < 40, mCG, np.nan), 1)
    return slope, ev, base


def gene_test(ct, srcs, label):
    slope, ev, base = covadj_slope(srcs)
    rc = reg.region_class.values; insec = reg.is_secretome.values; gene = reg.gene.values
    inreg = (rc == "promoter") | (rc == "proximal_peak")
    feat = np.column_stack([sf.gc.values, np.log1p(sf.cpg_density.values*1e4), np.log1p(reg.length.values),
                            np.log1p(tss), base])
    ok = np.isfinite(feat).all(1) & ev & inreg
    rna_ct = rna[rna.celltype == ct].set_index("gene"); froz_ct = set(froz[froz.celltype == ct].gene)
    out = []
    for setname, genes_in in [("frozen_hit", froz_ct), ("all_secretome", None)]:
        # gene-level oriented + features (mean over the gene's evaluable regulatory regions)
        def gene_vals(mask_genes_secretome):
            rows = {}
            for ri in np.where(ok & (insec == mask_genes_secretome))[0]:
                g = gene[ri]
                if mask_genes_secretome:
                    if g not in rna_ct.index: continue
                    if setname == "frozen_hit" and g not in genes_in: continue
                rows.setdefault(g, []).append(ri)
            return rows
        fg = gene_vals(True); bg = gene_vals(False)
        def collapse(rows, signed):
            G, OR, FT = [], [], []
            for g, idx in rows.items():
                idx = np.array(idx)
                if signed:
                    if g not in rna_ct.index: continue
                    s = np.sign(rna_ct.loc[g, "log2FoldChange"]);
                    s = s.iloc[0] if hasattr(s, "iloc") else s
                    if s == 0: continue
                    OR.append(-s * np.nanmean(slope[idx]))
                FT.append([np.nanmean(feat[idx, k]) for k in range(feat.shape[1])] + [np.log1p(len(idx))])
                G.append(g)
            return G, np.array(OR) if signed else None, np.array(FT)
        fG, fOR, fFT = collapse(fg, True); bG, _, bFT = collapse(bg, False)
        if len(fG) < 5 or len(bFT) < KNN+1:
            out.append(dict(celltype=label, set=setname, n_genes=len(fG), obs=np.nan, p_matched=np.nan, p_signperm=np.nan)); continue
        # background gene oriented needs a sign: use the gene's own rna sign in this celltype if available, else skip
        bg_or = {}
        for g, idx in [(g, np.array(i)) for g, i in bg.items()]:
            if g in rna_ct.index:
                s = np.sign(rna_ct.loc[g, "log2FoldChange"]); s = s.iloc[0] if hasattr(s, "iloc") else s
                if s != 0: bg_or[g] = -s * np.nanmean(slope[idx])
        bG2 = [g for g in bG if g in bg_or]; bOR = np.array([bg_or[g] for g in bG2])
        bFT2 = bFT[[bG.index(g) for g in bG2]]
        mu, sd = bFT2.mean(0), bFT2.std(0) + 1e-9
        tree = cKDTree((bFT2 - mu)/sd); _, nn = tree.query((fFT - mu)/sd, k=KNN)
        obs = float(np.mean(fOR))
        null = np.array([np.mean(bOR[nn[np.arange(len(fOR)), rng.integers(0, KNN, len(fOR))]]) for _ in range(NPERM)])
        p_match = (np.sum(null >= obs)+1)/(NPERM+1)
        sp = np.array([np.mean(fOR * rng.choice([-1, 1], len(fOR))) for _ in range(NPERM)])  # Rademacher sign-flip (two-sided)
        p_sign = (np.sum(np.abs(sp) >= abs(obs))+1)/(NPERM+1) if obs != 0 else 1.0
        out.append(dict(celltype=label, set=setname, n_genes=len(fG), obs=round(obs, 5),
                        p_matched=round(p_match, 4), p_signperm=round(p_sign, 4)))
    return out


res = []
for ct, srcs, lab in [("Astro", ["Astro1", "Astro2"], "Astro"), ("Micro", ["Micro1", "Micro2"], "Micro"),
                      ("OPC", ["OPC"], "OPC"), ("Oligo", ["Oligo"], "Oligo"),
                      ("Micro", ["Micro2"], "Micro2-only")]:
    res += gene_test(ct, srcs, lab)
df = pd.DataFrame(res); df.to_csv(f"{RES}/gene_level_robustness.csv", index=False)
pd.set_option("display.width", 200); print("=== GENE-LEVEL (covariate-adjusted slope; gene = unit; gene-level sign-perm two-sided) ===")
print(df.to_string(index=False))
print("\nInterpretation: Astro frozen_hit surviving gene-level = robust headline; broad sets attenuating = downgrade to region-level set effect;\nMicro2-only vs Micro = state-composition sensitivity.")
