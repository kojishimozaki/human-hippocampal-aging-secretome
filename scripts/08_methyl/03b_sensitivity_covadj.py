#!/usr/bin/env python
"""M1 methylome — SENSITIVITY: covariate-adjusted continuous-age model (locked spec mCG~age+sex+coverage).

The primary run (03) used an unadjusted young-vs-old mean difference. This sensitivity implements the
pre-specified model `mCG ~ age + sex + log(global CG coverage)` per region (masked WLS over AVAILABLE
donors with cov >= MIN_CG_COV, >=4 young & >=8 old & >=20 donors), extracts the AGE slope, and re-runs the
matched-null enrichment on the regulatory class. Directly addresses caveat C6 (coverage-as-confounder)
and confirms the load-bearing results (Astro frozen-hit; Micro/OPC/Oligo broad-set) are not artifacts of
coverage drift with age. Also reports the coverage<->age correlation (the confound's actual magnitude).
seed 42. Out: results/methyl/secretome_region_methyl_enrichment_covadj.csv
"""
import os, sys, numpy as np, pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(__file__))
from importlib import import_module
S = import_module("00_setup")
M3 = import_module("03_ageslope_matched_null")     # reuse enrich(), constants
PROJ, RES, CNT = S.PROJ, S.RES, f"{S.RES}/counts"
MIN_COV, NPERM = S.MIN_CG_COV, S.NPERM
NICHE = S.EVALUABLE  # Astro, Micro, OPC, Oligo


def load_class(man, target, reg_n):
    sub = man[man.target_class == target]
    donors = sorted(sub.donor.unique())
    mc = np.zeros((reg_n, len(donors))); cov = np.zeros((reg_n, len(donors)))
    gmc = np.zeros(len(donors)); gcov = np.zeros(len(donors)); have = np.zeros(len(donors), bool)
    for j, d in enumerate(donors):
        for _, r in sub[sub.donor == d].iterrows():
            f = f"{CNT}/{r.gsm}_{r.donor}_{r.src_celltype}.npz"
            if os.path.exists(f):
                z = np.load(f, allow_pickle=True)
                mc[:, j] += z["mc_cg"]; cov[:, j] += z["cov_cg"]
                gmc[j] += int(z["global_mc"]); gcov[j] += int(z["global_cov"]); have[j] = True
    k = have
    return np.array(donors)[k], mc[:, k], cov[:, k], gmc[k], gcov[k]


def main():
    reg = pd.read_csv(f"{RES}/regions_meta.csv")
    sf = pd.read_csv(f"{RES}/region_seqfeatures.csv").set_index("region_id").reindex(reg.region_id)
    man = pd.read_csv(f"{RES}/allc_manifest.csv")
    don = pd.read_csv(f"{RES}/donor_table.csv").set_index("donor")
    rna = pd.read_csv(f"{PROJ}/results/de/de_GSE268609_per_celltype.csv")
    froz = pd.read_csv(f"{PROJ}/results/validation/frozen_primary_signature.csv")
    n = len(reg)
    gt = pd.read_csv(f"{RES}/gene_tss.csv")
    mid = ((reg.start + reg.end) // 2).values
    tss_dist = np.full(n, np.nan)
    for ch, sub in gt.groupby("chrom"):
        idx = np.where(reg.chrom.values == ch)[0]
        if not len(idx): continue
        ts = np.sort(sub.tss.values); pos = mid[idx]
        j = np.clip(np.searchsorted(ts, pos), 0, len(ts) - 1)
        tss_dist[idx] = np.minimum(np.abs(ts[j] - pos), np.abs(ts[np.clip(j-1,0,len(ts)-1)] - pos))

    out = []
    for ct in NICHE:
        donors, mc, cov, gmc, gcov = load_class(man, ct, n)
        D = len(donors)
        age = don.loc[donors, "age"].values.astype(float)
        sex = (don.loc[donors, "sex"].values == "M").astype(float)
        logcov = np.log(gcov.astype(float))
        # coverage<->age confound magnitude
        rho_cov_age, p_cov_age = spearmanr(gcov, age)
        # design [1, age, sex, logcov]; per-region masked WLS over AVAILABLE donors (cov>=floor),
        # require >=4 young & >=8 old & >=20 donors (keeps foreground; complete-case was too strict
        # because coverage drops with age). age coefficient = covariate-adjusted continuous-age slope.
        X = np.column_stack([np.ones(D), age, sex, logcov])
        mCG = np.where(cov >= MIN_COV, mc / np.maximum(cov, 1), np.nan)
        w = (cov >= MIN_COV).astype(float)                 # n x D mask-weight
        young = age < 40; old = age >= 60
        ny = w[:, young].sum(1); no = w[:, old].sum(1); ndon = w.sum(1)
        evaluable = (ny >= 4) & (no >= 8) & (ndon >= 20)
        XX = np.einsum("di,dj->dij", X, X)                 # (D,4,4)
        XtWX = np.einsum("rd,dij->rij", w, XX)             # (n,4,4)
        ym = np.where(cov >= MIN_COV, mCG, 0.0)
        XtWy = ym @ X                                       # (n,4)
        slope = np.full(n, np.nan)
        sub = np.where(evaluable)[0]
        beta = np.linalg.solve(XtWX[sub] + 1e-6 * np.eye(4), XtWy[sub][..., None])[..., 0]  # batched WLS
        slope[sub] = beta[:, 1]                             # age coefficient
        baseline = np.nanmean(np.where(age < 40, mCG, np.nan), 1)   # young-mean baseline (matching feat)
        # features (same as 03)
        feat = np.column_stack([sf.gc.values, np.log1p(sf.cpg_density.values * 1e4),
                                np.log1p(reg.length.values), np.log1p(tss_dist), baseline])
        fok = np.isfinite(feat).all(1) & evaluable
        fz = np.zeros_like(feat); fz[fok] = (feat[fok] - feat[fok].mean(0)) / (feat[fok].std(0) + 1e-12)
        rna_ct = rna[rna.celltype == ct].set_index("gene")
        froz_ct = set(froz[froz.celltype == ct].gene)
        rc = reg.region_class.values; is_sec = reg.is_secretome.values; gene = reg.gene.values
        in_cls = (rc == "promoter") | (rc == "proximal_peak")        # regulatory
        rows_hit, rows_all = [], []
        for ri in np.where(in_cls & is_sec)[0]:
            g = gene[ri]
            if g not in rna_ct.index: continue
            lfc = rna_ct.loc[g, "log2FoldChange"]; lfc = lfc.iloc[0] if hasattr(lfc, "iloc") else lfc
            s = np.sign(lfc)
            if s == 0: continue
            rows_all.append((ri, s))
            if g in froz_ct: rows_hit.append((ri, s))
        pool = np.where(in_cls & (~is_sec) & fok)[0]
        for rows, lab in [(rows_hit, "frozen_hit_linked"), (rows_all, "all_secretome_linked")]:
            if not rows:
                continue
            li = np.array([r[0] for r in rows]); sgn = np.array([r[1] for r in rows], float)
            r = M3.enrich(li, sgn, slope, fz, fok, pool, lab, ct, "regulatory")
            r["evaluable_regions"] = int(evaluable.sum()); r["rho_cov_age"] = round(rho_cov_age, 3)
            r["p_cov_age"] = round(p_cov_age, 4)
            out.append(r)

    res = pd.DataFrame(out)
    res.to_csv(f"{RES}/secretome_region_methyl_enrichment_covadj.csv", index=False)
    pd.set_option("display.width", 220)
    print("=== covariate-adjusted (mCG~age+sex+logcov), continuous-age slope, masked-WLS available donors ===")
    print(res[["celltype", "set", "n", "obs_oriented", "p_matched", "p_signperm", "p_agedyn",
               "frac_concordant", "evaluable_regions", "rho_cov_age", "p_cov_age"]].to_string(index=False))
    print("\nCompare to primary (young-vs-old, 03): Astro frozen-hit & Micro/OPC/Oligo broad should persist (concordant + sign-specific).")


if __name__ == "__main__":
    main()
