#!/usr/bin/env python
"""M1 methylome — secretome regulatory-region mCG age-dynamics vs a feature-matched genomic null.

Pre-specified (RESOLUTION.md / plan M1; LOCKED before counting). Donor = unit; within-celltype; mCG.
Mirrors the A2 ATAC matched-null (scripts/04_atac/05*.py) so the two regulatory modalities are
methodologically consistent — same KDTree-NN matched null + 2000x permutation + sign-perm control.

Per evaluable niche celltype (Astro, Micro, Oligo, OPC; Endo excluded — too sparse):
  1. aggregate per-(donor x src) ALLC region counts -> donor x class (sum mc, cov over Astro1+Astro2 etc.)
  2. per region: mCG = mc/cov (>= MIN_CG_COV CpG basecalls/donor, else NA); YOUNG(<40,n<=10) vs OLD(>=60,n<=20)
     delta = mean_old - mean_young (>0 = hypermethylation with age); Welch t -> per-region p -> BH (bonus);
     baseline_mCG = young-group mean (a matching feature).
  3. matched-null enrichment, per region-class {promoter, gene_body, proximal_peak, regulatory=prom+peak}:
       oriented = (-delta) * sign(RNA log2FC of the linked gene in this celltype)   [>0 = hypo near UP-gene
         / hyper near DOWN-gene = RNA-concordant]; foreground = secretome regions, two sets
         (frozen_hit_linked: gene is a padj<0.1 hit in this celltype; all_secretome_linked).
       (i)  age-dynamism : mean(|delta| fg) vs matched non-secretome null
       (ii) direction    : mean(oriented fg) vs matched null given the SAME signs (p_matched)
       (iii) sign-perm    : shuffle signs among fg (negative control; p_signperm)
       match features (standardised, within class, non-secretome pool): gc, log cpg_density, log length,
         log tss_dist, baseline_mCG (+ mean_acc if results/methyl/region_accessibility.csv present).
         baseline_mCG is the within-methylation proxy for accessibility (mCG<->ATAC anti-correlated).
PASS (locked): direction-concordant (mean oriented>0) AND p_matched<0.05 in >=2 niche celltypes AND
  direction-specific (p_signperm<0.05). Per-region FDR survivors = bonus. Else honest null.
seed 42. Out: results/methyl/methyl_ageslope_<ct>.csv, results/methyl/secretome_region_methyl_enrichment.csv
"""
import os, glob, sys, warnings
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import t as tdist
warnings.filterwarnings("ignore", category=RuntimeWarning)

sys.path.insert(0, os.path.dirname(__file__))
from importlib import import_module
S = import_module("00_setup")  # CLASS_MAP, EVALUABLE, MIN_CG_COV, NPERM, paths
SEED = 42; rng = np.random.default_rng(SEED)
PROJ = S.PROJ; RES = S.RES; CNT = f"{S.RES}/counts"
MIN_COV = S.MIN_CG_COV
NPERM = S.NPERM
KNN = 50
MINY, MINO = 4, 8                       # min young / old donors (cov>=floor) for a region to be evaluable
CLASSES = ["promoter", "gene_body", "proximal_peak", "regulatory"]   # regulatory = promoter ∪ proximal_peak


def bh(p):
    p = np.asarray(p, float); ok = ~np.isnan(p); q = np.full(len(p), np.nan)
    pp = p[ok]; n = len(pp)
    if n == 0: return q
    o = np.argsort(pp); adj = np.empty(n)
    adj[o] = np.minimum.accumulate((pp[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    q[ok] = np.clip(adj, 0, 1); return q


def load_class_counts(man, target, reg_n):
    """Sum per-(donor x src) npz -> (mc, cov) arrays region x donor for one target class."""
    sub = man[man.target_class == target]
    donors = sorted(sub.donor.unique())
    mc = np.zeros((reg_n, len(donors))); cov = np.zeros((reg_n, len(donors)))
    have = np.zeros(len(donors), bool)
    for j, d in enumerate(donors):
        for _, r in sub[sub.donor == d].iterrows():
            f = f"{CNT}/{r.gsm}_{r.donor}_{r.src_celltype}.npz"
            if os.path.exists(f):
                z = np.load(f, allow_pickle=True)
                mc[:, j] += z["mc_cg"]; cov[:, j] += z["cov_cg"]; have[j] = True
    return np.array(donors)[have], mc[:, have], cov[:, have]


def yo_stats(mCG, young, old):
    Y, O = mCG[:, young], mCG[:, old]
    ny = (~np.isnan(Y)).sum(1); no = (~np.isnan(O)).sum(1)
    my = np.nanmean(Y, 1); mo = np.nanmean(O, 1)
    vy = np.nanvar(Y, 1, ddof=1); vo = np.nanvar(O, 1, ddof=1)
    delta = mo - my
    se = np.sqrt(vy / np.maximum(ny, 1) + vo / np.maximum(no, 1))
    with np.errstate(divide="ignore", invalid="ignore"):
        tval = delta / se
        df = se**4 / ((vy/np.maximum(ny,1))**2/np.maximum(ny-1,1) + (vo/np.maximum(no,1))**2/np.maximum(no-1,1))
        p = 2 * tdist.sf(np.abs(tval), np.clip(df, 1, None))
    evaluable = (ny >= MINY) & (no >= MINO)
    return delta, tval, p, my, evaluable


def enrich(li, sgn, delta, feat_z, fok, pool, label, ct, cls):
    """li=foreground region idx, sgn=RNA direction sign. Mirrors A2 enrich()+controls."""
    base = dict(celltype=ct, region_class=cls, set=label, n=len(li),
                obs_oriented=np.nan, p_matched=np.nan, p_signperm=np.nan,
                p_agedyn=np.nan, mag_ratio=np.nan, frac_concordant=np.nan)
    good = fok[li]; li, sgn = li[good], sgn[good]
    if len(li) < 5 or len(pool) < KNN + 1:
        base["n"] = len(li); return base
    oriented = (-delta[li]) * sgn                       # >0 = hypo near UP / hyper near DOWN = concordant
    obs = float(np.nanmean(oriented))
    tree = cKDTree(feat_z[pool])
    _, nn = tree.query(feat_z[li], k=KNN); nn_glob = pool[nn]            # (n_fg, KNN)
    # (ii) matched-null with same signs
    null = np.empty(NPERM)
    for j in range(NPERM):
        pick = nn_glob[np.arange(len(li)), rng.integers(0, KNN, len(li))]
        null[j] = np.nanmean((-delta[pick]) * sgn)
    p_matched = (np.sum(null >= obs) + 1) / (NPERM + 1)
    # (iii) sign-permutation negative control
    sp = np.array([np.nanmean((-delta[li]) * rng.permutation(sgn)) for _ in range(NPERM)])
    p_signperm = (np.sum(sp >= obs) + 1) / (NPERM + 1)
    # (i) age-dynamism: |delta| fg vs matched
    obs_mag = float(np.nanmean(np.abs(delta[li])))
    nullmag = np.empty(NPERM)
    for j in range(NPERM):
        pick = nn_glob[np.arange(len(li)), rng.integers(0, KNN, len(li))]
        nullmag[j] = np.nanmean(np.abs(delta[pick]))
    p_agedyn = (np.sum(nullmag >= obs_mag) + 1) / (NPERM + 1)
    mag_ratio = obs_mag / (np.nanmean(nullmag) + 1e-12)
    base.update(n=len(li), obs_oriented=round(obs, 5), p_matched=round(p_matched, 4),
                p_signperm=round(p_signperm, 4), p_agedyn=round(p_agedyn, 4),
                mag_ratio=round(mag_ratio, 3), frac_concordant=round(float(np.mean(oriented > 0)), 3))
    return base


def main():
    reg = pd.read_csv(f"{RES}/regions_meta.csv")
    sf = pd.read_csv(f"{RES}/region_seqfeatures.csv").set_index("region_id").reindex(reg.region_id)
    man = pd.read_csv(f"{RES}/allc_manifest.csv")
    don = pd.read_csv(f"{RES}/donor_table.csv").set_index("donor")
    rna = pd.read_csv(f"{PROJ}/results/de/de_GSE268609_per_celltype.csv")
    froz = pd.read_csv(f"{PROJ}/results/validation/frozen_primary_signature.csv")
    n = len(reg)

    # tss distance per region (nearest gene TSS to region midpoint; promoters ~0)
    gt = pd.read_csv(f"{RES}/gene_tss.csv")
    mid = ((reg.start + reg.end) // 2).values
    tss_dist = np.full(n, np.nan)
    for ch, sub in gt.groupby("chrom"):
        idx = np.where(reg.chrom.values == ch)[0]
        if not len(idx): continue
        ts = np.sort(sub.tss.values); pos = mid[idx]
        j = np.clip(np.searchsorted(ts, pos), 0, len(ts) - 1)
        d1 = np.abs(ts[j] - pos); d0 = np.abs(ts[np.clip(j - 1, 0, len(ts) - 1)] - pos)
        tss_dist[idx] = np.minimum(d0, d1)

    acc = None
    accf = f"{RES}/region_accessibility.csv"
    if os.path.exists(accf):
        acc = pd.read_csv(accf).set_index("region_id").reindex(reg.region_id)["mean_acc"].values

    enr_rows = []
    avail = sorted({c for c in S.EVALUABLE if (man.target_class == c).any()
                    and any(os.path.exists(f"{CNT}/{r.gsm}_{r.donor}_{r.src_celltype}.npz")
                            for _, r in man[man.target_class == c].iterrows())})
    print(f"celltypes with counts available: {avail}", flush=True)

    for ct in avail:
        donors, mc, cov = load_class_counts(man, ct, n)
        ages = don.loc[donors, "age"].values
        young = ages < 40; old = ages >= 60
        if young.sum() < MINY or old.sum() < MINO:
            print(f"  {ct}: donors young={young.sum()} old={old.sum()} -> skip"); continue
        mCG = np.where(cov >= MIN_COV, mc / np.maximum(cov, 1), np.nan)
        delta, tval, p, baseline, evaluable = yo_stats(mCG, young, old)
        padj = bh(np.where(evaluable, p, np.nan))
        # save per-region age table
        out = reg.copy()
        out["delta_mCG"] = delta; out["t"] = tval; out["p"] = p; out["padj"] = padj
        out["baseline_mCG"] = baseline; out["evaluable"] = evaluable
        out["n_donors"] = len(donors); out["gc"] = sf.gc.values; out["cpg_density"] = sf.cpg_density.values
        out.to_csv(f"{RES}/methyl_ageslope_{ct}.csv", index=False)
        print(f"  {ct}: donors={len(donors)} (Y{young.sum()}/O{old.sum()}) | evaluable regions={int(evaluable.sum())} "
              f"| per-region FDR<0.1={int((padj<0.1).sum())}", flush=True)

        # features for matching
        feat = np.column_stack([sf.gc.values, np.log1p(sf.cpg_density.values * 1e4),
                                np.log1p(reg.length.values), np.log1p(tss_dist), baseline])
        if acc is not None:
            feat = np.column_stack([feat, np.log1p(np.nan_to_num(acc))])
        fok = np.isfinite(feat).all(1) & evaluable
        fz = np.zeros_like(feat)
        fz[fok] = (feat[fok] - feat[fok].mean(0)) / (feat[fok].std(0) + 1e-12)

        rna_ct = rna[rna.celltype == ct].set_index("gene")
        froz_ct = set(froz[(froz.celltype == ct)].gene)

        # per region-class enrichment
        rc = reg.region_class.values; is_sec = reg.is_secretome.values; gene = reg.gene.values
        for cls in CLASSES:
            in_cls = (rc == "promoter") | (rc == "proximal_peak") if cls == "regulatory" else (rc == cls)
            # foreground sets
            rows_hit, rows_all = [], []
            fg_idx = np.where(in_cls & is_sec)[0]
            for ri in fg_idx:
                g = gene[ri]
                if g not in rna_ct.index: continue
                lfc = rna_ct.loc[g, "log2FoldChange"]
                lfc = lfc.iloc[0] if hasattr(lfc, "iloc") else lfc
                s = np.sign(lfc)
                if s == 0: continue
                rows_all.append((ri, s))
                if g in froz_ct: rows_hit.append((ri, s))
            pool = np.where(in_cls & (~is_sec) & fok)[0]
            for rows, lab in [(rows_hit, "frozen_hit_linked"), (rows_all, "all_secretome_linked")]:
                if not rows:
                    enr_rows.append(dict(celltype=ct, region_class=cls, set=lab, n=0, obs_oriented=np.nan,
                                         p_matched=np.nan, p_signperm=np.nan, p_agedyn=np.nan,
                                         mag_ratio=np.nan, frac_concordant=np.nan)); continue
                li = np.array([r[0] for r in rows]); sgn = np.array([r[1] for r in rows], float)
                enr_rows.append(enrich(li, sgn, delta, fz, fok, pool, lab, ct, cls))

    enr = pd.DataFrame(enr_rows)
    enr.to_csv(f"{RES}/secretome_region_methyl_enrichment.csv", index=False)
    if not len(enr):
        print("\n(no celltype counts available yet — pipeline validated, awaiting download)"); return
    print("\n===== M1 methylome enrichment — regulatory class (promoter ∪ proximal_peak) =====")
    print(enr[enr.region_class == "regulatory"].to_string(index=False))

    # Locked GO/NO-GO (RESOLUTION.md): per plan tests, on the REGULATORY class, per niche celltype —
    #   (i)  age-dynamism : broad all_secretome_linked, |delta| > matched null (p_agedyn<0.05)
    #   (ii) direction    : frozen_hit_linked, RNA-concordant above matched bg (obs>0 & p_matched<0.05)
    #   (iii) sign-specific: frozen_hit_linked dies under sign-permutation (p_signperm<0.05)
    # PASS = (ii) AND (iii) [headline frozen-hit direction] in >=2 niche glia, WITH (i) age-dynamism support.
    rh = enr[(enr["set"] == "frozen_hit_linked") & (enr.region_class == "regulatory")].set_index("celltype")
    ra = enr[(enr["set"] == "all_secretome_linked") & (enr.region_class == "regulatory")].set_index("celltype")
    print("\nper-celltype flags (regulatory):  i=age-dynamism(broad)  ii=concord(frozen-hit)  iii=sign-spec(frozen-hit)  +broad=concord(all-sec)")
    pass_ct = []
    for ct in rh.index:
        i = bool(ra.loc[ct, "p_agedyn"] < 0.05) if ct in ra.index else False
        ii = bool(rh.loc[ct, "obs_oriented"] > 0 and rh.loc[ct, "p_matched"] < 0.05)
        iii = bool(rh.loc[ct, "p_signperm"] < 0.05)
        broad = bool(ct in ra.index and ra.loc[ct, "obs_oriented"] > 0
                     and ra.loc[ct, "p_matched"] < 0.05 and ra.loc[ct, "p_signperm"] < 0.05)
        if ii and iii: pass_ct.append(ct)
        print(f"  {ct:6s} i={i!s:5} ii={ii!s:5} iii={iii!s:5} +broad={broad!s:5} "
              f"(hit n={int(rh.loc[ct,'n'])} oriented={rh.loc[ct,'obs_oriented']}; broad n={int(ra.loc[ct,'n']) if ct in ra.index else 0})")
    npass = len(pass_ct)
    print(f"\n>>> niche glia passing frozen-hit direction (ii & iii): {npass}  {pass_ct}")
    if npass >= 2:
        print(">>> M1 VERDICT: PASS — methylation supports the regulatory COORDINATION (not causation; mCG↔ATAC non-independent)")
    else:
        print(">>> M1 VERDICT: NO-GO / honest null on the frozen-hit direction test — report broad-set + age-dynamism results honestly")
    print("    (NOTE: single-celltype / partial runs are PRELIMINARY; final verdict needs all evaluable niche glia.)")
    print(f"wrote {RES}/secretome_region_methyl_enrichment.csv + methyl_ageslope_<ct>.csv")


if __name__ == "__main__":
    main()
