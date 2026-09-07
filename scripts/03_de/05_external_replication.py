#!/usr/bin/env python
"""A1 — test the FROZEN primary niche-secretome aging signature in an external cohort's
donor-pseudobulk DE. This *tests a fixed hypothesis*; it never re-discovers a signature.

Pre-specified in audit_log/2026-06-04_in_silico_hardening/RESOLUTION.md §A1 (locked before any
external DE was run). Donor-level only (the external DE is itself donor pseudobulk). seed 42.

Inputs
  results/validation/frozen_primary_signature.csv        60 hits: gene,celltype,log2FoldChange,direction,padj
  results/validation/frozen_primary_secretome_logfc.csv  all secretome×niche primary log2FC (fixed primary side)
  --external results/de/de_<COHORT>_per_celltype.csv      external donor-pseudobulk DE (same engine as primary)

Metrics (per niche celltype AND pooled), each with a permutation null (>=1000x):
  (i)   directional: external lfc of frozen UP hits > DOWN hits (one-sided MWU, plan-literal),
        plus a robust signed one-sample Wilcoxon on s_g = lfc_ext * sign(primary_dir); sign-flip null.
  (ii)  Spearman(primary lfc, external lfc) over shared secretome genes (one-sided rho>0); gene-shuffle null.
  (iii) 60-hit sign-concordance vs the cohort's OWN background (non-hit secretome concordance):
        binomial(hit k/n vs bg rate) + Fisher(hit vs non-hit); random-secretome-set null.
A metric PASSES if its parametric one-sided p<0.05 AND the permutation null_p<0.05.
OVERALL GO ('transcriptomic reproduction in a second cohort') if >=2 of (i pooled)/(ii >=2 celltypes)/(iii pooled)
PASS. Otherwise an honest null is reported (the protein-level external validation stays load-bearing).

Outputs
  results/validation/replication_<COHORT>.csv   tidy: cohort,scope,metric,statistic,p,null_p,n,pass,note
  figures/figure3/replication_<COHORT>.pdf       per-celltype primary-vs-external lfc scatter
"""
import argparse, os
import numpy as np
import pandas as pd
from scipy import stats

SEED = 42
NICHE = ["Astro", "Micro", "Endo", "OPC", "Oligo"]
NPERM = 2000


def onesided_spearman_gt0(x, y):
    """Spearman rho with one-sided p for the alternative rho>0."""
    if len(x) < 5:
        return np.nan, np.nan, len(x)
    rho, p2 = stats.spearmanr(x, y)
    if np.isnan(rho):
        return np.nan, np.nan, len(x)
    p1 = p2 / 2 if rho > 0 else 1 - p2 / 2
    return rho, p1, len(x)


def directional(merged, rng):
    """metric (i) on frozen hits merged with external lfc (col lfc_ext, dir col 'direction')."""
    up = merged.loc[merged.direction == "UP", "lfc_ext"].dropna().values
    dn = merged.loc[merged.direction == "DOWN", "lfc_ext"].dropna().values
    s = merged["lfc_ext"].values * np.where(merged.direction.values == "UP", 1.0, -1.0)
    s = s[~np.isnan(s)]
    out = {"n_up": len(up), "n_down": len(dn), "n": len(s), "mean_signed": float(np.mean(s)) if len(s) else np.nan}
    # plan-literal: one-sided MWU UP>DOWN
    if len(up) >= 1 and len(dn) >= 1:
        out["mwu_U"], out["mwu_p"] = stats.mannwhitneyu(up, dn, alternative="greater")
    else:
        out["mwu_U"], out["mwu_p"] = np.nan, np.nan
    # robust: signed one-sample Wilcoxon (handles single-direction celltypes)
    if len(s) >= 6 and np.any(s != 0):
        try:
            out["wilcox_p"] = stats.wilcoxon(s, alternative="greater").pvalue
        except ValueError:
            out["wilcox_p"] = np.nan
    else:
        out["wilcox_p"] = np.nan
    # sign-flip permutation null on mean_signed
    if len(s):
        flips = rng.choice([-1.0, 1.0], size=(NPERM, len(s)))
        null = (flips * s).mean(axis=1)
        out["null_p"] = float((null >= out["mean_signed"]).mean())
    else:
        out["null_p"] = np.nan
    return out


def concordance(hits_ext, nonhits_ext, rng):
    """metric (iii): hit sign-concordance vs the external cohort's non-hit secretome background."""
    hc = int(hits_ext["concord"].sum()); hn = int(hits_ext["concord"].count())
    bc = int(nonhits_ext["concord"].sum()); bn = int(nonhits_ext["concord"].count())
    hit_rate = hc / hn if hn else np.nan
    bg_rate = bc / bn if bn else np.nan
    out = {"n_hit": hn, "hit_concord": hit_rate, "n_bg": bn, "bg_concord": bg_rate}
    out["binom_p"] = stats.binomtest(hc, hn, bg_rate, alternative="greater").pvalue if hn and bn else np.nan
    if hn and bn:
        tbl = [[hc, hn - hc], [bc, bn - bc]]
        out["fisher_or"], out["fisher_p"] = stats.fisher_exact(tbl, alternative="greater")
    else:
        out["fisher_or"], out["fisher_p"] = np.nan, np.nan
    # null: random secretome sets of size n_hit drawn from the full secretome pool (hits+nonhits)
    pool = pd.concat([hits_ext["concord"], nonhits_ext["concord"]]).values
    if hn and len(pool) > hn:
        idx = np.array([rng.choice(len(pool), size=hn, replace=False) for _ in range(NPERM)])
        nulls = pool[idx].mean(axis=1)
        out["null_p"] = float((nulls >= hit_rate).mean())
    else:
        out["null_p"] = np.nan
    return out


def perm_spearman_null(x, y, obs_rho, rng):
    if len(x) < 5 or np.isnan(obs_rho):
        return np.nan
    nulls = np.empty(NPERM)
    yy = np.array(y)
    for i in range(NPERM):
        nulls[i] = stats.spearmanr(x, rng.permutation(yy))[0]
    return float((nulls >= obs_rho).mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--external", required=True, help="results/de/de_<COHORT>_per_celltype.csv")
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--proj", default=os.environ.get("PROJ", os.getcwd()))
    a = ap.parse_args()
    P = a.proj
    rng = np.random.default_rng(SEED)

    hits = pd.read_csv(f"{P}/results/validation/frozen_primary_signature.csv")
    full = pd.read_csv(f"{P}/results/validation/frozen_primary_secretome_logfc.csv")
    ext = pd.read_csv(a.external)[["gene", "celltype", "log2FoldChange", "padj"]].rename(
        columns={"log2FoldChange": "lfc_ext", "padj": "padj_ext"})

    rows = []
    scatter = {}
    ext_celltypes = set(ext.celltype.unique())

    def add(scope, metric, stat, p, null_p, n, passed, note=""):
        rows.append(dict(cohort=a.cohort, scope=scope, metric=metric, statistic=stat,
                         p=p, null_p=null_p, n=n, pass_=passed, note=note))

    # ---- per-celltype ----
    ii_pass_celltypes = 0
    for ct in NICHE:
        if ct not in ext_celltypes:
            add(ct, "ALL", np.nan, np.nan, np.nan, 0, False, "NOT_EVALUABLE (celltype absent in external DE)")
            continue
        h = hits[hits.celltype == ct].merge(ext[ext.celltype == ct], on=["gene", "celltype"], how="inner")
        f = full[full.celltype == ct].merge(ext[ext.celltype == ct], on=["gene", "celltype"], how="inner")
        f["concord"] = np.sign(f.log2FoldChange) == np.sign(f.lfc_ext)
        hitkeys = set(zip(h.gene, h.celltype))
        f["is_hit"] = [(g, c) in hitkeys for g, c in zip(f.gene, f.celltype)]
        h2 = f[f.is_hit].copy(); nh = f[~f.is_hit].copy()
        scatter[ct] = f
        # (i)
        di = directional(h.assign(direction=h.direction), rng)
        passed_i = (not np.isnan(di["mwu_p"]) and di["mwu_p"] < 0.05 and di["null_p"] < 0.05)
        add(ct, "(i)directional", di.get("mwu_U", np.nan), di["mwu_p"], di["null_p"], di["n"],
            passed_i, f"mean_signed={di['mean_signed']:.3f} wilcox_p={di['wilcox_p']:.3g} nUP={di['n_up']} nDOWN={di['n_down']}")
        # (ii)
        rho, p1, n = onesided_spearman_gt0(f.log2FoldChange.values, f.lfc_ext.values)
        np_ = perm_spearman_null(f.log2FoldChange.values, f.lfc_ext.values, rho, rng)
        passed_ii = (not np.isnan(p1) and rho > 0 and p1 < 0.05 and np_ < 0.05)
        if passed_ii:
            ii_pass_celltypes += 1
        add(ct, "(ii)spearman", rho, p1, np_, n, passed_ii, "shared secretome genes")
        # (iii)
        co = concordance(h2, nh, rng)
        passed_iii = (not np.isnan(co["binom_p"]) and co["binom_p"] < 0.05 and co["null_p"] < 0.05)
        add(ct, "(iii)concordance", co["hit_concord"], co["binom_p"], co["null_p"], co["n_hit"],
            passed_iii, f"bg={co['bg_concord']:.3f} fisherOR={co['fisher_or']:.2f} fisher_p={co['fisher_p']:.3g}")

    # ---- pooled across niche celltypes ----
    hp = hits.merge(ext, on=["gene", "celltype"], how="inner")
    fp = full.merge(ext, on=["gene", "celltype"], how="inner")
    fp["concord"] = np.sign(fp.log2FoldChange) == np.sign(fp.lfc_ext)
    hitkeys = set(zip(hits.gene, hits.celltype))
    fp["is_hit"] = [(g, c) in hitkeys for g, c in zip(fp.gene, fp.celltype)]
    di = directional(hp.assign(direction=hp.direction), rng)
    pooled_i = (not np.isnan(di["mwu_p"]) and di["mwu_p"] < 0.05 and di["null_p"] < 0.05)
    add("POOLED", "(i)directional", di.get("mwu_U", np.nan), di["mwu_p"], di["null_p"], di["n"], pooled_i,
        f"mean_signed={di['mean_signed']:.3f} wilcox_p={di['wilcox_p']:.3g} nUP={di['n_up']} nDOWN={di['n_down']}")
    rho, p1, n = onesided_spearman_gt0(fp.log2FoldChange.values, fp.lfc_ext.values)
    np_ = perm_spearman_null(fp.log2FoldChange.values, fp.lfc_ext.values, rho, rng)
    pooled_ii = (not np.isnan(p1) and rho > 0 and p1 < 0.05 and np_ < 0.05)
    add("POOLED", "(ii)spearman", rho, p1, np_, n, pooled_ii, "all secretome×niche")
    co = concordance(fp[fp.is_hit], fp[~fp.is_hit], rng)
    pooled_iii = (not np.isnan(co["binom_p"]) and co["binom_p"] < 0.05 and co["null_p"] < 0.05)
    add("POOLED", "(iii)concordance", co["hit_concord"], co["binom_p"], co["null_p"], co["n_hit"], pooled_iii,
        f"bg={co['bg_concord']:.3f} fisherOR={co['fisher_or']:.2f} fisher_p={co['fisher_p']:.3g}")

    res = pd.DataFrame(rows).rename(columns={"pass_": "pass"})
    out = f"{P}/results/validation/replication_{a.cohort}.csv"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    res.to_csv(out, index=False)

    # ---- verdict ----
    npass = sum([pooled_i, ii_pass_celltypes >= 2, pooled_iii])
    GO = npass >= 2
    print(f"\n===== A1 replication verdict: {a.cohort} =====")
    print(res.to_string(index=False))
    print(f"\n(i) pooled PASS={pooled_i} | (ii) celltypes passing={ii_pass_celltypes} (need>=2 ->{ii_pass_celltypes>=2}) | (iii) pooled PASS={pooled_iii}")
    print(f">>> {npass}/3 primary metrics PASS -> OVERALL {'GO: transcriptomic reproduction in a second cohort' if GO else 'NO-GO: reproduction NOT established (report honest null)'}")
    print(f"wrote {out}")

    # ---- scatter (publication-finalised; shared style in scripts/07_figures/_pubstyle.py) ----
    try:
        import sys
        os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_cache")
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        sys.path.insert(0, os.path.join(P, "scripts", "07_figures"))
        import _pubstyle as PS
        PS.apply()
        cts = [c for c in NICHE if c in scatter] + ["POOLED"]
        scatter["POOLED"] = fp
        ncol = 3; nrow = int(np.ceil(len(cts) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 3.3 * nrow), squeeze=False,
                                 layout="constrained")
        axes = axes.ravel()
        for k, ct in enumerate(cts):
            ax = axes[k]; d = scatter[ct]
            ax.axhline(0, c=PS.GRIDC, lw=0.8); ax.axvline(0, c=PS.GRIDC, lw=0.8)
            ax.scatter(d.log2FoldChange, d.lfc_ext, s=7, c=PS.NEUTRAL, alpha=0.5, lw=0,
                       rasterized=True, label="secretome")
            hh = d[d.is_hit]
            ax.scatter(hh.log2FoldChange, hh.lfc_ext, s=24, c=PS.SASP, lw=0, zorder=3, label="frozen hit")
            r, p1, nn = onesided_spearman_gt0(d.log2FoldChange.values, d.lfc_ext.values)
            ax.set_title(f"{ct}   ρ = {r:.2f}  ({PS.fmt_p(p1)}, n={nn})",
                         color=PS.CT.get(ct, PS.INK), fontsize=9)
            if k % ncol == 0:
                ax.set_ylabel(f"{a.cohort}\nlog$_2$FC")
            if k >= len(cts) - ncol:
                ax.set_xlabel("primary log$_2$FC (GSE268609)")
            PS.style_ax(ax)
        for k in range(len(cts), nrow * ncol):
            axes[k].axis("off")
        axes[0].legend(fontsize=7, loc="upper left")
        # B-4 (commit 4e3d27b): the claim label is "reproduction", not "replication" -- the
        # matched-background re-test shows the agreement is not demonstrably specific to the 60 hits.
        # The committed CSV filenames (replication_*.csv) keep the pre-registered metric name.
        fig.suptitle(f"Figure 3  ·  External transcriptomic reproduction of the frozen niche-secretome "
                     f"signature — {a.cohort}", fontsize=10.5, fontweight="bold")
        figout = f"{P}/figures/figure3/replication_{a.cohort}.pdf"
        os.makedirs(os.path.dirname(figout), exist_ok=True)
        fig.savefig(figout); plt.close(fig)
        print(f"wrote {figout}")
    except Exception as e:
        print(f"(scatter skipped: {type(e).__name__} {e})")


if __name__ == "__main__":
    main()
