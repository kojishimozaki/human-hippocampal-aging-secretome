#!/usr/bin/env python
"""Figure 5 (Visium spatial) — the assembled composite figure5.pdf + its standalone panel.

Left/mid: representative spatial maps (PROX1 = DG marker; Astro aging-UP signature), best-DG area.
Right: per-capture-area spatial Spearman(signature, PROX1) across ALL QC-passing areas → centred ~0
       ⇒ aging secretome signatures show NO DETECTABLE DG co-localization. A null association on sections
       that do not resolve the DG establishes neither DG exclusion nor uniformity. Avoids single-section artefact.

Uses ALL 34 QC-passing capture areas (no arbitrary truncation; "N of M used" logged). As of 2026-06-29
the main Figure 5 carries ONLY this Visium spatial-scope result: the former Fig 5A cross-cohort
boundary panel (GSE199243, signed reproducibility 0.02) was moved to the external-cohort replication-spectrum
supplement Fig S(E1) (scripts/07_figures/15_replication_spectrum_supp.py), so this script no longer assembles
a cross-cohort panel and the standalone panel is panel_A_visium_spatial.pdf. Computation unchanged from the
audited pipeline; only the rendering is publication-finalised (shared style in _pubstyle.py). Data-dependent
(needs the Visium matrices in processed/). seed 42.
"""
import os, sys
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")   # scanpy/numba + mpl reproducibility
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_cache")        # on hosts with a read-only/cacheless HOME
import scanpy as sc, pandas as pd, numpy as np, glob, scipy.io, scipy.stats as sst
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _pubstyle as PS

PS.apply()
np.random.seed(42)                                            # reproducible boxplot jitter
P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
D = f"{P}/processed/per_dataset/GSE264692_visium"; FIG = f"{P}/figures/figure5"
up = pd.read_csv(f"{P}/results/de/de_GSE268609_per_celltype.csv")


def sigs(A):
    out = {}
    for ct in ["Astro", "Micro", "Oligo"]:
        g = [x for x in up[(up.celltype == ct) & (up.padj < 0.2) & (up.log2FoldChange > 0)].gene
             if x in A.var_names]
        if len(g) >= 5:
            sc.tl.score_genes(A, g, score_name=f"{ct}_up"); out[f"{ct}_up"] = A.obs[f"{ct}_up"].values
    return out


def load_area(mx):
    pre = mx.replace("_matrix.mtx.gz", "")
    X = scipy.io.mmread(mx).tocsr().T
    gg = pd.read_csv(f"{pre}_features.tsv.gz", sep="\t", header=None)
    bc = pd.read_csv(f"{pre}_barcodes.tsv.gz", header=None)[0].tolist()
    a = sc.AnnData(X); a.var_names = gg[1].astype(str).values; a.obs_names = bc; a.var_names_make_unique()
    pos = pd.read_csv(f"{pre}-tissue_positions_list.csv.gz", header=None).set_index(0)
    a = a[a.obs_names.isin(pos.index)]; pp = pos.loc[a.obs_names]
    a = a[pp[1].values == 1].copy(); pp = pp[pp[1] == 1]
    a.obs["x"] = pp[5].values.astype(float); a.obs["y"] = -pp[4].values.astype(float)
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    return a


# ---- compute over ALL capture areas (QC inside the loop; number used is reported) ----
areas = sorted(glob.glob(f"{D}/*_matrix.mtx.gz"))
print(f"capture areas found: {len(areas)} (matrix files); QC-filtering for PROX1 + >=500 spots")
percorr = {"Astro_up": [], "Micro_up": [], "Oligo_up": []}
corr_rows, n_skip_qc, best = [], 0, None
for mx in areas:
    try:
        a = load_area(mx)
        if "PROX1" not in a.var_names or a.n_obs < 500:
            n_skip_qc += 1
            continue
        pr = np.asarray(a[:, "PROX1"].X.todense()).ravel()
        ss = sigs(a)
        area_id = os.path.basename(mx).replace("_matrix.mtx.gz", "")
        for k, v in ss.items():
            rho = sst.spearmanr(v, pr)[0]
            percorr[k].append(rho)
            corr_rows.append(dict(capture_area=area_id, signature=k, n_spots=int(a.n_obs),
                                  spearman_rho_vs_PROX1=rho))
        hi = pr >= np.quantile(pr, 0.95)
        conc = 1 - (a.obs[["x", "y"]].values[hi].std(0).mean() / (a.obs[["x", "y"]].values.std(0).mean() + 1e-9))
        if best is None or conc * np.sqrt(hi.sum()) > best[1]:
            best = (a, conc * np.sqrt(hi.sum()), pr, hi)
    except Exception as e:
        print("skip", os.path.basename(mx), str(e)[:50])
A, _, pr, hi = best
n_used = len(percorr["Astro_up"])
print(f"capture areas: {len(areas)} found, {n_used} used (PROX1 + >=500 spots); "
      f"{n_skip_qc} failed QC, {len(areas) - n_used - n_skip_qc} unreadable/no-position-list")
for k, v in percorr.items():
    print(f"  {k}: per-area Spearman(PROX1) median={np.median(v):+.3f} "
          f"(range {min(v):+.2f}..{max(v):+.2f}, n={len(v)})")

os.makedirs(f"{P}/results/figure5", exist_ok=True)
pd.DataFrame(corr_rows).to_csv(f"{P}/results/figure5/visium_spatial_correlation.csv", index=False)
_summ = [dict(signature=k, n_areas=len(v), median_rho=float(np.median(v)),
              median_abs_rho=float(np.median(np.abs(v))),
              min_rho=float(np.min(v)), max_rho=float(np.max(v))) for k, v in percorr.items()]
pd.DataFrame(_summ).to_csv(f"{P}/results/figure5/visium_spatial_correlation_summary.csv", index=False)
print("wrote results/figure5/visium_spatial_correlation{,_summary}.csv")

A.obs["Astro_up"] = sigs(A)["Astro_up"]


def panel_visium(cont):
    axes = cont.subplots(1, 3)
    for i, (col, ttl) in enumerate([("PROX1", "DG marker PROX1"),
                                    ("Astro_up", "Astrocyte aging-UP secretome")]):
        ax = axes[i]
        val = np.asarray(A[:, col].X.todense()).ravel() if col in A.var_names else A.obs[col].values
        sca = ax.scatter(A.obs.x, A.obs.y, c=val, s=12, cmap="viridis",
                         vmin=np.quantile(val, 0.02), vmax=np.quantile(val, 0.99), linewidths=0)
        ax.scatter(A.obs.x.values[hi], A.obs.y.values[hi], s=16, facecolors="none",
                   edgecolors="#e84a5f", lw=0.4, alpha=0.6)
        ax.set_title(ttl, fontsize=9); ax.set_aspect("equal"); ax.axis("off")
        cb = cont.colorbar(sca, ax=ax, fraction=0.045, pad=0.02); cb.ax.tick_params(labelsize=6)
        cb.outline.set_visible(False)
    ax = axes[2]
    labs = list(percorr.keys()); data = [percorr[k] for k in labs]
    cmap = {"Astro_up": PS.CT["Astro"], "Micro_up": PS.CT["Micro"], "Oligo_up": PS.CT["Oligo"]}
    bp = ax.boxplot(data, tick_labels=[l.replace("_up", " ↑") for l in labs], showfliers=False,
                    widths=0.6, medianprops=dict(color=PS.INK, lw=1.2), patch_artist=True)
    for patch, l in zip(bp["boxes"], labs):
        patch.set_facecolor(cmap[l]); patch.set_alpha(0.25); patch.set_edgecolor(cmap[l])
    for j, (v, l) in enumerate(zip(data, labs)):
        ax.scatter(np.full(len(v), j + 1) + np.random.uniform(-0.1, 0.1, len(v)), v,
                   s=14, color=cmap[l], alpha=0.65, lw=0)
    ax.axhline(0, c=PS.SUB, lw=0.8, ls=(0, (4, 3))); ax.set_ylim(-0.4, 0.4)
    ax.set_ylabel("spatial Spearman vs PROX1 (per capture area)", fontsize=8)
    ax.set_title(f"No detectable DG co-localization across {len(data[0])} areas\n(does not establish uniformity)", fontsize=9)
    PS.style_ax(ax, grid=True)
    return axes


# ---- standalone panel A (Visium scope; legend-referenced filename) ----
f = plt.figure(figsize=(14, 4.4), layout="constrained")
panel_visium(f)
f.suptitle("Visium: no detectable DG concentration of the aging secretome "
           "(red outline = DG / PROX1-high spots)", fontsize=10, fontweight="bold")
PS.save(f, f"{FIG}/panel_A_visium_spatial.pdf")

# ---- composite figure5.pdf (Visium-only scope figure; GSE199243 boundary moved to Fig S(E1)) ----
fig5 = plt.figure(figsize=(14, 4.7), layout="constrained")
fig5.suptitle("Figure 5  ·  The aging niche-secretome shows no detectable dentate-gyrus concentration",
              fontsize=10.5, fontweight="bold", ha="left", x=0.012)
panel_visium(fig5)
PS.save(fig5, f"{FIG}/figure5.pdf")
print("saved figures/figure5/panel_A_visium_spatial.pdf + composite figure5.pdf (Visium-only scope)")
