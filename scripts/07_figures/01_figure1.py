#!/usr/bin/env python
"""Figure 1 — GSE268609 atlas & normal-aging composition. READ-ONLY on the anchor.

A: UMAP × celltype_l1 | B: UMAP × Group | C: marker dot-plot (native, styled; colour = ABSOLUTE mean CP10K, not var-scaled -- B-1)
D: normal-aging (YA vs HA) donor-level composition + Mann-Whitney, BH-adjusted padj
   (the panel annotates the SAME adjusted p the text leads with: NSC padj=0.001, Micro padj=0.002)
S1: per-cell-type QC violins (nCount / nFeature / pct_mt, from counts)
S2: composition across the five deposited groups. SA = SuperAgers (aged >=80 with episodic memory at
    or above the norms for people in their 50s-60s), NOT "severe AD" -- so the groups are ordered
    cognitively-unimpaired first (YA, HA, SA) then impaired (MCI, AD), and the axis is not a severity
    ordering. Corrected 2026-08-26, audit finding F-8-001.

Standalone filenames are unchanged (umap_A_celltype.pdf, umap_B_group.pdf, dotplot__C_markers.pdf,
panel_D_composition_YAvsHA.pdf, supp_S1_QC.pdf, supp_S2_composition_allgroups.pdf); composite
figure1.pdf is added. Does NOT modify the anchor (audit fix 2026-06-03); requires the symbol +
UMAP prep steps upstream. Publication style in _pubstyle.py; cell-type colours match the rest of
the figure set. seed 42 (panel-D jitter).
"""
import os, sys
import scanpy as sc, pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, scipy.stats as sst
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _pubstyle as PS

PS.apply()
P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
FIG = f"{P}/figures/figure1"; os.makedirs(FIG, exist_ok=True)
os.makedirs(f"{P}/results/figure1", exist_ok=True)
rng = np.random.default_rng(42)

A = sc.read_h5ad(f"{P}/processed/per_dataset/GSE268609_anchor.h5ad")
if str(A.var_names[0]).startswith("ENSG"):
    sys.exit("anchor var_names are Ensembl IDs — run scripts/02_qc/05_build_268609_anchor.py first")
if "X_umap" not in A.obsm:
    sys.exit("anchor has no UMAP embedding — run scripts/02_qc/05b_anchor_umap.py first")

counts = A.layers["counts"]
A.X = counts.copy()
sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A); A.raw = A

CATS = list(A.obs["celltype_l1"].cat.categories) if hasattr(A.obs["celltype_l1"], "cat") \
    else sorted(A.obs["celltype_l1"].unique())
# cell-type palette: niche senders use the shared CT colours; remaining (neuronal/other) muted
_EXTRA = ["#b9c0cc", "#8c97a8", "#aab4c2", "#c9bfa6", "#9fb0c4", "#d4d9e1", "#7f8a9b"]
PAL = {}; ei = 0
for c in CATS:
    if c in PS.CT:
        PAL[c] = PS.CT[c]
    else:
        PAL[c] = _EXTRA[ei % len(_EXTRA)]; ei += 1

MARKERS = {
    "Astro": ["AQP4", "GFAP", "SLC1A2"], "Micro": ["CSF1R", "C1QB", "P2RY12"],
    "Oligo": ["PLP1", "MOG", "MOBP"], "OPC": ["PDGFRA", "CSPG4"],
    "Endo": ["CLDN5", "FLT1"], "Ependymal": ["FOXJ1", "PIFO"],
    "DG_GC": ["PROX1"], "CA_ExN": ["SLC17A7", "GRIK4"], "InN": ["GAD1", "GAD2"],
    "NSC": ["HOPX", "VIM"], "Neuroblast": ["DCX", "EOMES"], "Immature": ["CALB2", "STMN1"],
}
MARKERS = {k: [g for g in v if g in A.var_names] for k, v in MARKERS.items()}

# ---- normal-aging composition (computed once; used by panel D) ----------------
comp = A.obs.groupby(["donor_id", "celltype_l1"], observed=True).size().unstack(fill_value=0)
comp = comp.div(comp.sum(1), axis=0)
comp["Group"] = A.obs.groupby("donor_id", observed=True)["Group"].first()
norm = comp[comp.Group.isin(["YA", "HA"])]
cts = [c for c in comp.columns if c != "Group"]
res = []
for c in cts:
    ya = norm[norm.Group == "YA"][c]; ha = norm[norm.Group == "HA"][c]
    try:
        _, p = sst.mannwhitneyu(ha, ya)
    except ValueError:
        p = np.nan
    res.append({"celltype": c, "YA_mean": ya.mean(), "HA_mean": ha.mean(),
                "log2FC_HA_YA": np.log2((ha.mean() + 1e-4) / (ya.mean() + 1e-4)), "p": p})
rd = pd.DataFrame(res); rd["padj"] = sst.false_discovery_control(rd.p.fillna(1))
rd = rd.sort_values("p")
rd.to_csv(f"{P}/results/figure1/composition_YAvsHA.csv", index=False)
PADJ = dict(zip(rd.celltype, rd.padj))


# ---- panels ------------------------------------------------------------------
def _declutter_on_data_labels(ax, min_sep_frac=0.10, max_iter=800):
    """Nudge scanpy's on-data cluster labels apart so they stop overlapping.

    `legend_loc="on data"` drops each label at its cluster's median embedding position. Where
    two clusters sit on top of each other in the UMAP -- Immature inside DG_GC, Neuroblast
    inside Oligo -- the two labels land within a few points of each other and overprint, which
    is what a reader sees as a single unreadable smear (reported 2026-08-26).

    This does a small repulsion pass in axes coordinates: any two labels closer than
    `min_sep_frac` of the axes diagonal push each other apart along the line joining them, and
    the labels are then clamped back inside the axes. Only the LABEL positions move; no data
    point, colour or cluster assignment is touched.
    """
    texts = [t for t in ax.texts if t.get_text().strip()]
    if len(texts) < 2:
        return
    ax.figure.canvas.draw()                      # positions are only valid after a draw
    inv = ax.transData.inverted()
    pos = np.array([ax.transData.transform(t.get_position()) for t in texts], float)
    bb = ax.get_window_extent()
    min_sep = min_sep_frac * float(np.hypot(bb.width, bb.height))
    for _ in range(max_iter):
        moved = False
        for i in range(len(pos)):
            for j in range(i + 1, len(pos)):
                d = pos[j] - pos[i]
                dist = float(np.hypot(*d))
                if dist >= min_sep:
                    continue
                if dist < 1e-6:                  # exactly coincident: pick an arbitrary axis
                    d = np.array([1.0, 0.0]); dist = 1.0
                push = (min_sep - dist) / 2.0 * (d / dist)
                pos[i] -= push; pos[j] += push
                moved = True
        if not moved:
            break
    for t, p in zip(texts, pos):
        p = np.array([min(max(p[0], bb.x0 + 2), bb.x1 - 2),
                      min(max(p[1], bb.y0 + 2), bb.y1 - 2)])
        t.set_position(tuple(inv.transform(p)))


def panel_umap_celltype(cont):
    ax = cont.subplots()
    sc.pl.umap(A, color="celltype_l1", palette=[PAL[c] for c in CATS], ax=ax, show=False,
               legend_loc="on data", legend_fontsize=6, legend_fontoutline=2, frameon=False, title="")
    ax.set_title(f"Cell-type atlas (n = {A.n_obs:,} nuclei)", fontsize=10)
    _declutter_on_data_labels(ax)
    return ax


def panel_umap_group(cont):
    ax = cont.subplots()
    sc.pl.umap(A, color="Group", ax=ax, show=False, frameon=False, title="")
    ax.set_title("Diagnosis group", fontsize=10)
    return ax


def panel_dotplot(cont):
    """Native dot-plot: dot size = % cells expressing, colour = ABSOLUTE mean CP10K.

    Audit fix 2026-08-26 (finding F-4-001 / BLOCKER B-1). This panel used to var-scale each
    gene's mean to [0,1] across cell types (scanpy standard_scale='var'). That rescaling makes
    every gene look saturated wherever it is highest and pale everywhere else, which hides
    lineage-inappropriate expression: the NSC-labelled cluster carries PLP1 at 52% of mature
    oligodendrocytes' level and SLC1A2 at 51% of astrocytes', and var-scaling rendered both
    faint. Colour is now absolute mean CP10K on a log axis (ticks in real CP10K units) so a dot
    is comparable across cell types AND across genes, and the exclusion markers PLP1, MOBP and
    AQP4 sit in the same panel where that contamination is readable.
    """
    genes, groups, bounds = [], [], []
    pos = 0
    for grp, gg in MARKERS.items():
        if not gg:
            continue
        bounds.append((pos, pos + len(gg), grp)); genes += gg; pos += len(gg)
    rows = CATS
    sub = A[:, genes]
    Xln = np.asarray(sub.X.todense()) if hasattr(sub.X, "todense") else np.asarray(sub.X)
    Xct = np.asarray(sub.layers["counts"].todense()) if hasattr(sub.layers["counts"], "todense") \
        else np.asarray(sub.layers["counts"])
    ct_vec = A.obs["celltype_l1"].values
    frac = np.zeros((len(rows), len(genes))); mean = np.zeros((len(rows), len(genes)))
    for i, ctp in enumerate(rows):
        m = ct_vec == ctp
        frac[i] = (Xct[m] > 0).mean(0) * 100.0
        mean[i] = Xln[m].mean(0)
    # absolute mean CP10K per (cell type, gene), pseudobulk-style: summed counts / summed UMI.
    # NOT var-scaled -- see the docstring and finding F-4-001.
    tot_umi = np.asarray(A.layers["counts"].sum(1)).ravel()
    cp10k = np.zeros((len(rows), len(genes)))
    for i, ctp in enumerate(rows):
        m = ct_vec == ctp
        cp10k[i] = Xct[m].sum(0) / max(tot_umi[m].sum(), 1.0) * 1e4
    col = np.log10(cp10k + 1.0)                       # CP10K spans 0-42 with median 0.15
    vmax = float(np.log10(cp10k.max() + 1.0))
    ax = cont.subplots()
    xs, ys = np.meshgrid(np.arange(len(genes)), np.arange(len(rows)))
    sca = ax.scatter(xs.ravel(), ys.ravel(), s=(frac.ravel() / 100.0) * 90 + 1,
                     c=col.ravel(), cmap="Reds", vmin=0, vmax=vmax, edgecolors="none")
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=90, fontsize=7)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=8)
    for t, ctp in zip(ax.get_yticklabels(), rows):
        if ctp in PS.CT:
            t.set_color(PS.CT[ctp]); t.set_fontweight("bold")
    ax.set_xlim(-0.7, len(genes) - 0.3); ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.invert_yaxis(); ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    for (a, b, grp) in bounds:                                   # marker-group brackets
        ax.plot([a - 0.3, b - 0.7], [len(rows) - 0.4, len(rows) - 0.4], color=PS.SUB, lw=0.8,
                clip_on=False)
        ax.text((a + b - 1) / 2, len(rows) - 0.65, grp, ha="center", va="bottom", fontsize=7,
                color=PS.SUB)
    cb = cont.colorbar(sca, ax=ax, fraction=0.018, pad=0.01)
    ticks = [t for t in [0, 1, 3, 10, 30] if np.log10(t + 1.0) <= vmax + 1e-9]
    cb.set_ticks([np.log10(t + 1.0) for t in ticks])
    cb.set_ticklabels([str(t) for t in ticks])
    cb.set_label("mean CP10K\n(absolute, log axis)", fontsize=7); cb.ax.tick_params(labelsize=6)
    cb.outline.set_visible(False)
    handles = [Line2D([0], [0], marker="o", ls="", mfc=PS.SUB, mec="none",
                      ms=np.sqrt((f / 100) * 90 + 1), label=f"{f}%") for f in [20, 50, 100]]
    ax.legend(handles=handles, title="% cells", fontsize=6.5, title_fontsize=7,
              loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=3, handletextpad=0.1, columnspacing=0.8)
    return ax


def panel_composition(cont, plot_ct):
    nr = int(np.ceil(len(plot_ct) / 3))
    axes = cont.subplots(nr, 3)
    axes = np.atleast_1d(axes).ravel()
    for ax, c in zip(axes, plot_ct):
        col = PS.CT.get(c, PS.SUB)
        for i, g in enumerate(["YA", "HA"]):
            v = norm[norm.Group == g][c].values
            gc = PS.YA if g == "YA" else PS.HA
            ax.bar(i, v.mean(), 0.55, color=gc, alpha=0.28, lw=0)
            ax.scatter(np.full(len(v), i) + rng.uniform(-0.08, 0.08, len(v)), v, s=18, color=gc,
                       alpha=0.85, lw=0)
        padj = PADJ.get(c, np.nan)
        sig = "  *" if padj < 0.1 else ""
        ax.set_xticks([0, 1]); ax.set_xticklabels(["YA", "HA"], fontsize=8)
        ax.get_xticklabels()[0].set_color(PS.YA); ax.get_xticklabels()[1].set_color(PS.HA)
        ax.set_title(f"{c}{sig}\nadj. {PS.fmt_p(padj)}", fontsize=8.5,
                     color=col, fontweight="bold" if padj < 0.1 else "normal")
        ax.set_ylabel("fraction", fontsize=8); ax.set_ylim(bottom=0)
        PS.style_ax(ax, grid=True)
    for ax in axes[len(plot_ct):]:
        ax.axis("off")
    return axes


# ---- standalone panels (legend-referenced filenames) -------------------------
f = plt.figure(figsize=(5.2, 5.0), layout="constrained"); panel_umap_celltype(f)
PS.save(f, f"{FIG}/umap_A_celltype.pdf")
f = plt.figure(figsize=(5.6, 5.0), layout="constrained"); panel_umap_group(f)
PS.save(f, f"{FIG}/umap_B_group.pdf")
f = plt.figure(figsize=(9.5, 4.2), layout="constrained"); panel_dotplot(f)
PS.save(f, f"{FIG}/dotplot__C_markers.pdf")
PLOT_CT_FULL = ["NSC", "Astro", "Micro", "Oligo", "OPC", "Endo", "DG_GC", "CA_ExN", "InN"]
f = plt.figure(figsize=(9.5, 9.0), layout="constrained"); panel_composition(f, PLOT_CT_FULL)
PS.save(f, f"{FIG}/panel_D_composition_YAvsHA.pdf")

# ---- composite figure1.pdf ---------------------------------------------------
fig = plt.figure(figsize=(7.4, 9.4), layout="constrained")
fig.suptitle("Figure 1  ·  Human hippocampal multiome resolves the normal-aging niche",
             fontsize=11, fontweight="bold", ha="left", x=0.012)
rows = fig.subfigures(3, 1, height_ratios=[1.0, 0.78, 1.0])
toprow = rows[0].subfigures(1, 2, width_ratios=[1.0, 1.08])
panel_umap_celltype(toprow[0]); PS.letter(toprow[0], "A")
panel_umap_group(toprow[1]); PS.letter(toprow[1], "B")
panel_dotplot(rows[1]); PS.letter(rows[1], "C")
panel_composition(rows[2], ["NSC", "Micro", "Astro", "Oligo", "OPC", "Endo"]); PS.letter(rows[2], "D")
PS.save(fig, f"{FIG}/figure1.pdf")

# ---------- Supplementary S1: per-cell-type QC violins ----------
nCount = np.asarray(counts.sum(1)).ravel()
nFeature = np.asarray((counts > 0).sum(1)).ravel()
mt_idx = [A.var_names.get_loc(g) for g in A.var_names if str(g).startswith("MT-")]
mt_sum = np.asarray(counts[:, mt_idx].sum(1)).ravel() if mt_idx else np.zeros_like(nCount)
pct_mt = 100.0 * mt_sum / np.maximum(nCount, 1)
qc = pd.DataFrame({"celltype_l1": A.obs["celltype_l1"].values,
                   "nCount_RNA": nCount, "nFeature_RNA": nFeature, "pct_mt": pct_mt})
fig, axes = plt.subplots(3, 1, figsize=(11, 8.5), layout="constrained")
for ax, (col, lab) in zip(axes, [("nCount_RNA", "UMIs / nucleus"),
                                 ("nFeature_RNA", "genes / nucleus"), ("pct_mt", "% mitochondrial")]):
    data = [qc[qc.celltype_l1 == ct][col].values for ct in CATS]
    parts = ax.violinplot(data, showmedians=True, widths=0.85)
    for k, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(PAL[CATS[k]]); pc.set_alpha(0.6); pc.set_edgecolor("none")
    for key in ("cmedians", "cbars", "cmins", "cmaxes"):
        if key in parts:
            parts[key].set_color(PS.SUB); parts[key].set_linewidth(0.8)
    ax.set_xticks(range(1, len(CATS) + 1)); ax.set_xticklabels(CATS, rotation=45, ha="right")
    ax.set_ylabel(lab)
    if col != "pct_mt":
        ax.set_yscale("log")
    PS.style_ax(ax)
fig.suptitle(f"Figure S1  ·  Per-cell-type QC (GSE268609, n = {A.n_obs:,} nuclei)",
             fontsize=11, fontweight="bold")
PS.save(fig, f"{FIG}/supp_S1_QC.pdf")

# ---------- Supplementary S2: composition across all groups ----------
# SA = SuperAgers, a cognitive-resilience phenotype, NOT severe AD (F-8-001). Unimpaired groups
# (YA, HA, SA) first, then impaired (MCI, AD); SA is deliberately not adjacent to AD so that the
# panel order cannot be read as a severity axis.
GROUPS = ["YA", "HA", "SA", "MCI", "AD"]
allcomp = A.obs.groupby(["donor_id", "celltype_l1"], observed=True).size().unstack(fill_value=0)
allcomp = allcomp.div(allcomp.sum(1), axis=0)
allcomp["Group"] = A.obs.groupby("donor_id", observed=True)["Group"].first()
allcomp = allcomp[allcomp.Group.isin(GROUPS)]
cts_all = [c for c in allcomp.columns if c != "Group"]
gmean = allcomp.groupby("Group", observed=True)[cts_all].mean().reindex(GROUPS)
gsem = allcomp.groupby("Group", observed=True)[cts_all].sem().reindex(GROUPS)
gmean.to_csv(f"{P}/results/figure1/composition_allgroups.csv")
show_ct = [c for c in ["NSC", "Neuroblast", "Immature", "Micro", "Astro", "Oligo", "OPC", "Endo",
                       "InN", "DG_GC"] if c in cts_all]
ncol = 5; nrow = int(np.ceil(len(show_ct) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 2.9 * nrow), layout="constrained")
axes = np.atleast_1d(axes).ravel()
x = np.arange(len(GROUPS))
for ax, c in zip(axes, show_ct):
    col = PS.CT.get(c, PS.SUB)
    ax.errorbar(x, gmean[c].values, yerr=gsem[c].values, marker="o", capsize=3, color=col, lw=1.4)
    ax.set_xticks(x); ax.set_xticklabels(GROUPS, fontsize=8)
    ax.set_title(c, color=col, fontweight="bold", fontsize=9.5)
    ax.set_ylabel("fraction", fontsize=8); ax.set_ylim(bottom=0)
    PS.style_ax(ax, grid=True)
for ax in axes[len(show_ct):]:
    ax.axis("off")
fig.suptitle("Figure S2  ·  Composition across the five deposited groups "
             "(YA, HA, SA = SuperAgers, MCI, AD — not a severity ordering): "
             # Audit fix 2026-08-26. The old subtitle asserted two things the manuscript has since
             # withdrawn or downgraded: "neuroblast" is a depositor annotation this re-analysis does
             # not support (B-1 / F-4-001), and the microglial change is method-dependent -- it is
             # significant under the raw-proportion test plotted here but not under a centred
             # log-ratio (adjusted p = 0.48; P1-16 / F-5-003). The panel now describes what it plots.
             "relative composition by group; see Results for the compositional caveat",
             fontsize=10.5, fontweight="bold")
PS.save(fig, f"{FIG}/supp_S2_composition_allgroups.pdf")

print("\n=== Fig1D normal-aging composition (YA vs HA), BH-adjusted ===")
print(rd.round(4).to_string(index=False))
print("\nsaved panels A-D + composite figure1.pdf + S1 + S2; results/figure1/{composition_YAvsHA,composition_allgroups}.csv")
