#!/usr/bin/env python
"""Figure 2 (PRIMARY) — niche-secretome aging signature (GSE268609 HA vs YA).

A: per-niche-celltype volcano (secretome highlighted; SASP in red)
B: secretome-category x celltype mean-log2FC heatmap
C: SASP set directional enrichment (SenMayo + SASP Atlas) per celltype

Each panel is written standalone under its legend-referenced filename
(panel_A_volcano.pdf / panel_B_category_heatmap.pdf / panel_C_sasp_enrichment.pdf)
AND assembled into a single composite figure2.pdf with panel letters. The data /
statistics are unchanged from the audited pipeline (60 secretome hits at padj<0.1:
Astro 26, Micro 12, Endo 11, OPC 6, Oligo 5; SASP %up + directional p from
results/validation/sasp_overlap.csv); only the rendering is publication-finalised.
Shared look-and-feel lives in scripts/07_figures/_pubstyle.py.
"""
import os, sys
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _pubstyle as PS

PS.apply()
P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
FIG = f"{P}/figures/figure2"; os.makedirs(FIG, exist_ok=True)
NICHE = PS.NICHE

# ---- data (audited; unchanged) ------------------------------------------------
de = pd.read_csv(f"{P}/results/de/de_GSE268609_per_celltype.csv")
sec = pd.read_csv(f"{P}/refs/secretome_union.csv")
secset = set(sec.gene); cat = dict(zip(sec.gene, sec.category))
sasp = set(open(f"{P}/refs/senmayo.txt").read().split()) | \
       set(open(f"{P}/refs/sasp_atlas_core.txt").read().split())
ov = pd.read_csv(f"{P}/results/validation/sasp_overlap.csv")

CATLAB = {"cytokine_chemokine": "Cytokine / chemokine",
          "ecm_matrisome_like": "ECM (matrisome-like)",
          "growth_factor": "Growth factor", "neuropeptide": "Neuropeptide",
          "other_secreted": "Other secreted"}


# ---- helpers -----------------------------------------------------------------
def _place_labels(ax, df, n=6, fs=6.5):
    """Annotate the top-n secretome hits with a light vertical de-overlap + leader."""
    sub = df.nlargest(n, "nlp")
    if sub.empty:
        return
    ymin, ymax = ax.get_ylim()
    sep = (ymax - ymin) * 0.075
    rows = sub.sort_values("nlp", ascending=False)[["log2FoldChange", "nlp", "gene"]].values.tolist()
    ty_prev = None
    for r in rows:
        ty = r[1]
        if ty_prev is not None and ty > ty_prev - sep:
            ty = ty_prev - sep
        ty_prev = ty
        dx = 0.3 if r[0] >= 0 else -0.3
        ha = "left" if r[0] >= 0 else "right"
        ax.annotate(r[2], xy=(r[0], r[1]), xytext=(r[0] + dx, ty), ha=ha, va="center",
                    fontsize=fs, color=PS.INK,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color=PS.SUB, shrinkA=0, shrinkB=1))


# ---- panels (each takes a container = Figure or SubFigure) -------------------
def panel_volcano(cont):
    axes = cont.subplots(1, 5)
    for ax, ct in zip(axes, NICHE):
        d = de[de.celltype == ct].copy()
        d["nlp"] = -np.log10(d.padj.clip(lower=1e-300))
        d["issec"] = d.gene.isin(secset); d["issa"] = d.gene.isin(sasp)
        nhit = int((d.issec & (d.padj < 0.1)).sum())
        ax.scatter(d[~d.issec].log2FoldChange, d[~d.issec].nlp, s=3, c=PS.NEUTRAL,
                   alpha=0.35, lw=0, rasterized=True)
        ax.scatter(d[d.issec & ~d.issa].log2FoldChange, d[d.issec & ~d.issa].nlp, s=10,
                   c=PS.SECR, alpha=0.7, lw=0)
        ax.scatter(d[d.issa].log2FoldChange, d[d.issa].nlp, s=14, c=PS.SASP, alpha=0.9, lw=0)
        ax.axvline(0, c=PS.GRIDC, lw=0.8, zorder=0)
        ax.axhline(-np.log10(0.1), ls=(0, (4, 3)), c=PS.SUB, lw=0.7)
        ax.set_xlim(-5, 5)
        ax.set_ylim(0, max(2.5, float(d.nlp.max()) * 1.18))
        ax.set_title(f"{ct}  ({nhit})", color=PS.CT[ct], fontsize=10)
        ax.tick_params(length=3)
        if ct == "Astro":
            ax.set_ylabel("$-$log$_{10}$(adj. $p$)")
        else:
            ax.tick_params(labelleft=True)
        _place_labels(ax, d[d.issec & (d.padj < 0.1)])
    handles = [Line2D([0], [0], marker="o", ls="", mfc=PS.SECR, mec="none", ms=5, label="secretome"),
               Line2D([0], [0], marker="o", ls="", mfc=PS.SASP, mec="none", ms=5, label="SASP"),
               Line2D([0], [0], ls=(0, (4, 3)), c=PS.SUB, lw=0.8, label="adj. $p$ = 0.1")]
    axes[0].legend(handles=handles, loc="upper left", fontsize=6.5, handlelength=1.3)
    cont.supxlabel("log$_2$ fold-change (aging, HA $-$ YA)   ·   (n) = secretome hits at adj. $p$ < 0.1",
                   fontsize=8.5, y=0.0)
    return axes


def panel_heatmap(cont):
    d = de.copy(); d["cat"] = d.gene.map(cat)
    ds = d[d.gene.isin(secset) & d.padj.notna()]
    piv = (ds.groupby(["cat", "celltype"], observed=True).log2FoldChange.mean()
             .unstack().reindex(columns=NICHE))
    piv = piv.reindex([c for c in CATLAB if c in piv.index])
    ax = cont.subplots()
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-0.5, vmax=0.5, aspect="auto")
    ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns)
    for t, ct in zip(ax.get_xticklabels(), piv.columns):
        t.set_color(PS.CT[ct]); t.set_fontweight("bold")
    ax.set_yticks(range(len(piv.index))); ax.set_yticklabels([CATLAB[c] for c in piv.index])
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=7.5,
                        color="white" if abs(v) > 0.33 else PS.INK)
    cb = cont.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("mean log$_2$FC (HA $-$ YA)", fontsize=8); cb.ax.tick_params(labelsize=7)
    cb.outline.set_visible(False)
    ax.set_title("Secretome category × niche cell type")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    return ax, piv


def panel_sasp(cont):
    axes = cont.subplots(1, 2, sharey=True)
    titles = {"SenMayo": "SenMayo\n(transcriptomic SASP)",
              "SASP_Atlas": "SASP Atlas\n(mass-spec secretome)"}
    for ax, rs in zip(axes, ["SenMayo", "SASP_Atlas"]):
        sub = ov[ov.refset == rs].set_index("celltype").reindex(NICHE)
        y = -np.log10(sub.dir_p.values.astype(float))
        bars = ax.bar(range(len(NICHE)), y, width=0.66, color=[PS.CT[c] for c in NICHE])
        for b, p in zip(bars, sub.dir_p.values):
            b.set_alpha(0.95 if p < 0.05 else 0.4)
        for i, uf in enumerate(sub.up_frac.values):
            if not np.isnan(uf):
                ax.text(i, y[i] + 0.05, f"{uf*100:.0f}%↑", ha="center", va="bottom",
                        fontsize=7, color=PS.INK)
        ax.axhline(-np.log10(0.05), ls=(0, (4, 3)), c=PS.SUB, lw=0.7)
        ax.set_xticks(range(len(NICHE))); ax.set_xticklabels(NICHE, fontsize=8)
        for t, ct in zip(ax.get_xticklabels(), NICHE):
            t.set_color(PS.CT[ct])
        ax.set_title(titles[rs], fontsize=8.5)
        PS.style_ax(ax, grid=True)
    axes[0].set_ylabel("$-$log$_{10}$(directional $p$)")
    axes[1].text(len(NICHE) - 0.55, -np.log10(0.05), " p = 0.05", va="bottom", ha="right",
                 fontsize=6.5, color=PS.SUB)
    return axes


# ---- standalone panels (legend-referenced filenames, unchanged) --------------
f = plt.figure(figsize=(13, 2.9), layout="constrained"); panel_volcano(f)
PS.save(f, f"{FIG}/panel_A_volcano.pdf")

f = plt.figure(figsize=(5.0, 3.6), layout="constrained"); _, piv = panel_heatmap(f)
PS.save(f, f"{FIG}/panel_B_category_heatmap.pdf")

f = plt.figure(figsize=(6.8, 3.2), layout="constrained"); panel_sasp(f)
PS.save(f, f"{FIG}/panel_C_sasp_enrichment.pdf")

# ---- composite figure2.pdf ---------------------------------------------------
fig = plt.figure(figsize=(7.2, 6.6), layout="constrained")
fig.suptitle("Figure 2  ·  Aging reprograms the niche's secreted-protein output (primary result)",
             fontsize=11, fontweight="bold", ha="left", x=0.012)
sub = fig.subfigures(2, 1, height_ratios=[1.0, 1.22])
top = sub[0]; panel_volcano(top); PS.letter(top, "A")
bot = sub[1].subfigures(1, 2, width_ratios=[1.0, 1.18])
panel_heatmap(bot[0]); PS.letter(bot[0], "B")
panel_sasp(bot[1]); PS.letter(bot[1], "C")
PS.save(fig, f"{FIG}/figure2.pdf")

print("Fig2 secretome category × celltype (mean log2FC):\n", piv.round(2).to_string())
print("\nSASP %up / directional p (sasp_overlap.csv):")
print(ov[["refset", "celltype", "up_frac", "dir_p"]].to_string(index=False))
print("\nsaved panel_A_volcano.pdf, panel_B_category_heatmap.pdf, panel_C_sasp_enrichment.pdf, figure2.pdf")
