#!/usr/bin/env python
"""Figure 4 (SUPPORTING) — regulatory layer. The aging-secretome RNA change shows only a weak,
above-background direction-concordance with chromatin accessibility; the cis mechanism is unresolved.

A: secretome RNA log2FC vs linked-peak ATAC log2FC (direction-concordance; COL21A1 = lone ATAC-FDR hit)
B: chromVAR aging (HAvYA) nominal TF-motif activity, niche cells (ETS / Wnt-TCF / bHLH decline)
C: # FDR-sig TF motifs, aging (HAvYA) vs disease (ADvHA) — honest contrast
D: donor-level TF footprints (drawn by scripts/05_regulatory/05_footprint_figure.py)

Concise honesty stays on-figure ("supporting"; COL21A1-only; "exploratory" nominal points); the
long stat strings (binom/Fisher p, 0/4,971 & 0/369,393, matched-null enrichment) live in the legend
(figures/figure4/FIG4_legend_and_methods.md). Standalone panel filenames are unchanged; composite figure4.pdf
is added. Statistics are unchanged from the audited pipeline; only the rendering is finalised.
"""
import os, sys
import pandas as pd, numpy as np
from scipy.stats import binomtest, fisher_exact
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _pubstyle as PS

# Import draw_footprints() from the footprint script to embed panel D (donor-level TF
# footprints) in the composite. The module filename starts with a digit, so load via importlib.
import importlib.util as _ilu
_FP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "05_regulatory",
                        "05_footprint_figure.py")
_spec = _ilu.spec_from_file_location("footprint_figure", _FP_PATH)
_fpmod = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_fpmod)
draw_footprints = _fpmod.draw_footprints

NICHE = ["Astro", "Micro", "Oligo", "OPC", "Endo"]


# ---- panel A : RNA vs ATAC direction-concordance -----------------------------
def panel_concordance(cont, P):
    ac = pd.read_csv(f"{P}/results/de/atac_rna_concordance.csv")
    ac["conc"] = np.sign(ac.lfc_rna) == np.sign(ac.lfc_atac)
    rnasig = ac[ac.padj_rna < 0.1]; n = len(rnasig); k = int(rnasig.conc.sum())
    ns = ac[ac.padj_rna >= 0.1]; bg_r = ns.conc.mean()
    p_bg = binomtest(k, n, bg_r, alternative="greater").pvalue
    orr, p_fish = fisher_exact([[k, n - k], [int(ns.conc.sum()), len(ns) - int(ns.conc.sum())]],
                               alternative="greater")
    fdr = rnasig[rnasig.padj_atac_rnasig < 0.1]
    nom = ac[(ac.padj_rna < 0.1) & (ac.p_atac < 0.1)]
    ax = cont.subplots()
    ax.axhline(0, c=PS.GRIDC, lw=0.8); ax.axvline(0, c=PS.GRIDC, lw=0.8)
    ax.scatter(ns.lfc_rna, ns.lfc_atac, s=5, c=PS.NEUTRAL, alpha=0.3, lw=0, rasterized=True,
               label="RNA n.s.")
    ax.scatter(rnasig.lfc_rna, rnasig.lfc_atac, s=24, c=PS.SECR, alpha=0.7, lw=0,
               label=f"RNA-sig (n={n})")
    ax.scatter(nom.lfc_rna, nom.lfc_atac, s=46, facecolors="none", edgecolors=PS.SASP, lw=1.1,
               label=f"+nominal ATAC p<0.1 (n={len(nom)}, exploratory)")
    ax.scatter(fdr.lfc_rna, fdr.lfc_atac, s=95, marker="*", c=PS.INK, zorder=6,
               label=f"ATAC BH-FDR<0.1 within RNA-sig family (n={len(fdr)})")
    for _, r in fdr.iterrows():
        ax.annotate(r.gene, (r.lfc_rna, r.lfc_atac), xytext=(5, 3), textcoords="offset points",
                    fontsize=8, fontweight="bold")
    ax.set_xlabel("RNA log$_2$FC (HA $-$ YA)"); ax.set_ylabel("ATAC log$_2$FC (linked peaks, HA $-$ YA)")
    ax.set_title(f"{k}/{n} ({k/n:.0%}) concordant vs {bg_r*100:.1f}% background  ·  "
                 "COL21A1 is the only gene reaching ATAC FDR within the RNA-significant family; "
                   "0 of 369,393 genome-wide peak tests reach it",
                 fontsize=7.8, color=PS.SUB, fontweight="normal")
    cont.suptitle("Secretome RNA ↔ chromatin accessibility (HA vs YA) — supporting",
                  fontsize=10, fontweight="bold")
    ax.legend(fontsize=6.3, loc="lower right", handletextpad=0.4)
    PS.style_ax(ax)
    return ax, dict(n=n, k=k, bg=bg_r, p_bg=p_bg, orr=orr, p_fish=p_fish,
                    nom=len(nom), nom_conc=int(nom.conc.sum()), fdr=fdr.gene.tolist(), tot=len(ac))


# ---- panel B : chromVAR aging motif activity ---------------------------------
def panel_motifs(cont, P, n_tfs=20):
    cv = pd.read_csv(f"{P}/results/de/chromvar_motifs_strengthened.csv")
    aging = cv[cv.contrast == "HAvYA"].copy()
    top_tfs = aging.sort_values("p").drop_duplicates("tf").head(n_tfs).tf.tolist()
    piv = (aging[aging.tf.isin(top_tfs)].pivot_table(index="tf", columns="celltype", values="dz")
           .reindex(columns=NICHE))
    piv = piv.loc[piv.mean(1).sort_values().index]
    ax = cont.subplots()
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-0.3, vmax=0.3, aspect="auto")
    ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns, fontsize=8)
    for t, ct in zip(ax.get_xticklabels(), piv.columns):
        t.set_color(PS.CT[ct]); t.set_fontweight("bold")
    ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index, fontsize=7)
    cb = cont.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label("motif activity Δ (HA $-$ YA), nominal", fontsize=7.5); cb.ax.tick_params(labelsize=6.5)
    cb.outline.set_visible(False)
    ax.set_title("Aging TF-motif activity (chromVAR)\nETS / Wnt-TCF / bHLH decline", fontsize=9)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    return ax


# ---- panel C : aging vs disease motif counts ---------------------------------
def panel_aging_disease(cont, P):
    cv = pd.read_csv(f"{P}/results/de/chromvar_motifs_strengthened.csv")
    cnt = (cv[cv.padj < 0.1].groupby(["celltype", "contrast"], observed=True).size()
           .unstack(fill_value=0).reindex(NICHE).reindex(columns=["HAvYA", "ADvHA"]).fillna(0))
    # Audit fix 2026-08-26 (F-2-002 / F-3b-005): show the donor-centred AD count next to the raw
    # one, so the panel supports its own title instead of only showing the raw breadth. Source:
    # results/audit_reruns/chromvar_donor_centred_counts.csv (per-donor mean motif deviation
    # removed before the AD-vs-HA test). Per-donor deviations are not in the committed pipeline
    # outputs, which is why this comes from the audit table rather than being recomputed here.
    cc = pd.read_csv(f"{P}/results/audit_reruns/chromvar_donor_centred_counts.csv")
    cc = (cc[cc.contrast == "ADvHA"].set_index("celltype")["nsig_donor_centred"]
          .reindex(NICHE).fillna(0))
    cnt["ADvHA_centred"] = cc.values
    ax = cont.subplots()
    cnt.plot.bar(ax=ax, color={"HAvYA": "#6b8fb5", "ADvHA": PS.UP, "ADvHA_centred": "#e8a79c"},
                 width=0.80, legend=True)
    ax.set_ylabel("TF motifs at FDR < 0.1"); ax.set_xlabel("")
    ax.set_xticklabels(NICHE, rotation=0)
    for t, ct in zip(ax.get_xticklabels(), NICHE):
        t.set_color(PS.CT[ct])
    ax.legend(["normal aging (HA vs YA)", "disease (AD vs HA), raw",
               "disease (AD vs HA), donor-centred"], fontsize=6.5)
    # Audit fix 2026-08-26 (F-2-002 / F-3b-005). The old title read "AD = broad shift", but the
    # AD contrast's 1,935 FDR<0.1 motifs fall to 789 once each donor's genome-wide mean motif
    # deviation is centred out, and the same-direction share drops from 73-94% to 52-61%. The
    # breadth is largely a donor-level offset, so the title no longer asserts a broad AD shift.
    ax.set_title("Aging vs disease chromatin\naging: none; AD largely a donor-level offset", fontsize=9)
    PS.style_ax(ax, grid=True)
    return ax


def main():
    PS.apply()
    P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
    FIG = f"{P}/figures/figure4"; os.makedirs(FIG, exist_ok=True)

    # standalone panels (legend-referenced filenames)
    f = plt.figure(figsize=(6.6, 6.2), layout="constrained"); _, stat = panel_concordance(f, P)
    PS.save(f, f"{FIG}/panel_A_concordance.pdf")
    f = plt.figure(figsize=(5.0, 6.5), layout="constrained"); panel_motifs(f, P)
    PS.save(f, f"{FIG}/panel_B_aging_motifs.pdf")
    f = plt.figure(figsize=(6.2, 3.6), layout="constrained"); panel_aging_disease(f, P)
    PS.save(f, f"{FIG}/panel_C_aging_vs_disease.pdf")

    # composite figure4.pdf — top: A concordance | (B motifs over C aging-vs-disease);
    # bottom: D donor-level TF footprints full-width (embedded via draw_footprints).
    fig = plt.figure(figsize=(7.8, 9.6), layout="constrained")
    fig.suptitle("Figure 4  ·  Regulatory layer — chromatin coupling is directional, not locus-resolved (supporting)",
                 fontsize=10, fontweight="bold", ha="left", x=0.012)
    rows = fig.subfigures(2, 1, height_ratios=[5.0, 4.3])
    cols = rows[0].subfigures(1, 2, width_ratios=[1.12, 1.0])
    panel_concordance(cols[0], P); PS.letter(cols[0], "A")
    right = cols[1].subfigures(2, 1, height_ratios=[1.45, 1.0])
    panel_motifs(right[0], P, n_tfs=16); PS.letter(right[0], "B")
    panel_aging_disease(right[1], P); PS.letter(right[1], "C")
    draw_footprints(rows[1], P); PS.letter(rows[1], "D")
    rows[1].suptitle("Donor-level TF footprints, YA vs HA — no age-increase in SASP-TF / TEAD binding "
                     "(0/21 FDR<0.1; nominal trend = shallower with age)",
                     fontsize=8.5, fontweight="bold")
    PS.save(fig, f"{FIG}/figure4.pdf")

    print(f"Fig3A SUPPORTING: RNA-sig {stat['k']}/{stat['n']}={stat['k']/stat['n']:.0%} concordant "
          f"vs background {stat['bg']*100:.1f}%; enrichment binom p={stat['p_bg']:.3g}, "
          f"Fisher OR={stat['orr']:.2f} p={stat['p_fish']:.3g}; genome-wide ATAC FDR 0/{stat['tot']}; "
          f"nominal dual-sig {stat['nom_conc']}/{stat['nom']} (exploratory); ATAC-FDR survivors={stat['fdr']}")
    print("saved panel_A_concordance.pdf, panel_B_aging_motifs.pdf, panel_C_aging_vs_disease.pdf, "
          "panel_D_footprints.pdf, figure4.pdf")


if __name__ == "__main__":
    main()
