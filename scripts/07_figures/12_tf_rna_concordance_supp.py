#!/usr/bin/env python
"""Fig. S(R1) — Aging-dynamic motif transcription factors' own RNA does not track
motif-accessibility direction (SUPPORTING; honest negative completing the cis-mechanism
triangulation: motif enrichment / footprinting / TF expression are all null).

Reproducible from committed results/regulatory/motif_tf_rna_concordance_*.csv
(committed-data target — runs on a bare clone via `make figures`; no raw/processed needed).

A: HAvYA foreground TF RNA log2FC vs its motif chromVAR dz (concordance scatter; foreground 50%
   concordant vs matched empirical background 52.4%, stratified-permutation p = 0.68).
B: YA->HA->AD relative-expression trajectory of the foreground TFs (descriptive; cumulative
   pairwise log2FC; HAvYA primary + ADvHA projection — no new inference).
C: GSE278576 cross-cohort RNA-direction robustness of the foreground TFs (3/14 agree; rho ~ 0).
"""
import os, sys
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _pubstyle as PS

PROJ = os.environ.get("PROJ", os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
PFX = f"{PROJ}/results/regulatory/motif_tf_rna_concordance"
FIG = f"{PROJ}/figures/figureSR1"; os.makedirs(FIG, exist_ok=True)
NICHE = ["Astro", "Micro", "Oligo", "OPC", "Endo"]


def _letter(ax, s):
    ax.text(-0.12, 1.06, s, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")


def main():
    PS.apply()
    ha = pd.read_csv(f"{PFX}_HAvYA_primary_pairs.csv")
    ad = (pd.read_csv(f"{PFX}_ADvHA_projection_pairs.csv")[["gene", "celltype", "rna_lfc"]]
          .rename(columns={"rna_lfc": "rna_lfc_AD"}))
    enr = pd.read_csv(f"{PFX}_HAvYA_primary_enrichment.csv")
    rob = pd.read_csv(f"{PFX}_GSE278576_RNA_robustness.csv")
    prim = enr[enr.label.str.contains("PRIMARY")].iloc[0]

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), layout="constrained")

    # ---- A: HAvYA concordance scatter -------------------------------------------------
    ax = axes[0]
    t = ha[ha.testable & ha.concordant.notna()]
    for conc, col, lab in [(1.0, PS.GREEN, "concordant"), (0.0, PS.MUSTARD, "discordant")]:
        s = t[t.concordant == conc]
        ax.scatter(s.dz, s.rna_lfc, c=col, s=36, edgecolor="white", linewidth=0.5,
                   label=f"{lab} (n={len(s)})", zorder=3)
    ax.axhline(0, color="#9aa3b0", lw=0.7); ax.axvline(0, color="#9aa3b0", lw=0.7)
    ax.set_xlabel("motif activity Δ (chromVAR dz, HA $-$ YA)", fontsize=8)
    ax.set_ylabel("TF RNA log2FC (HA $-$ YA)", fontsize=8)
    ax.set_title(f"TF RNA vs its own motif (HAvYA)\n{int(prim.fg_concordant)}/{int(prim.fg_testable)} = "
                 f"{prim.fg_rate:.0%} concordant vs {prim.matched_bg_rate:.0%} matched bg "
                 f"(perm p = {prim.perm_p_onesided_greater:.2f})", fontsize=8)
    ax.legend(fontsize=6.5, loc="upper left", frameon=False)
    _letter(ax, "A")

    # ---- B: YA->HA->AD relative trajectory --------------------------------------------
    ax = axes[1]
    tr = ha[ha.testable].merge(ad, on=["gene", "celltype"], how="left")
    for _, r in tr.iterrows():
        xs, ys = [0, 1], [0.0, r.rna_lfc]
        if pd.notna(r.rna_lfc_AD):
            xs.append(2); ys.append(r.rna_lfc + r.rna_lfc_AD)
        ax.plot(xs, ys, "-", color=PS.CT.get(r.celltype, "#888"), lw=1.0, alpha=0.85, marker="o", ms=3)
    ax.axhline(0, color="#9aa3b0", lw=0.7)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["YA", "HA", "AD"]); ax.set_xlim(-0.2, 2.2)
    ax.set_ylabel("relative expression (cumulative log2FC)", fontsize=8)
    ax.set_title("Foreground-TF expression trajectory\n(descriptive; pairwise contrasts, no inference)", fontsize=8)
    ax.legend(handles=[Line2D([0], [0], color=PS.CT[c], lw=2, label=c) for c in NICHE],
              fontsize=6.5, frameon=False, loc="best", ncol=2)
    _letter(ax, "B")

    # ---- C: GSE278576 cross-cohort RNA-direction robustness ---------------------------
    ax = axes[2]
    rc = rob[rob.status == "TESTABLE"]
    for agree, col, mk in [(True, PS.GREEN, "o"), (False, PS.GREY, "x")]:
        s = rc[rc.direction_agree == agree]
        ax.scatter(s.lfc_268609, s.lfc_278576, c=col, s=38, marker=mk,
                   label=f"{'agree' if agree else 'disagree'} (n={len(s)})", zorder=3)
    lim = float(np.nanmax(np.abs(np.r_[rc.lfc_268609.values, rc.lfc_278576.values]))) * 1.15
    ax.plot([-lim, lim], [-lim, lim], "--", color="#cccccc", lw=0.7, zorder=1)
    ax.axhline(0, color="#9aa3b0", lw=0.6); ax.axvline(0, color="#9aa3b0", lw=0.6)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("TF RNA log2FC — GSE268609 (HA $-$ YA)", fontsize=8)
    ax.set_ylabel("TF RNA log2FC — GSE278576 (old $-$ young)", fontsize=8)
    ax.set_title(f"Cross-cohort RNA-direction robustness\n{int(rc.direction_agree.sum())}/{len(rc)} agree "
                 f"(limited / mixed)", fontsize=8)
    ax.legend(fontsize=6.5, frameon=False, loc="upper left")
    _letter(ax, "C")

    PS.save(fig, f"{FIG}/figureSR1.pdf", png=True)
    print("saved", f"{FIG}/figureSR1.pdf")


if __name__ == "__main__":
    main()
