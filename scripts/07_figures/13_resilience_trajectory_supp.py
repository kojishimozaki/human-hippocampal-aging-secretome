#!/usr/bin/env python
"""Fig. S(RT1) — The frozen aging niche-secretome projected into AD and resilience cohorts
(SUPPORTING / descriptive disease-projection arm + a reported negative).

Two pre-registered freeze->test projections of the frozen 60-hit niche aging-secretome signature:
  TRAJECTORY (GSE268609 = PRIMARY cohort, YA->HA->AD; WITHIN-COHORT, not independent replication): the frozen
    aging direction is recovered in 56/60 hits (gate A>=0.25; 0 reversed) = a self-consistency check under the
    3-group re-fit; by an effect-size rule (NOT an equivalence test) it is NOT broadly amplified in the AD group
    (category-maintained 54 / amplified 1 / AD-attenuated 1 / AD-divergent 0 of 56 eligible; only 5/56 reach
    AD-vs-HA padj<0.1; robust across 9 thresholds). Cross-sectional -- NO progression claim.
  CONSEQUENCE (GSE325391 RES vs SAD; separate cohort): the granule-cell receiver-competence/consequence leg shows
    no detectable resilience association (stage-averaged SAD-RES = -0.07 z, 95% CI [-0.335, 0.195], p = 0.54; NULL).

Reproducible from committed results/{trajectory,resilience}/*.tsv (committed-data target -- runs on a
bare clone via `make figures`; no raw/processed needed). seed 42 (cosmetic jitter only; no statistics here).

Panels:
  A schematic   freeze->test of the frozen signature onto the two cohorts (drawn).
  B resilience  GSE325391 stage-averaged receiver-competence score, RES vs SAD (NULL).
  C trajectory  GSE268609 module score YA->HA->[MCI context]->AD, UP & DOWN per niche cell type.
  D categories  per-hit healthy-aging->AD effect-size category counts (within cohort; not broadly AD-amplified).
  E celltype    trajectory category x niche cell type.
  F model       interpretive summary (drawn).
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _pubstyle as PS

np.random.seed(42)  # cosmetic strip-plot jitter only (rule 10)
PROJ = os.environ.get("PROJ", os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
TRJ = f"{PROJ}/results/trajectory"
RES = f"{PROJ}/results/resilience"
FIG = f"{PROJ}/figures/figureSRT1"
os.makedirs(FIG, exist_ok=True)

# trajectory categories: key, label, colour (verdict palette; maintained = dominant/retained)
CATS = [
    ("maintained", "maintained", PS.GREEN),
    ("amplified", "amplified", PS.UP),
    ("AD_attenuated", "AD-attenuated", PS.MUSTARD),
    ("AD_divergent", "AD-divergent", PS.GREY),
    ("weak_indeterminate_aging_projection", "weak/indet.", PS.NEUTRAL),
]
GROUPS = ["YA", "HA", "MCI", "AD"]


def _letter(ax, s):
    ax.text(-0.14, 1.07, s, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")


def _box(ax, x, y, w, h, text, fc, ec, fs=6.8):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0.012,rounding_size=0.03",
                                fc=fc, ec=ec, lw=1.1, zorder=2))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, zorder=3, color=PS.INK)


# ---- A: freeze->test schematic ---------------------------------------------------------
def panel_A(ax):
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    _box(ax, 0.5, 0.85, 0.94, 0.21,
         "Frozen 60-hit niche\naging-secretome\n(GSE268609 YA vs HA, padj < 0.1)", "#eef3ee", PS.GREEN, 7.0)
    ax.annotate("", xy=(0.27, 0.47), xytext=(0.40, 0.74),
                arrowprops=dict(arrowstyle="-|>", color=PS.INK, lw=1.1))
    ax.annotate("", xy=(0.73, 0.47), xytext=(0.60, 0.74),
                arrowprops=dict(arrowstyle="-|>", color=PS.INK, lw=1.1))
    ax.text(0.5, 0.62, "freeze → test", ha="center", fontsize=6.8, style="italic", color=PS.SUB)
    _box(ax, 0.27, 0.29, 0.46, 0.30,
         "GSE268609 (primary)\nYA→HA→AD\nwithin-cohort\ndisease projection", "#eef1f6", PS.YA, 6.6)
    _box(ax, 0.73, 0.29, 0.46, 0.30,
         "GSE325391 (separate)\nRES vs SAD\nreceiver\ncompetence", "#f6eef0", PS.HA, 6.6)
    _letter(ax, "A")


# ---- B: CONSEQUENCE null (GSE325391 receiver competence) -------------------------------
def panel_B(ax):
    d = pd.read_csv(f"{RES}/donor_stage_receptor_scores.tsv", sep="\t")
    pm = pd.read_csv(f"{RES}/primary_model_results.tsv", sep="\t").iloc[0]
    da = d.groupby(["donor", "group"], as_index=False)["receiver_score"].mean()  # stage-averaged per donor
    n = da.group.value_counts()
    assert int(n["RES"]) == 6 and int(n["SAD"]) == 7, dict(n)
    order = ["RES", "SAD"]
    col = {"RES": PS.GREEN, "SAD": PS.UP}
    for i, g in enumerate(order):
        s = da.loc[da.group == g, "receiver_score"].values
        x = np.random.normal(i, 0.055, len(s))
        ax.scatter(x, s, c=col[g], s=34, edgecolor="white", linewidth=0.5, zorder=3)
        ax.hlines(s.mean(), i - 0.22, i + 0.22, color=PS.INK, lw=1.6, zorder=4)
    ax.axhline(0, color="#9aa3b0", lw=0.7)
    ax.set_xlim(-0.5, 1.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"RES\n(n={int(n['RES'])})", f"SAD\n(n={int(n['SAD'])})"], fontsize=8)
    ax.set_ylabel("receiver-competence score (z)", fontsize=8)
    ax.set_title("GSE325391 granule-lineage receiver competence\n"
                 f"SAD $-$ RES = {pm.estimate:+.2f} z [{pm.ci_lo:.3f}, {pm.ci_hi:.3f}], "
                 f"{PS.fmt_p(pm.p_two_sided)} (n.s.)", fontsize=8)
    _letter(ax, "B")


# ---- C: aging-to-AD module trajectory (GSE268609) --------------------------------------
def panel_C(ax):
    ms = pd.read_csv(f"{TRJ}/frozen60_GSE268609_module_scores.tsv", sep="\t")
    xpos = {g: i for i, g in enumerate(GROUPS)}
    for ct in PS.NICHE:
        for mod, ls, mk, fill in [("UP", "-", "o", True), ("DOWN", "--", "s", False)]:
            sub = ms[(ms.celltype == ct) & (ms.module == mod)].set_index("group").reindex(GROUPS)
            ax.errorbar([xpos[g] for g in GROUPS], sub.mean_score.values, yerr=sub.ci95.values,
                        color=PS.CT[ct], ls=ls, lw=1.0, marker=mk, ms=3.4,
                        mfc=PS.CT[ct] if fill else "white", mec=PS.CT[ct],
                        capsize=1.4, elinewidth=0.6, alpha=0.9, zorder=3)
    ax.axvspan(xpos["MCI"] - 0.42, xpos["MCI"] + 0.42, color="#eef0f4", zorder=0)
    ax.text(xpos["MCI"], 0.99, "MCI\n(context)", transform=ax.get_xaxis_transform(),
            fontsize=6, ha="center", va="top", color=PS.SUB)
    ax.axhline(0, color="#9aa3b0", lw=0.7)
    ax.set_xlim(-0.3, 3.3)
    ax.set_xticks(range(4))
    ax.set_xticklabels(GROUPS, fontsize=8)
    ax.set_ylabel("module score (z)", fontsize=8)
    ax.set_title("GSE268609 aging-to-AD module trajectory\n(colours = niche cell type, as in Fig. 2)", fontsize=8)
    ax.legend(handles=[Line2D([0], [0], color=PS.INK, ls="-", marker="o", ms=3.4, label="UP module"),
                       Line2D([0], [0], color=PS.INK, ls="--", marker="s", mfc="white", ms=3.4, label="DOWN module")],
              fontsize=6, loc="upper right", frameon=False)
    _letter(ax, "C")


# ---- D: per-hit aging->AD trajectory-category counts -----------------------------------
def panel_D(ax):
    h = pd.read_csv(f"{TRJ}/frozen60_GSE268609_hit_trajectory.tsv", sep="\t")
    vc = h.trajectory_category.value_counts()
    labels = [lab for _, lab, _ in CATS]
    vals = [int(vc.get(key, 0)) for key, _, _ in CATS]
    cols = [c for _, _, c in CATS]
    # rule-2 guards: the verified frozen counts
    assert vals == [54, 1, 1, 0, 4], vals
    assert sum(vals) == 60 and (vals[0] + vals[1] + vals[2] + vals[3]) == 56
    y = np.arange(len(labels))[::-1]
    ax.barh(y, vals, color=cols, edgecolor="white", height=0.72)
    for yi, v in zip(y, vals):
        ax.text(v + 0.7, yi, str(v), va="center", fontsize=8, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel("frozen hits (n)", fontsize=8)
    ax.set_xlim(0, max(vals) * 1.2)
    ax.set_title("Per-hit healthy-aging→AD (within cohort)\n56/60 retain aging dir.; 54 maintained (effect-size)", fontsize=8)
    _letter(ax, "D")


# ---- E: trajectory category x niche cell type -----------------------------------------
def panel_E(ax):
    h = pd.read_csv(f"{TRJ}/frozen60_GSE268609_hit_trajectory.tsv", sep="\t")
    keys = [k for k, _, _ in CATS]
    labs = [lab for _, lab, _ in CATS]
    ct = (pd.crosstab(h.celltype, h.trajectory_category)
          .reindex(index=PS.NICHE, columns=keys).fillna(0).astype(int))
    assert ct.values.sum() == 60, ct
    vmax = ct.values.max()
    ax.imshow(ct.values, cmap="Greens", aspect="auto", vmin=0, vmax=vmax)
    for i in range(ct.shape[0]):
        for j in range(ct.shape[1]):
            v = ct.values[i, j]
            ax.text(j, i, str(v), ha="center", va="center", fontsize=7.5,
                    color="white" if v > vmax * 0.6 else PS.INK)
    ax.set_xticks(range(len(labs)))
    ax.set_xticklabels(labs, fontsize=6.3, rotation=35, ha="right")
    ax.set_yticks(range(len(PS.NICHE)))
    ax.set_yticklabels([PS.CT_FULL[c] for c in PS.NICHE], fontsize=7)
    ax.set_title("Trajectory category × niche cell type", fontsize=8)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    _letter(ax, "E")


# ---- F: interpretive model -------------------------------------------------------------
def panel_F(ax):
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.5, 0.95, "Interpretation", ha="center", fontsize=8.5, fontweight="bold", color=PS.INK)
    _box(ax, 0.5, 0.66, 0.96, 0.30,
         "OUTPUT (sender side; within primary cohort):\nthe frozen aging niche-secretome is largely\n"
         "retained, not broadly AD-amplified\n(54/56; effect-size, not equivalence)", "#eef3ee", PS.GREEN, 6.9)
    _box(ax, 0.5, 0.26, 0.96, 0.30,
         "RECEIVER side (boundary / supporting negative):\nreceiver-competence/consequence leg shows\n"
         "no detectable resilience association\n(RES vs SAD, p = 0.54)", "#f0f1f3", PS.GREY, 6.9)
    _letter(ax, "F")


def main():
    PS.apply()
    fig, axes = plt.subplots(2, 3, figsize=(11.6, 7.3), layout="constrained")
    panel_A(axes[0, 0])
    panel_B(axes[0, 1])
    panel_C(axes[0, 2])
    panel_D(axes[1, 0])
    panel_E(axes[1, 1])
    panel_F(axes[1, 2])
    PS.save(fig, f"{FIG}/figureSRT1.pdf", png=True)
    print("[Fig S(RT1)] wrote", f"{FIG}/figureSRT1.pdf  (+ .png)")


if __name__ == "__main__":
    main()
