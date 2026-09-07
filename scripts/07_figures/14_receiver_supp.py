#!/usr/bin/env python
"""Fig. S(Rcv) — Receiver compartments show no coordinated receptor up-regulation
(SUPPORTING; pre-registered negative completing the OUTPUT->receiver triangulation).

The 76 NicheNet receptors of the 22 aging-up secreted ligands, mapped onto donor-pseudobulk
receptor regulation across all 11 hippocampal compartments (5 niche senders also tested as
receivers; 3 mature-neuron types) plus, as an exploratory family, the deposit's three rare clusters
annotated NSC/Neuroblast/Immature. Those three labels are NOT supported by re-analysis (audit finding
F-4-001 / BLOCKER B-1: marker-floor re-assignment retains 0-0.9% of their nuclei under their own label
against 83-99.8% for controls), so the right-hand block is shown for transparency and is explicitly not
read as a neurogenic result. Corrected 2026-08-26. Primary test:
matched-null on mean signed log2FC, oriented aging-up (seed 42, 2000x). Pre-reg tags
prereg/receiver-map-v1 (+ v1.1). Result: 0/8 confirmatory pass (BH-q = 0.986); no shift survives
FDR in either direction (two-sided BH min-q = 0.12); the NSC-labelled cluster's repertoire is depth-explained, not biology.

Reproducible from committed results/receiver/*.csv (committed-data target — runs on a bare clone
via `make figures`; no raw/processed needed). seed 42 (no statistics here; CSVs are pre-computed).

A: receptor testability x regulation matrix — 76 receptors (rows) x 11 compartments (cols),
   coloured by aging log2FC where the receptor is testable (non-NA padj, after independent filtering), grey where not testable.
B: matched-null-on-lfc forest — observed mean signed receptor log2FC vs its matched-null per
   compartment, the pre-specified FAIL calls, and the NSC v1.1 depth/testability annotation.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _pubstyle as PS

np.random.seed(42)
PROJ = os.environ.get("PROJ", os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
R = f"{PROJ}/results/receiver"
FIG = f"{PROJ}/figures/figureS_receiver"; os.makedirs(FIG, exist_ok=True)

CONF = ["Astro", "Micro", "Oligo", "OPC", "Endo", "DG_GC", "CA_ExN", "InN"]
# Right-hand block: depositor annotations only. Re-analysis does not support these labels
# (F-4-001 / B-1), so they are exploratory and are NOT interpreted as neurogenic receivers.
NEUR = ["NSC", "Neuroblast", "Immature"]
ORDER = CONF + NEUR


def _letter(ax, s):
    ax.text(-0.13, 1.04, s, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")


def main():
    PS.apply()
    mat = pd.read_csv(f"{R}/detectability_matrix.csv")
    for c in ("in_table", "detectable"):
        mat[c] = mat[c].astype(str).str.lower().eq("true")
    univ = pd.read_csv(f"{R}/frozen_receptor_universe_76.csv")          # receptor, n_ligands_connecting
    mn = pd.read_csv(f"{R}/matched_null_per_compartment.csv")
    dt = pd.read_csv(f"{R}/depth_testability_v1_1.csv")

    # receptor order: by ligand connectivity then name
    nlig = dict(zip(univ.receptor, univ.n_ligands_connecting))
    receptors = sorted(nlig, key=lambda r: (-nlig.get(r, 0), r))

    # lfc matrix (detectable only) + detectability mask
    lfc = mat.pivot(index="receptor", columns="compartment", values="log2FoldChange").reindex(index=receptors, columns=ORDER)
    det = mat.pivot(index="receptor", columns="compartment", values="detectable").reindex(index=receptors, columns=ORDER).fillna(False)
    M = np.where(det.to_numpy(), lfc.to_numpy(), np.nan)

    fig = plt.figure(figsize=(13.2, 9.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.18, 1.0], wspace=0.34)

    # ---- A: detectability x regulation matrix ----------------------------------------
    axA = fig.add_subplot(gs[0, 0])
    cmap = plt.cm.RdBu_r.copy(); cmap.set_bad("#e9ecf1")
    im = axA.imshow(np.ma.masked_invalid(M), aspect="auto", cmap=cmap, vmin=-1, vmax=1)
    axA.set_xticks(range(len(ORDER))); axA.set_xticklabels(ORDER, rotation=45, ha="right", fontsize=7.5)
    axA.set_yticks(range(len(receptors))); axA.set_yticklabels(receptors, fontsize=4.3)
    axA.tick_params(length=0)
    sep = len(CONF) - 0.5
    axA.axvline(sep, color=PS.INK, lw=1.4)          # confirmatory | depositor-annotated rare clusters
    # group headers, so the panel cannot be read as testing a neurogenic lineage (B-1)
    axA.text((len(CONF) - 1) / 2, -1.6, "8 pre-specified confirmatory compartments",
             ha="center", va="bottom", fontsize=7.5, color=PS.INK, fontweight="bold")
    axA.text(len(CONF) + (len(NEUR) - 1) / 2, -1.6,
             "depositor-annotated rare clusters\n(exploratory; labels not supported — not a neurogenic result)",
             ha="center", va="bottom", fontsize=6.6, color=PS.SUB, style="italic")
    # testable (non-NA padj) count per compartment along the bottom
    for j, ct in enumerate(ORDER):
        n = int(det[ct].sum())
        axA.text(j, len(receptors) + 0.8, f"{n}", ha="center", va="top", fontsize=6.5, color=PS.INK)
    axA.text(-0.5, len(receptors) + 0.8, "n testable:", ha="right", va="top", fontsize=6.5, color=PS.SUB)
    cb = fig.colorbar(im, ax=axA, fraction=0.035, pad=0.02)
    cb.set_label("aging log2FC (HA $-$ YA), where testable", fontsize=7.5); cb.ax.tick_params(labelsize=6.5)
    axA.set_title("Receptor testability × aging regulation\n(76 ligand-mapped receptors × 8 confirmatory + 3 exploratory compartments; non-NA padj after independent filtering; grey = not testable)", fontsize=9, pad=26)
    _letter(axA, "A")

    # ---- B: matched-null-on-lfc forest -----------------------------------------------
    axB = fig.add_subplot(gs[0, 1])
    mn = mn.set_index("compartment")
    ys, ylab, dots, nulls, cols, qann = [], [], [], [], [], []
    y = 0
    for block in (CONF, NEUR):
        for ct in block:
            r = mn.loc[ct]
            obs = r["obs_mean_signed_lfc"]; nul = r["null_mean_signed_lfc"]
            if pd.isna(obs):
                continue
            ys.append(y); ylab.append(ct); dots.append(obs); nulls.append(nul)
            if ct in NEUR:
                cols.append(PS.CT["NSC"]); qann.append("exploratory")
            else:
                cols.append(PS.GREY); qann.append(f"q={r['bh_q_confirmatory']:.3f}")
            y -= 1
        y -= 0.7
    ys = np.array(ys)
    axB.axvline(0, color="#9aa3b0", lw=0.8, zorder=1)
    for yi, o, nu in zip(ys, dots, nulls):
        axB.plot([nu, o], [yi, yi], "-", color="#c9ced6", lw=1.2, zorder=2)
    axB.scatter(nulls, ys, marker="|", s=120, color="#7a8394", zorder=3, label="matched-null mean")
    axB.scatter(dots, ys, s=64, c=cols, edgecolor="white", linewidth=0.6, zorder=4, label="observed receptor-set mean")
    for yi, q in zip(ys, qann):
        axB.text(0.97, yi, q, transform=axB.get_yaxis_transform(), ha="right", va="center",
                 fontsize=6.5, color=(PS.SUB if "expl" in q else PS.INK))
    axB.set_yticks(ys); axB.set_yticklabels(ylab, fontsize=8)
    axB.set_ylim(min(ys) - 3.8, max(ys) + 2.6)                          # top headroom for legend; bottom margin for note
    axB.set_xlabel("mean signed receptor log2FC (oriented aging-up)", fontsize=8.5)
    axB.set_xlim(-0.4, 0.5)
    axB.set_title("Matched-null receptor remodeling: 8/8 confirmatory FAIL\n"
                  "directional BH-q = 0.986 · two-sided BH min-q = 0.12 (no FDR-surviving shift, either direction)",
                  fontsize=8.5)
    axB.legend(fontsize=6.8, loc="upper right", frameon=False, bbox_to_anchor=(1.0, 1.0))  # top headroom, clear of the bottom note
    # NSC depth annotation (v1.1) — alone in the empty bottom margin, below all dots
    nsc = dt.set_index("compartment").loc["NSC"]
    axB.text(0.015, 0.015,
             f"NSC-labelled cluster, depth-explained (v1.1): {int(nsc['n_in_table'])}/76 in DE table;\n"
             f"{int(nsc['n_detectable_in_table'])}/{int(nsc['n_in_table'])} testable after independent filtering "
             f"vs in-table null {nsc['matched_null_mean']:.0f}, p_low={nsc['p_low_depleted']:.2f}.\n"
             f"Neuroblast-/Immature-labelled: testability saturated — coverage descriptive only, not depletion.\n"
             f"These three labels are depositor annotations and are not supported by re-analysis (B-1).",
             transform=axB.transAxes, fontsize=6.3, color=PS.SUB, va="bottom", ha="left",
             bbox=dict(boxstyle="round,pad=0.35", fc="#f4f6f9", ec="#d7dce3", lw=0.6))
    _letter(axB, "B")

    PS.save(fig, f"{FIG}/figureS_receiver.pdf", png=True)
    print("saved", f"{FIG}/figureS_receiver.pdf")


if __name__ == "__main__":
    main()
