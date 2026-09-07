#!/usr/bin/env python
"""FP step 4: TF footprint meta-profiles, YA vs HA, niche cell types — Figure 4D.

Visualises the donor-level footprint result (PHASE_REPORT Layer C): no age-increase in
SASP-TF / TEAD footprints (0/21 FDR<0.1; nominal hits trend to SHALLOWER footprints with age).
Reads committed CSVs (results/regulatory/footprint_{metaprofiles,test}.csv). Exposes
draw_footprints() so the composite figure4.pdf can embed the same grid; run as __main__ to
write the standalone panel_D_footprints.pdf. Publication style in scripts/07_figures/_pubstyle.py.
"""
import os, sys
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                 "..", "07_figures")))
import _pubstyle as PS

CELLS = ["Astro", "Micro", "Endo"]
TFS = ["RELA", "IRF1", "STAT3", "CEBPB", "TEAD4"]


def _smooth(y, w=7):
    return np.convolve(y, np.ones(w) / w, mode="same")


def draw_footprints(cont, P):
    """3×5 footprint meta-profile grid (cell × TF), YA vs HA, into a Figure/SubFigure."""
    mp = pd.read_csv(f"{P}/results/regulatory/footprint_metaprofiles.csv")
    te = pd.read_csv(f"{P}/results/regulatory/footprint_test.csv").set_index(["celltype", "tf"])
    axes = cont.subplots(len(CELLS), len(TFS), sharex=True)
    for i, ct in enumerate(CELLS):
        for j, tf in enumerate(TFS):
            ax = axes[i, j]
            d = mp[(mp.celltype == ct) & (mp.tf == tf)]
            for grp, col in [("YA", PS.YA), ("HA", PS.HA)]:
                s = d[d.group == grp].sort_values("rel")
                if len(s):
                    ax.plot(s.rel, _smooth(s.cpm.values), color=col, lw=1.3, label=grp)
            ax.axvspan(-15, 15, color=PS.GRIDC, alpha=0.7, lw=0)
            try:
                r = te.loc[(ct, tf)]
                ttl = f"{ct} {tf}\nΔfp={r.fp_dz:+.3f}  {PS.fmt_p(r.fp_p)}"
            except KeyError:
                ttl = f"{ct} {tf}"
            ax.set_title(ttl, fontsize=7.5)
            ax.tick_params(labelsize=6, length=2)
            if i == len(CELLS) - 1:
                ax.set_xlabel("distance to motif (bp)", fontsize=7)
            if j == 0:
                ax.set_ylabel(f"{CELLS[i]}\nTn5 insertions (CPM)", fontsize=7.5, color=PS.CT[ct])
    handles = [Line2D([0], [0], color=PS.YA, lw=1.6, label="YA (young)"),
               Line2D([0], [0], color=PS.HA, lw=1.6, label="HA (aged)")]
    axes[0, 0].legend(handles=handles, fontsize=7, loc="upper right")
    return axes


def main():
    PS.apply()
    P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
    FIG = f"{P}/figures/figure4"; os.makedirs(FIG, exist_ok=True)
    fig = plt.figure(figsize=(4 * len(TFS), 3 * len(CELLS)), layout="constrained")
    draw_footprints(fig, P)
    fig.suptitle("TF footprints, YA vs HA, niche cells — no age-increase in SASP-TF / TEAD binding "
                 "(0/21 FDR<0.1; nominal trend = shallower footprints with age)",
                 fontsize=10, fontweight="bold")
    PS.save(fig, f"{FIG}/panel_D_footprints.pdf")
    print("saved figures/figure4/panel_D_footprints.pdf")


if __name__ == "__main__":
    main()
