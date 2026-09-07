#!/usr/bin/env python
"""Fig S(E1) — external-cohort reproduction spectrum of the frozen niche-secretome signature.

Naming note (2026-08-26): B-4 (commit 4e3d27b) replaced the "independent replication" claim label
with "reproduction" — the matched-background re-test shows the agreement is not demonstrably
specific to the 60 hits. The manuscript legend was changed then; this figure was not. Display text
here now says "reproduction"; the committed CSV filenames (replication_*.csv) keep the
pre-registered metric name on purpose, so provenance stays traceable.

Compact, committed-data supplement that contrasts ALL external cohorts on one descriptive
axis, so the main Figure 5 can carry only the Visium spatial scope and is
not visually over-weighted by a negative/boundary panel. This figure replaces the former
main "Fig 5A" cross-cohort panel (GSE199243), which is here demoted to a clearly-labelled
boundary comparison alongside the genuine GSE278576 replication and the GSE186538 sensitivity.

Panel A — cohort-level dashboard: pooled Spearman(primary secretome log2FC, external output),
          frozen-hit sign-concordance vs each cohort's own background, and the pre-specified
          replication pass-count (the five GSE278576/GSE186538 outputs only). GSE199243 is shown
          as a boundary comparison, NOT scored as a pre-specified replication.
Panel B — six pooled primary-vs-external scatter small-multiples (grey = shared secretome genes,
          accent = frozen primary hits); the glia-only GSE199243 panel is visually separated and
          labelled "boundary comparison, not primary replication".

Provenance (rule 7 — identical method/thresholds as the frozen-signature replication):
  * The five pre-specified cohorts' pooled rho / hit-concordance / background / pass-count are read
    from the committed results/validation/replication_<COHORT>.csv (the canonical A1 outputs that
    Fig 3 reports), NOT recomputed for the headline numbers.
  * GSE199243's analogous DESCRIPTIVE metrics are computed here from the same frozen tables, reusing
    the EXACT functions (onesided_spearman_gt0, concordance) imported from the primary replication
    engine scripts/03_de/05_external_replication.py — so the comparison is apples-to-apples — but it
    is labelled a boundary comparison and is never counted toward the pre-specified pass tally.
  * The per-gene scatter points for ALL six cohorts are the same full-secretome×niche merge used by
    the engine; for the five pre-specified cohorts the recomputed pooled rho is asserted to match the
    committed replication_*.csv value (self-check printed).

Committed-data (bare clone): reads results/validation/{frozen_primary_signature,
frozen_primary_secretome_logfc,replication_*}.csv + results/de/de_*_per_celltype.csv (all tracked).
Writes results/validation/replication_spectrum_summary.csv + figures/figureSE1/figureSE1.pdf. seed 42.
"""
import os, sys, importlib.util
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_cache")   # matplotlib cache on a read-only/cacheless HOME
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _pubstyle as PS

SEED = 42

# key, de table, replication_*.csv (None => compute descriptively), display, group, role, y-axis, scored
COHORTS = [
    ("GSE278576",        "de_GSE278576_per_celltype.csv",        "replication_GSE278576.csv",
     "GSE278576",             "GSE278576 reproduction family", "reproduction (3/3 pre-spec.)", "external log$_2$FC", True),
    ("GSE278576_sexadj", "de_GSE278576_sexadj_per_celltype.csv", "replication_GSE278576_sexadj.csv",
     "GSE278576 (sex-adj.)",  "GSE278576 reproduction family", "sensitivity",                  "external log$_2$FC", True),
    ("GSE278576_age",    "de_GSE278576_age_per_celltype.csv",    "replication_GSE278576_age.csv",
     "GSE278576 (cont. age)", "GSE278576 reproduction family", "sensitivity",                  "external log$_2$FC", True),
    ("GSE186538",        "de_GSE186538_per_celltype.csv",        "replication_GSE186538.csv",
     "GSE186538",             "GSE186538 leverage-sensitive",  "supportive (leverage-limited)","external log$_2$FC", True),
    ("GSE186538_drop79", "de_GSE186538_drop79_per_celltype.csv", "replication_GSE186538_drop79.csv",
     "GSE186538 (drop 79 y)", "GSE186538 leverage-sensitive",  "leverage sensitivity",         "external log$_2$FC", True),
    ("GSE199243v2",      "de_GSE199243v2_per_celltype.csv",      None,
     "GSE199243 (glia)",      "GSE199243 boundary comparison", "boundary — not scored",        "glia age slope",     False),
]
GROUP_COLOR = {
    "GSE278576 reproduction family": PS.GREEN,
    "GSE186538 leverage-sensitive": PS.MUSTARD,
    "GSE199243 boundary comparison": PS.GREY,
}


def load_engine(P):
    spec = importlib.util.spec_from_file_location(
        "ext_rep", f"{P}/scripts/03_de/05_external_replication.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def compute(P):
    """Return per-cohort dict with pooled rho, hit/bg concordance, pass-count, and scatter points."""
    ER = load_engine(P)
    full = pd.read_csv(f"{P}/results/validation/frozen_primary_secretome_logfc.csv")
    hits = pd.read_csv(f"{P}/results/validation/frozen_primary_signature.csv")
    hitkeys = set(zip(hits.gene, hits.celltype))
    out = []
    for key, de, repl, disp, group, role, ylab, scored in COHORTS:
        ext = pd.read_csv(f"{P}/results/de/{de}")[["gene", "celltype", "log2FoldChange", "padj"]].rename(
            columns={"log2FoldChange": "lfc_ext", "padj": "padj_ext"})
        fp = full.merge(ext, on=["gene", "celltype"], how="inner")
        fp["is_hit"] = [(g, c) in hitkeys for g, c in zip(fp.gene, fp.celltype)]
        rho_rec, _, n_shared = ER.onesided_spearman_gt0(fp.log2FoldChange.values, fp.lfc_ext.values)
        co = ER.concordance(fp[fp.is_hit].assign(concord=np.sign(fp[fp.is_hit].log2FoldChange) == np.sign(fp[fp.is_hit].lfc_ext)),
                            fp[~fp.is_hit].assign(concord=np.sign(fp[~fp.is_hit].log2FoldChange) == np.sign(fp[~fp.is_hit].lfc_ext)),
                            np.random.default_rng(SEED))
        if repl is not None:                      # canonical numbers from the committed A1 output
            rep = pd.read_csv(f"{P}/results/validation/replication_{key}.csv")
            pr = rep[rep.scope == "POOLED"].set_index("metric")
            rho = float(pr.loc["(ii)spearman", "statistic"])
            hitc = float(pr.loc["(iii)concordance", "statistic"])
            bg = float(pr.loc["(iii)concordance", "note"].split("bg=")[1].split()[0])
            npass = int(sum(str(pr.loc[m, "pass"]) == "True"
                            for m in ["(i)directional", "(ii)spearman", "(iii)concordance"]))
            assert abs(rho - rho_rec) < 1e-6, f"{key}: recomputed rho {rho_rec} != committed {rho}"
            passtxt = f"{npass}/3 pre-spec."
        else:                                     # GSE199243: descriptive, identical method, NOT scored
            rho, hitc, bg, npass = rho_rec, co["hit_concord"], co["bg_concord"], None
            passtxt = "boundary (not scored)"
        out.append(dict(key=key, display=disp, group=group, role=role, ylab=ylab, scored=scored,
                        pooled_rho=rho, n_shared=int(n_shared), hit_concord=hitc, n_hit=int(co["n_hit"]),
                        bg_concord=bg, metrics_passed=npass, passtxt=passtxt,
                        x=fp.log2FoldChange.values, y=fp.lfc_ext.values, hitmask=fp.is_hit.values))
    return out


def panel_dashboard(ax, rows):
    """Horizontal bar of pooled rho per cohort, colour-coded by group, with concordance + pass annotation."""
    n = len(rows)
    ypos = np.arange(n)[::-1]                      # GSE278576 at top
    ax.axvline(0, c=PS.INK, lw=0.8)
    for y, r in zip(ypos, rows):
        ax.barh(y, r["pooled_rho"], height=0.62, color=GROUP_COLOR[r["group"]],
                alpha=0.9 if r["scored"] else 0.7, edgecolor=PS.INK, linewidth=0.5,
                hatch="" if r["scored"] else "////")
        ax.text(r["pooled_rho"] + (0.006 if r["pooled_rho"] >= 0 else -0.006), y,
                f"{r['pooled_rho']:+.2f}", va="center",
                ha="left" if r["pooled_rho"] >= 0 else "right", fontsize=8, fontweight="bold")
        ax.text(0.30, y, f"hit-conc {r['hit_concord']:.2f} / bg {r['bg_concord']:.2f}    {r['passtxt']}",
                va="center", ha="left", fontsize=7.5, color=PS.SUB)
    ax.set_yticks(ypos); ax.set_yticklabels([r["display"] for r in rows], fontsize=8.5)
    ax.set_xlim(-0.05, 0.62); ax.set_xlabel("pooled Spearman ρ (primary secretome log$_2$FC vs external output)", fontsize=8.5)
    ax.set_title("Cohort-level summary — pooled effect-correlation, frozen-hit concordance, pre-specified pass count",
                 fontsize=9.5)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="x", length=3); ax.tick_params(axis="y", length=0)
    leg = [Patch(facecolor=GROUP_COLOR["GSE278576 reproduction family"], edgecolor=PS.INK, lw=0.5, label="GSE278576 reproduction family (Fig. 3)"),
           Patch(facecolor=GROUP_COLOR["GSE186538 leverage-sensitive"], edgecolor=PS.INK, lw=0.5, label="GSE186538 leverage-sensitive support"),
           Patch(facecolor=GROUP_COLOR["GSE199243 boundary comparison"], edgecolor=PS.INK, lw=0.5, hatch="////", label="GSE199243 boundary comparison (not scored)")]
    ax.legend(handles=leg, fontsize=7, loc="lower right", frameon=False)


def panel_scatter(cont, rows):
    """2×3 pooled primary-vs-external scatters; GSE199243 panel tinted + labelled as a boundary."""
    axes = cont.subplots(2, 3).ravel()
    xall = rows[0]["x"]
    xlo, xhi = np.quantile(xall, 0.01), np.quantile(xall, 0.99)
    xlim = max(abs(xlo), abs(xhi)) * 1.1
    for ax, r in zip(axes, rows):
        if not r["scored"]:
            ax.set_facecolor("#f3f1ee")            # tint the boundary-comparison panel
        ax.axhline(0, c=PS.GRIDC, lw=0.8); ax.axvline(0, c=PS.GRIDC, lw=0.8)
        ax.scatter(r["x"], r["y"], s=6, c=PS.NEUTRAL, alpha=0.45, lw=0, rasterized=True,
                   label="shared secretome genes")
        ax.scatter(r["x"][r["hitmask"]], r["y"][r["hitmask"]], s=22, c=PS.SASP, lw=0, zorder=3,
                   label="frozen primary hits")
        yr = np.quantile(np.abs(r["y"]), 0.99) * 1.15 + 1e-6
        ax.set_xlim(-xlim, xlim); ax.set_ylim(-yr, yr)
        title_c = GROUP_COLOR[r["group"]] if r["scored"] else PS.INK
        ax.set_title(f"{r['display']}\nρ = {r['pooled_rho']:+.2f}  ·  hit-conc {r['hit_concord']:.2f}",
                     fontsize=8.5, color=title_c)
        ax.set_xlabel("primary log$_2$FC (GSE268609)", fontsize=7.5)
        ax.set_ylabel(r["ylab"], fontsize=7.5)
        if not r["scored"]:
            ax.text(0.5, 0.97, "boundary comparison\n(not primary replication)", transform=ax.transAxes,
                    ha="center", va="top", fontsize=7.5, fontweight="bold", color=PS.GREY)
        PS.style_ax(ax)
    axes[0].legend(fontsize=6.5, loc="upper left", markerscale=1.2)
    return axes


def main():
    PS.apply()
    np.random.seed(SEED)
    P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
    FIG = f"{P}/figures/figureSE1"; os.makedirs(FIG, exist_ok=True)
    rows = compute(P)

    # ---- Source-Data summary (committed) ----
    summ = pd.DataFrame([{k: r[k] for k in ("key", "display", "group", "role", "scored",
                          "pooled_rho", "n_shared", "hit_concord", "n_hit", "bg_concord", "metrics_passed")}
                         for r in rows]).rename(columns={"key": "cohort"})
    os.makedirs(f"{P}/results/validation", exist_ok=True)
    summ.to_csv(f"{P}/results/validation/replication_spectrum_summary.csv", index=False)
    print("wrote results/validation/replication_spectrum_summary.csv")
    print(summ.to_string(index=False))

    # ---- composite figureSE1.pdf (A dashboard over B scatters) ----
    fig = plt.figure(figsize=(11.0, 9.2), layout="constrained")
    fig.suptitle("Figure S(E1)  ·  External-cohort reproduction spectrum of the frozen niche-secretome "
                 "signature (GSE278576 reproduction → GSE186538 sensitivity → GSE199243 boundary)",
                 fontsize=10.5, fontweight="bold", ha="left", x=0.012)
    rowsfig = fig.subfigures(2, 1, height_ratios=[0.78, 1.7])
    axA = rowsfig[0].subplots(1, 1); panel_dashboard(axA, rows); PS.letter(rowsfig[0], "A")
    panel_scatter(rowsfig[1], rows); PS.letter(rowsfig[1], "B")
    PS.save(fig, f"{FIG}/figureSE1.pdf")
    print("saved figures/figureSE1/figureSE1.pdf")


if __name__ == "__main__":
    main()
