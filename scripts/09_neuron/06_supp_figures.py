#!/usr/bin/env python
"""Figure S(N) — the aging niche-secretome is not reproduced in mature neurons (cell-type control).

A SUPPLEMENT that armors the GeroScience niche-secretome submission. In the SAME
GSE268609 donors, mature hippocampal neurons (DG_GC / CA_ExN / InN) do NOT
reproduce the niche aging-secretome — positive evidence for niche-enrichment, an
external-validity axis complementary to the spatial (no detectable DG concentration) and
cross-cohort axes. A descriptive add-on (D) records that InN nonetheless carries an
intrinsic, ambient-bounded inflammaging-like shift.

Panels (each writes its legend-referenced standalone PDF + assembles into the
composite figureS_neuron.pdf):
  A  pairwise secretome-log2FC Spearman ρ by pair-type (niche–niche / neuron–neuron /
     niche–neuron) vs the external reproduction benchmark ρ = 0.24
     -> results/validation/neuron_niche_L1_rho_aging.csv
  B  anchored frozen-60 concordance in neurons vs each neuron's OWN secretome
     background (the floor-cancelling test) -> neuron_niche_L1_anchored60.csv
  C  de-novo aging-secretome hit counts (adj. p < 0.1), niche vs neuron
     -> results/de/de_GSE268609_neuron_N0_summary.csv (secr_padj0_1)
  D  InN targeted niche-response sets (median log2FC; ambient-bounded, descriptive)
     -> results/validation/neuron_L2_targeted_response.csv

Every number is read live from the committed CSVs (no hard-coding). Donor-level
units; aging axis (HA vs YA) throughout. Committed-data target (bare-clone safe).
Caller owns the seed (rule 10).
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "07_figures"))
import _pubstyle as PS

PS.apply()
np.random.seed(42)                     # rule 10: one seed; jitter below uses default_rng(42)
P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
FIG = f"{P}/figures/figureS_neuron"; os.makedirs(FIG, exist_ok=True)
NPERM = 2000                            # permutation resolution -> null_p floor = 1/NPERM

# ---- committed inputs -------------------------------------------------------
rho  = pd.read_csv(f"{P}/results/validation/neuron_niche_L1_rho_aging.csv")
anc  = pd.read_csv(f"{P}/results/validation/neuron_niche_L1_anchored60.csv")
nsum = pd.read_csv(f"{P}/results/de/de_GSE268609_neuron_N0_summary.csv")
tr   = pd.read_csv(f"{P}/results/validation/neuron_L2_targeted_response.csv")

EXT_BENCH = 0.24                        # GSE278576 external-reproduction ρ (manuscript)


def _pfloor(p):
    """Journal p-string with a permutation floor: null_p==0 -> '< 1/NPERM' (= 5e-4)."""
    return rf"$p < 5\times10^{{-4}}$" if p <= 0 else PS.fmt_p(p)


# ---- A: ρ by pair-type ------------------------------------------------------
PAIR_ORDER = ["niche-niche", "neuron-neuron", "niche-neuron"]
PAIR_DISP  = {"niche-niche": "Niche–niche", "neuron-neuron": "Neuron–neuron",
              "niche-neuron": "Niche–neuron"}
PAIR_COL   = {"niche-niche": PS.CT["Astro"], "neuron-neuron": PS.CT["Oligo"],
              "niche-neuron": PS.INK}


def panel_rho(cont):
    ax = cont.subplots()
    rng = np.random.default_rng(42)
    meds = {}
    for i, pt in enumerate(PAIR_ORDER):
        v = rho.loc[rho.pair_type == pt, "rho"].values
        ax.scatter(np.full(len(v), i) + rng.uniform(-0.12, 0.12, len(v)), v,
                   s=26, color=PAIR_COL[pt], alpha=0.8, lw=0)
        m = float(np.median(v)); meds[pt] = m
        ax.plot([i - 0.26, i + 0.26], [m, m], color=PS.INK, lw=2.2)
        ax.text(i, m + 0.022, f"{m:.2f}", ha="center", va="bottom",
                fontsize=8.5, fontweight="bold")
    ax.axhline(EXT_BENCH, ls="--", lw=0.9, color=PS.SUB)
    ax.text(2.46, EXT_BENCH + 0.004, f"external reproduction ρ = {EXT_BENCH:.2f}",
            ha="right", va="bottom", fontsize=7, color=PS.SUB)
    ax.axhline(0, color=PS.INK, lw=0.6)
    ax.set_xlim(-0.5, 2.5)
    ax.set_xticks(range(3)); ax.set_xticklabels([PAIR_DISP[o] for o in PAIR_ORDER], fontsize=8)
    ax.set_ylabel("Spearman ρ (secretome log$_2$FC)")
    ax.set_title("Neurons resemble the niche less than\nniche types resemble each other",
                 fontsize=9.5)
    PS.style_ax(ax, grid=True)
    return meds


# ---- B: anchored frozen-60 vs own background --------------------------------
SCOPES = ["DG_GC", "CA_ExN", "InN", "POOLED_neuron"]
SCOPE_DISP = {"DG_GC": "DG_GC", "CA_ExN": "CA_ExN", "InN": "InN", "POOLED_neuron": "Pooled"}


def panel_anchored(cont):
    ax = cont.subplots()
    x = np.arange(len(SCOPES)); w = 0.38
    hit = [float(anc.loc[anc.scope == s, "hit_concord"].iloc[0]) for s in SCOPES]
    bg  = [float(anc.loc[anc.scope == s, "bg_concord"].iloc[0]) for s in SCOPES]
    above = [int(anc.loc[anc.scope == s, "concordant_above_bg"].iloc[0]) for s in SCOPES]
    ax.bar(x - w / 2, hit, w, color=PS.SUB, alpha=0.92, label="Frozen 60-hit signature")
    ax.bar(x + w / 2, bg, w, color=PS.NEUTRAL, alpha=0.95, label="Secretome background")
    ax.axhline(0.5, ls=":", lw=0.8, color=PS.INK)
    for i in range(len(SCOPES)):
        ax.text(i, max(hit[i], bg[i]) + 0.015, "n.s.", ha="center", va="bottom",
                fontsize=7.5, color=PS.SUB)
    ax.set_ylim(0, 0.82)
    ax.set_xticks(x); ax.set_xticklabels([SCOPE_DISP[s] for s in SCOPES], fontsize=8)
    ax.set_ylabel("Fraction concordant with\nniche aging direction")
    ax.set_title("Frozen niche signature not reproduced in\nneurons above their own background",
                 fontsize=9.5)
    ax.legend(loc="upper left", fontsize=7)
    ax.text(0.98, 0.03, f"{sum(above)}/{len(SCOPES)} scopes above background",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
            style="italic", color=PS.SUB)
    PS.style_ax(ax, grid=True)
    return hit, bg, above


# ---- C: de-novo secretome hit counts, niche vs neuron -----------------------
NICHE = ["Astro", "Micro", "Oligo", "OPC", "Endo"]
NEUR  = ["DG_GC", "CA_ExN", "InN"]


def panel_denovo(cont):
    ax = cont.subplots()
    s = nsum[nsum.contrast == "aging"].set_index("celltype")
    cells = NICHE + NEUR
    vals = [int(s.loc[c, "secr_padj0_1"]) for c in cells]
    cols = [PS.CT[c] for c in NICHE] + [PS.INK] * len(NEUR)
    x = np.arange(len(cells))
    ax.bar(x, vals, color=cols, width=0.72, alpha=0.92)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.25, str(v), ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    ax.axvline(len(NICHE) - 0.5, color=PS.SUB, lw=0.8, ls=":")
    ax.set_xticks(x); ax.set_xticklabels(cells, rotation=35, ha="right", fontsize=7.5)
    ax.set_ylabel("Secretome hits (adj. $p$ < 0.1)")
    ax.set_ylim(0, max(vals) * 1.18)
    ax.set_title("Neurons reprogram far fewer\nsecretome genes with aging", fontsize=9.5)
    ax.text((len(NICHE) - 1) / 2, max(vals) * 1.08, "Niche", ha="center",
            fontsize=8, color=PS.SUB, fontweight="bold")
    ax.text(len(NICHE) + (len(NEUR) - 1) / 2, max(vals) * 1.08, "Neuron", ha="center",
            fontsize=8, color=PS.SUB, fontweight="bold")
    PS.style_ax(ax, grid=True)
    return dict(zip(cells, vals))


# ---- D: InN targeted niche-response (descriptive; ambient-bounded) ----------
INN_SETS = ["IFN_alpha", "IFN_gamma", "NFkB_TNFA", "Complement", "Inflammatory", "SenMayo"]
SET_DISP = {"IFN_alpha": "IFN-α", "IFN_gamma": "IFN-γ", "NFkB_TNFA": "NF-κB/TNFα",
            "Complement": "Complement", "Inflammatory": "Inflammatory", "SenMayo": "SenMayo"}


def panel_inn(cont):
    ax = cont.subplots()
    d = tr[(tr.axis == "aging") & (tr.celltype == "InN")].set_index("set").loc[INN_SETS]
    y = np.arange(len(INN_SETS))[::-1]            # first row on top
    vals = d.median_lfc.values
    enr = d.up_enriched.values.astype(bool)
    cols = [PS.UP if e else PS.NEUTRAL for e in enr]
    ax.barh(y, vals, color=cols, alpha=0.92)
    for yi, v, p, e in zip(y, vals, d.null_p.values, enr):
        if e:
            ax.text(v + 0.006, yi, _pfloor(p), va="center", fontsize=6.8, color=PS.UP)
    ax.axvline(0, color=PS.INK, lw=0.6)
    ax.set_yticks(y); ax.set_yticklabels([SET_DISP[s] for s in INN_SETS], fontsize=8)
    ax.set_xlim(min(0, vals.min()) - 0.03, max(vals) + 0.11)
    ax.set_xlabel("InN median log$_2$FC (set, vs background)")
    ax.set_title("InN intrinsic inflammaging-like shift\n(set-level; ambient-bounded, descriptive)",
                 fontsize=9.5)
    ax.tick_params(length=3)
    return {s: (float(d.loc[s, "median_lfc"]), float(d.loc[s, "null_p"])) for s in INN_SETS}


# ---- standalone panels (legend-referenced filenames) ------------------------
f = plt.figure(figsize=(4.0, 3.4), layout="constrained"); meds = panel_rho(f)
PS.save(f, f"{FIG}/panel_A_rho_pairtype.pdf")
f = plt.figure(figsize=(4.4, 3.4), layout="constrained"); hit, bg, above = panel_anchored(f)
PS.save(f, f"{FIG}/panel_B_anchored60.pdf")
f = plt.figure(figsize=(4.4, 3.4), layout="constrained"); cnts = panel_denovo(f)
PS.save(f, f"{FIG}/panel_C_denovo_counts.pdf")
f = plt.figure(figsize=(4.2, 3.2), layout="constrained"); inn = panel_inn(f)
PS.save(f, f"{FIG}/panel_D_InN_response.pdf")

# ---- composite figureS_neuron.pdf -------------------------------------------
fig = plt.figure(figsize=(8.0, 7.0), layout="constrained")
fig.suptitle("Figure S(N)  ·  The aging niche-secretome is not "
             "reproduced in mature neurons",
             fontsize=10.5, fontweight="bold", ha="left", x=0.012)
sf = fig.subfigures(2, 2)
panel_rho(sf[0][0]);      PS.letter(sf[0][0], "A")
panel_anchored(sf[0][1]); PS.letter(sf[0][1], "B")
panel_denovo(sf[1][0]);   PS.letter(sf[1][0], "C")
panel_inn(sf[1][1]);      PS.letter(sf[1][1], "D")
PS.save(fig, f"{FIG}/figureS_neuron.pdf")   # PDF only (data-figure convention; no committed .png)

# ---- console trace (numbers must match manuscript) --------------------------
print("A  ρ medians by pair-type:", {PAIR_DISP[k]: round(v, 3) for k, v in meds.items()})
print("B  anchored-60 hit vs bg:",
      {SCOPE_DISP[s]: (round(hit[i], 3), round(bg[i], 3)) for i, s in enumerate(SCOPES)},
      "| scopes above background:", sum(above))
print("C  de-novo secretome hits (adj.p<0.1):", cnts)
print("D  InN set median log2FC / null_p:",
      {SET_DISP[k]: (round(v[0], 3), v[1]) for k, v in inn.items()})
print("saved -> figures/figureS_neuron/: panel_A_rho_pairtype.pdf, panel_B_anchored60.pdf, "
      "panel_C_denovo_counts.pdf, panel_D_InN_response.pdf, figureS_neuron.pdf (+ .png)")
