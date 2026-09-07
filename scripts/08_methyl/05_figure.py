#!/usr/bin/env python
"""M1 methylome — results figure (Fig S-M1). Descriptive display of mCG age-dynamics over
secretome regulatory regions. Reads committed result CSVs; regenerable.
Out: results/methyl/methyl_figure.{pdf,png} AND figures/figureM1/figure_SM1_methylation.pdf
(the submission path the manuscript bundle uses -- written here so the two cannot drift apart).

*** THE RNA-DIRECTED COORDINATION READING IS WITHDRAWN (manuscript Results / Fig S(M1) legend). ***
The only nominal positive — the astrocyte frozen-hit set (gene-level matched-null p = 0.018,
sign-permutation p = 0.043; results/methyl/gene_level_robustness.csv) — fails four pre-submission
controls: (i) it exists only when the two astrocyte substates are pooled (Astro1 alone p = 0.836,
Astro2 alone p = 0.291); (ii) orienting the same methylation data by another cell type's RNA does
not weaken it (OPC orientation p = 0.019 / 0.025); (iii) non-secretome RNA-significant genes match
it (p = 0.023, n = 155); (iv) the matched diagonal is indistinguishable from a 10,000x
scrambled-orientation null (empirical p = 0.333). The methylome layer is therefore reported as a
directionally NON-SPECIFIC mCG bias, and is NOT counted as support for the accessibility result.
Nothing on this figure may be labelled "coordination", "gene-robust" or "supporting positive".

Panels (honesty kept on-figure):
  A matched-null enrichment (regulatory, per cell type, frozen-hit vs broad); the overlay marks
    which bars are nominally positive at the gene level BEFORE the four controls — a bookkeeping
    label, not a robustness verdict;
  B Astro frozen-hit regions (RNA log2FC vs mCG age-change) — region-level UNADJUSTED reference,
    reported separately from the gene-level p-values because it is a covariate-unadjusted quantity
    over 92 correlated regions (manuscript Results);
  C cross-modal effect-size correlation (cross-cohort, NOT donor-paired; mCG<->ATAC inverse is
    expected by construction). Descriptive: rho is the result, and no p-value is shown because the
    rows are genes/peaks, not independent donors (CLAUDE.md rule 1);
  D genome-wide age-methylation (FDR<0.1; context only; Micro NOT composition-controlled).

Style: scripts/07_figures/_pubstyle.py.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                 "..", "07_figures")))
import _pubstyle as PS

PS.apply()
PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
RES = f"{PROJ}/results/methyl"
NICHE = ["Astro", "Micro", "OPC", "Oligo"]
FROZ_C, ALL_C = PS.UP, PS.SECR


# The two cells that are nominally positive at the gene level, and why neither is a positive
# result. Numbers traced to manuscript/manuscript.md Results (methylome paragraph); the four
# controls are recorded in commit 4e3d27b and audit_log/2026-08-26_manuscript_finalisation/.
CAVEAT = {
    ("Astro", "frozen_hit"): ("\u25b2",
        "Astro frozen-hit (gene-level p_matched 0.018 / p_signperm 0.043) FAILS all four "
        "pre-submission controls \u2014 substate pooling (Astro1 alone 0.836, Astro2 alone 0.291); "
        "orientation by another\n"
        "     cell type's RNA does not weaken it (OPC orientation 0.019 / 0.025); non-secretome "
        "RNA-significant genes match it (0.023, n = 155); indistinguishable from a\n"
        "     10,000\u00d7 scrambled-orientation null (empirical p = 0.333).  \u2192  the "
        "RNA-directed coordination reading is WITHDRAWN; this layer is not regulatory support."),
    ("OPC", "all_secretome"): ("\u25bc",
        "OPC broad set = 0.31 percentage points of methylation over the 51-year span \u2014 below "
        "what we consider biologically interpretable."),
}

enr = pd.read_csv(f"{RES}/secretome_region_methyl_enrichment.csv")
cm = pd.read_csv(f"{RES}/crossmodal_concordance.csv").set_index("celltype")
gl = pd.read_csv(f"{RES}/gene_level_robustness.csv")
glk = {(r.celltype, r["set"]): r for _, r in gl.iterrows()}
SETMAP = {"frozen_hit_linked": "frozen_hit", "all_secretome_linked": "all_secretome"}

fig, ax = plt.subplots(2, 2, figsize=(11.6, 9.8), layout="constrained")

# ---- A: matched-null enrichment (regulatory) ----
a = ax[0, 0]; reg = enr[enr.region_class == "regulatory"]
x = np.arange(len(NICHE)); wdt = 0.38
for k, (s, off, c) in enumerate([("frozen_hit_linked", -wdt / 2, FROZ_C),
                                 ("all_secretome_linked", wdt / 2, ALL_C)]):
    vals = []
    for ct in NICHE:
        r = reg[(reg.celltype == ct) & (reg.set == s)]
        vals.append(r.obs_oriented.values[0] if len(r) else np.nan)
    a.bar(x + off, vals, wdt, color=c, label=("frozen-hit" if k == 0 else "all-secretome"), alpha=0.9)
    for xi, v, ct in zip(x + off, vals, NICHE):
        if np.isnan(v):
            continue
        g = glk.get((ct, SETMAP[s]))
        # NOMINAL at the gene level BEFORE the four pre-submission controls. This is bookkeeping,
        # NOT a robustness verdict: both nominal cells carry a disqualifying caveat (see CAVEAT).
        nom = bool(g is not None and g.p_matched < 0.05 and g.p_signperm < 0.05)
        mark = CAVEAT.get((ct, SETMAP[s]), ("", ""))[0] if nom else ""
        a.text(xi, v + (0.0005 if v >= 0 else -0.0013),
               (f"nominal {mark}" if nom else "region-only"),
               ha="center", va="bottom" if v >= 0 else "top", fontsize=7,
               color=(PS.MUSTARD if nom else PS.GREY), fontweight="bold" if nom else "normal")
a.axhline(0, color=PS.INK, lw=0.6); a.set_xticks(x); a.set_xticklabels(NICHE)
for t, ct in zip(a.get_xticklabels(), NICHE):
    t.set_color(PS.CT[ct]); t.set_fontweight("bold")
a.set_ylabel("region-level oriented mCG effect\n(>0 = hypo near aging-UP / hyper near DOWN)", fontsize=9)
a.set_title("A   Secretome regulatory regions (region-level bars; descriptive)\n"
            "nominal = gene-level matched-null AND sign-perm p < 0.05 BEFORE the four\n"
            "controls (\u25b2\u25bc below) — bookkeeping, not a robustness verdict",
            fontsize=9, loc="left")
a.legend(fontsize=8, loc="upper right"); PS.style_ax(a, grid=True)
_lo, _hi = a.get_ylim(); a.set_ylim(_lo - 0.0013, _hi)

# ---- B: Astro frozen-hit biology ----
b = ax[0, 1]
al = pd.read_csv(f"{RES}/methyl_ageslope_Astro.csv")
rna = (pd.read_csv(f"{PROJ}/results/de/de_GSE268609_per_celltype.csv")
       .query("celltype=='Astro'").set_index("gene")["log2FoldChange"])
froz = set(pd.read_csv(f"{PROJ}/results/validation/frozen_primary_signature.csv")
           .query("celltype=='Astro'").gene)
d = al[al.evaluable & al.region_class.isin(["promoter", "proximal_peak"]) & al.is_secretome
       & al.gene.isin(froz)].copy()
d["rna"] = d.gene.map(rna); d = d.dropna(subset=["rna"]); d["conc"] = (-d.delta_mCG) * np.sign(d.rna) > 0
b.axhline(0, color=PS.GRIDC, lw=0.8); b.axvline(0, color=PS.GRIDC, lw=0.8)
b.scatter(d[d.conc].rna, d[d.conc].delta_mCG, s=20, c=PS.GREEN, alpha=0.85, lw=0, label="concordant")
b.scatter(d[~d.conc].rna, d[~d.conc].delta_mCG, s=20, c=PS.GREY, alpha=0.75, lw=0, label="discordant")
for _, rr in d.reindex(d.delta_mCG.abs().sort_values(ascending=False).index).head(6).iterrows():
    b.annotate(rr.gene, (rr.rna, rr.delta_mCG), xytext=(4, 2), textcoords="offset points", fontsize=7)
b.set_xlabel("RNA log$_2$FC (aging, Astro)", fontsize=9)
b.set_ylabel("Δ mCG (old − young), region-level", fontsize=9)
b.set_title(f"B   Astro frozen-hit regions — region-level UNADJUSTED reference "
            f"(n = {len(d)} regions, {100*d.conc.mean():.0f}% concordant)\n"
            "aging-UP→hypomethylation, aging-DOWN→hypermethylation. Covariate-unadjusted,\n"
            "and regions are correlated within genes — NOT on the same footing as A's p-values.",
            fontsize=9, loc="left")
b.legend(fontsize=8, loc="upper right"); PS.style_ax(b)

# ---- C: cross-modal ----
c = ax[1, 0]; x = np.arange(len(NICHE))
# CLAUDE.md rule 1: the rows here are genes/peaks, not independent donors, so the p-values that
# accompany these correlations are a product of the row count. rho (the effect size) is the result;
# the p-values are deliberately NOT plotted. n is annotated instead.
for k, (col, ncol, lab, cc) in enumerate([("rho_mCG_RNA", "n_genes", "mCG–RNA", "#7b5aa6"),
                                          ("rho_mCG_ATAC", "n_peaks", "mCG–ATAC", "#1f8a8a")]):
    vals = [cm.loc[ct, col] for ct in NICHE]; ns = [int(cm.loc[ct, ncol]) for ct in NICHE]
    c.bar(x + (k - .5) * wdt, vals, wdt, color=cc, label=lab, alpha=0.9)
    for xi, v, n in zip(x + (k - .5) * wdt, vals, ns):
        c.text(xi, v - 0.003, f"{n:,}", ha="center", va="top", fontsize=6, color=PS.SUB)
c.axhline(0, color=PS.INK, lw=0.6); c.set_xticks(x); c.set_xticklabels(NICHE)
for t, ct in zip(c.get_xticklabels(), NICHE):
    t.set_color(PS.CT[ct]); t.set_fontweight("bold")
c.set_ylabel("Spearman ρ", fontsize=9)
c.set_title("C   Cross-modal effect-size correlation (cross-cohort, NOT donor-paired)\n"
            "mCG[GSE299139] vs RNA/ATAC[GSE268609]; mCG↔ATAC inverse expected by construction.\n"
            "Descriptive: ρ is the result, no p shown — rows are genes/peaks (n under each\n"
            f"bar), not donors. max |ρ| = {cm[['rho_mCG_RNA','rho_mCG_ATAC']].abs().max().max():.2f}.",
            fontsize=9, loc="left")
c.legend(fontsize=8); PS.style_ax(c, grid=True)

# ---- D: genome-wide remodeling ----
dd = ax[1, 1]
fdr = []
for ct in NICHE:
    t = pd.read_csv(f"{RES}/methyl_ageslope_{ct}.csv", usecols=["padj"])
    fdr.append(int((t.padj < 0.1).sum()))
dd.bar(NICHE, np.maximum(fdr, 0.5), color=[PS.CT[ct] for ct in NICHE], alpha=0.9)
dd.set_yscale("log"); dd.set_ylabel("regions at BH-FDR < 0.1 (log)", fontsize=9)
for i, v in enumerate(fdr):
    dd.text(i, max(v, 0.5) * 1.15, f"{v:,}", ha="center", fontsize=9)
for t, ct in zip(dd.get_xticklabels(), NICHE):
    t.set_color(PS.CT[ct]); t.set_fontweight("bold")
dd.set_title("D   Genome-wide age-methylation (FDR<0.1; context only)\n"
             "Micro NOT composition-controlled (Micro1 absent in 6 old donors)", fontsize=9.5, loc="left")
PS.style_ax(dd)

fig.suptitle("Figure S-M1  ·  mCG age-dynamics over secretome regulatory regions "
             "(GSE299139 snm3C, 40 donors, hg38)  —  DESCRIPTIVE\n"
             "The RNA-directed coordination reading is WITHDRAWN: the one nominal positive "
             "(Astro frozen-hit) fails four controls;\nthe methylome layer is reported as a "
             "directionally non-specific mCG bias and is NOT counted as regulatory support.",
             fontsize=10, fontweight="bold")
# Figure-level footnote: why neither nominal cell in A is a positive result.
fig.text(0.012, 0.012,
         "\n".join(f"{m}  {txt}" for m, txt in CAVEAT.values()),
         fontsize=6.8, color=PS.SUB, va="bottom", ha="left", linespacing=1.45)
fig.get_layout_engine().set(rect=(0, 0.085, 1, 0.915))

fig.savefig(f"{RES}/methyl_figure.pdf"); fig.savefig(f"{RES}/methyl_figure.png", dpi=200)
# Submission path: figures/figureM1/ is what the manuscript bundle (pdfunite) pulls in. Writing it
# here instead of copying by hand is what stops the figure from drifting behind the manuscript --
# it had done exactly that, still claiming "coordination"/"gene-robust" after the reading was
# withdrawn (codex audit 2026-08-26).
FIGDIR = f"{PROJ}/figures/figureM1"
os.makedirs(FIGDIR, exist_ok=True)
fig.savefig(f"{FIGDIR}/figure_SM1_methylation.pdf")
plt.close(fig)
print(f"wrote {RES}/methyl_figure.pdf + .png and {FIGDIR}/figure_SM1_methylation.pdf "
      f"| FDR<0.1 counts:", dict(zip(NICHE, fdr)))
