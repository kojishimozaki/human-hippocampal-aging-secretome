#!/usr/bin/env python
"""Fig 2 supplement S3 — CSF aging proteome overlap BY DIRECTION (threshold-dependent, honest).

The aging-UP (SASP) niche secretome trends toward higher CSF age-estimates than the whole-proteome
background, but this is THRESHOLD-DEPENDENT and must be read with Results §3 and
`results/validation/proteome_validation.csv`:
  - frozen padj<0.1 60-hit signature (PRIMARY): aging-UP is a NON-significant directional trend
    (pooled n=25, 60% CSF-concordant, p = 0.21);
  - extended padj<0.2 secretome-DE set: significant (pooled n=42, 69% up, p = 0.001).
The aging-DOWN arm is NOT mirrored at either threshold. So CSF is SUGGESTIVE for the up/SASP arm,
not primary protein-level support. We emit BOTH thresholds so the dependence is explicit and the
`make figures` output cannot silently revert to the old over-claim. (Audit fix codex v2.6, 2026-06-05.)

Computation is unchanged from the audited pipeline; only the rendering is publication-finalised
(shared style in _pubstyle.py). Saves results/validation/csf_direction.csv (+ padj_thresh column)
and figures/figure2/supp_S3_csf_direction.pdf.
"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _pubstyle as PS

PS.apply()
np.random.seed(42)                                                   # reproducible jitter
P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
csf = pd.read_csv(f"{P}/refs/csf_aging_estimate.csv"); csf = csf.set_index(csf.columns[0])[csf.columns[1]]
sec = set(pd.read_csv(f"{P}/refs/secretome_union.csv").gene)
de = pd.read_csv(f"{P}/results/de/de_GSE268609_per_celltype.csv")
bg = csf.dropna()
NICHE = ["Astro", "Micro", "Oligo", "OPC", "Endo"]
THRESHOLDS = [0.1, 0.2]          # 0.1 = frozen 60-hit signature (PRIMARY); 0.2 = extended set


def ct_hits(ct, thr, sign):
    d = de[(de.celltype == ct) & de.padj.notna() & de.gene.isin(sec) & (de.padj < thr)]
    return csf.reindex(d[(d.log2FoldChange > 0) if sign > 0 else (d.log2FoldChange < 0)].gene).dropna()


def pooled_unique(thr, sign):
    """unique CSF-mapped secretome-DE genes of the given sign across niche cell types (POOLED_ALL,
    deduplicated — matches proteome_validation.csv so the % / p agree with the main text)."""
    genes = sum((list(ct_hits(ct, thr, sign).index) for ct in NICHE), [])
    s = csf.reindex(pd.Index(genes)).dropna()
    return s[~s.index.duplicated()]


rows, pool = [], {}
for thr in THRESHOLDS:
    up, dn = pooled_unique(thr, +1), pooled_unique(thr, -1)
    pool[thr] = (up, dn)
    for ct in NICHE:                                                 # per-cell-type (descriptive)
        u, v = ct_hits(ct, thr, +1), ct_hits(ct, thr, -1)
        rows.append(dict(padj_thresh=thr, celltype=ct,
                         n_up=len(u), up_medianCSF=round(float(u.median()), 4) if len(u) else np.nan,
                         pct_CSFup=round(100 * float((u > 0).mean()), 0) if len(u) else np.nan,
                         n_dn=len(v), dn_medianCSF=round(float(v.median()), 4) if len(v) else np.nan,
                         pct_CSFdn=round(100 * float((v < 0).mean()), 0) if len(v) else np.nan))
    rows.append(dict(padj_thresh=thr, celltype="POOLED_ALL",        # deduplicated pooled (matches proteome_validation.csv POOLED_ALL + main text)
                     n_up=len(up), up_medianCSF=round(float(up.median()), 4) if len(up) else np.nan,
                     pct_CSFup=round(100 * float((up > 0).mean()), 0) if len(up) else np.nan,
                     n_dn=len(dn), dn_medianCSF=round(float(dn.median()), 4) if len(dn) else np.nan,
                     pct_CSFdn=round(100 * float((dn < 0).mean()), 0) if len(dn) else np.nan))
res = pd.DataFrame(rows)
res.to_csv(f"{P}/results/validation/csf_direction.csv", index=False)


def pooled_p_up(thr):
    up = pool[thr][0]
    return mannwhitneyu(up, bg, alternative="greater")[1] if len(up) else np.nan


p_up_01, p_up_02 = pooled_p_up(0.1), pooled_p_up(0.2)
up01, dn01 = pool[0.1]
p_dn_01 = mannwhitneyu(dn01, bg, alternative="less")[1] if len(dn01) else np.nan

# ---- figure: lead with the FROZEN (primary) signature; threshold-dependence is the message ----
fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6), layout="constrained")

ax = axes[0]                                                         # boxplot at frozen padj<0.1
data = [bg.values, dn01.values, up01.values]
labs = [f"background\n(n={len(bg)})", f"aging-DOWN\n(n={len(dn01)})", f"aging-UP\n(n={len(up01)})"]
cols = [PS.GREY, PS.DN, PS.UP]
bp = ax.boxplot(data, tick_labels=labs, showfliers=False, widths=0.6,
                medianprops=dict(color=PS.INK, lw=1.2), patch_artist=True)
for patch, c in zip(bp["boxes"], cols):
    patch.set_facecolor(c); patch.set_alpha(0.22); patch.set_edgecolor(c)
for i, v in enumerate(data):
    ax.scatter(np.full(len(v), i + 1) + np.random.uniform(-0.12, 0.12, len(v)), v,
               s=9, c=cols[i], alpha=0.6, lw=0)
ax.axhline(0, c=PS.SUB, lw=0.7, ls=(0, (4, 3)))
ax.set_ylabel("CSF age estimate ( >0 = up with age )"); ax.set_ylim(-0.02, 0.025)
ax.set_title(f"Frozen padj<0.1 signature (primary)\naging-UP vs bg: {PS.fmt_p(p_up_01)} (n={len(up01)}, n.s.)",
             fontsize=9)
PS.style_ax(ax, grid=True)

ax = axes[1]                                                         # threshold sensitivity of pooled UP
pcts = [100 * float((up01 > 0).mean()), 100 * float((pool[0.2][0] > 0).mean())]
ns = [len(up01), len(pool[0.2][0])]
ps = [p_up_01, p_up_02]
xb = np.arange(2)
bars = ax.bar(xb, pcts, 0.55, color=[PS.GREY, PS.UP])
for b, al in zip(bars, [0.55, 0.92]):
    b.set_alpha(al)
ax.axhline(50, c=PS.SUB, ls=(0, (4, 3)), lw=0.7)
ax.set_xticks(xb); ax.set_xticklabels(["frozen\npadj<0.1\n(primary)", "extended\npadj<0.2"], fontsize=8)
ax.set_ylabel("% aging-UP genes concordant (CSF ↑ with age)"); ax.set_ylim(0, 100)
for i in range(2):
    sig = "n.s." if ps[i] >= 0.05 else "sig."
    ax.text(i, pcts[i] + 2, f"{pcts[i]:.0f}%\nn={ns[i]}\n{PS.fmt_p(ps[i])} ({sig})",
            ha="center", va="bottom", fontsize=7.5)
ax.set_title("aging-UP CSF concordance is threshold-dependent", fontsize=9)
PS.style_ax(ax, grid=True)

fig.suptitle("Fig S3  ·  CSF concordance with the aging-UP/SASP arm is threshold-dependent — "
             "suggestive, not primary protein-level support",
             fontsize=9.5, fontweight="bold")
PS.save(fig, f"{P}/figures/figure2/supp_S3_csf_direction.pdf")
print(res.to_string(index=False))
print(f"\nPOOLED aging-UP CSF-concordant: frozen padj<0.1 = {pcts[0]:.0f}% (n={ns[0]}, p={p_up_01:.2g}); "
      f"extended padj<0.2 = {pcts[1]:.0f}% (n={ns[1]}, p={p_up_02:.2g})")
print("saved figures/figure2/supp_S3_csf_direction.pdf + results/validation/csf_direction.csv")
