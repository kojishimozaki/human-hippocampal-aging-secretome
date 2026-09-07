#!/usr/bin/env python
"""Figure 6 — cognitive-resilience context (GSE325391 granule lineage).

A: # significant DE genes per contrast (2 for RES vs CTRL, 37 for RES vs SAD; counts, not equivalence)
B: RES-vs-SAD resilience signature (top genes; up-in-resilient red, down blue)
C: age by group (RES is older than SAD; age is not adjusted for)

Each panel writes its legend-referenced filename
(panel_A_de_counts.pdf / panel_B_resilience_signature.pdf / panel_C_age_by_group.pdf)
and assembles into composite figure6.pdf. Counts/genes/ages are computed from the
committed DE CSVs + GEO metadata (unchanged); only the rendering is finalised.
"""
import os, sys, re
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _pubstyle as PS

PS.apply()
P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
FIG = f"{P}/figures/figure6"; os.makedirs(FIG, exist_ok=True)

# group colours (resilient = preserved/green; disease = warm; tie to the rest of the set)
GCOL = {"CTRL": PS.GREY, "RES": PS.GREEN, "MAD": PS.MUSTARD, "SAD": PS.UP}

# ---- A: DE counts ------------------------------------------------------------
PAIRS = [("RESvCTRL", "RES vs CTRL", PS.GREEN),
         ("RESvSAD", "RES vs SAD", PS.MUSTARD),
         ("SADvCTRL", "SAD vs CTRL", PS.UP)]
cnt = {}
for f, lab, _ in PAIRS:
    d = pd.read_csv(f"{P}/results/de/de_GSE325391_{f}.csv")
    cnt[lab] = int((d.padj < 0.1).sum())


def panel_counts(cont):
    ax = cont.subplots()
    labs = [l for _, l, _ in PAIRS]; vals = [cnt[l] for l in labs]; cols = [c for _, _, c in PAIRS]
    bars = ax.bar(range(len(labs)), vals, color=cols, width=0.62, alpha=0.92)
    for i, v in enumerate(vals):
        ax.text(i, v + max(vals) * 0.02 + 0.3, str(v), ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs, fontsize=8)
    ax.set_ylabel("DE genes (adj. $p$ < 0.1)")
    ax.set_title("Differentially expressed genes per contrast")
    PS.style_ax(ax, grid=True)
    return ax


# ---- B: RES vs SAD signature -------------------------------------------------
def panel_signature(cont):
    rs = pd.read_csv(f"{P}/results/de/de_GSE325391_RESvSAD.csv")
    sig = rs[rs.padj < 0.1]
    top = sig.reindex(sig.log2FoldChange.abs().sort_values(ascending=False).index).head(20)
    ax = cont.subplots()
    colors = [PS.UP if x > 0 else PS.DN for x in top.log2FoldChange]
    ax.barh(range(len(top)), top.log2FoldChange, color=colors, alpha=0.92)
    ax.set_yticks(range(len(top))); ax.set_yticklabels(top.gene, fontsize=7.5)
    ax.invert_yaxis()
    ax.axvline(0, c=PS.INK, lw=0.6)
    ax.set_xlabel("log$_2$FC (RES $-$ SAD)")
    ax.set_title("Resilience signature\n(up in resilient = red)", fontsize=9.5)
    ax.tick_params(length=3)
    return ax


# ---- C: age by group ---------------------------------------------------------
def panel_age(cont):
    md = pd.read_csv(f"{P}/refs/metadata_raw.csv"); s = md[md.dataset == "GSE325391"]

    def kv(c):
        d = {}
        for p in str(c).split("|"):
            if ":" in p:
                k, v = p.split(":", 1); d[k.strip().lower()] = v.strip()
        return d
    rows = []
    for _, r in s.iterrows():
        d = kv(r.characteristics)
        a = re.search(r"\d+\.?\d*", d.get("age", "") or "")
        rows.append({"age": float(a.group()) if a else np.nan, "group": d.get("group")})
    ag = pd.DataFrame(rows).dropna()
    go = ["CTRL", "RES", "MAD", "SAD"]
    ax = cont.subplots()
    rng = np.random.default_rng(42)
    for i, g in enumerate(go):
        v = ag[ag.group == g].age.values
        ax.scatter(np.full(len(v), i) + rng.uniform(-0.09, 0.09, len(v)), v,
                   s=26, color=GCOL[g], alpha=0.8, lw=0)
        ax.plot([i - 0.22, i + 0.22], [v.mean(), v.mean()], color=PS.INK, lw=1.8)
    ax.set_xticks(range(len(go))); ax.set_xticklabels(go)
    for t, g in zip(ax.get_xticklabels(), go):
        t.set_color(GCOL[g]); t.set_fontweight("bold")
    ax.set_ylabel("Age (years)")
    rmean, smean = ag[ag.group == "RES"].age.mean(), ag[ag.group == "SAD"].age.mean()
    ax.set_title(f"Age by group: RES ({rmean:.0f} y) is older than SAD ({smean:.0f} y); age is not adjusted for", fontsize=8.5)
    PS.style_ax(ax, grid=True)
    return ax, ag


# ---- standalone panels (legend-referenced filenames) -------------------------
f = plt.figure(figsize=(4.0, 3.2), layout="constrained"); panel_counts(f)
PS.save(f, f"{FIG}/panel_A_de_counts.pdf")
f = plt.figure(figsize=(4.2, 5.2), layout="constrained"); panel_signature(f)
PS.save(f, f"{FIG}/panel_B_resilience_signature.pdf")
f = plt.figure(figsize=(4.2, 3.2), layout="constrained"); _, ag = panel_age(f)
PS.save(f, f"{FIG}/panel_C_age_by_group.pdf")

# ---- composite figure5.pdf ---------------------------------------------------
fig = plt.figure(figsize=(7.2, 5.0), layout="constrained")
# Descriptive boundary, not a demonstration that resilience preserves the transcriptome: the
# comparison rests on DE-gene counts, which are not calibrated at this n, and no test was pre-specified.
fig.suptitle("Figure 6  ·  Cognitive-resilience context — granule-lineage transcription by group (descriptive)",
             fontsize=11, fontweight="bold", ha="left", x=0.012)
cols = fig.subfigures(1, 2, width_ratios=[1.0, 0.92])
left = cols[0].subfigures(2, 1, height_ratios=[1.0, 1.0])
panel_counts(left[0]); PS.letter(left[0], "A")
panel_age(left[1]); PS.letter(left[1], "C")
panel_signature(cols[1]); PS.letter(cols[1], "B")
PS.save(fig, f"{FIG}/figure6.pdf")

print("Fig5 DE counts:", cnt)
print("group ages (mean):", ag.groupby("group").age.mean().round(1).to_dict())
print("saved panel_A_de_counts.pdf, panel_B_resilience_signature.pdf, panel_C_age_by_group.pdf, figure6.pdf")
