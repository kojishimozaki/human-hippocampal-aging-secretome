#!/usr/bin/env python
"""Move #4: external brain & CSF proteome validation of the niche aging secretome.
 - Wingo 2019 brain proteome: cognitive-trajectory beta (>0 stability, <0 decline).
 - Nat Aging 2025 CSF proteome: age estimate (>0 up with age).

Tests, for BOTH directions of the 268609 niche-secretome DE, whether the
external protein-level effect sizes are concordant:
  aging-UP   genes (log2FC>0, padj<thr)  -> Wingo expect beta<0 (decline);  CSF expect >0 (up with age)
  aging-DOWN genes (log2FC<0, padj<thr)  -> Wingo expect beta>0 (stability); CSF expect <0 (down with age)
Direction tested by one-sided Mann-Whitney of the mapped external effect sizes
vs the whole-proteome background.

Threshold (added 2026-06-05, audit response): the gene SET tested is defined by a padj
cut on the niche secretome DE. We emit BOTH thresholds so the dependence is explicit:
  padj<0.1 = PRIMARY = the frozen 60-hit signature (Results §2; frozen_primary_signature.csv)
  padj<0.2 = EXTENDED = a broader secretome-DE set (more mapped genes, more power, NOT frozen)
The manuscript's protein-level claims are anchored to the PRIMARY (frozen) signature, with
EXTENDED reported as a sensitivity. This matters for CSF: the CSF UP signal is significant
only in the EXTENDED set (p~0.001) and is a non-significant directional trend in the PRIMARY
frozen signature (n=25, 60% up, p~0.21); Wingo survives both (p=0.013 -> 0.023).

Design note (added 2026-06-03, audit response): the original script tested the
aging-UP direction only, but the manuscript claimed a bidirectional Wingo result
("aging-DOWN secretome tracks cognitive stability, median beta=+0.94"). That value
was the median of just 3 astrocyte DOWN genes and is NOT reproducible at the n>=4
inclusion threshold used for every other cell type. We therefore test BOTH
directions here and report the honest result. Two cell-type groupings are emitted:
  - per niche cell type (Astro/Micro/Oligo/OPC/Endo) when >=4 mapped genes, plus
    POOLED = concatenation of those qualifying cell types (UP rows reproduce the
    prior output exactly);
  - POOLED_ALL = all unique mapped secretome DE genes across the 5 niche cell types
    with NO per-cell-type gate, so the sparse aging-DOWN arm is still testable as an
    aggregate (this is what the bidirectional manuscript claim must rest on).

Saves results/validation/proteome_validation.csv
  cols: ref, direction, celltype, padj_thresh, n, median_effect, expected_dir_frac, p
"""
import os
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu

P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
NICHE = ["Astro", "Micro", "Oligo", "OPC", "Endo"]
PADJ_THRESHOLDS = [0.1, 0.2]    # 0.1 = PRIMARY frozen 60-hit signature; 0.2 = EXTENDED set (see header)

sec = set(pd.read_csv(f"{P}/refs/secretome_union.csv").gene)
de = pd.read_csv(f"{P}/results/de/de_GSE268609_per_celltype.csv")

wb = pd.read_csv(f"{P}/refs/wingo_cogtraj_beta.csv")
wb = wb.set_index(wb.columns[0])[wb.columns[1]]                       # gene -> cognitive-trajectory beta
csf = pd.read_csv(f"{P}/refs/csf_aging_estimate.csv")
csf = csf.set_index(csf.columns[0])[csf.columns[1]]                   # gene -> CSF age estimate

# (ref vector, alt for UP, alt for DOWN, note). For Wingo: UP genes expected beta<0
# ("less"), DOWN genes expected beta>0 ("greater"). For CSF: UP genes expected >0
# ("greater"), DOWN genes expected <0 ("less").
refs = {"Wingo_brain_cognition": (wb, "less", "greater", "beta<0=decline / beta>0=stability"),
        "CSF_aging":             (csf, "greater", "less", ">0=up-with-age / <0=down-with-age")}


def hits(ct, sign, padj):
    """secretome DE genes in cell type `ct` with the requested sign (1=UP, -1=DOWN), padj < `padj`."""
    d = de[(de.celltype == ct) & de.padj.notna() & de.gene.isin(sec)]
    mask = (d.log2FoldChange > 0) if sign > 0 else (d.log2FoldChange < 0)
    return d[(d.padj < padj) & mask].gene


def summarize(vec, vals, alt):
    frac = float((vals < 0).mean()) if alt == "less" else float((vals > 0).mean())
    p = mannwhitneyu(vals, vec.dropna(), alternative=alt)[1]
    return dict(n=len(vals), median_effect=round(float(vals.median()), 4),
                expected_dir_frac=round(frac, 2), p=p)


rows = []
for PADJ in PADJ_THRESHOLDS:
    for name, (vec, alt_up, alt_down, note) in refs.items():
        bg = vec.dropna()
        for direction, sign, alt in [("UP", 1, alt_up), ("DOWN", -1, alt_down)]:
            pooled = []
            for ct in NICHE:
                e = vec.reindex(hits(ct, sign, PADJ)).dropna()
                if len(e) < 4:
                    continue
                pooled += e.tolist()
                rows.append(dict(ref=name, direction=direction, celltype=ct,
                                 padj_thresh=PADJ, **summarize(vec, e, alt)))
            if pooled:                                               # POOLED of qualifying (>=4) cell types
                rows.append(dict(ref=name, direction=direction, celltype="POOLED",
                                 padj_thresh=PADJ, **summarize(vec, pd.Series(pooled), alt)))
            # POOLED_ALL: every unique mapped gene across the 5 niche cell types, no gate
            allg = pd.Index(sum((list(hits(ct, sign, PADJ)) for ct in NICHE), []))
            ea = vec.reindex(allg).dropna()
            ea = ea[~ea.index.duplicated()]
            if len(ea):
                rows.append(dict(ref=name, direction=direction, celltype="POOLED_ALL",
                                 padj_thresh=PADJ, **summarize(vec, ea, alt)))

res = pd.DataFrame(rows)[["ref", "direction", "celltype", "padj_thresh",
                          "n", "median_effect", "expected_dir_frac", "p"]]
os.makedirs(f"{P}/results/validation", exist_ok=True)
res.to_csv(f"{P}/results/validation/proteome_validation.csv", index=False)
print(res.to_string(index=False))
print("\nKEY: Wingo expected_dir = frac toward the concordant sign "
      "(UP->decline beta<0, DOWN->stability beta>0); CSF analogous (UP->up-with-age, DOWN->down).")
print("saved results/validation/proteome_validation.csv")
