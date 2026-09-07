#!/usr/bin/env python
"""[F-3b-002 / G-9] Rebuild the A2 foreground exactly as the published pipeline defines it.

WHY
---
`compute/pass3b_02_selection_circularity_null.py` (audit kit) builds its own foreground and
does NOT reproduce the published A2 observed statistic:

    celltype  published n / obs     pass3b_02 n / obs
    Astro        160 / +0.1834        106 / +0.1639
    Micro         54 / +0.2391         30 / +0.2926
    Oligo         30 / -0.0151         14 / +0.0131
    OPC           19 / +0.2863         55 / +0.0331     <-- the OPC arm of F-3b-002
    Endo          40 / +0.1118         57 / +0.0781

The peer session showed (in F-3b-012) that the missing ingredient is the "linked AND TESTED"
restriction: `scripts/04_atac/05_differential_accessibility.py:L183-L197` keeps a linked peak
only if it survived that cell type's testability filter (`pk in pos_in_ti.index`), and then
drops peaks whose matching features are NaN (`good = fok[li]`).

This script reproduces the published foreground from the committed DA tables, so the heavy
369,393 x 153,530 MTX is never touched. Gate G-9 (pre-registered in
audit_log/2026-08-26_remaining_still_open/PRESPEC.md) requires agreement with the published
`obs` to three decimals in all five cell types before any p-value built on it is used.

INPUTS   results/de/atac_DA_peaks_<celltype>.csv     (per-celltype tested peaks; committed)
         results/de/atac_linkedpeak_enrichment.csv   (published obs / n_linked)
         results/de/de_GSE268609_per_celltype.csv    (RNA log2FC + padj)
         refs/secretome_union.csv
         raw/GSE268609_lazarov2024/GSE268609_features.tsv.gz   (gene coords)
         processed/per_dataset/GSE268609_peaks_names.txt       (peak coords)
OUTPUT   results/audit_reruns/a2_foreground_reproduction.csv
RUNTIME  ~2 min (reads five ~70 MB CSVs).
"""
import gzip
import os

import numpy as np
import pandas as pd

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
PD = f"{PROJ}/processed/per_dataset"
OUT = f"{PROJ}/results/audit_reruns"
os.makedirs(OUT, exist_ok=True)
NICHE = ["Astro", "Micro", "Oligo", "OPC", "Endo"]
PROM = 2000                      # same gene-body +/- 2 kb window as the published script

# ---- peak coordinates, in the deposited row order -------------------------------------
peaks = np.array([l.strip() for l in open(f"{PD}/GSE268609_peaks_names.txt")])
parts = [pp.rsplit("-", 2) for pp in peaks]
pchr = np.array([x[0] for x in parts])
pstart = np.array([int(x[1]) for x in parts])
pend = np.array([int(x[2]) for x in parts])

# ---- secretome gene coordinates -------------------------------------------------------
sec = set(pd.read_csv(f"{PROJ}/refs/secretome_union.csv").iloc[:, 0].astype(str))
sec_coord = {}
with gzip.open(f"{PROJ}/raw/GSE268609_lazarov2024/GSE268609_features.tsv.gz", "rt") as f:
    for line in f:
        q = line.rstrip().split("\t")
        if len(q) >= 6 and q[2] == "Gene Expression" and q[1] in sec:
            try:
                sec_coord[q[1]] = (q[3], int(q[4]), int(q[5]))
            except ValueError:
                continue

# ---- secretome gene -> linked peak indices (identical to the published construction) ---
pkdf = pd.DataFrame({"chr": pchr, "start": pstart, "end": pend, "idx": np.arange(len(peaks))})
pk_by_chr = {c: g for c, g in pkdf.groupby("chr")}
g2peaks = {}
for gene, (ch, gs, ge) in sec_coord.items():
    if ch not in pk_by_chr:
        continue
    sub = pk_by_chr[ch]
    hit = sub[(sub.end >= gs - PROM) & (sub.start <= ge + PROM)].idx.values
    if len(hit):
        g2peaks[gene] = hit
print(f"secretome genes linked: {len(g2peaks)}", flush=True)

rna = pd.read_csv(f"{PROJ}/results/de/de_GSE268609_per_celltype.csv")
pub = pd.read_csv(f"{PROJ}/results/de/atac_linkedpeak_enrichment.csv")

rows = []
for celltype in NICHE:
    da = pd.read_csv(f"{PROJ}/results/de/atac_DA_peaks_{celltype}.csv")
    # `ti` is exactly the tested-peak index set the published script wrote out
    ti = da.peak_idx.values
    pos_in_ti = pd.Series(np.arange(len(ti)), index=ti)
    lfc_arr = da.lfc_atac.values
    padj_arr = da.padj.values
    # feature validity, same four columns and same NaN rule as the published `fok`
    feat = np.column_stack([da.gc.values, np.log1p(da.width.values),
                            np.log1p(da.mean_acc.values), np.log1p(da.tss_dist.values)])
    fok = ~np.isnan(feat).any(1)

    rna_ct = rna[rna.celltype == celltype].set_index("gene")
    rows_hit, rows_all = [], []
    for gene, pks in g2peaks.items():
        if gene not in rna_ct.index:
            continue
        rlfc = rna_ct.loc[gene, "log2FoldChange"]
        rpadj = rna_ct.loc[gene, "padj"]
        s = np.sign(rlfc)
        for pk in pks:
            if pk in pos_in_ti.index:
                loc = pos_in_ti[pk]
                rows_all.append((loc, s))
                if rpadj < 0.1:
                    rows_hit.append((loc, s))

    for rr, label in [(rows_hit, "frozen_hit_linked"), (rows_all, "all_secretome_linked")]:
        if len(rr) < 5:
            continue
        li = np.array([r[0] for r in rr])
        sgn = np.array([r[1] for r in rr], float)
        good = fok[li]
        li, sgn = li[good], sgn[good]
        oriented = lfc_arr[li] * sgn
        obs = float(np.mean(oriented))
        nfdr = int(np.sum((padj_arr[li] < 0.1) & (oriented > 0)))
        p = pub[(pub.celltype == celltype) & (pub["set"] == label)]
        pub_n = int(p.n_linked.iloc[0]) if len(p) else -1
        pub_obs = float(p.obs_mean_oriented.iloc[0]) if len(p) else np.nan
        rows.append(dict(celltype=celltype, set=label,
                         n_linked=len(li), published_n_linked=pub_n,
                         obs_mean_oriented=round(obs, 6), published_obs=round(pub_obs, 6),
                         n_genes=len({g for g, pk in g2peaks.items()}),  # placeholder, refined below
                         n_fdr_concordant=nfdr,
                         match_n=bool(len(li) == pub_n),
                         match_obs_3dp=bool(abs(obs - pub_obs) < 5e-4)))
        print(f"  {celltype:6s} {label:22s} n={len(li):5d} (published {pub_n:5d})  "
              f"obs={obs:+.4f} (published {pub_obs:+.4f})  "
              f"{'MATCH' if abs(obs - pub_obs) < 5e-4 else 'MISMATCH'}", flush=True)

R = pd.DataFrame(rows).drop(columns=["n_genes"])
R.to_csv(f"{OUT}/a2_foreground_reproduction.csv", index=False)
ok = bool(R.match_n.all() and R.match_obs_3dp.all())
print(f"\nG-9 reproduction gate: {'PASS' if ok else 'FAIL'} "
      f"({int(R.match_obs_3dp.sum())}/{len(R)} obs match to 3 dp, "
      f"{int(R.match_n.sum())}/{len(R)} n match)")
print(f"wrote {OUT}/a2_foreground_reproduction.csv")
