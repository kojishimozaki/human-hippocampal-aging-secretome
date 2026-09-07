#!/usr/bin/env python
"""A2 adversarial controls for the secretome-linked-peak matched-null enrichment (reads the cached
per-peak DA tables from 05 — no MTX reload). Verifies the directional enrichment is DIRECTION-SPECIFIC
and not an artifact of linked peaks merely having larger |effect|.

For each niche celltype and each linked set (frozen-hit / all-secretome):
  obs            = mean(lfc_atac * sign(linked gene's RNA lfc))           [concordant if >0]
  p_matched      = vs GC+width+acc+TSS-matched non-secretome-linked peaks given the SAME signs (reproduces 05)
  p_signperm     = NEGATIVE CONTROL: shuffle the sign labels among the linked peaks (H0: peak-lfc <-> gene-direction
                   pairing is random). A real directional coupling => obs >> sign-permuted null.
  mag_ratio      = mean(|z| linked) / mean(|z| matched)                   [undirected magnitude, descriptive]
seed 42.  Out: results/de/atac_linkedpeak_enrichment_controls.csv
"""
import os, gzip
import numpy as np, pandas as pd
from scipy.spatial import cKDTree

SEED = 42; rng = np.random.default_rng(SEED)
PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
PD = f"{PROJ}/processed/per_dataset"
NICHE = ["Astro", "Micro", "Oligo", "OPC", "Endo"]; PROM = 2000; NPERM = 2000; KNN = 50

sec = set(pd.read_csv(f"{PROJ}/refs/secretome_union.csv").gene)
peaks = np.array([l.strip() for l in open(f"{PD}/GSE268609_peaks_names.txt")])
parts = [pp.rsplit("-", 2) for pp in peaks]
pchr = np.array([x[0] for x in parts]); pstart = np.array([int(x[1]) for x in parts]); pend = np.array([int(x[2]) for x in parts])
pkdf = pd.DataFrame({"chr": pchr, "start": pstart, "end": pend, "idx": np.arange(len(peaks))})
pk_by_chr = {c: g for c, g in pkdf.groupby("chr")}

sec_coord = {}
with gzip.open(f"{PROJ}/raw/GSE268609_lazarov2024/GSE268609_features.tsv.gz", "rt") as f:
    for line in f:
        q = line.rstrip().split("\t")
        if len(q) >= 6 and q[2] == "Gene Expression" and q[1] in sec:
            try:
                sec_coord[q[1]] = (q[3], int(q[4]), int(q[5]))
            except ValueError:
                pass
g2peaks = {}
for gene, (ch, gs, ge) in sec_coord.items():
    if ch not in pk_by_chr:
        continue
    sub = pk_by_chr[ch]
    hit = sub[(sub.end >= gs - PROM) & (sub.start <= ge + PROM)].idx.values
    if len(hit):
        g2peaks[gene] = hit

rna = pd.read_csv(f"{PROJ}/results/de/de_GSE268609_per_celltype.csv")
rna = rna[rna.gene.isin(sec)]

out = []
for ct in NICHE:
    fp = f"{PROJ}/results/de/atac_DA_peaks_{ct}.csv"
    if not os.path.exists(fp):
        continue
    da = pd.read_csv(fp)
    pos = pd.Series(np.arange(len(da)), index=da.peak_idx.values)
    rna_ct = rna[rna.celltype == ct].set_index("gene")
    feat = np.column_stack([da.gc, np.log1p(da.width), np.log1p(da.mean_acc), np.log1p(da.tss_dist)])
    fok = ~np.isnan(feat).any(1)
    fz = (feat - np.nanmean(feat[fok], 0)) / np.nanstd(feat[fok], 0)
    lfc = da.lfc_atac.values; z = da.z.values
    pool = np.where(fok & (~da.is_sec_linked.values))[0]
    tree = cKDTree(fz[pool])

    rows_hit, rows_all = [], []
    for gene, pks in g2peaks.items():
        if gene not in rna_ct.index:
            continue
        s = np.sign(rna_ct.loc[gene, "log2FoldChange"]); rp = rna_ct.loc[gene, "padj"]
        for pk in pks:
            if pk in pos.index:
                loc = int(pos[pk]); rows_all.append((loc, s))
                if rp < 0.1:
                    rows_hit.append((loc, s))

    for rows, lab in [(rows_hit, "frozen_hit_linked"), (rows_all, "all_secretome_linked")]:
        if len(rows) < 5:
            out.append(dict(celltype=ct, set=lab, n=len(rows), obs=np.nan, p_matched=np.nan,
                            p_signperm=np.nan, mag_ratio=np.nan)); continue
        li = np.array([r[0] for r in rows]); sgn = np.array([r[1] for r in rows], float)
        good = fok[li]; li, sgn = li[good], sgn[good]
        oriented = lfc[li] * sgn; obs = float(np.mean(oriented))
        # matched-feature null (same signs, non-linked peers)
        _, nn = tree.query(fz[li], k=KNN); nn_glob = pool[nn]
        nm = np.array([np.mean(lfc[nn_glob[np.arange(len(li)), rng.integers(0, KNN, len(li))]] * sgn) for _ in range(NPERM)])
        p_matched = (np.sum(nm >= obs) + 1) / (NPERM + 1)
        # sign-permutation negative control (shuffle signs among linked peaks)
        sp = np.array([np.mean(lfc[li] * rng.permutation(sgn)) for _ in range(NPERM)])
        p_signperm = (np.sum(sp >= obs) + 1) / (NPERM + 1)
        # undirected magnitude vs matched
        matched_once = nn_glob[np.arange(len(li)), rng.integers(0, KNN, len(li))]
        mag_ratio = float(np.mean(np.abs(z[li])) / (np.mean(np.abs(z[matched_once])) + 1e-9))
        out.append(dict(celltype=ct, set=lab, n=len(li), obs=round(obs, 4),
                        p_matched=round(p_matched, 4), p_signperm=round(p_signperm, 4), mag_ratio=round(mag_ratio, 3)))

res = pd.DataFrame(out)
res.to_csv(f"{PROJ}/results/de/atac_linkedpeak_enrichment_controls.csv", index=False)
print(res.to_string(index=False))
print("\nInterpretation: a REAL directional coupling => obs>0 with small p_matched AND small p_signperm.")
print("If p_signperm is large while p_matched is small, the signal is magnitude (not direction) -> artifact.")
