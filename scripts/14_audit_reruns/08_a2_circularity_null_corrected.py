#!/usr/bin/env python
"""[F-3b-002] Selection-circularity null for A2, on the foreground the paper actually used.

WHY
---
`compute/pass3b_02_selection_circularity_null.py` reports p_label_permutation = 0.0299 (Astro),
0.0448 (Micro), 0.5323 (OPC). Those numbers rest on a foreground that does not reproduce the
published A2 statistic -- most severely for OPC, where it finds 55 peaks and obs = +0.0331
against the published 19 peaks and obs = +0.2863. Script 07 showed the missing ingredient is
the "linked AND TESTED (AND feature-valid)" restriction, and reproduced all ten published
values exactly. This script re-runs the circularity null on that corrected foreground.

WHAT THE NULL ASKS
------------------
The A2 foreground is chosen by RNA significance in the SAME donors whose ATAC is then measured.
Any donor-level nuisance shared by both modalities can therefore manufacture concordance. The
null permutes the donor group labels, re-selects the same NUMBER of genes by the permuted RNA
statistic, and recomputes the oriented ATAC mean under those same permuted labels. Selection
intensity is held constant; only the labels move.

Two observed statistics are reported, because they answer different questions:
  obs_published  the canonical statistic (real frozen hits, DESeq2 padj<0.1, published lfc_atac)
  obs_procedure  the same selection PROCEDURE the null uses (top-K by MWU p), under true labels
The null is procedure-matched, so `p_vs_procedure` is the internally consistent test and
`p_vs_published` says whether the published number sits outside the same null.

GATES (pre-registered in audit_log/2026-08-26_remaining_still_open/PRESPEC.md)
  G-9   script 07: foreground reproduces published obs in 5/5 cell types   [already PASS]
  G-9b  ATAC log2FC recomputed here from the linked-peak matrix must match the committed
        DA table's lfc_atac for the same peaks (Spearman >= 0.999 and max |diff| < 1e-6).
        Without this the null would be built on a different quantity from the observed.

INPUTS   processed/per_dataset/GSE268609_linkedpeaks.npz + _idx.npy   (5,888 linked peaks x cells)
         processed/per_dataset/GSE268609_metadata.csv                 (nCount_Peaks, SampleNumber, Group, Cluster)
         processed/per_dataset/GSE268609_anchor.h5ad                  (RNA counts, celltype_l1, donor_id)
         results/de/atac_DA_peaks_<celltype>.csv                      (tested peaks + published lfc_atac)
         results/de/de_GSE268609_per_celltype.csv, refs/secretome_union.csv
OUTPUT   results/audit_reruns/a2_circularity_null_corrected.csv
         results/audit_reruns/a2_circularity_null_corrected_draws.csv
RUNTIME  ~10 min (h5ad load dominates).
"""
import gzip
import os

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import norm, rankdata

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
PD = f"{PROJ}/processed/per_dataset"
OUT = f"{PROJ}/results/audit_reruns"
os.makedirs(OUT, exist_ok=True)
NICHE = ["Astro", "Micro", "Oligo", "OPC", "Endo"]
PROM = 2000
NPERM = int(os.environ.get("NPERM_LAB", "200"))
SEED = 42
rng = np.random.default_rng(SEED)

# ------------------------------------------------------------------ peak coords + links --
peaks = np.array([l.strip() for l in open(f"{PD}/GSE268609_peaks_names.txt")])
parts = [pp.rsplit("-", 2) for pp in peaks]
pchr = np.array([x[0] for x in parts])
pstart = np.array([int(x[1]) for x in parts])
pend = np.array([int(x[2]) for x in parts])

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

# ------------------------------------------------------------------ ATAC linked peaks ----
Xl = sp.load_npz(f"{PD}/GSE268609_linkedpeaks.npz").tocsr()          # linked peaks x cells
lp_idx = np.load(f"{PD}/GSE268609_linkedpeaks_idx.npy")              # -> global peak index
glob2lp = pd.Series(np.arange(len(lp_idx)), index=lp_idx)
# the matrix columns follow peaks_barcodes.txt; reindex the metadata onto that order
# exactly as scripts/04_atac/05_differential_accessibility.py:L132-L137 does.
CL = {"Astrocytes": "Astro", "Microglia": "Micro", "mOli": "Oligo",
      "OPCs": "OPC", "Endothelial": "Endo"}
barcodes = [l.strip() for l in open(f"{PD}/GSE268609_peaks_barcodes.txt")]
meta = pd.read_csv(f"{PD}/GSE268609_metadata.csv", index_col=0)
mr = meta.reindex(pd.Index(barcodes, dtype=str))
assert mr["Group"].notna().all(), "metadata does not cover every matrix barcode"
ncp = mr["nCount_Peaks"].values.astype(float)
donor_all = mr["SampleNumber"].astype(str).values
group_all = mr["Group"].astype(str).values
clus_all = np.array([CL.get(x, None) for x in mr["Cluster"].values], dtype=object)
assert Xl.shape[1] == len(mr), f"linkedpeaks cells {Xl.shape[1]} != metadata {len(mr)}"
print("cells per niche cluster:",
      {c: int((clus_all == c).sum()) for c in NICHE}, flush=True)

# ------------------------------------------------------------------ RNA ------------------
import scanpy as sc  # noqa: E402
A = sc.read_h5ad(f"{PD}/GSE268609_anchor.h5ad")
A = A[A.obs.Group.isin(["YA", "HA"])].copy()
sec_genes = [g for g in g2peaks if g in set(A.var_names.astype(str))]
gpos = {g: i for i, g in enumerate(A.var_names.astype(str))}
XR = A.layers["counts"].tocsc()[:, [gpos[g] for g in sec_genes]].tocsr()
rna_ct_all = A.obs.celltype_l1.astype(str).values
rna_don = A.obs.donor_id.astype(str).values
rna_grp = A.obs.Group.astype(str).values
print(f"RNA cells {A.n_obs}, secretome genes with links {len(sec_genes)}", flush=True)

rna_de = pd.read_csv(f"{PROJ}/results/de/de_GSE268609_per_celltype.csv")


def pseudobulk(Xcells_T, code, ndon):
    """Xcells_T: features x cells (csr). Returns donors x features dense."""
    OH = sp.csr_matrix((np.ones(code.shape[0]), (code, np.arange(code.shape[0]))),
                       shape=(ndon, code.shape[0]))
    return OH, np.asarray((OH @ Xcells_T.T).todense())


def mwu_p(mat, i1, i2):
    """Two-sided MWU normal approx, features in COLUMNS of `mat` (donors x features)."""
    sub = np.vstack([mat[i1], mat[i2]])
    n1, n2 = len(i1), len(i2)
    n = n1 + n2
    R = rankdata(sub, axis=0)
    U1 = R[:n1].sum(0) - n1 * (n1 + 1) / 2.0
    sd = np.sqrt(n1 * n2 * (n + 1) / 12.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (U1 - n1 * n2 / 2.0) / sd
    return np.where(np.isfinite(z), 2 * norm.sf(np.abs(z)), 1.0)


rows, draws = [], []
for celltype in NICHE:
    da = pd.read_csv(f"{PROJ}/results/de/atac_DA_peaks_{celltype}.csv")
    ti = da.peak_idx.values
    pos_in_ti = pd.Series(np.arange(len(ti)), index=ti)
    lfc_pub = da.lfc_atac.values
    feat = np.column_stack([da.gc.values, np.log1p(da.width.values),
                            np.log1p(da.mean_acc.values), np.log1p(da.tss_dist.values)])
    fok = ~np.isnan(feat).any(1)

    # ---- ATAC pseudobulk over the linked peaks, in this cell type ----------------------
    cm = (clus_all == celltype) & np.isin(group_all, ["YA", "HA"])
    if cm.sum() == 0:
        print(f"skip {celltype}: no cells"); continue
    donors = pd.unique(donor_all[cm])
    code_a = pd.Categorical(donor_all[cm], categories=donors).codes
    OH_a, pb_a = pseudobulk(Xl[:, cm], code_a, len(donors))          # donors x linked peaks
    depth_a = np.asarray(OH_a @ ncp[cm]).ravel() + 1.0
    cpm_a = pb_a / depth_a[:, None] * 1e6
    dgrp = np.array([group_all[cm][donor_all[cm] == d][0] for d in donors])
    i_ha = np.where(dgrp == "HA")[0]
    i_ya = np.where(dgrp == "YA")[0]
    if len(i_ha) < 4 or len(i_ya) < 4:
        print(f"skip {celltype}: YA={len(i_ya)} HA={len(i_ha)} (<4)"); continue

    def atac_lfc(lab):
        h = np.where(lab == "HA")[0]; y = np.where(lab == "YA")[0]
        return np.log2((cpm_a[h].mean(0) + 1.0) / (cpm_a[y].mean(0) + 1.0))   # per linked peak

    # ---- G-9b: does our recomputed lfc match the committed DA table? -------------------
    lfc_here = atac_lfc(dgrp)
    both = [(int(glob2lp[g]), int(pos_in_ti[g])) for g in ti if g in glob2lp.index]
    if len(both) < 50:
        print(f"{celltype}: too few shared peaks for G-9b ({len(both)})"); continue
    a = np.array([lfc_here[i] for i, _ in both])
    b = np.array([lfc_pub[j] for _, j in both])
    g9b_rho = float(pd.Series(a).corr(pd.Series(b), method="spearman"))
    g9b_max = float(np.max(np.abs(a - b)))
    g9b = bool(g9b_rho >= 0.999 and g9b_max < 1e-6)
    print(f"  {celltype:6s} G-9b: n={len(both)} rho={g9b_rho:.6f} max|diff|={g9b_max:.3g} "
          f"-> {'PASS' if g9b else 'FAIL'}", flush=True)

    # ---- RNA pseudobulk in this cell type ---------------------------------------------
    rm = (rna_ct_all == celltype) & np.isin(rna_grp, ["YA", "HA"])
    rdon = rna_don[rm]
    code_r = pd.Categorical(rdon, categories=donors).codes
    keep_r = code_r >= 0
    OH_r, pb_r = pseudobulk(XR[rm][keep_r].T.tocsr(), code_r[keep_r], len(donors))
    cpm_r = pb_r / (pb_r.sum(1, keepdims=True) + 1.0) * 1e6          # donors x genes

    K = int(rna_de[(rna_de.celltype == celltype) & (rna_de.gene.isin(sec))
                   & (rna_de.padj < 0.1)].gene.nunique())
    if K < 1:
        print(f"{celltype}: K=0"); continue

    def stat(lab):
        """top-K by permuted RNA p -> oriented mean of ATAC lfc over linked AND TESTED peaks."""
        h = np.where(lab == "HA")[0]; y = np.where(lab == "YA")[0]
        pr = mwu_p(cpm_r, h, y)
        lr = np.log2((cpm_r[h].mean(0) + 1) / (cpm_r[y].mean(0) + 1))
        la = atac_lfc(lab)
        sel = np.argsort(pr)[:K]
        vals, ng = [], 0
        for j in sel:
            g = sec_genes[j]
            used = False
            for pk in g2peaks[g]:
                if pk in pos_in_ti.index and fok[pos_in_ti[pk]] and pk in glob2lp.index:
                    vals.append(la[int(glob2lp[pk])] * np.sign(lr[j]))
                    used = True
            ng += int(used)
        if not vals:
            return np.nan, 0, 0
        return float(np.mean(vals)), ng, len(vals)

    obs_proc, ng_proc, np_proc = stat(dgrp)

    # canonical observed, straight from script 07's logic (real frozen hits, published lfc)
    rna_ct = rna_de[rna_de.celltype == celltype].set_index("gene")
    rr = []
    for gene, pks in g2peaks.items():
        # NB: `not (padj < 0.1)` and `padj >= 0.1` differ when padj is NaN. Script 07 (and the
        # published pipeline) use `if rpadj < 0.1`, which drops NaN. Keep that behaviour.
        if gene not in rna_ct.index or not (rna_ct.loc[gene, "padj"] < 0.1):
            continue
        s = np.sign(rna_ct.loc[gene, "log2FoldChange"])
        for pk in pks:
            if pk in pos_in_ti.index:
                rr.append((pos_in_ti[pk], s))
    li = np.array([r[0] for r in rr]); sgn = np.array([r[1] for r in rr], float)
    good = fok[li]
    obs_pub = float(np.mean(lfc_pub[li[good]] * sgn[good]))

    null = np.empty(NPERM)
    for j in range(NPERM):
        null[j] = stat(rng.permutation(dgrp))[0]
        draws.append(dict(celltype=celltype, draw=j, obs_star=null[j]))
    null = null[np.isfinite(null)]
    p_proc = (np.sum(null >= obs_proc) + 1) / (len(null) + 1)
    p_pub = (np.sum(null >= obs_pub) + 1) / (len(null) + 1)

    rows.append(dict(celltype=celltype, K_selected=K, g9b_pass=g9b, g9b_rho=round(g9b_rho, 6),
                     g9b_max_abs_diff=g9b_max,
                     n_genes_procedure=ng_proc, n_peaks_procedure=np_proc,
                     n_peaks_published=int(good.sum()),
                     obs_procedure=round(obs_proc, 4), obs_published=round(obs_pub, 4),
                     null_mean=round(float(null.mean()), 4), null_sd=round(float(null.std()), 4),
                     null_p95=round(float(np.quantile(null, 0.95)), 4),
                     p_vs_procedure=round(p_proc, 4), p_vs_published=round(p_pub, 4),
                     n_draws=len(null)))
    print(f"  {celltype:6s} K={K:3d} obs_proc={obs_proc:+.4f} obs_pub={obs_pub:+.4f} | "
          f"null mean={null.mean():+.4f} sd={null.std():.4f} p95={np.quantile(null, 0.95):+.4f} "
          f"-> p_proc={p_proc:.4g} p_pub={p_pub:.4g}", flush=True)

pd.DataFrame(rows).to_csv(f"{OUT}/a2_circularity_null_corrected.csv", index=False)
pd.DataFrame(draws).to_csv(f"{OUT}/a2_circularity_null_corrected_draws.csv", index=False)
print(f"\nwrote {OUT}/a2_circularity_null_corrected.csv")
