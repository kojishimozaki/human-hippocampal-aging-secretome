#!/usr/bin/env python
"""A2 — genome-wide donor-pseudobulk differential accessibility (DA) + secretome-linked-peak
matched-null enrichment. GSE268609 multiome, HA vs YA. Pre-specified: RESOLUTION.md §A2.

Completes the accessibility-matched null that results/regulatory/PHASE_REPORT.md deferred: tests
whether secretome-linked peaks are enriched for DA in the RNA-concordant direction (open near RNA-UP
genes / close near RNA-DOWN) vs a background matched on GC + width + mean-accessibility + distance-to-TSS.

Method (donor = unit; no cell-level p):
  per niche celltype:
    1. donor pseudobulk over ALL 369,393 peaks (sum counts per donor over that celltype's cells)
    2. depth-normalise: CPM = peak_counts / (per-donor total fragments-in-peaks, nCount_Peaks) * 1e6
    3. per-peak vectorised Mann-Whitney (HA donors vs YA donors) -> signed z, lfc, BH-FDR genome-wide
    4. matched-null enrichment of the secretome-linked peak set, oriented by each linked gene's RNA
       direction: observed mean(oriented lfc) vs feature-matched non-secretome-linked peaks given the
       same direction labels (KDTree NN on standardised [gc, log width, log mean-acc, log TSS-dist]).
GO (per celltype): concordant shift (mean oriented lfc>0) AND matched-null p<0.05 AND >=3 linked peaks
  at DA FDR<0.1 with RNA-concordant direction. Else ATAC stays SUPPORTING (genome-wide per-peak family
  is expected underpowered at n=8/9).

seed 42. Out: results/de/atac_DA_peaks_<celltype>.csv, results/de/atac_linkedpeak_enrichment.csv.
"""
import os, gzip
import numpy as np, pandas as pd, scipy.io, scipy.sparse as sp
from scipy.stats import norm, rankdata
from scipy.spatial import cKDTree

SEED = 42
PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
PD = f"{PROJ}/processed/per_dataset"
CL = {"Astrocytes": "Astro", "Microglia": "Micro", "mOli": "Oligo", "OPCs": "OPC", "Endothelial": "Endo"}
NICHE = ["Astro", "Micro", "Oligo", "OPC", "Endo"]
PROM = 2000
NPERM = 2000
KNN = 50
rng = np.random.default_rng(SEED)


def bh(p):
    p = np.asarray(p, float); n = len(p)
    if n == 0:
        return p
    o = np.argsort(p); adj = np.empty(n)
    adj[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


def vec_mwu(cpm, ha, ya):
    """Vectorised two-sided Mann-Whitney per peak (rows), normal approx with tie correction.
    cpm: peaks x donors (YA+HA only). Returns signed z (sign=HA-YA mean diff), two-sided p, lfc."""
    sel = ha | ya
    sub = cpm[:, sel]
    ha_local = ha[sel]
    n1 = int(ha_local.sum()); n = int(sel.sum()); n2 = n - n1
    R = rankdata(sub, axis=1)                                     # avg ranks (ties handled)
    U1 = R[:, ha_local].sum(1) - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    # vectorised tie term: per row, sum over tie-groups of (t^3 - t). n donors -> <=n groups.
    ss = np.sort(sub, axis=1)
    gid = np.zeros_like(ss, dtype=np.int64)
    gid[:, 1:] = np.cumsum(np.diff(ss, axis=1) != 0, axis=1)
    tie_term = np.zeros(sub.shape[0])
    for g in range(n):
        cgt = (gid == g).sum(1).astype(float)
        tie_term += cgt ** 3 - cgt
    sd = np.sqrt(n1 * n2 / 12.0 * ((n + 1) - tie_term / (n * (n - 1))))
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (U1 - mu) / sd
    p = np.where(np.isfinite(z), 2 * norm.sf(np.abs(z)), 1.0)
    meanHA = cpm[:, ha].mean(1); meanYA = cpm[:, ya].mean(1)
    lfc = np.log2((meanHA + 1.0) / (meanYA + 1.0))
    z_signed = np.where(np.isfinite(z), np.sign(meanHA - meanYA) * np.abs(z), 0.0)
    return z_signed, p, lfc, meanHA, meanYA


# ---- peak coords, gc/width ----
peaks = np.array([l.strip() for l in open(f"{PD}/GSE268609_peaks_names.txt")])
parts = [pp.rsplit("-", 2) for pp in peaks]
pchr = np.array([x[0] for x in parts]); pstart = np.array([int(x[1]) for x in parts]); pend = np.array([int(x[2]) for x in parts])
pmid = (pstart + pend) // 2
gcw = pd.read_csv(f"{PD}/GSE268609_peak_gc.csv")
assert (gcw.peak.values == peaks).all(), "peak order mismatch gc vs names"
gc = gcw.gc.values; width = gcw.width.values.astype(float)

# ---- gene coords (all genes -> TSS; secretome subset -> linking) ----
sec = set(pd.read_csv(f"{PROJ}/refs/secretome_union.csv").gene)
allg_chr, allg_tss = [], []
sec_coord = {}
with gzip.open(f"{PROJ}/raw/GSE268609_lazarov2024/GSE268609_features.tsv.gz", "rt") as f:
    for line in f:
        q = line.rstrip().split("\t")
        if len(q) >= 6 and q[2] == "Gene Expression":
            try:
                ch, s, e = q[3], int(q[4]), int(q[5])
            except ValueError:
                continue
            allg_chr.append(ch); allg_tss.append(s)          # TSS approx = gene start (no strand in features)
            if q[1] in sec:
                sec_coord[q[1]] = (ch, s, e)
allg = pd.DataFrame({"chr": allg_chr, "tss": allg_tss})

# distance to nearest TSS per peak (per chr, searchsorted)
tss_dist = np.full(len(peaks), np.nan)
for ch, sub in allg.groupby("chr"):
    idx = np.where(pchr == ch)[0]
    if not len(idx):
        continue
    ts = np.sort(sub.tss.values)
    pos = pmid[idx]
    j = np.searchsorted(ts, pos)
    j = np.clip(j, 0, len(ts) - 1)
    d1 = np.abs(ts[j] - pos); d0 = np.abs(ts[np.clip(j - 1, 0, len(ts) - 1)] - pos)
    tss_dist[idx] = np.minimum(d0, d1)

# secretome gene -> linked peak indices (gene body +/- 2kb)
pkdf = pd.DataFrame({"chr": pchr, "start": pstart, "end": pend, "idx": np.arange(len(peaks))})
pk_by_chr = {c: g for c, g in pkdf.groupby("chr")}
g2peaks = {}
linked_any = set()
for gene, (ch, gs, ge) in sec_coord.items():
    if ch not in pk_by_chr:
        continue
    sub = pk_by_chr[ch]
    hit = sub[(sub.end >= gs - PROM) & (sub.start <= ge + PROM)].idx.values
    if len(hit):
        g2peaks[gene] = hit
        linked_any.update(hit.tolist())
linked_any = np.array(sorted(linked_any))
print(f"secretome genes linked: {len(g2peaks)}; linked peaks: {len(linked_any)}", flush=True)

# ---- cell annotations (aligned to MTX columns) ----
barcodes = [l.strip() for l in open(f"{PD}/GSE268609_peaks_barcodes.txt")]
m = pd.read_csv(f"{PD}/GSE268609_metadata.csv", index_col=0); m.index = m.index.astype(str)
mr = m.reindex(pd.Index(barcodes, dtype=str))
donor = mr["SampleNumber"].values.astype(object)
group = mr["Group"].values.astype(object)
ct = np.array([CL.get(x, None) for x in mr["Cluster"].values], dtype=object)
ncp = mr["nCount_Peaks"].values.astype(float)

# ---- RNA direction (primary DE) ----
rna = pd.read_csv(f"{PROJ}/results/de/de_GSE268609_per_celltype.csv")
rna = rna[rna.gene.isin(sec)]

# ---- load peak MTX (heavy) ----
print("loading peak MTX (369,393 x 153,530; 636M nnz) ...", flush=True)
X = scipy.io.mmread(f"{PD}/GSE268609_peaks_counts.mtx").tocsc().astype(np.float32)  # CSC: fast column (cell) slicing
print("loaded:", X.shape, flush=True)

enr_rows = []
for celltype in NICHE:
    cells = (ct == celltype) & np.isin(group, ["YA", "HA"])
    if cells.sum() == 0:
        continue
    Xc = X[:, cells]
    dc = donor[cells]; gc_grp = group[cells]; ncp_c = ncp[cells]
    donors = pd.unique(dc)
    code = pd.Categorical(dc, categories=donors).codes
    OH = sp.csr_matrix((np.ones(len(dc)), (code, np.arange(len(dc)))), shape=(len(donors), len(dc)))
    pb = np.asarray((Xc @ OH.T).todense())                       # peaks x donors (raw counts)
    depth = np.asarray(OH @ ncp_c).ravel() + 1.0                 # per-donor fragments-in-peaks
    dgrp = np.array([gc_grp[dc == d][0] for d in donors])
    ha = dgrp == "HA"; ya = dgrp == "YA"
    if ha.sum() < 4 or ya.sum() < 4:
        print(f"  skip {celltype}: YA={ya.sum()} HA={ha.sum()} (<4)"); continue
    cpm = pb / depth * 1e6                                        # peaks x donors, depth-normalised
    mean_acc = cpm.mean(1)
    # test peaks: accessible in >=3 donors and nonzero variance
    detected = (pb > 0).sum(1)
    testable = (detected >= 3) & (cpm.std(1) > 0)
    ti = np.where(testable)[0]
    z, p, lfc, mHA, mYA = vec_mwu(cpm[ti], ha, ya)
    padj = bh(p)
    da = pd.DataFrame({
        "peak": peaks[ti], "chr": pchr[ti], "start": pstart[ti], "end": pend[ti],
        "z": z, "p": p, "padj": padj, "lfc_atac": lfc, "mean_acc": mean_acc[ti],
        "gc": gc[ti], "width": width[ti], "tss_dist": tss_dist[ti],
        "is_sec_linked": np.isin(ti, linked_any), "peak_idx": ti})
    da.to_csv(f"{PROJ}/results/de/atac_DA_peaks_{celltype}.csv", index=False)
    nsig = int((da.padj < 0.1).sum())
    print(f"  {celltype}: tested {len(ti)} peaks; DA FDR<0.1 = {nsig}", flush=True)

    # ---- matched-null enrichment (oriented by RNA direction) ----
    rna_ct = rna[rna.celltype == celltype].set_index("gene")
    # build linked-peak -> (rna lfc sign) for this celltype, two sets: frozen-hit (padj_rna<0.1) and all-secretome
    rows_hit, rows_all = [], []
    pos_in_ti = pd.Series(np.arange(len(ti)), index=ti)
    for gene, pks in g2peaks.items():
        if gene not in rna_ct.index:
            continue
        rlfc = rna_ct.loc[gene, "log2FoldChange"]; rpadj = rna_ct.loc[gene, "padj"]
        s = np.sign(rlfc)
        for pk in pks:
            if pk in pos_in_ti.index:
                loc = pos_in_ti[pk]
                rows_all.append((loc, s))
                if rpadj < 0.1:
                    rows_hit.append((loc, s))

    # features for matching (standardised), defined over tested peaks
    feat = np.column_stack([da.gc.values, np.log1p(da.width.values), np.log1p(da.mean_acc.values), np.log1p(da.tss_dist.values)])
    fok = ~np.isnan(feat).any(1)
    fz = (feat - np.nanmean(feat[fok], 0)) / np.nanstd(feat[fok], 0)
    lfc_arr = da.lfc_atac.values; padj_arr = da.padj.values

    def enrich(rows, label):
        if len(rows) < 5:
            return dict(celltype=celltype, set=label, n_linked=len(rows), obs_mean_oriented=np.nan,
                        null_mean=np.nan, p=np.nan, n_fdr_concordant=0, pass_=False)
        li = np.array([r[0] for r in rows]); sgn = np.array([r[1] for r in rows], float)
        good = fok[li]
        li, sgn = li[good], sgn[good]
        oriented = lfc_arr[li] * sgn
        obs = float(np.mean(oriented))
        pool = np.where(fok & (~da.is_sec_linked.values))[0]       # feature-valid, non-secretome-linked
        tree = cKDTree(fz[pool])
        _, nn = tree.query(fz[li], k=KNN)
        nn_glob = pool[nn]                                          # (n_linked, KNN)
        null = np.empty(NPERM)
        for j in range(NPERM):
            pick = nn_glob[np.arange(len(li)), rng.integers(0, KNN, len(li))]
            null[j] = np.mean(lfc_arr[pick] * sgn)
        pval = (np.sum(null >= obs) + 1) / (NPERM + 1)
        nfdr = int(np.sum((padj_arr[li] < 0.1) & (oriented > 0)))
        passed = bool(obs > 0 and pval < 0.05 and nfdr >= 3)
        return dict(celltype=celltype, set=label, n_linked=len(li), obs_mean_oriented=obs,
                    null_mean=float(null.mean()), p=pval, n_fdr_concordant=nfdr, pass_=passed)

    for rows, lab in [(rows_hit, "frozen_hit_linked"), (rows_all, "all_secretome_linked")]:
        r = enrich(rows, lab); enr_rows.append(r)
        print(f"    [{lab}] n={r['n_linked']} obs_oriented_lfc={r['obs_mean_oriented']:.4f} "
              f"null={r['null_mean']:.4f} p={r['p']:.4g} nFDR_concord={r['n_fdr_concordant']} PASS={r['pass_']}", flush=True)

del X
enr = pd.DataFrame(enr_rows).rename(columns={"pass_": "pass"})
enr.to_csv(f"{PROJ}/results/de/atac_linkedpeak_enrichment.csv", index=False)
print("\n===== A2 enrichment summary =====")
print(enr.to_string(index=False))
go = enr[(enr["set"] == "frozen_hit_linked")]
overall = bool(go["pass"].any())
print(f"\n>>> A2 OVERALL (frozen-hit-linked, any niche celltype PASS): "
      f"{'GO: strong-ish ATAC' if overall else 'NO-GO: ATAC stays SUPPORTING (report honest set-level result)'}")
print(f"wrote results/de/atac_linkedpeak_enrichment.csv + atac_DA_peaks_<celltype>.csv")
