#!/usr/bin/env python
"""#3 ATAC regulatory layer: do secretome RNA aging changes (YA->HA) have concordant
chromatin-accessibility changes? GSE268609 same-nucleus multiome.

1. link peaks -> secretome genes (peak overlaps gene body +/- 2kb promoter)
2. donor pseudobulk linked-peak accessibility per niche celltype (YA+HA)
3. per gene: summed linked-peak accessibility, CPM-normalised by each donor's TOTAL
   fragments-in-peaks (metadata nCount_Peaks; NOT by the linked-peak subset, which is a biased
   ~0.5% of all peaks) -> Mann-Whitney HA vs YA. The old linked-peak-only normalisation is kept
   as a robustness column (audit v2.4 fix: per-donor ATAC depth spans ~35x, so the denominator
   is not negligible; the direction-concordance result is reported under the proper depth norm).
4. concordance of ATAC direction with RNA log2FC (HA vs YA)

The (slow) 8.9G peak MTX load is cached to a tiny npz of just the linked peaks.
Output: results/de/atac_rna_concordance.csv
"""
import os, gzip
import numpy as np, pandas as pd, scipy.io, scipy.sparse as sp, scipy.stats as sst

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
PD = f"{PROJ}/processed/per_dataset"
sec = set(pd.read_csv(f"{PROJ}/refs/secretome_union.csv").gene)


def bh(p):
    """Benjamini-Hochberg FDR (no statsmodels dependency)."""
    p = np.asarray(p, float); n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    adj = np.empty(n)
    adj[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)

# secretome gene coords from features.tsv (ENSG, symbol, type, chr, start, end)
gcoord = {}
with gzip.open(f"{PROJ}/raw/GSE268609_lazarov2024/GSE268609_features.tsv.gz", "rt") as f:
    for line in f:
        p = line.rstrip().split("\t")
        if len(p) >= 6 and p[2] == "Gene Expression" and p[1] in sec:
            try:
                gcoord[p[1]] = (p[3], int(p[4]), int(p[5]))
            except ValueError:
                pass
print("secretome genes with coords:", len(gcoord), flush=True)

peaks = [l.strip() for l in open(f"{PD}/GSE268609_peaks_names.txt")]
pk = pd.DataFrame([pp.rsplit("-", 2) for pp in peaks], columns=["chr", "start", "end"])
pk["start"] = pk.start.astype(int); pk["end"] = pk.end.astype(int); pk["idx"] = range(len(pk))
pk_by_chr = {c: g for c, g in pk.groupby("chr")}

PROM = 2000
links = []
for gene, (ch, gs, ge) in gcoord.items():
    if ch not in pk_by_chr:
        continue
    sub = pk_by_chr[ch]
    hit = sub[(sub.end >= gs - PROM) & (sub.start <= ge + PROM)]
    links += [(i, gene) for i in hit.idx]
linkdf = pd.DataFrame(links, columns=["peak_idx", "gene"])
linked_peaks = np.array(sorted(linkdf.peak_idx.unique()))
print(f"linked peaks: {len(linked_peaks)} for {linkdf.gene.nunique()} genes", flush=True)

cache = f"{PD}/GSE268609_linkedpeaks.npz"
cache_idx = f"{PD}/GSE268609_linkedpeaks_idx.npy"
if os.path.exists(cache) and os.path.exists(cache_idx) and np.array_equal(np.load(cache_idx), linked_peaks):
    print("loading cached linked-peak matrix...", flush=True)
    Xl = sp.load_npz(cache)
else:
    print("loading peaks MTX (slow)...", flush=True)
    X = scipy.io.mmread(f"{PD}/GSE268609_peaks_counts.mtx").tocsr()  # peaks x cells
    print("peaks x cells:", X.shape, flush=True)
    Xl = X[linked_peaks].tocsr()
    del X
    sp.save_npz(cache, Xl); np.save(cache_idx, linked_peaks)
    print("cached linked-peak matrix", Xl.shape, flush=True)

barcodes = [l.strip() for l in open(f"{PD}/GSE268609_peaks_barcodes.txt")]
m = pd.read_csv(f"{PD}/GSE268609_metadata.csv", index_col=0)
m.index = m.index.astype(str)
mr = m.reindex(pd.Index(barcodes, dtype=str))
CL = {"Astrocytes": "Astro", "Microglia": "Micro", "mOli": "Oligo", "OPCs": "OPC", "Endothelial": "Endo"}
donor = np.asarray(mr["SampleNumber"].values, dtype=object)
group = np.asarray(mr["Group"].values, dtype=object)
ct = np.array([CL.get(x, None) for x in mr["Cluster"].values], dtype=object)
ncp_all = np.asarray(mr["nCount_Peaks"].values, dtype=float)   # per-cell total fragments in ALL peaks (ATAC depth)

keep = np.array([g in ("YA", "HA") for g in group]) & np.array([c is not None for c in ct])
Xl = Xl[:, keep]; donor = donor[keep]; group = group[keep]; ct = ct[keep]; ncp = ncp_all[keep]
print("YA+HA niche cells:", Xl.shape[1], flush=True)

peak_local = pd.Series(range(len(linked_peaks)), index=linked_peaks)
g2lp = {g: [peak_local[i] for i in grp.peak_idx if i in peak_local.index]
        for g, grp in linkdf.groupby("gene")}

rna = pd.read_csv(f"{PROJ}/results/de/de_GSE268609_per_celltype.csv")
rna = rna[rna.gene.isin(sec)]
out = []
for celltype in ["Astro", "Micro", "Oligo", "OPC", "Endo"]:
    cm = np.asarray(ct == celltype, dtype=bool)
    Xc = Xl[:, cm]; dc = donor[cm]; gc = group[cm]
    donors = pd.unique(dc)
    row = pd.Categorical(dc, categories=donors).codes
    OH = sp.csr_matrix((np.ones(len(dc)), (row, np.arange(len(dc)))), shape=(len(donors), len(dc)))
    pb = (Xc @ OH.T).toarray()        # linked_peaks x donor
    ncp_c = ncp[cm]
    depth = np.asarray(OH @ ncp_c).ravel() + 1.0   # per-donor TOTAL fragments-in-peaks (proper ATAC depth)
    tot_linked = pb.sum(0) + 1.0                   # OLD denominator (linked peaks only) -- robustness check
    dgrp = np.array([gc[dc == d][0] for d in donors])
    ya = dgrp == "YA"; ha = dgrp == "HA"
    if ya.sum() < 3 or ha.sum() < 3:
        print(f"  skip {celltype}: YA={ya.sum()} HA={ha.sum()}"); continue
    rna_ct = rna[rna.celltype == celltype].set_index("gene")
    for g, lps in g2lp.items():
        if not lps or g not in rna_ct.index:
            continue
        acc = pb[lps].sum(0) / depth * 1e6          # CPM-in-peaks, all-peak depth (value used downstream)
        acc_lnk = pb[lps].sum(0) / tot_linked * 1e4 # linked-peak-only norm -- robustness only
        if acc.sum() == 0:
            continue
        try:
            _, p = sst.mannwhitneyu(acc[ha], acc[ya], alternative="two-sided")
        except ValueError:
            continue
        lfc_atac = np.log2((acc[ha].mean() + 1) / (acc[ya].mean() + 1))
        lfc_lnk = np.log2((acc_lnk[ha].mean() + 1) / (acc_lnk[ya].mean() + 1))
        lfc_rna = rna_ct.loc[g, "log2FoldChange"]
        out.append({"gene": g, "celltype": celltype, "n_peaks": len(lps),
                    "lfc_atac": lfc_atac, "p_atac": p, "lfc_rna": lfc_rna,
                    "padj_rna": rna_ct.loc[g, "padj"],
                    "concord": np.sign(lfc_atac) == np.sign(lfc_rna),
                    "lfc_atac_linkednorm": lfc_lnk,
                    "concord_linkednorm": np.sign(lfc_lnk) == np.sign(lfc_rna)})

res = pd.DataFrame(out)
# Multiple-testing correction on the ATAC side (audit fix 2026-06-03). Previously only the
# raw Mann-Whitney p_atac was reported, and Fig3 called nominal p_atac<0.1 "ATAC-significant"
# -> overstated. We now carry:
#   padj_atac        = BH-FDR over ALL gene-celltype ATAC tests (the honest genome-wide family)
#   padj_atac_rnasig = BH-FDR computed WITHIN the RNA-significant family only (the defensible,
#                      pre-conditioned subset: ask the ATAC question only where RNA already moved)
if len(res):
    res["padj_atac"] = bh(res.p_atac.values)
    res["padj_atac_rnasig"] = np.nan
    rmask = res.padj_rna < 0.1
    res.loc[rmask, "padj_atac_rnasig"] = bh(res.loc[rmask, "p_atac"].values)
res.to_csv(f"{PROJ}/results/de/atac_rna_concordance.csv", index=False)

print(f"\n=== ATAC-RNA concordance (secretome, HA vs YA) ===")
print("total gene-celltype ATAC tests:", len(res))
if len(res):
    from scipy.stats import binomtest, fisher_exact
    print(f"[a] ATAC BH-FDR<0.1 over all {len(res)} tests: {(res.padj_atac < 0.1).sum()} "
          f"(the genome-wide ATAC family is underpowered -> ~0 survive)")
    sig = res[res.padj_rna < 0.1]
    nconc = int(sig.concord.sum()); n = len(sig)
    # Direction-concordance is 57.6% (NOT 50%) even among RNA-non-significant genes -- there is a
    # global RNA<->accessibility directional coupling. The defensible LEAD statistic is therefore
    # the ENRICHMENT of RNA-significant concordance ABOVE that empirical background, not a binomial
    # against 0.5 (which overstates the effect: p=0.003 vs the honest ~0.06). [audit v2.2/v2.4]
    # (Fig 3 itself is SUPPORTING, not the paper headline -- this is the lead stat *within* Fig 3.)
    bg = res[res.padj_rna >= 0.1]; bg_c = int(bg.concord.sum()); bg_n = len(bg)
    bg_r = bg_c / bg_n if bg_n else float("nan")
    p_half = binomtest(nconc, n, 0.5, alternative="greater").pvalue if n else float("nan")
    p_bg = binomtest(nconc, n, bg_r, alternative="greater").pvalue if n else float("nan")
    orr, p_fish = fisher_exact([[nconc, n - nconc], [bg_c, bg_n - bg_c]], alternative="greater")
    print(f"[b] LEAD stat (Fig 3 is SUPPORTING, not the paper headline): "
          f"RNA-sig {nconc}/{n} = {nconc / n:.0%} concordant")
    print(f"    empirical background concordance: all tests {res.concord.mean():.1%}, RNA-non-sig {bg_r:.1%}")
    print(f"    ENRICHMENT above background: binom vs {bg_r:.3f} p={p_bg:.3g}; "
          f"Fisher OR={orr:.2f} p={p_fish:.3g}   (naive vs 0.5 p={p_half:.3g} = OVERSTATED)")
    atac_fdr = sig[sig.padj_atac_rnasig < 0.1]
    print(f"[c] ATAC BH-FDR<0.1 within the RNA-sig family: {len(atac_fdr)} "
          f"(concordant {int(atac_fdr.concord.sum())}/{len(atac_fdr)})")
    if len(atac_fdr):
        print(atac_fdr.sort_values("padj_atac_rnasig")
              [["gene", "celltype", "lfc_rna", "lfc_atac", "p_atac", "padj_atac_rnasig", "concord"]].to_string(index=False))
    rob = res[res.padj_rna < 0.1]
    print(f"[e] NORMALISATION ROBUSTNESS (audit v2.4): RNA-sig direction-concordance is "
          f"{int(rob.concord.sum())}/{len(rob)}={rob.concord.mean():.0%} under all-peak-depth (depth-norm) "
          f"vs {int(rob.concord_linkednorm.sum())}/{len(rob)}={rob.concord_linkednorm.mean():.0%} under the "
          f"old linked-peak-only norm; per-pair sign agrees "
          f"{int((rob.concord == rob.concord_linkednorm).sum())}/{len(rob)}")
    both = res[(res.padj_rna < 0.1) & (res.p_atac < 0.1)]
    print(f"[d] EXPLORATORY (nominal dual-sig, padj_rna<0.1 & p_atac<0.1): n={len(both)}, "
          f"concordant {int(both.concord.sum())} ({both.concord.mean():.0%})")
    print("\ntop RNA-sig + nominal-ATAC-concordant secretome hits (exploratory labels for Fig3):")
    print(both[both.concord].sort_values("padj_rna").head(20)
          [["gene", "celltype", "lfc_rna", "lfc_atac", "padj_rna", "p_atac"]].to_string(index=False))
