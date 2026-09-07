#!/usr/bin/env python
"""#2 refinement: validate proximity peak-gene links by MULTIOME co-accessibility.
Same-nucleus 268609: link a peak to a secretome gene only if its accessibility
co-varies with that gene's RNA across SAME-CELLTYPE metacells (Signac-LinkPeaks-style),
instead of naive +/-2 kb proximity. The correlation is computed within each niche celltype
(NOT pooled across cell types -- that would let between-celltype mean differences pose as
regulation; audit fix 2026-06-03), and a Fig3 dual-sig hit counts as validated only if it
carries a validated link in its OWN celltype.
"""
import os
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")   # scanpy/numba + matplotlib reproducibility
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_cache")        # on hosts with a read-only/cacheless HOME
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, scipy.stats as sst, gzip, warnings
warnings.filterwarnings("ignore")
P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project"); PD = f"{P}/processed/per_dataset"

# linked peaks (peaks x cells) cache + recompute peak->gene proximity map
Xl = sp.load_npz(f"{PD}/GSE268609_linkedpeaks.npz").tocsr()
linked_peaks = np.load(f"{PD}/GSE268609_linkedpeaks_idx.npy")
peaks = [l.strip() for l in open(f"{PD}/GSE268609_peaks_names.txt")]
bc = [l.strip() for l in open(f"{PD}/GSE268609_peaks_barcodes.txt")]
sec = set(pd.read_csv(f"{P}/refs/secretome_union.csv").gene)
gcoord = {}
with gzip.open(f"{P}/raw/GSE268609_lazarov2024/GSE268609_features.tsv.gz", "rt") as f:
    for line in f:
        p = line.rstrip().split("\t")
        if len(p) >= 6 and p[2] == "Gene Expression" and p[1] in sec:
            try: gcoord[p[1]] = (p[3], int(p[4]), int(p[5]))
            except ValueError: pass
pk = pd.DataFrame([peaks[i].rsplit("-", 2) for i in linked_peaks], columns=["chr", "s", "e"])
pk["s"] = pk.s.astype(int); pk["e"] = pk.e.astype(int); pk["local"] = range(len(linked_peaks))
links = []
for g, (ch, gs, ge) in gcoord.items():
    sub = pk[(pk.chr == ch) & (pk.e >= gs - 2000) & (pk.s <= ge + 2000)]
    links += [(int(i), g) for i in sub.local]
linkdf = pd.DataFrame(links, columns=["plocal", "gene"])
print("proximity links:", len(linkdf), "genes:", linkdf.gene.nunique())

# anchor RNA (cells x genes, symbols), align to ATAC cells
A = sc.read_h5ad(f"{PD}/GSE268609_anchor.h5ad")
A.X = A.layers["counts"].copy(); sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
bc_idx = pd.Series(range(len(bc)), index=bc)
common = [b for b in A.obs_names if b in bc_idx.index]
A = A[common].copy()
atac_cols = bc_idx.loc[common].values
Xa = Xl[:, atac_cols]                      # peaks x common-cells (accessibility)

# metacells: random bins of ~80 cells WITHIN each celltype; remember each metacell's celltype
rng = np.random.default_rng(42)
ct = A.obs["celltype_l1"].astype(str).values
mc = np.empty(len(ct), dtype=object); mc_ct = {}; k = 0
for c in pd.unique(ct):
    idx = np.where(ct == c)[0]; rng.shuffle(idx)
    for j in range(0, len(idx), 80):
        lab = f"mc{k}"; mc[idx[j:j + 80]] = lab; mc_ct[lab] = c; k += 1
mcs = pd.unique(mc)
oh = sp.csr_matrix((np.ones(len(mc)), (pd.Categorical(mc, categories=mcs).codes, np.arange(len(mc)))),
                   shape=(len(mcs), len(mc)))
sizes = np.asarray(oh.sum(1)).ravel()
peak_mc = (Xa @ oh.T).toarray() / sizes              # peaks x metacell mean access
rna_mc = (oh @ A.X).toarray() / sizes[:, None]       # metacell x genes mean RNA
gidx = pd.Series(range(A.n_vars), index=A.var_names)
mc_ct_vec = np.array([mc_ct[m] for m in mcs])        # celltype of each metacell
NICHE = ["Astro", "Micro", "Oligo", "OPC", "Endo"]
print(f"metacells: {len(mcs)} total; per niche celltype: "
      f"{ {c: int((mc_ct_vec == c).sum()) for c in NICHE} }")

# Validate each proximity link by peak<->RNA co-variation across metacells, computed
# WITHIN each niche celltype (audit fix 2026-06-03). The previous version pooled metacells
# across ALL cell types, so a peak merely "open in astrocytes" correlated with a gene merely
# "high in astrocytes" -> between-celltype mean differences masqueraded as peak-gene
# regulation, inflating the validation rate. Within-celltype Spearman tests genuine
# co-variation among same-type metacells; a link is then validated FOR A SPECIFIC celltype.
MIN_MC = 10  # need enough same-type metacells for a meaningful within-type correlation
val = []
for _, r in linkdf.iterrows():
    g = r.gene
    if g not in gidx.index:
        continue
    pa_all = peak_mc[r.plocal]; ge_all = rna_mc[:, gidx[g]]
    for c in NICHE:
        sel = mc_ct_vec == c
        if sel.sum() < MIN_MC:
            continue
        pa = pa_all[sel]; ge = ge_all[sel]
        if pa.std() == 0 or ge.std() == 0:
            continue
        rho, pv = sst.spearmanr(pa, ge)
        val.append((int(r.plocal), g, c, int(sel.sum()), rho, pv))
vdf = pd.DataFrame(val, columns=["plocal", "gene", "celltype", "n_mc", "rho", "p"])
vdf["validated"] = (vdf.rho > 0.1) & (vdf.p < 0.01)
print(f"within-celltype link tests: {len(vdf)}, validated (rho>0.1,p<0.01): "
      f"{vdf.validated.sum()} ({100 * vdf.validated.mean():.0f}%); median rho={vdf.rho.median():+.3f}")
vdf.to_csv(f"{P}/results/de/coaccess_links.csv", index=False)

# (gene,celltype) pairs with >=1 within-celltype validated link
val_gc = set(map(tuple, vdf[vdf.validated][["gene", "celltype"]].to_numpy()))
val_genes_any = set(vdf[vdf.validated].gene)
print(f"(gene,celltype) pairs with >=1 within-celltype validated link: {len(val_gc)}; "
      f"distinct genes (any celltype): {len(val_genes_any)}")
# how the Fig3 nominal dual-sig hits fare, REQUIRING the validated link in the SAME celltype
ac = pd.read_csv(f"{P}/results/de/atac_rna_concordance.csv")
both = ac[(ac.padj_rna < 0.1) & (ac.p_atac < 0.1)][["gene", "celltype"]]
match_same = [(g, c) for g, c in both.to_numpy() if (g, c) in val_gc]
match_any = [(g, c) for g, c in both.to_numpy() if g in val_genes_any]
print(f"\nFig3 nominal dual-sig (gene,celltype) pairs: {len(both)}")
print(f"  with SAME-celltype validated co-accessible link: {len(match_same)}/{len(both)} "
      f"(e.g. {match_same[:8]})")
print(f"  with a validated link in ANY celltype (weaker bar): {len(match_any)}/{len(both)}")
