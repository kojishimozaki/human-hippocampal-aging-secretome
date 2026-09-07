#!/usr/bin/env python
"""M1 methylome — cross-modal coordination (SECONDARY, descriptive; plan M1 step 5).

Within each niche celltype, does secretome-region mCG age-change co-vary with (a) RNA log2FC and
(b) ATAC age-change? Framed as COORDINATION, not causation (methylation may be a consequence; and
mCG<->ATAC are intrinsically anti-correlated, so (b) is partly definitional). Donor-level inputs:
  - mCG age-change: results/methyl/methyl_ageslope_<ct>.csv (delta_mCG = old-young, evaluable regions)
  - RNA log2FC: results/de/de_GSE268609_per_celltype.csv (secretome genes)
  - ATAC age-change: results/de/atac_DA_peaks_<ct>.csv (lfc_atac per peak; peak_idx)
Expected sign: mCG vs RNA NEGATIVE (hypo near aging-UP); mCG vs ATAC NEGATIVE (mCG<->accessibility).
Out: results/methyl/crossmodal_concordance.csv
"""
import os, numpy as np, pandas as pd
from scipy.stats import spearmanr

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
RES = f"{PROJ}/results/methyl"
NICHE = ["Astro", "Micro", "OPC", "Oligo"]
sec = set(pd.read_csv(f"{PROJ}/refs/secretome_union.csv").gene)
rna_all = pd.read_csv(f"{PROJ}/results/de/de_GSE268609_per_celltype.csv")

rows = []
for ct in NICHE:
    a = pd.read_csv(f"{RES}/methyl_ageslope_{ct}.csv")
    a = a[a.evaluable & a.is_secretome].copy()
    rna = rna_all[rna_all.celltype == ct].set_index("gene")["log2FoldChange"]

    # (a) gene-level mCG (regulatory = promoter ∪ proximal_peak) vs RNA log2FC
    regu = a[a.region_class.isin(["promoter", "proximal_peak"])]
    gmcg = regu.groupby("gene").delta_mCG.mean()
    g = pd.DataFrame({"dmCG": gmcg}).join(rna.rename("rna")).dropna()
    rho_rna, p_rna = spearmanr(g.dmCG, g.rna) if len(g) >= 10 else (np.nan, np.nan)

    # (b) peak-level mCG vs ATAC lfc (secretome-linked proximal peaks)
    pk = a[a.region_class == "proximal_peak"].copy()
    pk["pidx"] = pk.region_id.str.replace("peak:", "", regex=False).astype(int)
    atac = pd.read_csv(f"{PROJ}/results/de/atac_DA_peaks_{ct}.csv")[["peak_idx", "lfc_atac"]]
    m = pk.merge(atac, left_on="pidx", right_on="peak_idx").dropna(subset=["delta_mCG", "lfc_atac"])
    rho_atac, p_atac = spearmanr(m.delta_mCG, m.lfc_atac) if len(m) >= 10 else (np.nan, np.nan)

    rows.append(dict(celltype=ct, n_genes=len(g), rho_mCG_RNA=round(rho_rna, 4), p_mCG_RNA=p_rna,
                     n_peaks=len(m), rho_mCG_ATAC=round(rho_atac, 4), p_mCG_ATAC=p_atac))

res = pd.DataFrame(rows)
res.to_csv(f"{RES}/crossmodal_concordance.csv", index=False)
pd.set_option("display.width", 200)
print(res.to_string(index=False))
print("\nExpected (coordination): rho_mCG_RNA<0 (hypo near aging-UP), rho_mCG_ATAC<0 (mCG vs accessibility).")
