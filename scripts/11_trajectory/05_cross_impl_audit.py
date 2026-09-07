#!/usr/bin/env python
"""Phase B (TRAJECTORY) B7 — R-DESeq2 vs PyDESeq2 raw-LFC cross-implementation audit (frozen 60, descriptive).
Recomputes YA-reference single-fit raw (MLE) LFC in PyDESeq2 and compares to R DESeq2 results() raw LFC.
Supporting implementation sensitivity ONLY — R DESeq2+ashr is primary; nothing excluded/changed. env: bio.
"""
import os, numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
PROJ=os.environ.get("PROJ",os.getcwd()); OUT=f"{PROJ}/results/trajectory"; NICHE=["Astro","Micro","Endo","OPC","Oligo"]
R=pd.read_csv(f"{OUT}/_frozen60_contrasts_R.tsv",sep="\t").set_index(["gene","celltype"])
rows=[]
for ct in NICHE:
    cts=pd.read_csv(f"{OUT}/pseudobulk_counts_{ct}.tsv.gz",sep="\t",index_col=0)   # genes x donors
    md=pd.read_csv(f"{OUT}/pseudobulk_metadata_{ct}.tsv",sep="\t",index_col=0)
    md.index=md.index.astype(str); cts.columns=cts.columns.astype(str)
    meta=md.copy(); meta["grp"]=pd.Categorical(meta.grp,categories=["YA","HA","AD"])
    dds=DeseqDataSet(counts=cts.T.loc[meta.index].astype(int),metadata=meta,design="~arm + grp",quiet=True); dds.deseq2()
    lf={}
    for nm,con in [("HA_YA",["grp","HA","YA"]),("AD_HA",["grp","AD","HA"]),("AD_YA",["grp","AD","YA"])]:
        st=DeseqStats(dds,contrast=con,quiet=True); st.summary(); lf[nm]=st.results_df["log2FoldChange"]
    fz=R.index[R.index.get_level_values("celltype")==ct]
    for (g,c) in fz:
        if g not in lf["HA_YA"].index: continue
        rows.append(dict(gene=g,celltype=ct,
            R_HA_YA=R.loc[(g,c),"HA_YA_LFC_raw"],py_HA_YA=lf["HA_YA"][g],
            R_AD_HA=R.loc[(g,c),"AD_HA_LFC_raw"],py_AD_HA=lf["AD_HA"][g],
            R_AD_YA=R.loc[(g,c),"AD_YA_LFC_raw"],py_AD_YA=lf["AD_YA"][g]))
d=pd.DataFrame(rows); d.to_csv(f"{OUT}/r_vs_pydeseq2_raw_lfc_audit.tsv",sep="\t",index=False)
print("=== R DESeq2 vs PyDESeq2 raw-LFC (frozen 60; supporting, descriptive) ===")
for nm in ["HA_YA","AD_HA","AD_YA"]:
    r=d[f"R_{nm}"].values; p=d[f"py_{nm}"].values
    print(f"  {nm}: Pearson={pearsonr(r,p)[0]:.4f} Spearman={spearmanr(r,p)[0]:.4f} "
          f"median|Δ|={np.median(np.abs(r-p)):.3f} max|Δ|={np.max(np.abs(r-p)):.3f} "
          f"dir-agree={np.mean(np.sign(r)==np.sign(p))*100:.0f}%")
print(f"\nn frozen hits compared: {len(d)}")
