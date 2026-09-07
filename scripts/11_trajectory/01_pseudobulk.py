#!/usr/bin/env python
"""Phase B (TRAJECTORY) B1 — fresh donor x celltype pseudobulk from the GSE268609 anchor.
Pre-reg v1 + amendment v1.3. YA/HA/AD donors (>=10 cells) for the single R DESeq2 fit; gene filter >=10 total
count (pre-registered). Also a frozen-hit log2-CPM matrix for ALL donors incl. MCI (module context only).
env: bio. No outcome computed here (counts only).
"""
import os, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
PROJ=os.environ.get("PROJ",os.getcwd()); NICHE=["Astro","Micro","Endo","OPC","Oligo"]
OUT=f"{PROJ}/results/trajectory"; os.makedirs(OUT,exist_ok=True)
STEP0={"Astro":(8,9,10),"Micro":(7,9,10),"Endo":(7,9,9),"OPC":(8,9,10),"Oligo":(8,9,10)}
sig=pd.read_csv(f"{PROJ}/results/validation/frozen_primary_signature.csv")
A=ad.read_h5ad(f"{PROJ}/processed/per_dataset/GSE268609_anchor.h5ad",backed='r'); obs=A.obs
ccol="celltype_l1"; dcol="donor_id"; gcol="Group"
def arm(x):
    try: return "whole" if int(str(x).split("_")[0])<=14 else "DG"
    except: return "NA"
for ct in NICHE:
    cells=obs.index[(obs[ccol].astype(str)==ct)&(obs[gcol].astype(str).isin(["YA","HA","AD","MCI"]))]
    sub=A[cells].to_memory(); X=sub.X; X=X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    don=sub.obs[dcol].astype(str).values; grp=sub.obs[gcol].astype(str).values
    alld=sorted(set(don)); keep=[d for d in alld if (don==d).sum()>=10]
    pd.DataFrame({"donor":alld,"n_cells":[(don==d).sum() for d in alld],
                  "group":[grp[don==d][0] for d in alld],
                  "included":[d in keep for d in alld],"reason":["" if d in keep else "<10 cells" for d in alld]}
                 ).to_csv(f"{OUT}/inclusion_log_{ct}.tsv",sep="\t",index=False)
    # pseudobulk sum per donor
    pb=pd.DataFrame(np.vstack([np.asarray(X[don==d].sum(0)).ravel() for d in keep]),
                    index=keep,columns=sub.var_names.astype(str)).astype(int)
    g=pd.Series({d:grp[don==d][0] for d in keep})
    # DE pseudobulk: YA/HA/AD only, gene filter >=10 total count
    de_d=[d for d in keep if g[d] in ("YA","HA","AD")]
    cnt=tuple(int((g[de_d]==x).sum()) for x in ["YA","HA","AD"])
    assert cnt==STEP0[ct], f"{ct} donor inclusion {cnt} != Step-0 {STEP0[ct]}"
    pde=pb.loc[de_d]; pde=pde.loc[:,pde.sum(0)>=10]
    pde.T.to_csv(f"{OUT}/pseudobulk_counts_{ct}.tsv.gz",sep="\t")          # genes x donors
    pd.DataFrame({"grp":g[de_d].values,"arm":[arm(d) for d in de_d]},index=de_d).to_csv(f"{OUT}/pseudobulk_metadata_{ct}.tsv",sep="\t")
    # module expression: frozen hits of this celltype, log2 CPM, ALL donors incl MCI
    fz=[gn for gn in sig[sig.celltype==ct].gene if gn in pb.columns]
    cpm=pb[fz].div(pb.sum(1),axis=0)*1e6; lcpm=np.log2(cpm+1)
    lcpm["grp"]=g.values; lcpm["arm"]=[arm(d) for d in keep]
    lcpm.to_csv(f"{OUT}/module_expression_{ct}.tsv",sep="\t")
    print(f"  {ct}: DE donors {cnt} ({len(de_d)}), genes {pde.shape[1]}; module hits {len(fz)}; +MCI donors {(g=='MCI').sum()}")
print("[B1] pseudobulk + module-expression written (YA/HA/AD for DE; +MCI for module context).")
