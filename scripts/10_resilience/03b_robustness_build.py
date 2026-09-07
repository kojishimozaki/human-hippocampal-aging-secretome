#!/usr/bin/env python
# CONSEQUENCE (Phase A) — STEP 03b (part 1/2): build robustness-variant score tables.
# Provenance: VERBATIM from origin session 027c5fdb (transcript line 849, python part); feeds 03b_robustness_fit.R.
# env: bio.  DATA-DEPENDENT -> run from $PROJ.  Out: results/resilience/{_robustness_scores_2stage.tsv,
#   _robustness_scores_3stage.tsv, _robustness_rawpooled.tsv}.
import pandas as pd, numpy as np, anndata as ad, scipy.sparse as sp
PRIMARY=["DiffN","MatN_SGCZ_NTF3","MatN_SGCZ_CHRM3"]
frozen=pd.read_csv("results/resilience/frozen_receiver_receptor_set.csv").receptor.tolist()
edges=pd.read_csv("results/resilience/frozen_ligand_receptor_edges.csv")
meta=pd.read_csv("processed/per_dataset/GSE325391_adultgc_metadata.csv",index_col=0).drop_duplicates('sample').set_index('sample')
runmap=meta['Run'].astype(str)
A=ad.read_h5ad("processed/per_dataset/GSE325391_resilience.h5ad", backed='r'); obs=A.obs
mask=(obs['group'].astype(str).isin(["RES","SAD"]))&(obs['cell_type'].astype(str).isin(PRIMARY))
sub=A[mask.values].to_memory(); X=sub.X; X=X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
vn=pd.Index(sub.var_names.astype(str)); ridx=[vn.get_loc(r) for r in frozen]
tot=np.asarray(X.sum(1)).ravel(); d=sub.obs['donor_id'].astype(str).values; ct=sub.obs['cell_type'].astype(str).values
gmap=sub.obs.drop_duplicates('donor_id').set_index(sub.obs.drop_duplicates('donor_id').donor_id.astype(str))['group'].astype(str)

# per donor x subtype logCPM (>=20)
idx=[]; mat=[]
for don in sorted(set(d)):
    for s in PRIMARY:
        cm=(d==don)&(ct==s)
        if cm.sum()<20: continue
        rc=np.asarray(X[cm][:,ridx].sum(0)).ravel(); lib=tot[cm].sum()
        idx.append((don,s)); mat.append(np.log2(rc/lib*1e6+1))
L=pd.DataFrame(mat,index=pd.MultiIndex.from_tuples(idx,names=["donor","subtype"]),columns=frozen)
Z=(L-L.mean(0))/L.std(0,ddof=1)   # per-receptor z across 39 samples

def two_stage(scoreseries):
    ss=scoreseries.reset_index(); ss.columns=["donor","subtype","score"]
    imm=ss[ss.subtype=="DiffN"][["donor","score"]].rename(columns={"score":"receiver_score"}); imm["maturation_stage"]="immature"
    mat=ss[ss.subtype.isin(["MatN_SGCZ_NTF3","MatN_SGCZ_CHRM3"])].groupby("donor").score.mean().reset_index().rename(columns={"score":"receiver_score"}); mat["maturation_stage"]="mature"
    f=pd.concat([imm,mat],ignore_index=True); f["group"]=f.donor.map(gmap); f["Run"]=f.donor.map(runmap)
    return f

variants={}
variants["primary"]=two_stage(Z.mean(1))
# ligand-balanced: mean over ligands of mean(its frozen receptors' z)
ligrec={lg:[r for r in g if r in frozen] for lg,g in edges.groupby("ligand").receptor.apply(list).items()}
lb=pd.DataFrame({lg:Z[rs].mean(1) for lg,rs in ligrec.items() if rs}).mean(1)
variants["ligand_balanced"]=two_stage(lb)
# drop ECM + growth-factor-connected receptors
ECM_GF=set(["TNC","COL4A2","COL8A1","COL12A1","COL21A1","EFEMP1","CCN2","SERPINE1","SPON1","FAP","KITLG"])
drop_rec=set(edges[edges.ligand.isin(ECM_GF)].receptor) - set(edges[~edges.ligand.isin(ECM_GF)].receptor)
keep=[r for r in frozen if r not in drop_rec]
variants["drop_ECM_GF_specific"]=two_stage(Z[keep].mean(1))
# leave-one-receptor-out -> range of primary estimates (built as many variants)
loo={f"LOO_{r}":two_stage(Z[[c for c in frozen if c!=r]].mean(1)) for r in frozen}
# write long table
rows=[]
for v,f in {**variants,**loo}.items():
    ff=f.copy(); ff["variant"]=v; rows.append(ff)
long=pd.concat(rows,ignore_index=True)[["variant","donor","group","Run","maturation_stage","receiver_score"]]
long.to_csv("results/resilience/_robustness_scores_2stage.tsv",sep="\t",index=False)
# 3-stage table (subtype = stage)
s3=Z.mean(1).reset_index(); s3.columns=["donor","maturation_stage","receiver_score"]; s3["group"]=s3.donor.map(gmap); s3["Run"]=s3.donor.map(runmap)
s3.to_csv("results/resilience/_robustness_scores_3stage.tsv",sep="\t",index=False)
# raw-count pooled mature
rows=[]
for don in sorted(set(d)):
    cmI=(d==don)&(ct=="DiffN"); cmM=(d==don)&(np.isin(ct,["MatN_SGCZ_NTF3","MatN_SGCZ_CHRM3"]))
    for nm,cm in [("immature",cmI),("mature",cmM)]:
        if cm.sum()<20: continue
        rc=np.asarray(X[cm][:,ridx].sum(0)).ravel(); lib=tot[cm].sum()
        rows.append({"donor":don,"maturation_stage":nm,"logcpm":np.log2(rc/lib*1e6+1)})
RP=pd.DataFrame(rows); M=np.vstack(RP.logcpm.values); Zr=(M-M.mean(0))/M.std(0,ddof=1)
RP["receiver_score"]=Zr.mean(1); RP["group"]=RP.donor.map(gmap); RP["Run"]=RP.donor.map(runmap)
RP[["donor","group","Run","maturation_stage","receiver_score"]].to_csv("results/resilience/_robustness_rawpooled.tsv",sep="\t",index=False)
print(f"[A5] built variants: primary, ligand_balanced, drop_ECM_GF_specific({len(drop_rec)} recs dropped), {len(loo)} LOO-receptor, 3-stage, raw-pooled")
