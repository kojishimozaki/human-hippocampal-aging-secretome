#!/usr/bin/env python
# CONSEQUENCE (Phase A) — STEP 02: per-donor x maturation-stage receiver-competence score.
# Provenance: recovered VERBATIM from origin session 027c5fdb (transcript line 822); produced the committed
#   donor_stage_receptor_scores.tsv (26 rows; RES 6 / SAD 7) byte-for-byte (agent-verified).
# env: bio.  DATA-DEPENDENT (GSE325391 h5ad + metadata) -> run from $PROJ.
# Out: results/resilience/{donor_stage_receptor_scores.tsv, pseudobulk_cell_counts.csv, pseudobulk_inclusion_log.csv}
import pandas as pd, numpy as np, anndata as ad, scipy.sparse as sp
PRIMARY=["DiffN","MatN_SGCZ_NTF3","MatN_SGCZ_CHRM3"]
frozen=pd.read_csv("results/resilience/frozen_receiver_receptor_set.csv").receptor.tolist()
A=ad.read_h5ad("processed/per_dataset/GSE325391_resilience.h5ad", backed='r')
obs=A.obs
mask=(obs['group'].astype(str).isin(["RES","SAD"])) & (obs['cell_type'].astype(str).isin(PRIMARY))
sub=A[mask.values].to_memory()
X=sub.X; X=X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
vn=pd.Index(sub.var_names.astype(str)); ridx=[vn.get_loc(r) for r in frozen]
tot=np.asarray(X.sum(1)).ravel()
d=sub.obs['donor_id'].astype(str).values; ct=sub.obs['cell_type'].astype(str).values
gmap=sub.obs.drop_duplicates('donor_id').set_index(sub.obs.drop_duplicates('donor_id').donor_id.astype(str))['group'].astype(str)
# Run from metadata
meta=pd.read_csv("processed/per_dataset/GSE325391_adultgc_metadata.csv",index_col=0).drop_duplicates('sample').set_index('sample')
runmap=meta['Run'].astype(str)

# pseudobulk logCPM per donor x original subtype (>=20 cells)
cc=[]; inc=[]; pbrows=[]
for don in sorted(set(d)):
    for s in PRIMARY:
        cm=(d==don)&(ct==s); n=int(cm.sum())
        cc.append({"donor":don,"group":gmap[don],"subtype":s,"n_cells":n})
        if n<20: inc.append({"donor":don,"subtype":s,"n_cells":n,"included":False,"reason":"<20 cells"}); continue
        inc.append({"donor":don,"subtype":s,"n_cells":n,"included":True,"reason":""})
        rc=np.asarray(X[cm][:,ridx].sum(0)).ravel(); lib=tot[cm].sum()
        cpm=rc/lib*1e6; pbrows.append(pd.Series(np.log2(cpm+1.0),index=frozen,name=(don,s)))
pd.DataFrame(cc).to_csv("results/resilience/pseudobulk_cell_counts.csv",index=False)
pd.DataFrame(inc).to_csv("results/resilience/pseudobulk_inclusion_log.csv",index=False)
print("inclusion: RES/SAD donor x subtype included =",sum(r['included'] for r in inc),"/",len(inc),
      "| group-wise excluded:", pd.DataFrame(inc).query("included==False").shape[0])

logcpm=pd.DataFrame(pbrows)  # rows=(donor,subtype), cols=receptors
# per-receptor z across all RES+SAD donor x subtype samples
z=(logcpm-logcpm.mean(0))/logcpm.std(0,ddof=1)
subtype_score=z.mean(1)  # mean over 56 receptors -> per donor x subtype score
ss=subtype_score.reset_index(); ss.columns=["donor","subtype","score"]
# 2-level maturation: immature=DiffN ; mature = equal-weight mean(NTF3,CHRM3)
imm=ss[ss.subtype=="DiffN"][["donor","score"]].rename(columns={"score":"receiver_score"}); imm["maturation_stage"]="immature"
mat=ss[ss.subtype.isin(["MatN_SGCZ_NTF3","MatN_SGCZ_CHRM3"])].groupby("donor").score.mean().reset_index().rename(columns={"score":"receiver_score"}); mat["maturation_stage"]="mature"
final=pd.concat([imm,mat],ignore_index=True)
final["group"]=final.donor.map(gmap); final["Run"]=final.donor.map(runmap)
final=final[["donor","group","Run","maturation_stage","receiver_score"]].sort_values(["group","donor","maturation_stage"])
final.to_csv("results/resilience/donor_stage_receptor_scores.tsv",sep="\t",index=False)
print(f"[A3] final analysis table: {final.shape[0]} rows (donor x stage); donors RES {(final.group=='RES').sum()//2} SAD {(final.group=='SAD').sum()//2}")
print("saved donor_stage_receptor_scores.tsv (z-scored 56-receptor composite; subtype-balanced mature)")
