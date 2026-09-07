#!/usr/bin/env python
# CONSEQUENCE (Phase A) — STEP 03c (part 1/2): build leave-one-ligand-family-out score tables.
# Provenance: VERBATIM from origin session 027c5fdb (transcript line 883, python part); feeds 03c_robustness_extended_fit.R.
# env: bio.  DATA-DEPENDENT -> run from $PROJ.  Out: results/resilience/{_famdrop_2stage.tsv, _famdrop_log.tsv}.
import pandas as pd, numpy as np, anndata as ad, scipy.sparse as sp
PRIMARY=["DiffN","MatN_SGCZ_NTF3","MatN_SGCZ_CHRM3"]
frozen=pd.read_csv("results/resilience/frozen_receiver_receptor_set.csv").receptor.tolist()
edges=pd.read_csv("results/resilience/frozen_ligand_receptor_edges.csv")
meta=pd.read_csv("processed/per_dataset/GSE325391_adultgc_metadata.csv",index_col=0).drop_duplicates('sample').set_index('sample'); runmap=meta['Run'].astype(str)
A=ad.read_h5ad("processed/per_dataset/GSE325391_resilience.h5ad",backed='r'); obs=A.obs
mask=(obs['group'].astype(str).isin(["RES","SAD"]))&(obs['cell_type'].astype(str).isin(PRIMARY))
sub=A[mask.values].to_memory(); X=sub.X; X=X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
vn=pd.Index(sub.var_names.astype(str)); ridx=[vn.get_loc(r) for r in frozen]; tot=np.asarray(X.sum(1)).ravel()
d=sub.obs['donor_id'].astype(str).values; ct=sub.obs['cell_type'].astype(str).values
gmap=sub.obs.drop_duplicates('donor_id').set_index(sub.obs.drop_duplicates('donor_id').donor_id.astype(str))['group'].astype(str)
idx=[]; mat=[]
for don in sorted(set(d)):
    for s in PRIMARY:
        cm=(d==don)&(ct==s)
        if cm.sum()<20: continue
        rc=np.asarray(X[cm][:,ridx].sum(0)).ravel(); idx.append((don,s)); mat.append(np.log2(rc/tot[cm].sum()*1e6+1))
L=pd.DataFrame(mat,index=pd.MultiIndex.from_tuples(idx,names=["donor","subtype"]),columns=frozen)
Z=(L-L.mean(0))/L.std(0,ddof=1)
def two_stage(sc,label):
    ss=sc.reset_index(); ss.columns=["donor","subtype","score"]
    imm=ss[ss.subtype=="DiffN"][["donor","score"]].rename(columns={"score":"receiver_score"}); imm["maturation_stage"]="immature"
    mt=ss[ss.subtype!="DiffN"].groupby("donor").score.mean().reset_index().rename(columns={"score":"receiver_score"}); mt["maturation_stage"]="mature"
    f=pd.concat([imm,mt],ignore_index=True); f["group"]=f.donor.map(gmap); f["Run"]=f.donor.map(runmap); f["variant"]=label; return f
FAM={"cytokine":["CXCL2","IL15"],"collagen":["COL4A2","COL8A1","COL12A1","COL21A1"],
     "matricellular":["TNC","EFEMP1","CCN2","SPON1","SERPINE1"],"apolipoprotein":["APOC1","APOD","APOE"],
     "other_growth":["KITLG","SEMA3D","DKK2","FAP","LGALS9","CYTL1","AZGP1"]}
rows=[]; fam_drop_log=[]
for fam,ligs in FAM.items():
    rec_fam=set(edges[edges.ligand.isin(ligs)].receptor); rec_other=set(edges[~edges.ligand.isin(ligs)].receptor)
    drop=sorted((rec_fam-rec_other)&set(frozen)); keep=[r for r in frozen if r not in drop]
    fam_drop_log.append({"family":fam,"n_exclusive_receptors_dropped":len(drop),"dropped":";".join(drop),"n_remaining":len(keep)})
    rows.append(two_stage(Z[keep].mean(1), f"drop_family_{fam}"))
pd.concat(rows,ignore_index=True).to_csv("results/resilience/_famdrop_2stage.tsv",sep="\t",index=False)
pd.DataFrame(fam_drop_log).to_csv("results/resilience/_famdrop_log.tsv",sep="\t",index=False)
print("[A5+] family-exclusive receptor drops:")
print(pd.DataFrame(fam_drop_log)[["family","n_exclusive_receptors_dropped","n_remaining"]].to_string(index=False))
