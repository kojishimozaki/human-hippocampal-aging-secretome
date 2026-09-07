#!/usr/bin/env python
# CONSEQUENCE (Phase A) — STEP 04b (part 1/2): competitive-Wilcoxon downstream enrichment (top-20/50/100).
# Provenance: VERBATIM from origin session 027c5fdb (transcript line 896, python part). Run AFTER 04a_nichenet_targets.R.
# env: bio.  DATA-DEPENDENT (GSE325391 h5ad + results/de/de_GSE325391_RESvSAD.csv) -> run from $PROJ.
# Out: results/resilience/{nichenet_top20_results.tsv, nichenet_top50_results.tsv, nichenet_top100_results.tsv,
#   _ctrl_detectable_genes.txt, _sadminusres_ranking.tsv, _nichenet_wilcox_stability.tsv}.  Result: 0/22 ligands BH<0.05.
import pandas as pd, numpy as np, anndata as ad, scipy.sparse as sp
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
PRIMARY=["DiffN","MatN_SGCZ_NTF3","MatN_SGCZ_CHRM3"]
de=pd.read_csv("results/de/de_GSE325391_RESvSAD.csv").dropna(subset=["stat"]).drop_duplicates("gene")
de["sad_minus_res"]=-de["stat"]; rankv=de.set_index("gene")["sad_minus_res"]
rankv.to_frame().to_csv("results/resilience/_sadminusres_ranking.tsv",sep="\t")   # for fgsea
# CTRL-detectable set
A=ad.read_h5ad("processed/per_dataset/GSE325391_resilience.h5ad",backed='r'); obs=A.obs
cm=(obs['group'].astype(str)=="CTRL")&(obs['cell_type'].astype(str).isin(PRIMARY))
sub=A[cm.values].to_memory(); X=sub.X; X=X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
vn=np.array(sub.var_names.astype(str)); tot=np.asarray(X.sum(1)).ravel()
d=sub.obs['donor_id'].astype(str).values; ct=sub.obs['cell_type'].astype(str).values
det=np.zeros(len(vn),dtype=int)
for don in sorted(set(d)):
    anysub=np.zeros(len(vn),dtype=bool)
    for s in PRIMARY:
        m=(d==don)&(ct==s)
        if m.sum()<20: continue
        anysub |= (np.asarray(X[m].sum(0)).ravel()/tot[m].sum()*1e6>=1)
    det+=anysub.astype(int)
detectable=set(vn[det>=3]); pd.Series(sorted(detectable)).to_csv("results/resilience/_ctrl_detectable_genes.txt",index=False,header=False)
T=pd.read_csv("results/resilience/_nichenet_targets_top100_ranked.tsv",sep="\t")
def comp(tset):
    tset=[g for g in set(tset) if g in rankv.index]
    if len(tset)<5: return (np.nan,np.nan,len(tset))
    ins=rankv[tset]; out=rankv[~rankv.index.isin(tset)]
    U,p=mannwhitneyu(ins,out,alternative="two-sided")
    return (float(ins.median()-out.median()),float(p),len(tset))
stab=[]
for N in [20,50,100]:
    tN=T[T['rank']<=N]; tN=tN[tN.target.isin(detectable)]
    union=sorted(set(tN.target))
    eu,pu,nu=comp(union)
    rows=[{"ligand":"UNION","n_targets_detectable":nu,"median_diff_SADminusRES":eu,"p":pu,"padj_BH":np.nan}]
    pl=[]
    for lg,g in tN.groupby("ligand"):
        e,p,n=comp(list(g.target)); pl.append({"ligand":lg,"n_targets_detectable":n,"median_diff_SADminusRES":e,"p":p})
    pl=pd.DataFrame(pl).dropna(subset=["p"]); pl["padj_BH"]=multipletests(pl.p,method="fdr_bh")[1]
    out=pd.concat([pd.DataFrame(rows),pl.sort_values("p")],ignore_index=True)
    out.to_csv(f"results/resilience/nichenet_top{N}_results.tsv",sep="\t",index=False)
    stab.append({"topN":N,"union_n":nu,"union_median_diff":eu,"union_p":pu,
                 "n_ligand_BH<0.05":int((pl.padj_BH<0.05).sum()),"min_ligand_padj":float(pl.padj_BH.min()),
                 "union_direction":"SAD-high" if eu>0 else "RES-high"})
sdf=pd.DataFrame(stab)
print("=== A6 threshold stability (competitive Wilcoxon) ==="); print(sdf.to_string(index=False))
sdf.to_csv("results/resilience/_nichenet_wilcox_stability.tsv",sep="\t",index=False)
