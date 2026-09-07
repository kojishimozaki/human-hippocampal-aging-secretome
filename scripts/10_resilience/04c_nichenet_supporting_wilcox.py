#!/usr/bin/env python
# CONSEQUENCE (Phase A) — STEP 04c (part 2/2): competitive-Wilcoxon on the top-50 union (supporting result).
# Provenance: VERBATIM from origin session 027c5fdb (transcript line 853, python part). Run AFTER 04c_nichenet_supporting_targets.R.
#   See the CAVEAT in 04c_nichenet_supporting_targets.R about the committed nichenet_supporting_results.tsv.
# env: bio.  DATA-DEPENDENT (GSE325391 h5ad + results/de/de_GSE325391_RESvSAD.csv) -> run from $PROJ.
# Out: results/resilience/nichenet_supporting_results.tsv.  Result: 0/22 ligands BH<0.05.
import pandas as pd, numpy as np, anndata as ad, scipy.sparse as sp
from scipy.stats import mannwhitneyu
# SAD-RES ranking = -(RES-SAD stat). de stat>0 = up in RES -> SAD-RES = -stat
de=pd.read_csv("results/de/de_GSE325391_RESvSAD.csv")
de=de.dropna(subset=["stat"]).drop_duplicates("gene")
de["sad_minus_res"]=-de["stat"]
rankv=de.set_index("gene")["sad_minus_res"]

# CTRL-granule detectable gene set (same rule as receptors, all genes)
PRIMARY=["DiffN","MatN_SGCZ_NTF3","MatN_SGCZ_CHRM3"]
A=ad.read_h5ad("processed/per_dataset/GSE325391_resilience.h5ad",backed='r'); obs=A.obs
cm=(obs['group'].astype(str)=="CTRL")&(obs['cell_type'].astype(str).isin(PRIMARY))
sub=A[cm.values].to_memory(); X=sub.X; X=X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
vn=np.array(sub.var_names.astype(str)); tot=np.asarray(X.sum(1)).ravel()
d=sub.obs['donor_id'].astype(str).values; ct=sub.obs['cell_type'].astype(str).values
det=np.zeros(len(vn),dtype=int);
import collections
perdon=collections.defaultdict(lambda: np.zeros(len(vn),dtype=bool))
for don in sorted(set(d)):
    anysub=np.zeros(len(vn),dtype=bool)
    for s in PRIMARY:
        m=(d==don)&(ct==s)
        if m.sum()<20: continue
        cpm=np.asarray(X[m].sum(0)).ravel()/tot[m].sum()*1e6
        anysub |= (cpm>=1)
    det += anysub.astype(int)
detectable=set(vn[det>=3])
print("CTRL-detectable genes:", len(detectable))

targets=pd.read_csv("results/resilience/_nichenet_targets_top50.tsv",sep="\t")
targets=targets[targets.target.isin(detectable)]   # restrict to CTRL-detectable
union=sorted(set(targets.target)&set(rankv.index))
def comp_enrich(tset):
    tset=[g for g in tset if g in rankv.index]
    if len(tset)<5: return (np.nan,np.nan,len(tset))
    inset=rankv[tset]; out=rankv[~rankv.index.isin(tset)]
    U,p=mannwhitneyu(inset,out,alternative="two-sided")
    # direction: median(in) vs median(out); positive = SAD-high
    return (float(inset.median()-out.median()), float(p), len(tset))
es,ep,en=comp_enrich(union)
print(f"[A6] UNION target enrichment on SAD-RES ranking: n={en}, median diff={es:+.3f}, p={ep:.3f} (positive=SAD-high)")
# per-ligand
rows=[]
for lg,g in targets.groupby("ligand"):
    e,p,n=comp_enrich(list(g.target)); rows.append({"ligand":lg,"n_targets_detectable":n,"median_diff_SADminusRES":e,"p":p})
pl=pd.DataFrame(rows).dropna()
from statsmodels.stats.multitest import multipletests
pl["padj_BH"]=multipletests(pl.p,method="fdr_bh")[1]
pl=pl.sort_values("p")
out=pd.concat([pd.DataFrame([{"ligand":"UNION","n_targets_detectable":en,"median_diff_SADminusRES":es,"p":ep,"padj_BH":np.nan}]),pl],ignore_index=True)
out.to_csv("results/resilience/nichenet_supporting_results.tsv",sep="\t",index=False)
print(f"per-ligand: {(pl.padj_BH<0.05).sum()}/{len(pl)} ligands BH<0.05 | min padj={pl.padj_BH.min():.3f}")
print(out.head(8).to_string(index=False))
