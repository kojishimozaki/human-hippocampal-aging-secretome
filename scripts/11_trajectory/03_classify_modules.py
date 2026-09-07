#!/usr/bin/env python
"""Phase B (TRAJECTORY) B6-B7 — projection gate + 4-category classification (ashr T_class), threshold
sensitivity (9 combos), shrinkage sensitivities, and UP/DOWN module trajectories (YA/HA/AD + MCI context).
Pre-reg v1 §B + amendments v1.3/v1.4. Classification uses ashr LFC ONLY; T_class = A + D. env: bio. seed 42.
"""
import os, json, numpy as np, pandas as pd
PROJ=os.environ.get("PROJ",os.getcwd()); OUT=f"{PROJ}/results/trajectory"; np.random.seed(42)
NICHE=["Astro","Micro","Endo","OPC","Oligo"]
h=pd.read_csv(f"{OUT}/_frozen60_contrasts_R.tsv",sep="\t")
s=np.sign(h.frozen_log2FC)
# ashr-based classification effects (primary)
h["A"]=s*h.HA_YA_LFC_ashr; h["D"]=s*h.AD_HA_LFC_ashr; h["T_class"]=h.A+h.D
h["T_direct_shrunk"]=s*h.AD_YA_LFC_ashr; h["T_discrepancy"]=h.T_direct_shrunk-h.T_class
# raw-based (supporting only)
h["A_raw"]=s*h.HA_YA_LFC_raw; h["D_raw"]=s*h.AD_HA_LFC_raw
h["T_raw_direct"]=s*h.AD_YA_LFC_raw; h["T_raw_sum"]=h.A_raw+h.D_raw; h["raw_identity_error"]=h.T_raw_direct-h.T_raw_sum

def classify(A,D,T,rel=0.5,abso=0.25):
    if not np.isfinite(A): return ("not_testable","")
    if A<0: return ("aging_direction_not_recapitulated","")
    if A<abso: return ("weak_indeterminate_aging_projection","")
    d=max(abso,rel*A)
    if D>=d: return ("amplified","")
    if D<=-d:
        if T>0: return ("AD_attenuated","")
        return ("AD_divergent", "reverted_toward_YA" if abs(T)<0.25 else "direction_reversed")
    return ("maintained","")

# primary classification (rel=0.5, abs=0.25)
cls=[classify(a,d,t) for a,d,t in zip(h.A,h.D,h.T_class)]
h["delta_final"]=[max(0.25,0.5*a) if (np.isfinite(a) and a>=0.25) else np.nan for a in h.A]
h["aging_projection_status"]=np.where(h.A<0,"not_recapitulated",np.where(h.A<0.25,"weak_indeterminate","eligible"))
h["trajectory_category"]=[c[0] for c in cls]; h["trajectory_subflag"]=[c[1] for c in cls]
# raw-LFC category (supporting sensitivity)
h["trajectory_category_rawLFC"]=[classify(a,d,t)[0] for a,d,t in zip(h.A_raw,h.D_raw,h.A_raw+h.D_raw)]
h["category_margin_D_minus_delta"]=h.D-h.delta_final
cols=["gene","celltype","frozen_direction","frozen_log2FC","baseMean","betaConv",
 "HA_YA_LFC_coef","AD_HA_LFC_coef","AD_YA_LFC_coef","coefficient_identity_error",
 "HA_YA_LFC_raw","AD_HA_LFC_raw","AD_YA_LFC_raw","results_identity_error",
 "HA_YA_LFC_ashr","AD_HA_LFC_ashr","AD_YA_LFC_ashr",
 "A_raw","D_raw","T_raw_sum","T_raw_direct","raw_identity_error",
 "A","D","T_class","T_direct_shrunk","T_discrepancy","delta_final","category_margin_D_minus_delta",
 "aging_projection_status","trajectory_category","trajectory_subflag","trajectory_category_rawLFC",
 "pvalue_HA_YA","padj_HA_YA","pvalue_AD_HA","padj_AD_HA","pvalue_AD_YA","padj_AD_YA"]
h[cols].to_csv(f"{OUT}/frozen60_GSE268609_hit_trajectory.tsv",sep="\t",index=False)

# B7 threshold sensitivity (9 combos); main = 0.50 + 0.25
rows=[]
for rel in [0.33,0.50,0.67]:
    for abso in [0.20,0.25,0.30]:
        cc=pd.Series([classify(a,d,t,rel,abso)[0] for a,d,t in zip(h.A,h.D,h.T_class)]).value_counts()
        rows.append(dict(relative=rel,absolute=abso,is_main=(rel==0.50 and abso==0.25),
            amplified=int(cc.get("amplified",0)),maintained=int(cc.get("maintained",0)),
            AD_attenuated=int(cc.get("AD_attenuated",0)),AD_divergent=int(cc.get("AD_divergent",0)),
            not_recapitulated=int(cc.get("aging_direction_not_recapitulated",0)),
            weak_indeterminate=int(cc.get("weak_indeterminate_aging_projection",0))))
pd.DataFrame(rows).to_csv(f"{OUT}/frozen60_GSE268609_threshold_sensitivity.tsv",sep="\t",index=False)

# modules: UP/DOWN per celltype, donor-level z-scored log2CPM, group summary incl MCI
sig=pd.read_csv(f"{PROJ}/results/validation/frozen_primary_signature.csv")
mod_rows=[]
for ct in NICHE:
    me=pd.read_csv(f"{OUT}/module_expression_{ct}.tsv",sep="\t",index_col=0)
    grp=me["grp"]; expr=me.drop(columns=["grp","arm"])
    for direction in ["UP","DOWN"]:
        genes=[g for g in sig[(sig.celltype==ct)&(sig.direction==direction)].gene if g in expr.columns]
        if not genes: continue
        z=(expr[genes]-expr[genes].mean(0))/expr[genes].std(0,ddof=1)
        score=z.mean(1)
        for g in ["YA","HA","MCI","AD"]:
            v=score[grp==g].values
            if len(v)==0: continue
            ci=1.96*v.std(ddof=1)/np.sqrt(len(v)) if len(v)>1 else np.nan
            mod_rows.append(dict(celltype=ct,module=direction,group=g,n_donors=len(v),
                mean_score=float(v.mean()),ci95=float(ci) if np.isfinite(ci) else np.nan,
                primary_inference=(g in ["YA","HA","AD"])))
pd.DataFrame(mod_rows).to_csv(f"{OUT}/frozen60_GSE268609_module_scores.tsv",sep="\t",index=False)

# classification manifest
vc=h.trajectory_category.value_counts()
conc_raw_ashr=float((h.trajectory_category==h.trajectory_category_rawLFC).mean())
manifest=dict(amendments=["v1","v1.1","v1.2(superseded)","v1.3","v1.4"],
  classification_effect="ashr posterior-mean LFC; T_class=A+D; raw-LFC category = supporting sensitivity",
  primary_thresholds={"relative":0.5,"absolute":0.25},
  category_counts_primary={k:int(v) for k,v in vc.items()},
  raw_vs_ashr_category_concordance=conc_raw_ashr,
  T_discrepancy_abs_median=float(h.T_discrepancy.abs().median()),
  T_discrepancy_abs_max=float(h.T_discrepancy.abs().max()),
  n_hits=int(len(h)), n_eligible=int((h.aging_projection_status=="eligible").sum()))
json.dump(manifest,open(f"{OUT}/frozen60_GSE268609_classification_manifest.json","w"),indent=2)
print("=== B6 primary classification (rel0.5/abs0.25) ===")
print(vc.to_string())
print(f"\neligible (A>=0.25): {int((h.aging_projection_status=='eligible').sum())}/60 | "
      f"not-recapitulated: {int((h.aging_projection_status=='not_recapitulated').sum())} | "
      f"weak: {int((h.aging_projection_status=='weak_indeterminate').sum())}")
print(f"raw-vs-ashr category concordance: {conc_raw_ashr:.2f} | T_discrepancy |median|={h.T_discrepancy.abs().median():.3f} max={h.T_discrepancy.abs().max():.3f}")
print("\n[B6-B7] hit_trajectory + threshold_sensitivity + module_scores + manifest saved.")
