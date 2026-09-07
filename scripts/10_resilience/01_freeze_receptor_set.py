#!/usr/bin/env python
# CONSEQUENCE (Phase A) — STEP 01: build + freeze the receiver receptor set (76 -> 56 via CTRL-only detectability).
# Provenance: recovered VERBATIM from origin session 027c5fdb (transcript line 809); produced the committed
#   frozen receptor set (tag freeze/consequence-receptor-set-v1, be864a3) byte-for-byte (agent-verified).
# env: bio.  DATA-DEPENDENT (reads processed/per_dataset/GSE325391_resilience.h5ad) -> run from $PROJ.
# Out: results/resilience/{frozen_receiver_receptor_set.csv, frozen_ligand_receptor_edges.csv,
#      excluded_receptors_with_reasons.tsv, frozen_receiver_set_manifest.json}
import pandas as pd, numpy as np, anndata as ad, json, datetime, subprocess
PRIMARY=["DiffN","MatN_SGCZ_NTF3","MatN_SGCZ_CHRM3"]
recs=[l.strip() for l in open("results/step0/_receptors_prefilter.txt")]
edges=pd.read_csv("results/step0/consequence_ligand_receptor_edge_audit.tsv",sep="\t")
ligs=[l.strip() for l in open("results/step0/_nichenet_ligands.txt")]

A=ad.read_h5ad("processed/per_dataset/GSE325391_resilience.h5ad", backed='r')
obs=A.obs
ctrl_mask=(obs['group'].astype(str)=="CTRL") & (obs['cell_type'].astype(str).isin(PRIMARY))
sub=A[ctrl_mask.values].to_memory()
X=sub.X.tocsr() if hasattr(sub.X,"tocsr") else sub.X
import scipy.sparse as sp
if not sp.issparse(X): X=sp.csr_matrix(X)
vn=pd.Index(sub.var_names.astype(str))
ridx=[vn.get_loc(r) for r in recs]   # all 76 present (verified A0)
tot=np.asarray(X.sum(1)).ravel()     # per-cell library size (this matrix)
d=sub.obs['donor_id'].astype(str).values; ct=sub.obs['cell_type'].astype(str).values

# CTRL pseudobulk CPM per donor x subtype for the 76 receptors
rows=[]
ctrl_donors=sorted(set(d))
for don in ctrl_donors:
    for s in PRIMARY:
        cm=(d==don)&(ct==s)
        if cm.sum()<20: continue   # min-cell consistent with primary
        rc=np.asarray(X[cm][:,ridx].sum(0)).ravel()
        lib=tot[cm].sum()
        cpm=rc/lib*1e6
        rows.append(pd.DataFrame({"donor":don,"subtype":s,"receptor":recs,"cpm":cpm,"n_cells":cm.sum()}))
pb=pd.concat(rows,ignore_index=True)
pb["detectable"]=pb.cpm>=1.0
# adopt rule: detectable in >=3 of 6 CTRL donors in >=1 subtype
det=pb[pb.detectable].groupby(["receptor","subtype"]).donor.nunique().reset_index(name="n_ctrl_donors_detect")
best=det.groupby("receptor").n_ctrl_donors_detect.max().reindex(recs).fillna(0).astype(int)
adopt=best>=3
frozen=sorted(best[adopt].index)
excluded=sorted(best[~adopt].index)

# save artifacts
fr=pd.DataFrame({"receptor":frozen})
fr["n_ligands_connecting"]=fr.receptor.map(edges.groupby("receptor").ligand.nunique())
fr.to_csv("results/resilience/frozen_receiver_receptor_set.csv",index=False)
edges[edges.receptor.isin(frozen)].to_csv("results/resilience/frozen_ligand_receptor_edges.csv",index=False)
exdf=pd.DataFrame({"receptor":excluded})
exdf["max_ctrl_donors_detectable"]=exdf.receptor.map(best)
exdf["reason"]="not detectable (CPM>=1) in >=3/6 CTRL donors in any primary subtype (DiffN/NTF3/CHRM3)"
exdf.to_csv("results/resilience/excluded_receptors_with_reasons.tsv",sep="\t",index=False)

chk=dict(l.split("\t")[0:2][::-1] if False else (l.split("\t")[1],l.split("\t")[0]) for l in open("preregistration/manifests/frozen_input_checksums_v1.tsv").read().splitlines()[1:])
manifest=dict(
  arm="CONSEQUENCE_A1_frozen_receiver_set", prereg_tag="prereg/secretome-resilience-trajectory-v1",
  prereg_commit="5b29c74b91ba68886a239b066da10919899dfbf7",
  freeze_timestamp="2026-06-20", git_commit=subprocess.run(["git","rev-parse","HEAD"],capture_output=True,text=True).stdout.strip(),
  nichenet_lr_network="refs/nichenet/lr_network_human_21122021.rds",
  nichenet_lr_sha256="47c971d2fbba4ecd0ba7485d1846a74054432a0d1a97ea3e6a79ae27d0da8094",
  frozen_signature_sha256="956e8743cae6e31c9f6b1c6cdd5eb040b81734c56dcd416f9b27ff415765f5ab",
  gse325391_h5ad_sha256="ea814d1753bb6b865d97bb47f68d4b55262d19857d1df24ae21e3dc75d60b350",
  symbol_mapping="direct NicheNet symbol match; 76/76 receptors present in GSE325391 feature space (exact)",
  ctrl_only_filter="CPM>=1 in >=3 of 6 CTRL donors in >=1 of {DiffN, MatN_SGCZ_NTF3, MatN_SGCZ_CHRM3} (>=20 cells/donor-subtype)",
  filter_note="'>=1 SGCZ maturation subtype' operationalised as >=1 of the 3 primary score subtypes (DiffN immature + the two MatN_SGCZ mature), to match the primary score's stage coverage; recorded as an implementation clarification, not a design change",
  n_up_ligands=len(ligs), n_receptors_prefilter=len(recs), n_receptors_frozen=len(frozen), n_excluded=len(excluded),
  score_formula="per-receptor pseudobulk logCPM -> per-receptor z across RES+SAD donor x subtype -> mean over frozen receptors -> one score per donor x subtype",
  ligands=ligs, frozen_receptors=frozen, excluded_receptors=excluded)
json.dump(manifest, open("results/resilience/frozen_receiver_set_manifest.json","w"), indent=2)
print(f"[A1] frozen receptor set: {len(frozen)}/{len(recs)} adopted (excluded {len(excluded)})")
print("frozen receptors:", ", ".join(frozen))
print("excluded:", ", ".join(excluded))
