#!/usr/bin/env python
"""Build GSE278576 aging human hippocampus snRNA AnnData (RAW counts) for A1-PRIMARY external
replication. Independent lab/donors vs primary GSE268609.

Dataset: GSE278576 (PMID 39463924, "Epigenetic and 3D genome reprogramming during the aging of
human hippocampus", 2024; verified via !Series_pubmed_id).

The deposited filtered Seurat object holds ONLY an SCT assay (corrected, non-integer counts) — wrong
input for DESeq2. So we assemble RAW UMI counts from the 40 per-donor CellRanger matrices
(`GSE278576_hcXXXX_raw_feature_bc_matrix.h5`) and subset each to that donor's FILTERED barcodes taken
from the deposited cell-metadata tsv (which also carries the author annotations Age/Gender/subclass).
This mirrors the primary pipeline (raw counts -> donor pseudobulk -> pyDESeq2).

Out: processed/per_dataset/GSE278576_aging.h5ad  (X=raw counts, layers['counts'],
     obs: donor_id, age_years, sex, age_group, subclass, celltype_l1, grp[YA/HA/mid]).
"""
import os, glob
import numpy as np, pandas as pd, scanpy as sc, anndata as ad

P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
D = f"{P}/raw/GSE278576_aging_hippo_multiome"
PD = f"{D}/per_donor"
MD = f"{D}/GSE278576_hippocampus_RNA_seurat_object_filtered_cells_metadata.tsv.gz"
NICHE_MAP = {"Astro": "Astro", "Microglia": "Micro", "OPC": "OPC", "Oligo": "Oligo", "Endo": "Endo"}

md = pd.read_csv(MD, sep="\t", low_memory=False)
bc_col = "bacrode" if "bacrode" in md.columns else md.columns[0]
md = md.rename(columns={bc_col: "barcode"})
# recover the bare CellRanger barcode = last "_"-delimited token. Metadata barcodes are
# "<donor>_<AAAC...-1>" but some donors carry a deep-reseq infix "<donor>_deep_<AAAC...-1>",
# so a fixed "<donor>_" strip fails for them; the bare barcode is always the final token.
md["cell_bc"] = md["barcode"].str.rsplit("_", n=1).str[-1]
donors = sorted(md["orig.ident"].astype(str).unique())
print(f"metadata: {len(md)} filtered cells, {len(donors)} donors")

missing = [d for d in donors if not os.path.exists(f"{PD}/GSE278576_{d}_raw_feature_bc_matrix.h5")]
if missing:
    raise SystemExit(f"MISSING per-donor matrices ({len(missing)}): {missing[:6]}... — wait for download")

parts = []
for d in donors:
    a = sc.read_10x_h5(f"{PD}/GSE278576_{d}_raw_feature_bc_matrix.h5")
    a.var_names_make_unique()
    sub = md[md["orig.ident"].astype(str) == d]
    bc2full = dict(zip(sub["cell_bc"], sub["barcode"]))   # bare barcode -> full metadata barcode
    a = a[a.obs_names.isin(bc2full)].copy()
    a.obs_names = [bc2full[b] for b in a.obs_names]       # -> full metadata barcode (for join)
    parts.append(a)
    print(f"  {d}: kept {a.n_obs}/{len(sub)} cells")

A = ad.concat(parts, join="outer", fill_value=0, merge="same")
A.X = A.X.tocsr()
# annotate from metadata (barcode-aligned)
m = md.set_index("barcode").loc[A.obs_names]
A.obs["donor_id"] = m["orig.ident"].astype(str).values
A.obs["age_years"] = pd.to_numeric(m["Age"], errors="coerce").values
A.obs["sex"] = m["Gender"].astype(str).values
A.obs["age_group"] = m["age_group"].astype(str).values
A.obs["subclass"] = m["subclass"].astype(str).values
A.obs["celltype_l1"] = A.obs["subclass"].map(NICHE_MAP).fillna("other").astype("category")
A.obs["grp"] = np.where(A.obs.age_years < 40, "YA", np.where(A.obs.age_years >= 60, "HA", "mid"))
A.layers["counts"] = A.X.copy()

os.makedirs(f"{P}/processed/per_dataset", exist_ok=True)
A.write(f"{P}/processed/per_dataset/GSE278576_aging.h5ad")

print(f"\nassembled {A.n_obs} cells x {A.n_vars} genes (RAW counts; X.max={A.X.max():.0f}, integer={np.allclose(A.X.data, np.round(A.X.data))})")
print("grp -> n_donors:", A.obs.groupby("grp")["donor_id"].nunique().to_dict())
sub = A.obs[A.obs.celltype_l1 != "other"]
print("\nniche celltype x grp (donors with >=10 cells):")
for ct in ["Astro", "Micro", "Endo", "OPC", "Oligo"]:
    s = sub[sub.celltype_l1 == ct]
    pc = s.groupby(["grp", "donor_id"]).size()
    ya = (pc.loc["YA"] >= 10).sum() if "YA" in pc.index.get_level_values(0) else 0
    ha = (pc.loc["HA"] >= 10).sum() if "HA" in pc.index.get_level_values(0) else 0
    print(f"  {ct:6s} cells={len(s):7d}  YA_donors>=10={ya:2d}  HA_donors>=10={ha:2d}")
fz = pd.read_csv(f"{P}/results/validation/frozen_primary_signature.csv")
pres = fz.gene.isin(set(A.var_names))
print(f"\nfrozen genes present: {pres.sum()}/{len(fz)} hit rows ({fz.loc[pres,'gene'].nunique()}/{fz.gene.nunique()} unique)")
print("wrote processed/per_dataset/GSE278576_aging.h5ad")
