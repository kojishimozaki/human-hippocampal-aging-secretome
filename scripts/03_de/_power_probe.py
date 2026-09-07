#!/usr/bin/env python
"""Probe: donors/group with >=10 cells per cell type in GSE268609 anchor.
Read-only obs probe to assess neuron-vs-niche DE feasibility. Backed mode (no X load)."""
import os, anndata as ad, pandas as pd
PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
A = ad.read_h5ad(f"{PROJ}/processed/per_dataset/GSE268609_anchor.h5ad", backed="r")
obs = A.obs.copy()
print("total nuclei:", len(obs))
print("\n=== celltype_l1 totals ===")
print(obs["celltype_l1"].value_counts().to_string())
print("\n=== donors per Group ===")
print(obs.groupby("Group", observed=True)["donor_id"].nunique().to_string())
print("\n=== donors/group with >=10 cells, per cell type ===")
groups = ["YA", "HA", "MCI", "AD", "SA"]
rows = []
for ct in obs["celltype_l1"].cat.categories if hasattr(obs["celltype_l1"], "cat") else sorted(obs["celltype_l1"].unique()):
    sub = obs[obs["celltype_l1"].astype(str) == ct]
    rec = {"celltype": ct, "n_cells": len(sub)}
    for g in groups:
        sg = sub[sub["Group"].astype(str) == g]
        nd = (sg.groupby("donor_id", observed=True).size() >= 10).sum()
        rec[g] = int(nd)
    rows.append(rec)
tab = pd.DataFrame(rows).set_index("celltype")
print(tab.to_string())
print("\nLegend: cell counts of donors (>=10 cells) per group. DE needs >=4/group.")
print("Normal-aging neuron contrast = YA vs HA; AD neuron contrast = AD vs HA.")
