#!/usr/bin/env python
"""Data-intrinsic tissue-prep check (read-only, backed obs).
DG-microdissection lacks CA pyramidal neurons; whole-hippocampus has both.
So per-donor coexistence of DG_GC and CA_ExN reveals the prep without metadata.
Decisive for the niche-neuron N0 gate (docs/NEURON_NICHE_PLAN.md)."""
import os, anndata as ad, pandas as pd, numpy as np
PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
A = ad.read_h5ad(f"{PROJ}/processed/per_dataset/GSE268609_anchor.h5ad", backed="r")
obs = A.obs.copy()
neur = ["DG_GC", "CA_ExN", "InN"]
# per-donor cell counts for each neuron type
ct = (obs.assign(n=1).pivot_table(index=["donor_id", "Group"], columns="celltype_l1",
        values="n", aggfunc="sum", observed=True, fill_value=0))
for c in neur + ["Astro", "Micro"]:
    if c not in ct.columns: ct[c] = 0
ct = ct.reset_index()
ct["tot"] = obs.groupby("donor_id", observed=True).size().reindex(ct["donor_id"]).values
# infer prep: whole-hippo if CA_ExN well represented; DG-enriched if CA_ExN ~absent
ct["CA_frac"] = ct["CA_ExN"] / ct["tot"]
ct["DGGC_frac"] = ct["DG_GC"] / ct["tot"]
ct["prep_guess"] = np.where(ct["CA_ExN"] >= 10, "has_CA(whole-hippo-like)", "CA-depleted(DG-microdissect?)")
print("=== per-donor neuron cell counts (sorted by Group, CA_ExN) ===")
cols = ["donor_id", "Group", "tot", "DG_GC", "CA_ExN", "InN", "Astro", "Micro", "CA_frac", "prep_guess"]
show = ct[cols].sort_values(["Group", "CA_ExN"]).reset_index(drop=True)
pd.set_option("display.width", 200, "display.max_rows", 60, "display.float_format", lambda x: f"{x:.3f}")
print(show.to_string(index=False))
print("\n=== how many donors have BOTH DG_GC>=10 AND CA_ExN>=10, per Group? ===")
both = ct.assign(both=(ct["DG_GC"] >= 10) & (ct["CA_ExN"] >= 10),
                 dg_only=(ct["DG_GC"] >= 10) & (ct["CA_ExN"] < 10))
print(both.groupby("Group", observed=True)[["both", "dg_only"]].sum().to_string())
print("\nIf 'both' ~= each group's donor count and 'dg_only' ~0 => DG_GC & CA_ExN coexist")
print("within donors => NO DG-vs-CA tissue-prep confound in the anchor.")
