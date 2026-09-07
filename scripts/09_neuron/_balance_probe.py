#!/usr/bin/env python
"""Read-only N0 checkpoint probe: per-celltype arm x group donor balance + ~arm+grp estimability.

Backs the pre-DE checkpoint for the neuron-niche arm (docs/NEURON_N0_HANDOFF.md). NO DE is run
and the 9.3 GB anchor X is NOT loaded (backed='r', obs only) -- same cheap pattern as
scripts/03_de/_power_probe.py / _tissueprep_probe.py.

For every (celltype x contrast) it reports, among donors with >=10 cells of that celltype:
  - donors per group, split by preparation arm (WH = whole-hippocampus orig.ident<=14;
    DG = DG-microdissected orig.ident>=15);
  - whether design ~arm+grp is estimable, using the SAME rule the validated niche check
    (scripts/03_de/06_arm_robustness.py:72) used: arm has 2 levels AND every arm x group
    cell is populated (>=1). If not estimable -> the DE will fall back to ~grp for that
    celltype x contrast (and say so).

Contrasts (reference pivot = HA): YA-vs-HA (aging), AD-vs-HA, MCI-vs-HA.
Cell types: neurons DG_GC/CA_ExN/InN (even-handed), niche Astro/Micro/Oligo/OPC/Endo (symmetry),
immature NSC/Neuroblast/Immature (YA-vs-HA exploratory only; AD infeasible n=1/0).
"""
import os
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_cache")
import numpy as np, pandas as pd, anndata as ad

np.random.seed(42)
P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")

NEURON   = ["DG_GC", "CA_ExN", "InN"]
NICHE    = ["Astro", "Micro", "Oligo", "OPC", "Endo"]
IMMATURE = ["NSC", "Neuroblast", "Immature"]
MIN_CELLS = 10

# ---- (A) donor (SampleNumber) -> arm, from orig.ident (<=14 WH, >=15 DG) -------------------
meta = pd.read_csv(f"{P}/processed/per_dataset/GSE268609_metadata.csv", index_col=0)
donors = (meta[["orig.ident", "SampleNumber", "Group"]].drop_duplicates()
          .rename(columns={"orig.ident": "oid"}))
donors["SampleNumber"] = donors.SampleNumber.astype(str)
donors["arm"] = np.where(donors.oid.astype(int) <= 14, "WH", "DG")
# rigor: each donor must map to exactly ONE arm (no donor straddling both preps)
straddle = donors.groupby("SampleNumber").arm.nunique()
assert straddle.max() == 1, f"donors straddling arms: {straddle[straddle>1].index.tolist()}"
arm_map = dict(zip(donors.SampleNumber, donors.arm))
print("=== donor-level arm x Group (all donors, from orig.ident) ===")
dd = donors.drop_duplicates("SampleNumber")
print(pd.crosstab(dd.Group, dd.arm).to_string())
print(f"(total donors: {dd.SampleNumber.nunique()})")

# ---- (B) backed obs read of the anchor (no X) ----------------------------------------------
A = ad.read_h5ad(f"{P}/processed/per_dataset/GSE268609_anchor.h5ad", backed="r")
obs = A.obs.copy()
obs["donor_id"] = obs["donor_id"].astype(str)
obs["arm"] = obs["donor_id"].map(arm_map)
obs["ct"] = obs["celltype_l1"].astype(str)
print(f"\nanchor obs: {len(obs)} nuclei; donors {obs.donor_id.nunique()}; "
      f"arm NA cells: {int(obs.arm.isna().sum())}")
print("\n=== celltype_l1 totals (label sanity-check) ===")
print(obs["ct"].value_counts().to_string())

# per (celltype, donor) cell count, with that donor's group + arm
cd = (obs.groupby(["ct", "donor_id"], observed=True).size().rename("n_cells").reset_index())
gmap = dd.set_index("SampleNumber").Group.to_dict()
cd["Group"] = cd.donor_id.map(gmap)
cd["arm"] = cd.donor_id.map(arm_map)
cd = cd[cd.n_cells >= MIN_CELLS]   # >=10 cells/donor for that celltype

CONTRASTS = [("YA-vs-HA", "YA", "HA"), ("AD-vs-HA", "AD", "HA"), ("MCI-vs-HA", "MCI", "HA")]

def block(label, celltypes, contrasts):
    print(f"\n########## {label} ##########")
    rows = []
    for ct in celltypes:
        sub = cd[cd.ct == ct]
        for cname, test, ref in contrasts:
            d = sub[sub.Group.isin([test, ref])]
            bal = pd.crosstab(d.Group, d.arm)
            for g in (test, ref):
                if g not in bal.index: bal.loc[g] = 0
            for a in ("WH", "DG"):
                if a not in bal.columns: bal[a] = 0
            bal = bal.loc[[test, ref], ["WH", "DG"]]
            n_test = int(bal.loc[test].sum()); n_ref = int(bal.loc[ref].sum())
            arm_levels = int((bal.sum(0) > 0).sum())
            estimable = (arm_levels >= 2) and (int(bal.values.min()) >= 1)
            powered = (n_test >= 4) and (n_ref >= 4)
            rows.append(dict(
                celltype=ct, contrast=cname,
                n_test=n_test, test_WH=int(bal.loc[test, "WH"]), test_DG=int(bal.loc[test, "DG"]),
                n_ref=n_ref, ref_WH=int(bal.loc[ref, "WH"]), ref_DG=int(bal.loc[ref, "DG"]),
                powered_4pg=powered, arm_estimable=estimable,
                design=("~arm+grp" if estimable else "~grp (arm not estimable)")))
    t = pd.DataFrame(rows)
    pd.set_option("display.width", 220, "display.max_rows", 80)
    print(t.to_string(index=False))
    return t

t_neuron = block("NEURONS (even-handed: DG_GC / CA_ExN / InN) -- all 3 contrasts", NEURON, CONTRASTS)
t_niche  = block("NICHE (symmetry: Astro / Micro / Oligo / OPC / Endo) -- all 3 contrasts", NICHE, CONTRASTS)
t_imm    = block("IMMATURE (NSC / Neuroblast / Immature) -- YA-vs-HA exploratory ONLY",
                 IMMATURE, [("YA-vs-HA", "YA", "HA")])

print("\n=== checkpoint summary ===")
allt = pd.concat([t_neuron, t_niche, t_imm], ignore_index=True)
ne = allt[~allt.arm_estimable]
print(f"celltype x contrast cells evaluated: {len(allt)}")
print(f"  ~arm+grp estimable: {int(allt.arm_estimable.sum())}/{len(allt)}")
print(f"  powered (>=4/group): {int(allt.powered_4pg.sum())}/{len(allt)}")
if len(ne):
    print("  NOT estimable under ~arm+grp (will fall back to ~grp):")
    print(ne[["celltype", "contrast", "n_test", "test_WH", "test_DG", "n_ref", "ref_WH", "ref_DG"]].to_string(index=False))
else:
    print("  (all evaluated cells estimable under ~arm+grp)")
