#!/usr/bin/env python
"""N2 (Layer 2) intrinsic — (N2.2) unbiased pathway GSEA + (N2.3) targeted niche-response readout.

Pre-registered in audit_log/2026-06-18_neuron_niche/RESOLUTION.md. Operates on the N0 genome-wide neuron
DE (no anchor load). seed 42.

N2.3 targeted (FAST, written first): one-sided Mann-Whitney of neuron log2FC, set vs non-set, within the
tested-gene universe, + 2,000x random-set permutation null; per (neuron CT x axis); direction = UP-with-aging
(response to a niche SASP). Sets = refs/neuron_response_sets.gmt. DESCRIPTIVE/correlative, no causal language.

N2.2 GSEA (slower): gseapy prerank on the DESeq2 stat-ranked gene list per (neuron CT x axis) vs cached
Hallmark_2020 / Reactome_2022 / GO_BP_2023; NPERM=1000; seed 42. Descriptive pathway-naming of the
donor-level DE ranking (NOT a donor-level test). AD axis is PMI/arm-confounded.

Outputs: results/validation/neuron_L2_targeted_response.csv, results/validation/neuron_L2_gsea.csv
"""
import os
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache"); os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_cache")
import numpy as np, pandas as pd
from scipy import stats
import warnings; warnings.filterwarnings("ignore")

P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
NEURON = ["DG_GC", "CA_ExN", "InN"]
NPERM = 2000
rng = np.random.default_rng(42)

aging = pd.read_csv(f"{P}/results/de/de_GSE268609_neuron_aging_per_celltype.csv")
advha = pd.read_csv(f"{P}/results/de/de_GSE268609_neuron_ADvHA_per_celltype.csv")
AXES = [("aging", aging), ("ADvHA", advha)]

def read_gmt(path):
    s = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3:
            s[p[0]] = [g for g in p[2:] if g]
    return s
targeted = read_gmt(f"{P}/refs/neuron_response_sets.gmt")
print(f"targeted sets: {[(k, len(v)) for k, v in targeted.items()]}")

# ===== N2.3 targeted niche-response readout (fast) ==============================================
print("\n=== N2.3 targeted niche-response directional readout (UP-with-aging vs background) ===")
trows = []
for axis, de in AXES:
    for ct in NEURON:
        d = de[de.celltype == ct].dropna(subset=["log2FoldChange"])
        universe = dict(zip(d.gene, d.log2FoldChange))
        uni_vals = np.array(list(universe.values()))
        for sname, sgenes in targeted.items():
            sset = set(sgenes)
            inset = np.array([universe[g] for g in sset if g in universe])
            n_in = len(inset)
            if n_in < 5:
                trows.append(dict(axis=axis, celltype=ct, set=sname, n_set=n_in, median_lfc=np.nan,
                                  mean_lfc=np.nan, frac_up=np.nan, mwu_p=np.nan, null_p=np.nan,
                                  up_enriched=False, note="n_set<5 NOT_EVALUABLE")); continue
            bg = np.array([v for g, v in universe.items() if g not in sset])
            U, mwu_p = stats.mannwhitneyu(inset, bg, alternative="greater")
            obs = float(inset.mean())
            null = np.array([rng.choice(uni_vals, size=n_in, replace=False).mean() for _ in range(NPERM)])
            null_p = float((null >= obs).mean())
            trows.append(dict(axis=axis, celltype=ct, set=sname, n_set=n_in,
                              median_lfc=round(float(np.median(inset)), 3), mean_lfc=round(obs, 3),
                              frac_up=round(float((inset > 0).mean()), 3), mwu_p=mwu_p, null_p=null_p,
                              up_enriched=bool(mwu_p < 0.05 and null_p < 0.05), note=""))
targ = pd.DataFrame(trows)
targ.to_csv(f"{P}/results/validation/neuron_L2_targeted_response.csv", index=False)
for axis, _ in AXES:
    print(f"\n--- {axis}: up_enriched sets (mwu_p<0.05 AND null_p<0.05) ---")
    s = targ[(targ.axis == axis) & (targ.up_enriched)]
    if len(s):
        print(s[["celltype", "set", "n_set", "median_lfc", "frac_up", "mwu_p", "null_p"]].to_string(index=False))
    else:
        print("  (none — no targeted niche-response set is UP-enriched above background)")
print("\nfull targeted table head (aging):")
print(targ[targ.axis == "aging"][["celltype", "set", "n_set", "median_lfc", "frac_up", "mwu_p", "null_p"]].to_string(index=False))
print("wrote results/validation/neuron_L2_targeted_response.csv")

# ===== N2.2 unbiased pathway GSEA (prerank) =====================================================
print("\n=== N2.2 GSEA prerank (descriptive pathway-naming of the donor-level DE ranking) ===")
import gseapy
LIBS = {"Hallmark": f"{P}/refs/enrichr/MSigDB_Hallmark_2020.gmt",
        "Reactome": f"{P}/refs/enrichr/Reactome_2022.gmt",
        "GO_BP": f"{P}/refs/enrichr/GO_Biological_Process_2023.gmt"}
def col(df, *cands):
    for c in cands:
        if c in df.columns: return c
    raise KeyError(cands)
grows = []
for axis, de in AXES:
    for ct in NEURON:
        d = de[de.celltype == ct].dropna(subset=["stat"])[["gene", "stat"]].drop_duplicates("gene")
        rnk = d.sort_values("stat", ascending=False).reset_index(drop=True)
        for lname, gmt in LIBS.items():
            try:
                pre = gseapy.prerank(rnk=rnk, gene_sets=gmt, permutation_num=1000, seed=42,
                                     min_size=15, max_size=500, threads=4, outdir=None, no_plot=True, verbose=False)
                r = pre.res2d.copy()
                tc = col(r, "Term"); nesc = col(r, "NES"); fdrc = col(r, "FDR q-val")
                nomc = col(r, "NOM p-val"); ledc = col(r, "Lead_genes", "Leading_edge")
                r = r.rename(columns={tc: "Term", nesc: "NES", fdrc: "FDR_q", nomc: "NOM_p", ledc: "lead_genes"})
                r["NES"] = pd.to_numeric(r.NES, errors="coerce"); r["FDR_q"] = pd.to_numeric(r.FDR_q, errors="coerce")
                r["axis"] = axis; r["celltype"] = ct; r["library"] = lname
                grows.append(r[["axis", "celltype", "library", "Term", "NES", "NOM_p", "FDR_q", "lead_genes"]])
                nsig = int((r.FDR_q < 0.25).sum())
                print(f"  {axis}/{ct}/{lname}: {len(r)} sets, FDR<0.25: {nsig}")
            except Exception as e:
                print(f"  FAIL {axis}/{ct}/{lname}: {type(e).__name__} {str(e)[:120]}")
gsea = pd.concat(grows, ignore_index=True) if grows else pd.DataFrame()
gsea.to_csv(f"{P}/results/validation/neuron_L2_gsea.csv", index=False)
print(f"\nwrote results/validation/neuron_L2_gsea.csv: {len(gsea)} rows")
# top enrichments per (axis, celltype) at FDR<0.25
if len(gsea):
    print("\n=== top GSEA enrichments (FDR<0.25, |NES| desc, max 6 per CT x axis) ===")
    for axis, _ in AXES:
        for ct in NEURON:
            s = gsea[(gsea.axis == axis) & (gsea.celltype == ct) & (gsea.FDR_q < 0.25)].copy()
            if not len(s):
                print(f"\n{axis}/{ct}: (no set FDR<0.25)"); continue
            s = s.reindex(s.NES.abs().sort_values(ascending=False).index).head(6)
            print(f"\n{axis}/{ct}:")
            for r in s.itertuples():
                print(f"   NES={r.NES:+.2f} FDR={r.FDR_q:.3f}  {r.Term[:64]}  [{r.library}]")
print("\nDONE (N2.2 + N2.3).")
