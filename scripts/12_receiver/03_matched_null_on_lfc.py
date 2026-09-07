#!/usr/bin/env python
"""receiver-map STEP 03 — PRIMARY: matched-null on log2FC (ATAC/methyl analogue).

Pre-registration: docs/RECEIVER_MAP_PLAN.md (sec 5,7,9), tag prereg/receiver-map-v1.  env: bio.
seed=42, 2000 permutations. Operates on the locked ~arm+grp DE table; NO recompute.

Per compartment C:
  foreground = detectable frozen-76 receptors (in-table, non-NA padj -> have lfc).
  background = non-receptor detectable genes in C, EXCLUDING all 76, matched on log10(baseMean+1) decile.
  primary statistic = mean SIGNED receptor lfc, oriented age-UP (all 22 ligands are age-UP -> the
    coordination hypothesis is receptors age-UP -> a single expected direction; the oriented statistic
    is simply mean(signed lfc) tested one-sided positive). one-sided p_up + two-sided p reported.
  secondary = mean |lfc| (age dynamism), two-sided. DESCRIPTIVE ONLY (never an activation claim).
  not-testable if < 3 detectable receptors.
GO bar (confirmatory family of 8): >=3 receptors AND mean signed lfc>0 AND BH-q(p_up)<0.1 AND LOO-robust.
Neurogenic (NSC/Neuroblast/Immature): exploratory; nominal p_up + all-compartment BH-q reference only.
"""
import os, json, numpy as np, pandas as pd

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
os.chdir(PROJ)
OUT = "results/receiver"
DE = "results/de/de_GSE268609_neuron_aging_per_celltype.csv"
NPERM = 2000; SEED = 42; NBIN = 10; MIN_REC = 3

CONFIRMATORY = ["Astro", "Micro", "Oligo", "OPC", "Endo", "DG_GC", "CA_ExN", "InN"]
NEUROGENIC = ["NSC", "Neuroblast", "Immature"]
COMPARTMENTS = CONFIRMATORY + NEUROGENIC

rec76set = set(pd.read_csv(f"{OUT}/_universe76_from_manifest.csv")["receptor"])
de = pd.read_csv(DE, usecols=["gene", "baseMean", "log2FoldChange", "padj", "celltype"])
rng = np.random.default_rng(SEED)


def bh(pvals):
    p = np.asarray(pvals, float); n = len(p); order = np.argsort(p)
    q = np.empty(n); prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        prev = min(prev, p[i] * n / (rank + 1))
        q[i] = prev
    return q


def matched_null(fg_lfc, fg_bm, pool_lfc, pool_bm):
    pool_logbm = np.log10(pool_bm + 1.0)
    edges = np.quantile(pool_logbm, np.linspace(0, 1, NBIN + 1)); edges[0] = -np.inf; edges[-1] = np.inf
    pool_bin = np.clip(np.digitize(pool_logbm, edges[1:-1], right=False), 0, NBIN - 1)
    bin_to_idx = {b: np.where(pool_bin == b)[0] for b in range(NBIN)}
    fg_bin = np.clip(np.digitize(np.log10(fg_bm + 1.0), edges[1:-1], right=False), 0, NBIN - 1)
    counts = {b: int((fg_bin == b).sum()) for b in range(NBIN) if (fg_bin == b).sum() > 0}
    ns = np.empty(NPERM); na = np.empty(NPERM)
    for p in range(NPERM):
        picks = []
        for b, k in counts.items():
            idx = bin_to_idx[b]
            picks.append(rng.choice(idx, size=k, replace=(len(idx) < k)))
        sel = pool_lfc[np.concatenate(picks)]
        ns[p] = sel.mean(); na[p] = np.abs(sel).mean()
    return ns, na


res = []
for ct in COMPARTMENTS:
    d = de[de.celltype == ct]
    det = d[d["padj"].notna()].copy()                       # detectable (have lfc)
    fg = det[det.gene.isin(rec76set)]
    n_fg = len(fg)
    base = dict(compartment=ct, family=("confirmatory" if ct in CONFIRMATORY else "neurogenic"),
                n_detectable_receptors=n_fg)
    if n_fg < MIN_REC:
        res.append({**base, "call": "NOT_TESTABLE", "obs_mean_signed_lfc": None})
        continue
    pool = det[~det.gene.isin(rec76set)]
    fg_lfc = fg["log2FoldChange"].to_numpy(); fg_bm = fg["baseMean"].to_numpy()
    obs_signed = float(fg_lfc.mean()); obs_abs = float(np.abs(fg_lfc).mean())
    ns, na = matched_null(fg_lfc, fg_bm, pool["log2FoldChange"].to_numpy(), pool["baseMean"].to_numpy())
    p_up = (1 + int((ns >= obs_signed).sum())) / (NPERM + 1)
    p_dn = (1 + int((ns <= obs_signed).sum())) / (NPERM + 1)
    p_two = min(1.0, 2 * min(p_up, p_dn))
    pa_hi = (1 + int((na >= obs_abs).sum())) / (NPERM + 1)
    pa_lo = (1 + int((na <= obs_abs).sum())) / (NPERM + 1)
    p_dyn = min(1.0, 2 * min(pa_hi, pa_lo))
    res.append({**base, "call": "pending",
                "obs_mean_signed_lfc": round(obs_signed, 4),
                "null_mean_signed_lfc": round(float(ns.mean()), 4),
                "effect_signed": round(obs_signed - float(ns.mean()), 4),
                "p_up_oneside": round(p_up, 4), "p_two_signed": round(p_two, 4),
                "obs_mean_abs_lfc": round(obs_abs, 4), "null_mean_abs_lfc": round(float(na.mean()), 4),
                "p_dynamism_two": round(p_dyn, 4)})

R = pd.DataFrame(res)

# BH across the 8 confirmatory testable compartments (on p_up)
conf = R[(R.family == "confirmatory") & (R.call != "NOT_TESTABLE")].copy()
R["bh_q_confirmatory"] = np.nan
if len(conf):
    q = bh(conf["p_up_oneside"].to_numpy())
    R.loc[conf.index, "bh_q_confirmatory"] = np.round(q, 4)

# all-compartment BH-q (reference for neurogenic)
testable = R[R.call != "NOT_TESTABLE"].copy()
R["bh_q_all_ref"] = np.nan
if len(testable):
    R.loc[testable.index, "bh_q_all_ref"] = np.round(bh(testable["p_up_oneside"].to_numpy()), 4)

# GO calls (LOO still pending for any confirmatory passer)
def call_row(r):
    if r["call"] == "NOT_TESTABLE":
        return "NOT_TESTABLE"
    if r["family"] == "confirmatory":
        passed = (r["n_detectable_receptors"] >= MIN_REC and r["obs_mean_signed_lfc"] > 0
                  and pd.notna(r["bh_q_confirmatory"]) and r["bh_q_confirmatory"] < 0.1)
        return "PASS_pending_LOO" if passed else "FAIL"
    else:
        # neurogenic: exploratory; never PASS on nominal alone (needs depth-null + donor-composite, sec 9.2)
        return "EXPLORATORY"

R["call"] = R.apply(call_row, axis=1)
R.to_csv(f"{OUT}/matched_null_per_compartment.csv", index=False)

print("STEP 03 matched-null-on-lfc (PRIMARY; seed=42, 2000x, oriented age-UP)")
cols = ["compartment", "family", "n_detectable_receptors", "obs_mean_signed_lfc", "effect_signed",
        "p_up_oneside", "p_two_signed", "bh_q_confirmatory", "bh_q_all_ref", "p_dynamism_two", "call"]
print(R[cols].to_string(index=False))
n_pass = int((R.call == "PASS_pending_LOO").sum())
print(f"\nconfirmatory passers (pending LOO): {n_pass}")
print(f"wrote: {OUT}/matched_null_per_compartment.csv")
if n_pass == 0:
    print("\nNO-GO (pre-authorized): The frozen age-UP ligand set is not accompanied by coordinated "
          "age-UP receptor regulation in candidate receiver compartments.")
