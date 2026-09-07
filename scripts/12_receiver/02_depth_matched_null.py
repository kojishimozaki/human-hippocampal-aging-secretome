#!/usr/bin/env python
"""receiver-map STEP 02 — depth/expression-matched detectability null (neurogenic repertoire control).

Pre-registration: docs/RECEIVER_MAP_PLAN.md (sec 6), tag prereg/receiver-map-v1.  env: bio.
seed=42, 2000 permutations.

Question (sec 6): is the receptor set detected LESS than genes of comparable expression depth in the
same compartment? Guards the NSC low-repertoire observation against the depth/rarity confound.

Statistic: observed = # of the frozen 76 detectable (in-table & non-NA padj) in compartment C.
Null: for each of the 76 receptors, draw one matched gene from C's full in-table gene pool, matched by
log10(baseMean+1) decile (receptors absent from the table -> baseMean 0 -> lowest decile; this is
conservative, slightly OVER-stating receptor expression). Count how many drawn genes are detectable.
2000x. Empirical p_low = receptors LESS detectable than matched (the depletion direction). An unmatched
random-gene null is also reported as a reference. Selection is blind to receptor age-regulation.
"""
import os, json, numpy as np, pandas as pd

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
os.chdir(PROJ)
OUT = "results/receiver"
DE = "results/de/de_GSE268609_neuron_aging_per_celltype.csv"
NPERM = 2000; SEED = 42; NBIN = 10

CONFIRMATORY = ["Astro", "Micro", "Oligo", "OPC", "Endo", "DG_GC", "CA_ExN", "InN"]
NEUROGENIC = ["NSC", "Neuroblast", "Immature"]
COMPARTMENTS = CONFIRMATORY + NEUROGENIC

rec76 = pd.read_csv(f"{OUT}/_universe76_from_manifest.csv")["receptor"].tolist()
rec76set = set(rec76)
de = pd.read_csv(DE, usecols=["gene", "baseMean", "padj", "celltype"])
rng = np.random.default_rng(SEED)


def decile_bins(logbm_pool):
    # bin edges from the pool distribution; digitize returns 0..NBIN
    qs = np.quantile(logbm_pool, np.linspace(0, 1, NBIN + 1))
    qs[0] = -np.inf; qs[-1] = np.inf
    return qs


def assign_bin(vals, edges):
    return np.clip(np.digitize(vals, edges[1:-1], right=False), 0, NBIN - 1)


res = []
for ct in COMPARTMENTS:
    d = de[de.celltype == ct].copy()
    d["detect"] = d["padj"].notna().to_numpy()
    pool = d[~d.gene.isin(rec76set)].copy()          # background = non-receptor in-table genes
    pool_logbm = np.log10(pool["baseMean"].to_numpy() + 1.0)
    pool_detect = pool["detect"].to_numpy()
    edges = decile_bins(pool_logbm)
    pool_bin = assign_bin(pool_logbm, edges)
    bin_to_idx = {b: np.where(pool_bin == b)[0] for b in range(NBIN)}

    # receptor baseMean (0 if absent) + observed detectability
    drow = d.set_index("gene")
    rec_bm = np.array([float(drow.loc[r, "baseMean"]) if r in drow.index else 0.0 for r in rec76])
    rec_det = np.array([bool(drow.loc[r, "detect"]) if r in drow.index else False for r in rec76])
    obs_detect = int(rec_det.sum())
    rec_bin = assign_bin(np.log10(rec_bm + 1.0), edges)
    bin_counts = {b: int((rec_bin == b).sum()) for b in range(NBIN) if (rec_bin == b).sum() > 0}

    # matched null
    null_detect = np.empty(NPERM, dtype=int)
    for p in range(NPERM):
        tot = 0
        for b, k in bin_counts.items():
            idx = bin_to_idx[b]
            if len(idx) >= k:
                pick = rng.choice(idx, size=k, replace=False)
            else:  # tiny bin: nearest non-empty fallback (should not occur; pools are large)
                pick = rng.choice(np.where(pool_bin >= 0)[0], size=k, replace=False)
            tot += int(pool_detect[pick].sum())
        null_detect[p] = tot

    # unmatched reference null
    n76 = 76
    unm = np.array([int(pool_detect[rng.choice(len(pool_detect), size=n76, replace=False)].sum())
                    for _ in range(NPERM)])

    p_low = (1 + int((null_detect <= obs_detect).sum())) / (NPERM + 1)
    p_high = (1 + int((null_detect >= obs_detect).sum())) / (NPERM + 1)
    p_two = min(1.0, 2 * min(p_low, p_high))
    res.append(dict(
        compartment=ct, family=("confirmatory" if ct in CONFIRMATORY else "neurogenic"),
        obs_detectable=obs_detect, n_universe=76,
        matched_null_mean=round(float(null_detect.mean()), 2),
        matched_null_sd=round(float(null_detect.std()), 2),
        effect_obs_minus_null=round(obs_detect - float(null_detect.mean()), 2),
        p_low_depleted=round(p_low, 4), p_two=round(p_two, 4),
        unmatched_null_mean=round(float(unm.mean()), 2),
    ))

R = pd.DataFrame(res)
R.to_csv(f"{OUT}/depth_matched_null.csv", index=False)
print("STEP 02 depth-matched detectability null (seed=42, 2000x)")
print(R.to_string(index=False))
print("\nNEUROGENIC focus:")
print(R[R.family == "neurogenic"].to_string(index=False))
print(f"\nwrote: {OUT}/depth_matched_null.csv")
print("interpretation: p_low_depleted small => receptors detected LESS than depth-matched genes "
      "(receiver-intrinsic depletion); p_low ~0.5 => depth-explained (bounded caveat, NOT biology).")
