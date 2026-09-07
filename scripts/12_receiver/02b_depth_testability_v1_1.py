#!/usr/bin/env python
"""receiver-map STEP 02b (v1.1) — decomposed depth/testability analysis.

Pre-registration: docs/RECEIVER_MAP_PLAN.md AMENDMENT v1.1, tag prereg/receiver-map-v1.1.  env: bio.
seed=42, 2000 permutations. Refines STEP 02 only; STEP 03 (primary) is unchanged.

Decomposition per compartment (replaces the conflated v1 "detectable/76 vs in-dds pool" statistic):
  A. Coverage / table inclusion : n_in_table / 76  (descriptive; rare low-depth compartments include
     fewer genes overall -> contextualised by n_genes_in_table).
  B. In-table detectability      : among in-table receptors, n_detectable / n_in_table  (direction-blind).
  C. Expression/depth-matched null WITHIN the in-table universe : in-table receptors vs baseMean-matched
     in-table NON-receptor genes. Absent receptors are EXCLUDED from the null (handled in A). For compartments
     with SATURATED detectability (background-pool detectable fraction >= 0.98) layer C has no resolving power
     and is reported "not_interpretable_saturated" rather than a misleading p-value.
v1 result (results/receiver/depth_matched_null.csv) is preserved unchanged as the original record.
"""
import os, numpy as np, pandas as pd

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
os.chdir(PROJ)
OUT = "results/receiver"
DE = "results/de/de_GSE268609_neuron_aging_per_celltype.csv"
NPERM = 2000; SEED = 42; NBIN = 10; SAT = 0.98

CONFIRMATORY = ["Astro", "Micro", "Oligo", "OPC", "Endo", "DG_GC", "CA_ExN", "InN"]
NEUROGENIC = ["NSC", "Neuroblast", "Immature"]
COMPARTMENTS = CONFIRMATORY + NEUROGENIC

rec76set = set(pd.read_csv(f"{OUT}/_universe76_from_manifest.csv")["receptor"])
de = pd.read_csv(DE, usecols=["gene", "baseMean", "padj", "celltype"])
rng = np.random.default_rng(SEED)

res = []
for ct in COMPARTMENTS:
    d = de[de.celltype == ct].copy()
    d["detect"] = d["padj"].notna().to_numpy()
    n_genes = len(d)
    rec = d[d.gene.isin(rec76set)]
    pool = d[~d.gene.isin(rec76set)]                       # in-table non-receptor genes
    pool_detect_frac = float(pool["detect"].mean())

    # A. coverage
    n_in_table = len(rec)
    # B. in-table detectability
    n_detect_intab = int(rec["detect"].sum())
    in_table_detect_frac = (n_detect_intab / n_in_table) if n_in_table else np.nan

    row = dict(compartment=ct, family=("confirmatory" if ct in CONFIRMATORY else "neurogenic"),
               n_genes_in_table=n_genes, global_detectable_frac=round(float(d["detect"].mean()), 4),
               # A
               n_in_table=n_in_table, coverage_in_table_of_76=round(n_in_table / 76, 4),
               n_absent_from_table=76 - n_in_table,
               # B
               n_detectable_in_table=n_detect_intab,
               in_table_detectable_frac=round(in_table_detect_frac, 4) if n_in_table else None,
               pool_detectable_frac=round(pool_detect_frac, 4))

    # C. matched-null WITHIN the in-table universe (only if not saturated and enough eligible receptors)
    if pool_detect_frac >= SAT:
        row.update(layerC="not_interpretable_saturated", obs_detect_intab=n_detect_intab,
                   matched_null_mean=None, effect_obs_minus_null=None, p_low_depleted=None)
    elif n_in_table < 3:
        row.update(layerC="not_testable_few_receptors", obs_detect_intab=n_detect_intab,
                   matched_null_mean=None, effect_obs_minus_null=None, p_low_depleted=None)
    else:
        pool_logbm = np.log10(pool["baseMean"].to_numpy() + 1.0)
        pool_detect = pool["detect"].to_numpy()
        edges = np.quantile(pool_logbm, np.linspace(0, 1, NBIN + 1)); edges[0] = -np.inf; edges[-1] = np.inf
        pool_bin = np.clip(np.digitize(pool_logbm, edges[1:-1], right=False), 0, NBIN - 1)
        bin_to_idx = {b: np.where(pool_bin == b)[0] for b in range(NBIN)}
        rec_logbm = np.log10(rec["baseMean"].to_numpy() + 1.0)
        rec_bin = np.clip(np.digitize(rec_logbm, edges[1:-1], right=False), 0, NBIN - 1)
        counts = {b: int((rec_bin == b).sum()) for b in range(NBIN) if (rec_bin == b).sum() > 0}
        null = np.empty(NPERM, dtype=int)
        for p in range(NPERM):
            tot = 0
            for b, k in counts.items():
                idx = bin_to_idx[b]
                pick = rng.choice(idx, size=k, replace=(len(idx) < k))
                tot += int(pool_detect[pick].sum())
            null[p] = tot
        obs = n_detect_intab
        p_low = (1 + int((null <= obs).sum())) / (NPERM + 1)
        row.update(layerC="interpretable", obs_detect_intab=obs,
                   matched_null_mean=round(float(null.mean()), 2),
                   effect_obs_minus_null=round(obs - float(null.mean()), 2),
                   p_low_depleted=round(p_low, 4))
    res.append(row)

R = pd.DataFrame(res)
R.to_csv(f"{OUT}/depth_testability_v1_1.csv", index=False)

show = ["compartment", "family", "n_in_table", "coverage_in_table_of_76", "n_detectable_in_table",
        "in_table_detectable_frac", "pool_detectable_frac", "layerC", "obs_detect_intab",
        "matched_null_mean", "effect_obs_minus_null", "p_low_depleted"]
print("STEP 02b (v1.1) decomposed depth/testability  [A coverage | B in-table detect | C in-table matched-null]")
print(R[show].to_string(index=False))
print("\nNEUROGENIC focus:")
print(R[R.family == "neurogenic"][show].to_string(index=False))
print(f"\nwrote: {OUT}/depth_testability_v1_1.csv")
print("rule: layer C 'not_interpretable_saturated' when in-table background detectable frac >= 0.98 "
      "(no resolving power); p_low_depleted small => receptors detected LESS than baseMean-matched "
      "in-table genes (receiver-intrinsic); p_low ~0.5 => depth/expression-explained.")
