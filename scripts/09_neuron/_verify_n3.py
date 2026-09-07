#!/usr/bin/env python
"""Read-only check of L3-B: is the ligand-activity DISCRIMINATIVE for the niche candidate ligands,
or does the perm null pass ~everything (the weak-calibration artifact the plan warned about)?
Compare candidate-ligand pass-rate vs ALL-ligand pass-rate per (receiver, response). Also confirm the
receptor expression/age-regulation summary."""
import os
import numpy as np, pandas as pd
P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
la = pd.read_csv(f"{P}/results/regulatory/neuron_nichenet_ligand_activity.csv")
print("=== ligand-activity DISCRIMINATION: candidate vs ALL-ligand perm_padj<0.10 pass rate ===")
print(f"(non-discriminative if frac_pass_ALL is high and ~= frac_pass_candidate; NA pearson = zero-variance)")
for (rc, resp), g in la.groupby(["receiver", "response"]):
    nNA = int(g.pearson.isna().sum())
    allp = g.perm_padj < 0.10
    cand = g[g.is_candidate]
    candp = cand.perm_padj < 0.10
    print(f"\n-- {rc} / {resp} --  (n_ligands={len(g)}, NA_pearson={nNA})")
    print(f"   ALL ligands  perm_padj<0.10: {int(allp.sum())}/{len(g)} ({100*allp.mean():.0f}%)")
    print(f"   CANDIDATES   perm_padj<0.10: {int(candp.sum())}/{len(cand)} ({100*candp.mean():.0f}%)")
    # enrichment of candidates among passers (Fisher-ish): is candidate pass-rate > all pass-rate?
    base = allp.mean()
    print(f"   -> candidate pass-rate {100*candp.mean():.0f}% vs background {100*base:.0f}%  "
          f"=> {'DISCRIMINATIVE (candidates enriched)' if candp.mean() > base + 0.15 else 'NON-discriminative (candidates ~= background)'}")

print("\n=== receptors: expression + age-regulation (the ambient-vs-response discriminator) ===")
rec = pd.read_csv(f"{P}/results/regulatory/neuron_nichenet_receptors.csv")
for rc, g in rec.groupby("receiver"):
    ur = g.drop_duplicates("receptor")
    n_expr = int(ur.expressed.sum()); n_areg = int(ur.age_regulated.sum())
    # any receptor with the most extreme (smallest padj) regulation, even if n.s.
    ex = ur[ur.expressed].copy()
    best = ex.reindex(ex.padj.fillna(1).sort_values().index).head(3)
    print(f"\n-- {rc}: {len(ur)} unique receptors | expressed {n_expr} | age-regulated(padj<0.1) {n_areg}")
    print("   most-regulated expressed receptors (smallest padj, even if n.s.):")
    for r in best.itertuples():
        print(f"     {r.receptor}: lfc={r.lfc:+.2f} padj={r.padj}")
