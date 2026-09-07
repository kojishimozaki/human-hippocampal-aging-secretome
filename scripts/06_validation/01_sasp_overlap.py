#!/usr/bin/env python
"""Move #1 (orthogonal validation): is the 268609 niche aging-secretome signature
a bona-fide SASP? Tests overlap/directional concordance with
  (a) SenMayo (Saul 2022, transcriptomic curated SASP; refs/senmayo.txt)
  (b) SASP Atlas (Basisty 2020, mass-spec secreted senescence proteins; refs/sasp_atlas_core.txt)
Output: results/validation/sasp_overlap.csv
"""
import os
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu, fisher_exact

P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
sets = {
    "SenMayo": set(open(f"{P}/refs/senmayo.txt").read().split()),
    "SASP_Atlas": set(open(f"{P}/refs/sasp_atlas_core.txt").read().split()),
}
de = pd.read_csv(f"{P}/results/de/de_GSE268609_per_celltype.csv")   # HA vs YA, genome-wide, symbols
secr = set(pd.read_csv(f"{P}/refs/secretome_union.csv").gene)       # the analysis universe

# Audit fix (2026-06-03): the SASP claim is "the niche *secretome* aging signature is a
# SASP", so the correct background for both the directional and the enrichment test is the
# SECRETOME universe (what Methods states), not all genome-wide tested genes. The genome-wide
# enrichment (OR up to 13) is inflated because secretome genes are intrinsically SASP-like.
# We therefore LEAD with the directional concordance (which survives restriction) and report
# the Fisher enrichment for BOTH universes so the universe-dependence is explicit. See
# Limitations: enrichment is positioned descriptively, not as a load-bearing enrichment claim.

def fisher_enrich(d, gs, universe_genes):
    """OR/p for SASP-set over-representation among sig-UP genes, within `universe_genes`."""
    d = d[d.gene.isin(universe_genes)]
    tested = set(d.gene); smt = gs & tested
    sigup = set(d[(d.padj < 0.1) & (d.log2FoldChange > 0)].gene)
    n11 = len(smt & sigup); n10 = len(smt - sigup); n01 = len(sigup - smt)
    n00 = len(tested) - n11 - n10 - n01
    orr, p = fisher_exact([[n11, n10], [n01, n00]], alternative="greater")
    return n11, orr, p

rows = []
for sname, gs in sets.items():
    for ct in ["Astro", "Micro", "Oligo", "OPC", "Endo"]:
        d = de[(de.celltype == ct) & de.padj.notna()].copy()
        # PRIMARY directional test (Methods-matching): within the SECRETOME universe,
        # are SASP-set secretome genes more up-regulated than non-SASP secretome genes?
        dsec = d[d.gene.isin(secr)]
        sm_s, rest_s = dsec[dsec.gene.isin(gs)], dsec[~dsec.gene.isin(gs)]
        if len(sm_s) < 5:
            continue
        _, dir_p = mannwhitneyu(sm_s.log2FoldChange, rest_s.log2FoldChange, alternative="greater")
        up_frac = float((sm_s.log2FoldChange > 0).mean())
        # reference directional test on the genome-wide background (for transparency)
        sm_g, rest_g = d[d.gene.isin(gs)], d[~d.gene.isin(gs)]
        _, dir_p_genome = mannwhitneyu(sm_g.log2FoldChange, rest_g.log2FoldChange, alternative="greater")
        # Fisher enrichment under BOTH universes
        n11_s, or_s, p_s = fisher_enrich(d, gs, secr)
        n11_g, or_g, p_g = fisher_enrich(d, gs, set(de.gene))
        rows.append(dict(refset=sname, celltype=ct, n_secr=len(sm_s), up_frac=round(up_frac, 3),
                         median_lfc=round(float(sm_s.log2FoldChange.median()), 3),
                         dir_p=dir_p, dir_p_genome=dir_p_genome,
                         sigUP_overlap_secr=n11_s, OR_secr=round(or_s, 2), p_secr=p_s,
                         sigUP_overlap_genome=n11_g, OR_genome=round(or_g, 2), p_genome=p_g))
res = pd.DataFrame(rows)
os.makedirs(f"{P}/results/validation", exist_ok=True)
res.to_csv(f"{P}/results/validation/sasp_overlap.csv", index=False)
pd.set_option("display.width", 200)
print(res.to_string(index=False))
print("\nLEAD (Methods-matching) = directional `dir_p` + `up_frac` on the SECRETOME universe.")
print("Enrichment: `OR_secr`/`p_secr` (correct universe) vs `OR_genome`/`p_genome` (inflated).")
print("wrote results/validation/sasp_overlap.csv")
