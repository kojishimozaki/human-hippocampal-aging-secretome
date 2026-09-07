#!/usr/bin/env python
"""Read-only verification of the N0 DE outputs (no recompute, no anchor).

(1) DETERMINISTIC orientation/sign cross-check: my niche YA-vs-HA log2FC (de_..._aging) must
    reproduce the frozen ~arm+grp refit (results/de/aging_arm_sensitivity_perhit.csv arm_lfc),
    because both are the SAME design (~arm+grp), SAME contrast ["grp","HA","YA"], SAME YA/HA
    donors. Exact match (to 3-dp rounding) validates the engine AND the aging sign convention.
(2) Overlap of my niche aging secretome padj<0.1 hits with the frozen 60-hit signature.
(3) AD-vs-HA magnitude sanity (the large hit counts): padj distribution + log2FC range +
    top |log2FC| hits, to rule out a degenerate model and frame the magnitude honestly.
"""
import os
import numpy as np, pandas as pd
P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
NICHE = ["Astro", "Micro", "Oligo", "OPC", "Endo"]

aging = pd.read_csv(f"{P}/results/de/de_GSE268609_neuron_aging_per_celltype.csv")
advha = pd.read_csv(f"{P}/results/de/de_GSE268609_neuron_ADvHA_per_celltype.csv")
perhit = pd.read_csv(f"{P}/results/de/aging_arm_sensitivity_perhit.csv")   # frozen ~arm+grp refit
frozen = pd.read_csv(f"{P}/results/validation/frozen_primary_signature.csv")
sec = set(pd.read_csv(f"{P}/refs/secretome_union.csv").gene)

print("=== (1) niche aging log2FC vs frozen ~arm+grp refit (arm_lfc) — must match ===")
m = aging.merge(perhit[["celltype", "gene", "arm_lfc"]], on=["celltype", "gene"], how="inner")
m = m[m.celltype.isin(NICHE)]
m["absdiff"] = (m.log2FoldChange - m.arm_lfc).abs()
m["sign_match"] = np.sign(m.log2FoldChange) == np.sign(m.arm_lfc)
print(f"frozen hits matched: {len(m)} (perhit rows: {len(perhit)})")
print(f"max |log2FC - arm_lfc|: {m.absdiff.max():.5f}  (expect <0.005 = 3-dp rounding of arm_lfc)")
print(f"mean |diff|: {m.absdiff.mean():.5f} | sign match: {int(m.sign_match.sum())}/{len(m)} | "
      f"Pearson r: {m.log2FoldChange.corr(m.arm_lfc):.5f}")
worst = m.nlargest(3, "absdiff")[["celltype", "gene", "log2FoldChange", "arm_lfc", "absdiff"]]
print("largest discrepancies:\n" + worst.to_string(index=False))
VERDICT1 = "PASS" if (m.absdiff.max() < 0.005 and m.sign_match.all()) else "FAIL — INVESTIGATE"
print(f">>> orientation/engine check: {VERDICT1}")

print("\n=== (2) my niche aging secretome padj<0.1 hits vs the frozen 60-hit signature ===")
my_niche_hits = aging[(aging.celltype.isin(NICHE)) & (aging.padj < 0.1) & (aging.gene.isin(sec))]
fz = frozen[["celltype", "gene"]].drop_duplicates()
mh = my_niche_hits[["celltype", "gene"]].drop_duplicates()
inter = mh.merge(fz, on=["celltype", "gene"], how="inner")
print(f"frozen signature size: {len(fz)} | my niche aging secretome padj<0.1: {len(mh)}")
print(f"overlap (same celltype+gene): {len(inter)}  "
      f"(arm-robustness reported 42/60 frozen hits retain padj<0.1 under ~arm+grp)")
print(f"my hits NOT in frozen (new under ~arm+grp): {len(mh) - len(inter)} | "
      f"frozen hits I lose at padj<0.1: {len(fz) - len(inter)}")

print("\n=== (3) AD-vs-HA magnitude sanity (large counts) ===")
for ct in ["DG_GC", "Astro", "Endo"]:
    d = advha[advha.celltype == ct]
    nn = d.padj.notna()
    sig = d[d.padj < 0.1]
    print(f"\n-- {ct}: genes={len(d)}, padj notna={int(nn.sum())}, padj<0.1={len(sig)} "
          f"({100*len(sig)/max(1,int(nn.sum())):.0f}% of tested)")
    print(f"   padj quantiles (notna): " +
          ", ".join(f"{q}={d.padj.quantile(q):.3g}" for q in [0.01, 0.05, 0.25, 0.5]))
    print(f"   |log2FC| all: median={d.log2FoldChange.abs().median():.2f}, "
          f"max={d.log2FoldChange.abs().max():.2f}; among padj<0.1: "
          f"median |lfc|={sig.log2FoldChange.abs().median():.2f}")
    top = sig.reindex(sig.log2FoldChange.abs().sort_values(ascending=False).index).head(5)
    print("   top |log2FC| hits (gene/lfc/baseMean/padj):")
    for r in top.itertuples():
        flag = " [secretome]" if r.gene in sec else ""
        print(f"     {r.gene:>10}  lfc={r.log2FoldChange:+.2f}  base={r.baseMean:.1f}  padj={r.padj:.2g}{flag}")

print("\n=== AD-vs-HA arm imbalance reminder ===")
print("AD 6WH/4DG vs HA 4WH/5DG: arm imbalanced on this axis (modeled via ~arm+grp). "
      "PMI gap (AD 5.65 vs HA 7.67 h) is NOT adjusted (SOFT PMI unreliable for AD; 04 caveat).")
