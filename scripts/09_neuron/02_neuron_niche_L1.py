#!/usr/bin/env python
"""N1 (Layer 1) — is the niche aging-secretome niche-specific? niche <-> neuron secretome comparison.

Pre-registered in audit_log/2026-06-18_neuron_niche/RESOLUTION.md (frozen before this ran). Reuses the
EXACT statistical machinery of scripts/03_de/05_external_replication.py (imported, not reimplemented) so
the niche<->neuron test uses the same permutation null / concordance design as the validated niche
external replication. Donor-level effect sizes only; effect size (rho / concordance) is the result, the
gene-count Spearman p is de-emphasized, the 2,000x permutation null gives significance. seed 42.

Inputs (N0 outputs, ~arm+grp, same donors for niche AND neuron):
  results/de/de_GSE268609_neuron_aging_per_celltype.csv    YA-vs-HA (niche+neuron)
  results/de/de_GSE268609_neuron_ADvHA_per_celltype.csv    AD-vs-HA (niche+neuron)
  results/validation/frozen_primary_signature.csv          frozen 60 niche hits (anchored test)
  results/validation/frozen_primary_secretome_logfc.csv    frozen full niche secretome lfc (background)
  refs/secretome_union.csv                                  secretome universe + category

Metrics: (1) de-novo neuron secretome hits; (2) all-pairwise rho over shared secretome (niche-niche /
neuron-neuron / niche-neuron), both axes; (3) anchored frozen-60 concordance vs neuron own background
(aging); (4) category-level secretome log2FC.

Outputs: results/validation/neuron_niche_L1_{denovo_hits,rho_aging,rho_ADvHA,anchored60,category}.csv
"""
import os, itertools, importlib.util
import numpy as np, pandas as pd
from scipy import stats

P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
NICHE  = ["Astro", "Micro", "Oligo", "OPC", "Endo"]
NEURON = ["DG_GC", "CA_ExN", "InN"]
ALL8   = NICHE + NEURON

# --- import 05's exact functions (single source of truth for the stats) ---
spec = importlib.util.spec_from_file_location("rep05", f"{P}/scripts/03_de/05_external_replication.py")
rep05 = importlib.util.module_from_spec(spec); spec.loader.exec_module(rep05)
onesided_spearman_gt0 = rep05.onesided_spearman_gt0
perm_spearman_null = rep05.perm_spearman_null
concordance = rep05.concordance
NPERM = rep05.NPERM
rng = np.random.default_rng(rep05.SEED)
print(f"reusing 05_external_replication: NPERM={NPERM}, seed={rep05.SEED}")

# --- data ---
secdf = pd.read_csv(f"{P}/refs/secretome_union.csv")
sec = set(secdf.gene); cat_map = dict(zip(secdf.gene, secdf.category))
aging = pd.read_csv(f"{P}/results/de/de_GSE268609_neuron_aging_per_celltype.csv")
advha = pd.read_csv(f"{P}/results/de/de_GSE268609_neuron_ADvHA_per_celltype.csv")
frozen_hits = pd.read_csv(f"{P}/results/validation/frozen_primary_signature.csv")     # 60 hits
frozen_full = pd.read_csv(f"{P}/results/validation/frozen_primary_secretome_logfc.csv")
print(f"secretome categories: {secdf.category.value_counts().to_dict()}")

def sec_lfc(de, ct):
    """secretome-only (gene -> log2FC) for one celltype, non-null lfc."""
    d = de[(de.celltype == ct) & (de.gene.isin(sec))][["gene", "log2FoldChange"]].dropna()
    return d

# ===== (1) de-novo neuron secretome program =====================================================
print("\n=== (1) de-novo neuron secretome hits (no niche anchoring) ===")
denovo = []
for axis, de in [("aging", aging), ("ADvHA", advha)]:
    for ct in NEURON:
        d = de[(de.celltype == ct) & (de.gene.isin(sec))].copy()
        for thr in (0.1, 0.2):
            h = d[d.padj < thr]
            for r in h.itertuples():
                denovo.append(dict(axis=axis, celltype=ct, gene=r.gene, log2FoldChange=round(r.log2FoldChange, 3),
                                   padj=r.padj, direction=("UP" if r.log2FoldChange > 0 else "DOWN"),
                                   threshold=f"padj<{thr}", category=cat_map.get(r.gene, "NA")))
den = pd.DataFrame(denovo).drop_duplicates(subset=["axis", "celltype", "gene", "threshold"])
# collapse: a gene passing 0.1 also passes 0.2; tag the tightest threshold
den = den.sort_values(["axis", "celltype", "threshold"]).drop_duplicates(["axis", "celltype", "gene"], keep="first")
den.to_csv(f"{P}/results/validation/neuron_niche_L1_denovo_hits.csv", index=False)
for axis in ("aging", "ADvHA"):
    for ct in NEURON:
        sub = den[(den.axis == axis) & (den.celltype == ct)]
        n01 = (sub.threshold == "padj<0.1").sum(); n02 = len(sub)
        genes = ", ".join(sub[sub.threshold == "padj<0.1"].gene.tolist()[:12]) or "(none)"
        print(f"  {axis:6s}/{ct:7s}: padj<0.1={n01}, <0.2={n02} | hits: {genes}")

# ===== (2) all-pairwise rho over shared secretome ===============================================
def pairwise_rho(de, axis):
    rows = []
    lfcs = {ct: sec_lfc(de, ct) for ct in ALL8}
    for a, b in itertools.combinations(ALL8, 2):
        m = lfcs[a].merge(lfcs[b], on="gene", suffixes=("_a", "_b"))
        if len(m) < 5:
            continue
        rho, p1, n = onesided_spearman_gt0(m.log2FoldChange_a.values, m.log2FoldChange_b.values)
        nullp = perm_spearman_null(m.log2FoldChange_a.values, m.log2FoldChange_b.values, rho, rng)
        pt = ("niche-niche" if a in NICHE and b in NICHE else
              "neuron-neuron" if a in NEURON and b in NEURON else "niche-neuron")
        rows.append(dict(axis=axis, a=a, b=b, pair_type=pt, rho=round(rho, 3),
                         p_genecount=p1, null_p=nullp, n_shared=n,
                         null_sig=(nullp < 0.05 if pd.notna(nullp) else False)))
    return pd.DataFrame(rows)

print("\n=== (2) pairwise Spearman rho over shared secretome (effect size is the result) ===")
rho_tabs = {}
for axis, de in [("aging", aging), ("ADvHA", advha)]:
    t = pairwise_rho(de, axis)
    t.to_csv(f"{P}/results/validation/neuron_niche_L1_rho_{axis}.csv", index=False)
    rho_tabs[axis] = t
    print(f"\n--- {axis}: rho by pair type (median [min,max]; n pairs; n null_p<0.05) ---")
    for pt in ["niche-niche", "neuron-neuron", "niche-neuron"]:
        s = t[t.pair_type == pt]
        if not len(s): continue
        print(f"  {pt:14s}: rho median={s.rho.median():+.3f} [{s.rho.min():+.3f},{s.rho.max():+.3f}] "
              f"| n={len(s)} | null_p<0.05: {int(s.null_sig.sum())}/{len(s)}")
    print("  niche-neuron detail (rho | null_p | n):")
    nn = t[t.pair_type == "niche-neuron"].sort_values("rho", ascending=False)
    for r in nn.itertuples():
        print(f"     {r.a:6s} x {r.b:7s}: rho={r.rho:+.3f}  null_p={r.null_p:.3f}  n={r.n_shared}")

# ===== (3) anchored frozen-60 concordance in neurons (aging) ====================================
print("\n=== (3) anchored frozen-60 niche-hit concordance in neurons vs neuron's own background ===")
hit_keys = set(zip(frozen_hits.gene, frozen_hits.celltype))          # (gene, niche_ct) that are hits
anch_rows = []
pooled_hits = []; pooled_nh = []
for ct in NEURON:
    nl = sec_lfc(aging, ct)
    nlmap = dict(zip(nl.gene, nl.log2FoldChange))
    # hits: 60 frozen niche hits whose gene is testable in this neuron
    hc = []
    for r in frozen_hits.itertuples():
        if r.gene in nlmap:
            hc.append(np.sign(nlmap[r.gene]) == (1 if r.direction == "UP" else -1))
    # background: frozen FULL niche secretome (gene, niche_ct) NON-hits, gene testable in this neuron
    bc = []
    for r in frozen_full.itertuples():
        if (r.gene, r.celltype) in hit_keys:
            continue
        if r.gene in nlmap:
            bc.append(np.sign(nlmap[r.gene]) == np.sign(r.log2FoldChange))
    hits_ext = pd.DataFrame({"concord": pd.Series(hc, dtype=bool)})
    nonh_ext = pd.DataFrame({"concord": pd.Series(bc, dtype=bool)})
    pooled_hits.append(hits_ext); pooled_nh.append(nonh_ext)
    co = concordance(hits_ext, nonh_ext, rng)
    passed = (pd.notna(co["binom_p"]) and co["binom_p"] < 0.05 and pd.notna(co["null_p"]) and co["null_p"] < 0.05)
    anch_rows.append(dict(scope=ct, n_hit=co["n_hit"], hit_concord=round(co["hit_concord"], 3),
                          n_bg=co["n_bg"], bg_concord=round(co["bg_concord"], 3),
                          binom_p=co["binom_p"], fisher_or=round(co["fisher_or"], 2), fisher_p=co["fisher_p"],
                          null_p=co["null_p"], concordant_above_bg=passed))
# pooled over the 3 neurons
ph = pd.concat(pooled_hits, ignore_index=True); pn = pd.concat(pooled_nh, ignore_index=True)
co = concordance(ph, pn, rng)
passed = (pd.notna(co["binom_p"]) and co["binom_p"] < 0.05 and pd.notna(co["null_p"]) and co["null_p"] < 0.05)
anch_rows.append(dict(scope="POOLED_neuron", n_hit=co["n_hit"], hit_concord=round(co["hit_concord"], 3),
                      n_bg=co["n_bg"], bg_concord=round(co["bg_concord"], 3),
                      binom_p=co["binom_p"], fisher_or=round(co["fisher_or"], 2), fisher_p=co["fisher_p"],
                      null_p=co["null_p"], concordant_above_bg=passed))
anch = pd.DataFrame(anch_rows)
anch.to_csv(f"{P}/results/validation/neuron_niche_L1_anchored60.csv", index=False)
print(anch.to_string(index=False))

# ===== (4) category-level secretome log2FC ======================================================
print("\n=== (4) category-level secretome log2FC (neuropeptide + growth-factor = candidate neuron modules) ===")
cat_rows = []
for axis, de in [("aging", aging), ("ADvHA", advha)]:
    d = de[de.gene.isin(sec)].copy()
    d["category"] = d.gene.map(cat_map)
    for ct in ALL8:
        dc = d[d.celltype == ct]
        for cat, g in dc.groupby("category"):
            cat_rows.append(dict(axis=axis, celltype=ct, category=cat, n=len(g),
                                 mean_lfc=round(float(g.log2FoldChange.mean()), 3),
                                 median_lfc=round(float(g.log2FoldChange.median()), 3),
                                 n_up_padj02=int(((g.padj < 0.2) & (g.log2FoldChange > 0)).sum()),
                                 n_dn_padj02=int(((g.padj < 0.2) & (g.log2FoldChange < 0)).sum())))
cat = pd.DataFrame(cat_rows)
cat.to_csv(f"{P}/results/validation/neuron_niche_L1_category.csv", index=False)
print("aging axis, neuropeptide + growth_factor (mean log2FC by celltype):")
piv = (cat[(cat.axis == "aging") & (cat.category.isin(["neuropeptide", "growth_factor"]))]
       .pivot_table(index="celltype", columns="category", values="mean_lfc"))
if len(piv): print(piv.reindex(ALL8).to_string())

# ===== pre-registered GO/NO-GO readout (aging axis, primary) ====================================
print("\n=== PRE-REGISTERED READOUT (aging axis) ===")
ta = rho_tabs["aging"]; nn = ta[ta.pair_type == "niche-neuron"]
nnn = ta[ta.pair_type == "niche-niche"]
n_nn_sig = int(nn.null_sig.sum())
anch_neuron = anch[anch.scope.isin(NEURON + ["POOLED_neuron"])]
n_anch_pass = int(anch_neuron.concordant_above_bg.sum())
print(f"niche-neuron rho: median={nn.rho.median():+.3f}, max={nn.rho.max():+.3f}, "
      f"null_p<0.05 in {n_nn_sig}/{len(nn)} pairs")
print(f"niche-niche rho (calibration anchor): median={nnn.rho.median():+.3f} "
      f"(external GSE278576 reproduction benchmark rho~0.24)")
print(f"anchored frozen-60 concordant-above-background: {n_anch_pass}/{len(anch_neuron)} neuron scopes")
verdict = ("NICHE-SPECIFIC (neuron does NOT mirror the niche secretome)"
           if (n_nn_sig == 0 and n_anch_pass == 0) else
           "SHARED/GLOBAL (neuron mirrors the niche)" if (n_nn_sig >= len(nn) // 2 and n_anch_pass >= 2)
           else "MIXED — read pair-by-pair (report honestly, no forcing)")
print(f">>> pre-registered verdict: {verdict}")
print("\nDONE (N1 L1). STOP per protocol — await PI go for N2.")
