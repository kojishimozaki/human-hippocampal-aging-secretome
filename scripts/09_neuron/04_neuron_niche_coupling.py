#!/usr/bin/env python
"""N2.4 — donor-level niche<->neuron coupling (the matched-design strength: niche+neuron, SAME donors).

Pre-registered in audit_log/2026-06-18_neuron_niche/RESOLUTION.md. donor = unit (pseudoreplication-safe):
one score per donor per cell-pool, correlated ACROSS donors. Deterministic pseudobulk z-score signature
(NO sc.tl.score_genes hidden RNG). seed 42. Coordination, NOT causation (no causal verbs in the write-up).

Per donor (contrast donors only): donor-pseudobulk CPM(log1p) per gene, z-scored across the contrast's
donors, then mean-z over a gene set:
  niche_SASP   = frozen niche aging-UP secretome genes, scored in the donor's NICHE cells (output signal)
  neuron_resp  = inflammatory niche-response composite (NFkB/IFN/Inflammatory/Complement), in NEURON cells
  neuron_ISR / neuron_senes / neuron_recept = the other targeted sets, in NEURON cells
  neuron_secr  = secretome_union, in NEURON cells
Spearman across donors: niche_SASP <-> each neuron score. aging primary (YA+HA); AD secondary (AD+HA).
Caveats reported: small n; both scores co-vary with age -> raw + within-HA-only (n=9) + age-trend note.

Output: results/validation/neuron_L2_donor_coupling.csv
"""
import os
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache"); os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_cache")
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from scipy import stats
import warnings; warnings.filterwarnings("ignore")

np.random.seed(42)
P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
NICHE  = ["Astro", "Micro", "Oligo", "OPC", "Endo"]
NEURON = ["DG_GC", "CA_ExN", "InN"]

def read_gmt(path):
    s = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3: s[p[0]] = [g for g in p[2:] if g]
    return s
targeted = read_gmt(f"{P}/refs/neuron_response_sets.gmt")
sec = list(pd.read_csv(f"{P}/refs/secretome_union.csv").gene)
frozen = pd.read_csv(f"{P}/results/validation/frozen_primary_signature.csv")
niche_up = sorted(set(frozen[frozen.direction == "UP"].gene))          # niche SASP-output markers
resp_composite = sorted(set(targeted["NFkB_TNFA"]) | set(targeted["IFN_gamma"]) | set(targeted["IFN_alpha"])
                        | set(targeted["Inflammatory"]) | set(targeted["Complement"]))
SETS_NEURON = {"neuron_resp": resp_composite, "neuron_ISR": targeted["ISR_ATF4"],
               "neuron_senes": targeted["Senescence_core"], "neuron_recept": targeted["Cyt_Compl_Receptors"],
               "neuron_secr": sec}
print(f"niche_SASP genes (frozen UP): {len(niche_up)} | neuron_resp composite: {len(resp_composite)}")

print("loading anchor ...")
A = sc.read_h5ad(f"{P}/processed/per_dataset/GSE268609_anchor.h5ad")
A.obs["donor_id"] = A.obs["donor_id"].astype(str); A.obs["ct"] = A.obs["celltype_l1"].astype(str)
genes = A.var_names.astype(str)
gidx = {g: i for i, g in enumerate(genes)}

def donor_pool_logcpm(adata, donors, celltypes):
    """donor x gene log1p-CPM for the pooled cells of `celltypes`, one row per donor (>=10 pooled cells)."""
    sub = adata[(adata.obs.ct.isin(celltypes)) & (adata.obs.donor_id.isin(donors))]
    X = sub.layers["counts"]; X = X.tocsr() if sp.issparse(X) else X
    don = sub.obs.donor_id.values
    rows = {}
    for d in donors:
        m = don == d
        if m.sum() < 10:  # need a minimum pooled-cell count per donor
            continue
        v = np.asarray(X[m].sum(0)).ravel().astype(float)
        tot = v.sum()
        if tot <= 0: continue
        rows[d] = np.log1p(v / tot * 1e6)
    M = pd.DataFrame(rows, index=genes).T  # donors x genes
    return M

def sig_score(M, setgenes):
    """mean z (across donors) over set genes present & variable -> one score per donor."""
    g = [x for x in setgenes if x in M.columns]
    sub = M[g]
    sd = sub.std(0); sub = sub.loc[:, sd > 0]
    if sub.shape[1] < 3:
        return pd.Series(np.nan, index=M.index), 0
    z = (sub - sub.mean(0)) / sub.std(0)
    return z.mean(1), sub.shape[1]

def couple(adata, donors_meta, axis, ha_label="HA"):
    """donors_meta: index=donor, col 'grp'. Compute scores + Spearman niche_SASP <-> neuron sets."""
    donors = list(donors_meta.index)
    Mn = donor_pool_logcpm(adata, donors, NICHE)
    Me = donor_pool_logcpm(adata, donors, NEURON)
    common = Mn.index.intersection(Me.index)
    Mn, Me = Mn.loc[common], Me.loc[common]
    grp = donors_meta.grp.reindex(common)
    niche_s, n_ng = sig_score(Mn, niche_up)
    out = {}
    rows = []
    for sname, sg in SETS_NEURON.items():
        neur_s, n_sg = sig_score(Me, sg)
        df = pd.DataFrame({"niche": niche_s, "neur": neur_s, "grp": grp}).dropna()
        if len(df) < 5:
            continue
        rho, p = stats.spearmanr(df.niche, df.neur)
        # within-HA-only (residual coupling beyond the YA/HA age jump)
        h = df[df.grp == ha_label]
        rho_h, p_h = (stats.spearmanr(h.niche, h.neur) if len(h) >= 5 else (np.nan, np.nan))
        rows.append(dict(axis=axis, neuron_set=sname, n_donor=len(df), n_niche_genes=n_ng, n_set_genes=n_sg,
                         rho_all=round(float(rho), 3), p_all=float(p),
                         n_HA=int((df.grp == ha_label).sum()), rho_HAonly=round(float(rho_h), 3) if pd.notna(rho_h) else np.nan,
                         p_HAonly=float(p_h) if pd.notna(p_h) else np.nan))
    return rows

# aging (YA+HA) primary; AD (AD+HA) secondary
meta = (A.obs[["donor_id", "Group"]].drop_duplicates().set_index("donor_id").rename(columns={"Group": "grp"}))
meta["grp"] = meta.grp.astype(str)
all_rows = []
for axis, grps in [("aging", ["YA", "HA"]), ("ADvHA", ["AD", "HA"])]:
    dm = meta[meta.grp.isin(grps)]
    print(f"\n=== coupling {axis}: donors {dict(dm.grp.value_counts())} ===")
    all_rows += couple(A, dm, axis)
cp = pd.DataFrame(all_rows)
cp.to_csv(f"{P}/results/validation/neuron_L2_donor_coupling.csv", index=False)
print("\n=== N2.4 donor-level niche_SASP <-> neuron coupling (Spearman across donors) ===")
print(cp.to_string(index=False))
print("\nNOTE: rho_all conflates with the YA/HA age exposure (both scores track age); rho_HAonly tests")
print("residual within-aged coupling. Coordination, not causation. n is small -> suggestive at best.")
print("\nDONE (N2.4).")
