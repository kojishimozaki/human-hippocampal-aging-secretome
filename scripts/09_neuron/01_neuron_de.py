#!/usr/bin/env python
"""N0 neuron-niche differential expression (GSE268609) — donor pseudobulk, ~arm+grp.

Foundation for the neuron-niche arm (docs/NEURON_NICHE_PLAN.md, N0). Donor-level pyDESeq2 on
the matched niche + neuron (+ immature) cell types in the SAME donors, three contrasts, with the
N0-gate-validated preparation-arm batch covariate.

Design / contrasts (reference pivot = HA; see RESOLUTION 2026-06-18 sign-convention decision):
  aging    YA-vs-HA : contrast ["grp","HA","YA"]  (+log2FC = aging-UP; matches the FROZEN niche)
  ADvHA    AD-vs-HA : contrast ["grp","AD","HA"]  (+log2FC = AD-UP vs aged-healthy)
  MCIvHA   MCI-vs-HA: contrast ["grp","MCI","HA"] (+log2FC = MCI-UP vs aged-healthy, exploratory)
Design = ~arm+grp where estimable (arm 2 levels AND every arm x group cell >=1, the validated niche
rule 06_arm_robustness.py:72); else fall back to ~grp for that celltype x contrast (flagged).

Cell types: neurons DG_GC/CA_ExN/InN (even-handed, NO DG privilege) + niche Astro/Micro/Oligo/OPC/Endo
(symmetry) on all 3 contrasts; immature NSC/Neuroblast/Immature on YA-vs-HA ONLY (exploratory; AD
infeasible n=1/0). SA excluded throughout.

Engine reuse (no modification of the shared scripts/03_de/01_pseudobulk_de.py engine):
  - pseudobulk = sum raw counts (layers["counts"]) per donor   [06_arm_robustness.py]
  - gene filter sum>=10 & nonzero>=max(3,0.5n); >=10 cells/donor; >=4 donors/group; >=50 genes
  - one row per (gene,celltype): keep higher-baseMean duplicate (read-through symbols)

Sensitivities — YA-vs-HA ONLY (04_aging_robustness.py pattern):
  ~pmi+arm+grp (PMI donor-mapped via GEO SOFT/orig.ident, clean only for YA/HA) + leave-one-donor-out.
  Disease axes are NOT PMI-adjusted (04 caveat: SOFT pmi/diagnosis inconsistent for some AD/MCI GSMs).

Outputs (results/de/):
  de_GSE268609_neuron_aging_per_celltype.csv      (YA-vs-HA; niche+neuron+immature)
  de_GSE268609_neuron_ADvHA_per_celltype.csv      (AD-vs-HA; niche+neuron)
  de_GSE268609_neuron_MCIvHA_per_celltype.csv     (MCI-vs-HA; niche+neuron)
    cols: gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,celltype,contrast,n_donor,n_wh,n_dg,
          design,arm_estimable,exploratory
  de_GSE268609_neuron_aging_pmi_sensitivity.csv   (~arm+grp vs ~pmi+arm+grp, per celltype)
  de_GSE268609_neuron_aging_loo.csv               (LOO sign/sig retention, per sig secretome hit)
seed=42; donor is the unit; exact p + effect + n.
"""
import os, re, gzip
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_cache")
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import warnings; warnings.filterwarnings("ignore")

np.random.seed(42)
P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")

NEURON   = ["DG_GC", "CA_ExN", "InN"]
NICHE    = ["Astro", "Micro", "Oligo", "OPC", "Endo"]
IMMATURE = ["NSC", "Neuroblast", "Immature"]
MIN_CELLS, MIN_DONORS, MIN_GENES = 10, 4, 50

# contrast = ["grp", test, ref]; immature only on the aging axis (exploratory)
CONTRASTS = [
    dict(name="aging",  test="HA",  ref="YA", celltypes=NEURON + NICHE + IMMATURE, out="aging"),
    dict(name="ADvHA",  test="AD",  ref="HA", celltypes=NEURON + NICHE,            out="ADvHA"),
    dict(name="MCIvHA", test="MCI", ref="HA", celltypes=NEURON + NICHE,            out="MCIvHA"),
]

# ---- (A) donor -> arm (orig.ident<=14 WH / >=15 DG), keyed donor_id=SampleNumber --------------
meta = pd.read_csv(f"{P}/processed/per_dataset/GSE268609_metadata.csv", index_col=0)
donors = (meta[["orig.ident", "SampleNumber", "Group"]].drop_duplicates()
          .rename(columns={"orig.ident": "oid"}))
donors["SampleNumber"] = donors.SampleNumber.astype(str)
donors["arm"] = np.where(donors.oid.astype(int) <= 14, "WH", "DG")
assert donors.groupby("SampleNumber").arm.nunique().max() == 1, "a donor straddles two prep arms"
arm_map = dict(zip(donors.SampleNumber, donors.arm))

# ---- (B) donor -> PMI via GEO SOFT (orig.ident), clean ONLY for YA/HA (04 caveat) -------------
recs = []; cur = None
with gzip.open(f"{P}/raw/refs_soft/GSE268609_family.soft.gz", "rt", errors="ignore") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("^SAMPLE"):
            if cur: recs.append(cur)
            cur = {"title": None, "chars": []}
        elif cur is not None:
            if line.startswith("!Sample_title"):
                cur["title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1"):
                cur["chars"].append(line.split("=", 1)[1].strip())
    if cur: recs.append(cur)

def fld(ch, k):
    for c in ch:
        m = re.match(rf"\s*{k}\s*:\s*(.+)", c, re.I)
        if m: return m.group(1).strip()
    return None

soft = pd.DataFrame([dict(
    oid=int(re.match(r"\s*(\d+)", r["title"]).group(1)) if re.match(r"\s*(\d+)", r["title"] or "") else None,
    pmi=pd.to_numeric(fld(r["chars"], "pmi"), errors="coerce")) for r in recs])
_inc = soft.groupby("oid").pmi.nunique()
_bad = _inc[_inc > 1].index.tolist()
if _bad:
    print(f"[caveat] {len(_bad)} orig.ident with inconsistent SOFT pmi across RNA+ATAC GSM: {_bad} "
          f"(disease axes are NOT PMI-adjusted; YA/HA verified clean)")
soft_d = soft.groupby("oid").pmi.first().reset_index()
donors = donors.merge(soft_d, on="oid", how="left")
pmi_map = dict(zip(donors.SampleNumber, donors.pmi))   # SampleNumber(str) -> PMI

# ---- (C) load anchor once; restrict to analysed groups (drop SA); attach arm -------------------
print("loading anchor ...")
A = sc.read_h5ad(f"{P}/processed/per_dataset/GSE268609_anchor.h5ad")
A = A[A.obs.Group.isin(["YA", "HA", "MCI", "AD"])].copy()   # SA excluded from primary
A.obs["donor_id"] = A.obs["donor_id"].astype(str)
A.obs["arm"] = A.obs["donor_id"].map(arm_map)
A.obs["ct"] = A.obs["celltype_l1"].astype(str)
assert A.obs.arm.notna().all(), "unmapped arm after SA removal"
print(f"cells {A.n_obs}; donors {A.obs.donor_id.nunique()}; groups {sorted(A.obs.Group.unique())}")


def pseudobulk(sub):
    X = sub.layers["counts"]; X = X.tocsr() if sp.issparse(X) else X
    don = sub.obs["donor_id"].astype(str).values; uniq = pd.unique(don)
    mat = np.vstack([np.asarray(X[don == d].sum(0)).ravel() for d in uniq])
    return pd.DataFrame(mat, index=uniq, columns=sub.var_names.astype(str))


def gene_keep(pb):
    return (pb.sum(0) >= 10) & (pb.astype(bool).sum(0) >= max(3, int(0.5 * len(pb))))


def estimable_design(md):
    """~arm+grp iff arm has 2 levels AND every arm x grp cell is populated (06_arm_robustness rule)."""
    bal = pd.crosstab(md.grp, md.arm)
    est = (md.arm.nunique() >= 2) and (bal.size == 4) and (int(bal.values.min()) >= 1)
    return ("~arm + grp" if est else "~grp"), bool(est)


def fit(pb, md, design, contrast):
    sub = pb.loc[md.index, gene_keep(pb.loc[md.index])]
    if sub.shape[1] < MIN_GENES:
        return None, sub.shape[1]
    dds = DeseqDataSet(counts=sub.astype(int), metadata=md, design=design, quiet=True)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=contrast, quiet=True); ds.summary()
    return ds.results_df, sub.shape[1]


def build_md(sub, extra_pmi=False):
    pb = pseudobulk(sub)
    md = pd.DataFrame(index=pb.index)
    md["grp"] = sub.obs.groupby("donor_id", observed=True).Group.first().astype(str).reindex(pb.index).values
    md["arm"] = sub.obs.groupby("donor_id", observed=True).arm.first().astype(str).reindex(pb.index).values
    if extra_pmi:
        md["pmi"] = [pmi_map.get(d, np.nan) for d in pb.index]
    return pb, md


# ---- (D) main DE: per contrast x celltype -----------------------------------------------------
def run_contrast(cfg):
    test, ref = cfg["test"], cfg["ref"]
    contrast = ["grp", test, ref]
    rows = []; summary = []
    for ct in cfg["celltypes"]:
        sub = A[(A.obs.ct == ct) & (A.obs.Group.isin([test, ref]))].copy()
        nc = sub.obs.groupby("donor_id", observed=True).size()
        sub = sub[sub.obs.donor_id.isin(nc[nc >= MIN_CELLS].index)].copy()
        pb, md = build_md(sub)
        ntest = int((md.grp == test).sum()); nref = int((md.grp == ref).sum())
        if min(ntest, nref) < MIN_DONORS:
            print(f"  skip {cfg['name']}/{ct}: per-group {test}={ntest} {ref}={nref} (<{MIN_DONORS}/grp)")
            continue
        design, est = estimable_design(md)
        try:
            res, ngenes = fit(pb, md, design, contrast)
        except Exception as e:   # honest loud failure; NO method substitution
            print(f"  FAIL {cfg['name']}/{ct}: {type(e).__name__} {str(e)[:160]}")
            continue
        if res is None:
            print(f"  skip {cfg['name']}/{ct}: too few genes ({ngenes})")
            continue
        nwh = int((md.arm == "WH").sum()); ndg = int((md.arm == "DG").sum())
        exploratory = ct in IMMATURE
        r = res.copy()
        r["gene"] = r.index; r["celltype"] = ct; r["contrast"] = cfg["name"]
        r["n_donor"] = len(md); r["n_wh"] = nwh; r["n_dg"] = ndg
        r["design"] = design; r["arm_estimable"] = est; r["exploratory"] = exploratory
        rows.append(r.reset_index(drop=True))
        s10 = int((r.padj < 0.1).sum()); s20 = int((r.padj < 0.2).sum())
        summary.append(dict(contrast=cfg["name"], celltype=ct, n_donor=len(md),
                            test=f"{test}:{ntest}", ref=f"{ref}:{nref}", wh=nwh, dg=ndg,
                            design=design, n_genes=ngenes, padj_lt0_1=s10, padj_lt0_2=s20,
                            exploratory=exploratory))
        print(f"  {cfg['name']}/{ct}: n={len(md)} ({test} {ntest}/{ref} {nref}; WH {nwh}/DG {ndg}), "
              f"genes={ngenes}, design='{design}', padj<0.1={s10}, padj<0.2={s20}"
              + ("  [exploratory]" if exploratory else ""))
    de = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if len(de):   # one row per (gene,celltype): keep higher-baseMean copy (read-through symbols)
        before = len(de)
        de = (de.sort_values("baseMean", ascending=False)
                .drop_duplicates(subset=["gene", "celltype"], keep="first").reset_index(drop=True))
        if len(de) < before:
            print(f"    collapsed {before - len(de)} duplicate (gene,celltype) rows")
    cols = ["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "celltype",
            "contrast", "n_donor", "n_wh", "n_dg", "design", "arm_estimable", "exploratory"]
    if len(de): de = de[cols]
    out = f"{P}/results/de/de_GSE268609_neuron_{cfg['out']}_per_celltype.csv"
    de.to_csv(out, index=False)
    print(f"  wrote {out}: {len(de)} rows, {de.celltype.nunique() if len(de) else 0} celltypes")
    return de, pd.DataFrame(summary)


print("\n=== N0 differential expression (donor pseudobulk, ~arm+grp) ===")
all_summ = []
de_by = {}
for cfg in CONTRASTS:
    print(f"\n## contrast {cfg['name']} (test={cfg['test']} ref={cfg['ref']}) ##")
    de, summ = run_contrast(cfg)
    de_by[cfg["name"]] = de
    all_summ.append(summ)

# ---- (E) YA-vs-HA sensitivities: ~pmi+arm+grp + LOO -------------------------------------------
print("\n=== YA-vs-HA sensitivities (~pmi+arm+grp + LOO) ===")
sec = set(pd.read_csv(f"{P}/refs/secretome_union.csv").gene)
contrast_aging = ["grp", "HA", "YA"]
sens_rows = []; loo_rows = []
aging_cts = [c for c in NEURON + NICHE + IMMATURE]
for ct in aging_cts:
    sub = A[(A.obs.ct == ct) & (A.obs.Group.isin(["YA", "HA"]))].copy()
    nc = sub.obs.groupby("donor_id", observed=True).size()
    sub = sub[sub.obs.donor_id.isin(nc[nc >= MIN_CELLS].index)].copy()
    pb, md = build_md(sub, extra_pmi=True)
    if min((md.grp == "YA").sum(), (md.grp == "HA").sum()) < MIN_DONORS:
        continue
    if md.pmi.isna().any():           # none expected for YA/HA; drop to be safe + report
        miss = md.index[md.pmi.isna()].tolist()
        print(f"  [pmi] {ct}: dropping {len(miss)} donor(s) lacking SOFT PMI: {miss}")
        keep = md.pmi.notna(); md = md[keep]; pb = pb.loc[md.index]
    design, est = estimable_design(md)
    pmi_design = ("~pmi + arm + grp" if est else "~pmi + grp")
    try:
        r_pri, _ = fit(pb, md, design, contrast_aging)                   # primary ~arm+grp (or ~grp)
        r_pmi, _ = fit(pb, md, pmi_design, contrast_aging)               # + PMI
    except Exception as e:   # honest loud failure; NO method substitution
        print(f"  FAIL pmi-sens/{ct}: {type(e).__name__} {str(e)[:160]}")
        continue
    if r_pri is None or r_pmi is None:
        continue
    common = r_pri.index.intersection(r_pmi.index)
    cs = pd.DataFrame({"lfc": r_pri.loc[common, "log2FoldChange"], "padj": r_pri.loc[common, "padj"],
                       "lfc_pmi": r_pmi.loc[common, "log2FoldChange"], "padj_pmi": r_pmi.loc[common, "padj"]})
    cs_sec = cs[cs.index.isin(sec)]
    sig = cs_sec[cs_sec.padj < 0.1]
    sens_rows.append(dict(celltype=ct, design=design, pmi_design=pmi_design, n_donor=len(md),
        exploratory=(ct in IMMATURE), n_sig_secr=len(sig),
        lfc_corr=round(float(cs_sec.lfc.corr(cs_sec.lfc_pmi)), 3) if len(cs_sec) > 1 else np.nan,
        sign_retained=int((np.sign(sig.lfc) == np.sign(sig.lfc_pmi)).sum()),
        padj01_retained=int((sig.padj_pmi < 0.1).sum())))
    for g in sig.index:
        loo_rows.append(dict(celltype=ct, gene=g, in_secretome=True,
            full_lfc=round(float(sig.loc[g, "lfc"]), 3), full_padj=float(sig.loc[g, "padj"]),
            pmiadj_lfc=round(float(sig.loc[g, "lfc_pmi"]), 3), pmiadj_padj=float(sig.loc[g, "padj_pmi"]),
            exploratory=(ct in IMMATURE)))

# leave-one-donor-out on the primary aging design (per sig secretome hit)
loo = pd.DataFrame(loo_rows)
if len(loo):
    loo = loo.set_index(["celltype", "gene"])
    loo["loo_n"] = 0; loo["loo_sign_stable_frac"] = np.nan; loo["loo_max_padj"] = np.nan; loo["loo_grp_fallback"] = 0
    for ct in loo.index.get_level_values(0).unique():
        sub = A[(A.obs.ct == ct) & (A.obs.Group.isin(["YA", "HA"]))].copy()
        nc = sub.obs.groupby("donor_id", observed=True).size()
        sub = sub[sub.obs.donor_id.isin(nc[nc >= MIN_CELLS].index)].copy()
        pb, md = build_md(sub)
        genes = [g for (c, g) in loo.index if c == ct]
        full_sign = {g: np.sign(loo.loc[(ct, g), "full_lfc"]) for g in genes}
        signs = {g: [] for g in genes}; padjs = {g: [] for g in genes}; n_loo = 0; nfb = 0
        for d in md.index:
            md_d = md.drop(index=d)
            if min((md_d.grp == "YA").sum(), (md_d.grp == "HA").sum()) < 2:
                continue
            design_d, est_d = estimable_design(md_d)
            if not est_d: nfb += 1
            try:
                r, _ = fit(pb.drop(index=d), md_d, design_d, contrast_aging)
            except Exception:
                continue
            if r is None: continue
            n_loo += 1
            for g in genes:
                if g in r.index:
                    signs[g].append(np.sign(r.loc[g, "log2FoldChange"])); padjs[g].append(r.loc[g, "padj"])
        for g in genes:
            loo.loc[(ct, g), "loo_n"] = n_loo
            loo.loc[(ct, g), "loo_grp_fallback"] = nfb
            loo.loc[(ct, g), "loo_sign_stable_frac"] = (
                round(np.mean([s == full_sign[g] for s in signs[g]]), 3) if signs[g] else np.nan)
            loo.loc[(ct, g), "loo_max_padj"] = (
                round(float(np.nanmax(padjs[g])), 4) if padjs[g] and not all(pd.isna(padjs[g])) else np.nan)
    loo = loo.reset_index()
loo.to_csv(f"{P}/results/de/de_GSE268609_neuron_aging_loo.csv", index=False)
sens = pd.DataFrame(sens_rows)
sens.to_csv(f"{P}/results/de/de_GSE268609_neuron_aging_pmi_sensitivity.csv", index=False)

print("\n--- PMI-adjusted sensitivity (~arm+grp vs ~pmi+arm+grp), aging axis ---")
if len(sens): print(sens.to_string(index=False))
if len(loo):
    ndir = int((loo.loo_sign_stable_frac == 1).sum()); nsig = int((loo.loo_max_padj <= 0.1).sum())
    print(f"\n--- LOO (honest) --- {len(loo)} sig secretome hits | dir stable ALL refits: {ndir}/{len(loo)} | "
          f"also padj<=0.1 ALL refits: {nsig}/{len(loo)}")

# ---- (F) consolidated summary + secretome intersect (orientation only) ------------------------
print("\n=== N0 SUMMARY: per (celltype x contrast) hits + secretome intersect (orientation only) ===")
summ = pd.concat(all_summ, ignore_index=True)
def sec_hits(name, ct, thr):
    de = de_by[name]
    if not len(de): return 0
    d = de[(de.celltype == ct)]
    return int(((d.padj < thr) & (d.gene.isin(sec))).sum())
summ["secr_padj0_1"] = [sec_hits(r.contrast, r.celltype, 0.1) for r in summ.itertuples()]
summ["secr_padj0_2"] = [sec_hits(r.contrast, r.celltype, 0.2) for r in summ.itertuples()]
pd.set_option("display.width", 240, "display.max_rows", 80)
print(summ.to_string(index=False))
summ.to_csv(f"{P}/results/de/de_GSE268609_neuron_N0_summary.csv", index=False)
print("\nwrote results/de/de_GSE268609_neuron_N0_summary.csv")
print("\nDONE (N0 DE). STOP per protocol — await PI go for N1.")
