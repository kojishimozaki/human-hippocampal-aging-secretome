#!/usr/bin/env python
"""Aging-secretome PMI handling + robustness for GSE268609 YA-vs-HA.

Audit v2.1 correction (2026-06-04). The previous version wrongly claimed per-donor PMI was
unmappable and reported a PSEUDOREPLICATED PMI p-value. Both are fixed here:

 * PMI/age ARE donor-mappable for the YA/HA primary contrast. The GEO SOFT `!Sample_title`
   integer equals the Seurat `orig.ident`; each orig.ident has exactly 2 GSM (RNA+ATAC) with
   identical pmi/age, and SOFT diagnosis matches the Seurat Group for **YA 8/8 and HA 9/9**.
   (Only AD has 2 orig.ident with a SOFT-MCI/Seurat-AD label discrepancy; irrelevant to YA/HA.)
   So we map orig.ident -> SOFT pmi/age, key to the donor (Seurat SampleNumber), and:

 * group PMI/age is summarised at the DONOR level (one row per donor), NOT per GSM. The earlier
   GSM-level summary (YA n=16 / HA n=18) double-counted the RNA+ATAC libraries -> pseudoreplication,
   which is what produced the bogus p=0.011. Donor-level (YA n=8 / HA n=9): same means
   (5.84 / 7.67 h) but two-sided Mann-Whitney p = 0.081.

 * Because PMI is mappable, we now fit the PMI-ADJUSTED sensitivity model `~pmi + grp` and compare
   it to the primary `~grp` (log2FC correlation + significant-hit sign/significance retention).

 * Leave-one-donor-out is reported honestly: direction is stable in all refits, but only a subset
   keep padj<=0.1 (per-gene significance is donor-sensitive at n=8-9).

YA (age 21-38) vs HA (age 60-93) remains a young-vs-old contrast: age is the exposure, not an
adjustable covariate. Sex is unavailable in the deposit (sex='U').

Outputs:
  refs/gse268609_group_pmi_age.csv      DONOR-level group PMI/age summary (GEO SOFT via orig.ident)
  results/de/aging_pmi_sensitivity.csv  per niche celltype: ~grp vs ~pmi+grp, secretome sig hits
  results/de/aging_loo_robustness.csv   per significant secretome gene: LOO sign + significance
"""
import os, re, gzip
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")   # scanpy/numba + matplotlib reproducibility
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_cache")        # on hosts with a read-only/cacheless HOME
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from scipy.stats import mannwhitneyu
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import warnings; warnings.filterwarnings("ignore")

P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
NICHE = ["Astro", "Micro", "Oligo", "OPC", "Endo"]

# ---- (A) donor-level PMI/age via orig.ident -> GEO SOFT ----------------------------------
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
    diagnosis=fld(r["chars"], "diagnosis"),
    age=pd.to_numeric(fld(r["chars"], "age[ _]at_death"), errors="coerce"),
    pmi=pd.to_numeric(fld(r["chars"], "pmi"), errors="coerce")) for r in recs])
# collapse the 2 GSM (RNA+ATAC) per orig.ident to one donor row
soft_d = soft.groupby("oid").agg(diagnosis=("diagnosis", "first"), age=("age", "first"),
                                 pmi=("pmi", "first")).reset_index()
# Data-quality caveat: the group PMI/age summary below is verified clean ONLY for YA/HA
# (orig.ident -> SOFT is 8/8 + 9/9 group-concordant). AD has 2 orig.ident SOFT-labelled MCI,
# and some disease-group orig.ident disagree on pmi/age across their RNA+ATAC GSM. Do NOT reuse
# this CSV for a disease-axis PMI adjustment without resolving those per-donor.
_inc = soft.groupby("oid").agg(npmi=("pmi", "nunique"), nage=("age", "nunique"), ndiag=("diagnosis", "nunique"))
_bad = _inc[(_inc.npmi > 1) | (_inc.nage > 1) | (_inc.ndiag > 1)]
if len(_bad):
    print(f"[caveat] {len(_bad)} orig.ident with inconsistent SOFT pmi/age/diagnosis across RNA+ATAC: {_bad.index.tolist()}")

meta = pd.read_csv(f"{P}/processed/per_dataset/GSE268609_metadata.csv", index_col=0)
donors = (meta[["orig.ident", "SampleNumber", "Group"]].drop_duplicates()
          .rename(columns={"orig.ident": "oid"}))
donors = donors.merge(soft_d, on="oid", how="left")
donors["SampleNumber"] = donors.SampleNumber.astype(str)
# transparency: SOFT-diagnosis vs Seurat-Group concordance
print("=== SOFT diagnosis (by orig.ident) vs Seurat Group ===")
for grp in ["YA", "HA", "AD", "MCI", "SA"]:
    s = donors[donors.Group == grp]
    print(f"  {grp}: n={len(s)} match={int((s.Group == s.diagnosis).sum())}/{len(s)} "
          f"SOFT={s.diagnosis.value_counts(dropna=False).to_dict()}")

# DONOR-level group PMI/age summary (overwrite the old GSM-level file)
grp = donors.groupby("Group").agg(
    n_donor=("SampleNumber", "nunique"), pmi_mean=("pmi", "mean"), pmi_sd=("pmi", "std"),
    pmi_min=("pmi", "min"), pmi_max=("pmi", "max"),
    age_mean=("age", "mean"), age_min=("age", "min"), age_max=("age", "max")).round(2)
grp.to_csv(f"{P}/refs/gse268609_group_pmi_age.csv")
print("\n=== DONOR-level group PMI / age (GEO SOFT via orig.ident) ===")
print(grp.to_string())
ya = donors[donors.Group == "YA"].pmi.dropna(); ha = donors[donors.Group == "HA"].pmi.dropna()
pmi_p = mannwhitneyu(ya, ha, alternative="two-sided").pvalue
print(f"\nYA vs HA donor-level PMI: YA n={len(ya)} mean={ya.mean():.2f}; HA n={len(ha)} mean={ha.mean():.2f}; "
      f"two-sided Mann-Whitney p={pmi_p:.4f}  (prior GSM-level n=16/18 gave a pseudoreplicated p=0.011)")
pmi_map = dict(zip(donors.SampleNumber, donors.pmi))   # SampleNumber(str) -> PMI

# ---- (B) load anchor, YA/HA ----------------------------------------------------------------
sec = set(pd.read_csv(f"{P}/refs/secretome_union.csv").gene)
A = sc.read_h5ad(f"{P}/processed/per_dataset/GSE268609_anchor.h5ad")
A.X = A.layers["counts"]
A = A[A.obs.Group.isin(["YA", "HA"])].copy()
print(f"\nYA+HA cells: {A.n_obs}; donors: {A.obs.donor_id.nunique()}")

def pseudobulk(sub):
    X = sub.layers["counts"]; X = X.tocsr() if sp.issparse(X) else X
    don = sub.obs["donor_id"].astype(str).values; uniq = pd.unique(don)
    mat = np.vstack([np.asarray(X[don == d].sum(0)).ravel() for d in uniq])
    return pd.DataFrame(mat, index=uniq, columns=sub.var_names.astype(str))

def fit(pb, meta, design, drop=None):
    donors_ = [d for d in pb.index if drop is None or d != drop]
    md = meta.loc[donors_]
    if (md.grp == "YA").sum() < 2 or (md.grp == "HA").sum() < 2:
        return None
    sub = pb.loc[donors_]
    keep = (sub.sum(0) >= 10) & (sub.astype(bool).sum(0) >= max(3, int(0.5 * len(sub))))
    sub = sub.loc[:, keep]
    dds = DeseqDataSet(counts=sub.astype(int), metadata=md, design=design, quiet=True)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=["grp", "HA", "YA"], quiet=True); ds.summary()
    return ds.results_df

# ---- (C) PMI-adjusted sensitivity (~grp vs ~pmi+grp) + LOO --------------------------------
sens_rows = []; loo_rows = []
for ct in NICHE:
    sub = A[A.obs.celltype_l1.astype(str) == ct].copy()
    nc = sub.obs.groupby("donor_id", observed=True).size()
    sub = sub[sub.obs.donor_id.isin(nc[nc >= 10].index)].copy()
    pb = pseudobulk(sub)
    md = pd.DataFrame(index=pb.index)
    md["grp"] = sub.obs.groupby("donor_id", observed=True).Group.first().astype(str).reindex(pb.index).values
    md["pmi"] = [pmi_map.get(d, np.nan) for d in pb.index]
    if md.pmi.isna().any():                       # drop donors lacking PMI (none expected for YA/HA)
        ok = md.pmi.notna(); pb = pb.loc[md.index[ok]]; md = md[ok]

    r_grp = fit(pb, md, "~grp")
    r_adj = fit(pb, md, "~pmi + grp")
    if r_grp is None or r_adj is None:
        print(f"skip {ct}"); continue
    common = r_grp.index.intersection(r_adj.index)
    cmp = pd.DataFrame({"lfc_grp": r_grp.loc[common, "log2FoldChange"], "padj_grp": r_grp.loc[common, "padj"],
                        "lfc_pmiadj": r_adj.loc[common, "log2FoldChange"], "padj_pmiadj": r_adj.loc[common, "padj"]})
    cs = cmp[cmp.index.isin(sec)]
    sig = cs[cs.padj_grp < 0.1]
    sens_rows.append(dict(celltype=ct, n_donor=len(md), n_sig_secr=len(sig),
        lfc_corr=round(float(cs.lfc_grp.corr(cs.lfc_pmiadj)), 3),
        sig_sign_retained=int((np.sign(sig.lfc_grp) == np.sign(sig.lfc_pmiadj)).sum()),
        sig_still_padj01=int((sig.padj_pmiadj < 0.1).sum())))
    for g in sig.index:
        loo_rows.append(dict(celltype=ct, gene=g, full_lfc=round(float(sig.loc[g, "lfc_grp"]), 3),
                             full_padj=float(sig.loc[g, "padj_grp"]),
                             pmiadj_lfc=round(float(sig.loc[g, "lfc_pmiadj"]), 3),
                             pmiadj_padj=float(sig.loc[g, "padj_pmiadj"])))

# leave-one-donor-out on ~grp (direction + significance retention)
loo = pd.DataFrame(loo_rows).set_index(["celltype", "gene"])
loo["loo_n"] = 0; loo["loo_sign_stable_frac"] = np.nan; loo["loo_max_padj"] = np.nan
for ct in NICHE:
    sub = A[A.obs.celltype_l1.astype(str) == ct].copy()
    nc = sub.obs.groupby("donor_id", observed=True).size()
    sub = sub[sub.obs.donor_id.isin(nc[nc >= 10].index)].copy()
    pb = pseudobulk(sub)
    md = pd.DataFrame(index=pb.index)
    md["grp"] = sub.obs.groupby("donor_id", observed=True).Group.first().astype(str).reindex(pb.index).values
    genes = [g for (c, g) in loo.index if c == ct]
    if not genes: continue
    signs = {g: [] for g in genes}; padjs = {g: [] for g in genes}; n_loo = 0
    full_sign = {g: np.sign(loo.loc[(ct, g), "full_lfc"]) for g in genes}
    for d in pb.index:
        r = fit(pb, md, "~grp", drop=d)
        if r is None: continue
        n_loo += 1
        for g in genes:
            if g in r.index:
                signs[g].append(np.sign(r.loc[g, "log2FoldChange"])); padjs[g].append(r.loc[g, "padj"])
    for g in genes:
        loo.loc[(ct, g), "loo_n"] = n_loo
        loo.loc[(ct, g), "loo_sign_stable_frac"] = round(np.mean([s == full_sign[g] for s in signs[g]]), 3) if signs[g] else np.nan
        loo.loc[(ct, g), "loo_max_padj"] = round(float(np.nanmax(padjs[g])), 4) if padjs[g] else np.nan
loo = loo.reset_index()
loo.to_csv(f"{P}/results/de/aging_loo_robustness.csv", index=False)

sens = pd.DataFrame(sens_rows)
sens.to_csv(f"{P}/results/de/aging_pmi_sensitivity.csv", index=False)
print("\n=== PMI-adjusted sensitivity (~grp vs ~pmi+grp), niche celltypes ===")
print(sens.to_string(index=False))
ndir = int((loo.loo_sign_stable_frac == 1).sum()); nsig = int((loo.loo_max_padj <= 0.1).sum()); n = len(loo)
print(f"\n=== LOO (honest) === {n} sig secretome hits | direction stable in ALL refits: {ndir}/{n} | "
      f"also padj<=0.1 in ALL refits: {nsig}/{n}")
tot = sens.n_sig_secr.sum()
print(f"\n=== PMI-adjustment summary === across niche: log2FC(~grp vs ~pmi+grp) corr ~"
      f"{sens.lfc_corr.min():.2f}-{sens.lfc_corr.max():.2f}; "
      f"sign retained {sens.sig_sign_retained.sum()}/{tot}; padj<0.1 retained {sens.sig_still_padj01.sum()}/{tot}")
print("\nwrote refs/gse268609_group_pmi_age.csv, results/de/aging_pmi_sensitivity.csv, results/de/aging_loo_robustness.csv")
