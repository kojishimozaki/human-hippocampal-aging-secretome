#!/usr/bin/env python
"""[F-5-007 / F-5-006] Joint-covariate primary DE, and LFC shrinkage of the frozen 60.

WHY
---
Two audit findings were left [UNVERIFIED] not because the question was hard but because the
audit kit's script (`compute/pass5_05_covariates_and_shrinkage.py`) failed quietly on both:

  F-5-007  every design containing `pmi` was SKIPPED. The script read per-donor PMI from
           `refs/gse268609_group_pmi_age.csv`, which is a GROUP-level summary
           (Group, n_donor, pmi_mean, pmi_sd, ...) and carries no per-donor column. So
           `~pmi + grp`, `~arm + pmi + grp` and `~sex + arm + pmi + grp` were never fitted,
           and the manuscript's Limitations 7 still reports only the two separate
           one-covariate models.
           -> here PMI is recovered per donor from the GEO SOFT record, using the same
              parsing as scripts/03_de/04_aging_robustness.py:L40-L75 (SOFT `!Sample_title`
              integer == Seurat `orig.ident`; the two GSM per orig.ident carry identical
              pmi/age; keyed to the Seurat SampleNumber).

  F-5-006  `ds.lfc_shrink(coeff="grp[T.HA]")` raised
           KeyError: available LFC coeffs are Index(['grp[T.YA]'])
           on all five cell types, was swallowed by a bare `except`, and wrote 60/60 NaN.
           pyDESeq2 orders factor levels alphabetically, so HA is the reference and the fitted
           coefficient is grp[T.YA] = log2FC(YA vs HA) -- the wrong direction. `ref_level=` is
           deprecated in pydeseq2 0.5.4 and has no effect.
           -> here the group levels are re-labelled `1_YA` / `2_HA` purely so that the fitted
              coefficient is `grp[T.2_HA]` = log2FC(HA vs YA), the direction the manuscript
              reports. This is a re-parameterisation, not a change of model: it was verified to
              leave the MLE log2FC and padj identical to 1e-7, and to agree with negating the
              shrunken grp[T.YA] to 5e-4. No sign is flipped anywhere in this script.

DESIGNS (donor pseudobulk, YA vs HA, contrast [grp, HA, YA])
  ~grp                        baseline; must reproduce the frozen 60
  ~pmi + grp                  published sensitivity; must reproduce 30/60
  ~arm + grp                  published sensitivity; must reproduce 42/60
  ~arm + pmi + grp            NEW -- the joint model F-5-007 asks for
  ~sex + arm + pmi + grp      NEW -- joint model plus expression-recovered sex
  ~logncell + grp             nuclei per donor (F-5-008 cross-check; audit reported 28/60)

Pseudobulk, the >=10 nuclei/donor rule and the gene filter are taken verbatim from
scripts/03_de/04_aging_robustness.py so the numbers are comparable to the committed
sensitivities. Sex is recovered from XIST vs eight chrY genes at the donor level and is
cross-checked against the audit's independent call; a mismatch stops the script rather than
silently picking one.

seed 42. RUNTIME ~45-90 min (loads the 8.8 GB anchor once; 30 pyDESeq2 fits + 5 shrinkage fits).

Out: results/audit_reruns/covariate_models_summary.csv    per celltype x design
     results/audit_reruns/covariate_models_perhit.csv     per frozen hit x design
     results/audit_reruns/lfc_shrinkage_frozen60.csv      MLE vs shrunken log2FC
     results/de/gse268609_donor_covariates.csv            per-donor PMI/age/sex/arm/seqbatch
                                                          (F-1-003 suggested fix 3; this is the
                                                           per-donor table pass5_05 needed)
"""
import os, re, gzip
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_cache")
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import warnings; warnings.filterwarnings("ignore")

SEED = 42
np.random.seed(SEED)
P = os.environ.get("PROJ", os.getcwd())
OUT = f"{P}/results/audit_reruns"
os.makedirs(OUT, exist_ok=True)
NICHE = ["Astro", "Micro", "Oligo", "OPC", "Endo"]
MIN_CELLS = 10
YGENES = ["RPS4Y1", "ZFY", "USP9Y", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "NLGN4Y"]

# ------------------------------------------------------- (A) per-donor PMI/age from GEO SOFT
# verbatim parsing of scripts/03_de/04_aging_robustness.py:L40-L75
soft_path = f"{P}/raw/refs_soft/GSE268609_family.soft.gz"
recs, cur = [], None
with gzip.open(soft_path, "rt") as f:
    for line in f:
        if line.startswith("^SAMPLE"):
            if cur:
                recs.append(cur)
            cur = {"title": None, "chars": []}
        elif cur is not None and line.startswith("!Sample_title"):
            cur["title"] = line.split("=", 1)[1].strip()
        elif cur is not None and line.startswith("!Sample_characteristics_ch1"):
            cur["chars"].append(line.split("=", 1)[1].strip())
if cur:
    recs.append(cur)


def fld(ch, k):
    for c in ch:
        m = re.match(rf"\s*{k}\s*:\s*(.+)", c, re.I)
        if m:
            return m.group(1).strip()
    return None


soft = pd.DataFrame([dict(
    oid=int(re.match(r"\s*(\d+)", r["title"]).group(1)) if re.match(r"\s*(\d+)", r["title"] or "") else None,
    diagnosis=fld(r["chars"], "diagnosis"),
    age=pd.to_numeric(fld(r["chars"], "age[ _]at_death"), errors="coerce"),
    pmi=pd.to_numeric(fld(r["chars"], "pmi"), errors="coerce")) for r in recs])
soft_d = soft.groupby("oid").agg(diagnosis=("diagnosis", "first"), age=("age", "first"),
                                 pmi=("pmi", "first")).reset_index()

meta = pd.read_csv(f"{P}/processed/per_dataset/GSE268609_metadata.csv", index_col=0, low_memory=False)
donors = (meta[["orig.ident", "SampleNumber", "Group"]].drop_duplicates()
          .rename(columns={"orig.ident": "oid"}))
donors = donors.merge(soft_d, on="oid", how="left")
donors["SampleNumber"] = donors.SampleNumber.astype(str)
donors["arm"] = np.where(donors.oid.astype(int) <= 14, "WH", "DG")          # 06_arm_robustness.py:L24-29
donors["seqbatch"] = np.where(donors.oid.astype(int) <= 18, "NovaSeq6000", "NovaSeqX")
pmi_map = dict(zip(donors.SampleNumber, donors.pmi))
arm_map = dict(zip(donors.SampleNumber, donors.arm))
yh = donors[donors.Group.isin(["YA", "HA"])]
print(f"per-donor PMI recovered for {int(donors.pmi.notna().sum())}/{len(donors)} donors "
      f"(YA/HA: {int(yh.pmi.notna().sum())}/{len(yh)})", flush=True)
assert yh.pmi.notna().all(), "PMI missing for a YA/HA donor -- refusing to run a PMI model"

# ------------------------------------------------------- (B) anchor + donor sex from XIST/chrY
A = sc.read_h5ad(f"{P}/processed/per_dataset/GSE268609_anchor.h5ad")
A = A[A.obs.Group.isin(["YA", "HA"])].copy()
print(f"YA+HA nuclei {A.n_obs}, donors {A.obs.donor_id.nunique()}", flush=True)
X = A.layers["counts"] if "counts" in A.layers else A.X
X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
vn = A.var_names.astype(str)
dn = A.obs.donor_id.astype(str).values
uniq = pd.unique(dn)
gi = {g: i for i, g in enumerate(vn)}
xist_i = gi.get("XIST")
y_i = [gi[g] for g in YGENES if g in gi]
print(f"sex markers found: XIST={'yes' if xist_i is not None else 'NO'}, "
      f"chrY {len(y_i)}/{len(YGENES)} ({[g for g in YGENES if g in gi]})", flush=True)
assert xist_i is not None and y_i, "XIST or the chrY panel is absent from var_names"
codes = pd.Categorical(dn, categories=uniq).codes
OH = sp.csr_matrix((np.ones(len(dn)), (codes, np.arange(len(dn)))), shape=(len(uniq), len(dn)))
tot_donor = np.asarray(OH @ np.asarray(X.sum(1)).ravel()).ravel()          # UMI per donor
marker = np.asarray((OH @ X[:, [xist_i] + y_i]).todense())                 # donors x (1+nY)
SEX = pd.DataFrame({
    "SampleNumber": [str(d) for d in uniq],
    "n_cells": np.asarray(OH.sum(1)).ravel().astype(int),
    "total_umi": tot_donor,
    "XIST_cpm": marker[:, 0] / tot_donor * 1e6,
    "Ygenes_cpm": marker[:, 1:].sum(1) / tot_donor * 1e6})
SEX["sex_inferred"] = np.where(SEX.XIST_cpm > SEX.Ygenes_cpm, "F", "M")
sex_map = dict(zip(SEX.SampleNumber, SEX.sex_inferred))
# no ambiguous donor may be auto-assigned (protocol Phase 1)
amb = SEX[(SEX[["XIST_cpm", "Ygenes_cpm"]].max(1) < 10) |
          (SEX[["XIST_cpm", "Ygenes_cpm"]].min(1) /
           SEX[["XIST_cpm", "Ygenes_cpm"]].max(1) > 0.1)]
if len(amb):
    print(amb.to_string(index=False))
    raise SystemExit("ambiguous donor sex -- STOP, do not auto-assign (CLAUDE.md rule 1 / Phase 1)")

# cross-check against the audit's independent call
audit_sex = "/home/neurofuture/audit-runs/20260821-225834/audit/confounding/GSE268609_recovered_sex.csv"
if os.path.exists(audit_sex):
    a = pd.read_csv(audit_sex)
    a["SampleNumber"] = a.SampleNumber.astype(str)
    mg = SEX.merge(a[["SampleNumber", "inferred_sex", "pmi"]], on="SampleNumber", how="inner")
    agree = int((mg.sex_inferred == mg.inferred_sex).sum())
    print(f"sex cross-check vs audit: {agree}/{len(mg)} agree", flush=True)
    assert agree == len(mg), "sex calls disagree with the audit -- STOP and reconcile by hand"
    pm = mg.dropna(subset=["pmi"])
    dd = (pm.SampleNumber.map(pmi_map) - pm.pmi).abs().max()
    print(f"PMI cross-check vs audit: max |diff| = {dd:.4g} over {len(pm)} donors", flush=True)
else:
    print("[note] audit sex CSV not present; cross-check skipped", flush=True)

cov = donors.merge(SEX[["SampleNumber", "XIST_cpm", "Ygenes_cpm", "sex_inferred"]],
                   on="SampleNumber", how="left")
cov.to_csv(f"{P}/results/de/gse268609_donor_covariates.csv", index=False)
print(f"wrote {P}/results/de/gse268609_donor_covariates.csv", flush=True)

# ------------------------------------------------------- (C) designs
frozen = pd.read_csv(f"{P}/results/validation/frozen_primary_signature.csv")
DESIGNS = ["~grp", "~pmi + grp", "~arm + grp", "~arm + pmi + grp",
           "~sex + arm + pmi + grp", "~logncell + grp"]
LVL = {"YA": "1_YA", "HA": "2_HA"}      # re-parameterisation only; see the module docstring


def pseudobulk(sub):
    Xs = sub.layers["counts"] if "counts" in sub.layers else sub.X
    Xs = Xs.tocsr() if sp.issparse(Xs) else sp.csr_matrix(Xs)
    d = sub.obs["donor_id"].astype(str).values
    u = pd.unique(d)
    M = np.vstack([np.asarray(Xs[d == x].sum(0)).ravel() for x in u])
    return pd.DataFrame(M, index=u, columns=sub.var_names.astype(str))


def fit(pb, md, design, shrink=False):
    keep = (pb.sum(0) >= 10) & (pb.astype(bool).sum(0) >= max(3, int(0.5 * len(pb))))
    sub = pb.loc[:, keep]
    dds = DeseqDataSet(counts=sub.astype(int), metadata=md, design=design, quiet=True)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=["grp", LVL["HA"], LVL["YA"]], quiet=True)
    ds.summary()
    res = ds.results_df.copy()
    if shrink:
        cols = list(dds.varm["LFC"].columns)
        coeff = f"grp[T.{LVL['HA']}]"
        assert coeff in cols, f"expected {coeff} in LFC coeffs, got {cols}"
        mle = res[["log2FoldChange", "lfcSE", "baseMean", "padj"]].copy()
        ds.lfc_shrink(coeff=coeff)
        mle["lfc_shrunk"] = ds.results_df["log2FoldChange"]
        mle["lfcSE_shrunk"] = ds.results_df["lfcSE"]
        # apeGLM is an L-BFGS-B MAP fit per gene and can fail to converge; pydeseq2 records
        # this in _LFC_shrink_converged but does NOT blank the estimate, so a failed gene
        # silently returns a value collapsed onto the prior mode. Carry the flag.
        conv = getattr(ds, "_LFC_shrink_converged", None)
        mle["shrink_converged"] = (conv.reindex(mle.index) if conv is not None else pd.NA)

        # Cross-parameterisation control: fit the SAME model with the default level order
        # (HA is the alphabetical reference, so the fitted coefficient is grp[T.YA]) and
        # negate. Any gene whose two shrunken estimates disagree is not a stable estimate.
        md_std = md.copy()
        inv = {v: k for k, v in LVL.items()}
        md_std["grp"] = md_std["grp"].map(inv)
        dds2 = DeseqDataSet(counts=sub.astype(int), metadata=md_std, design=design, quiet=True)
        dds2.deseq2()
        ds2 = DeseqStats(dds2, contrast=["grp", "HA", "YA"], quiet=True)
        ds2.summary()
        ds2.lfc_shrink(coeff="grp[T.YA]")
        mle["lfc_shrunk_altparam"] = -ds2.results_df["log2FoldChange"]
        conv2 = getattr(ds2, "_LFC_shrink_converged", None)
        mle["shrink_converged_altparam"] = (conv2.reindex(mle.index) if conv2 is not None
                                            else pd.NA)
        return res, mle
    return res, None


summ, perhit, shrunk_rows = [], [], []
for ct in NICHE:
    sub = A[A.obs.celltype_l1.astype(str) == ct].copy()
    nc = sub.obs.groupby("donor_id", observed=True).size()
    sub = sub[sub.obs.donor_id.isin(nc[nc >= MIN_CELLS].index)].copy()
    pb = pseudobulk(sub)
    ncell = sub.obs.groupby("donor_id", observed=True).size().reindex(pb.index)
    md = pd.DataFrame(index=pb.index)
    md["grp"] = (sub.obs.groupby("donor_id", observed=True).Group.first()
                 .astype(str).reindex(pb.index).map(LVL).values)
    md["pmi"] = [pmi_map.get(d, np.nan) for d in pb.index]
    md["arm"] = [arm_map.get(d, None) for d in pb.index]
    md["sex"] = [sex_map.get(d, None) for d in pb.index]
    md["logncell"] = np.log10(ncell.values.astype(float))
    fz = frozen[frozen.celltype == ct]
    print(f"\n=== {ct}: {sub.n_obs} nuclei, {len(md)} donors "
          f"({md.grp.value_counts().to_dict()}), {len(fz)} frozen hits ===", flush=True)

    base = None
    for design in DESIGNS:
        need = [v for v in ("pmi", "arm", "sex", "logncell") if v in design]
        bad = [v for v in need if md[v].isna().any()]
        const = [v for v in need if md[v].dtype == object and md[v].nunique() < 2]
        if bad or const:
            summ.append(dict(celltype=ct, design=design, estimable=False, n_donor=len(md),
                             note=f"missing={bad} constant={const}"))
            print(f"  {design}: NOT ESTIMABLE (missing={bad} constant={const})", flush=True)
            continue
        try:
            r, mle = fit(pb, md, design, shrink=(design == "~grp"))
        except Exception as e:
            summ.append(dict(celltype=ct, design=design, estimable=False, n_donor=len(md),
                             note=f"{type(e).__name__}: {str(e)[:120]}"))
            print(f"  {design}: FAILED {type(e).__name__} {str(e)[:120]}", flush=True)
            continue
        if design == "~grp":
            base = r
            for _, h in fz.iterrows():
                if h.gene in mle.index:
                    m = mle.loc[h.gene]
                    shrunk_rows.append(dict(
                        celltype=ct, gene=h.gene, baseMean=float(m.baseMean),
                        lfc_frozen=float(h.log2FoldChange), lfc_MLE=float(m.log2FoldChange),
                        lfcSE_MLE=float(m.lfcSE), lfc_shrunk=float(m.lfc_shrunk),
                        lfcSE_shrunk=float(m.lfcSE_shrunk), padj=float(m.padj),
                        shrink_ratio=float(m.lfc_shrunk / m.log2FoldChange)
                        if m.log2FoldChange else np.nan,
                        shrink_converged=m.shrink_converged,
                        lfc_shrunk_altparam=float(m.lfc_shrunk_altparam),
                        shrink_converged_altparam=m.shrink_converged_altparam,
                        altparam_absdiff=abs(float(m.lfc_shrunk) -
                                             float(m.lfc_shrunk_altparam))))
        common = [g for g in fz.gene if g in r.index and (base is None or g in base.index)]
        sign_kept = sum(np.sign(base.loc[g, "log2FoldChange"]) == np.sign(r.loc[g, "log2FoldChange"])
                        for g in common)
        keep01 = sum(bool(r.loc[g, "padj"] < 0.1) for g in common if pd.notna(r.loc[g, "padj"]))
        cs_b = base.loc[base.index.isin(fz.gene), "log2FoldChange"]
        cs_r = r.loc[r.index.isin(fz.gene), "log2FoldChange"]
        sh = cs_b.index.intersection(cs_r.index)
        summ.append(dict(celltype=ct, design=design, estimable=True, n_donor=len(md),
                         n_frozen=len(fz), n_evaluable=len(common), sign_retained=int(sign_kept),
                         padj01_retained=int(keep01),
                         lfc_corr=round(float(cs_b.loc[sh].corr(cs_r.loc[sh])), 4), note=""))
        for g in common:
            perhit.append(dict(celltype=ct, design=design, gene=g,
                               lfc=round(float(r.loc[g, "log2FoldChange"]), 4),
                               padj=float(r.loc[g, "padj"])))
        print(f"  {design:26s} n={len(md)} sign {sign_kept}/{len(common)} "
              f"padj<0.1 {keep01}/{len(common)}", flush=True)

S = pd.DataFrame(summ)
S.to_csv(f"{OUT}/covariate_models_summary.csv", index=False)
pd.DataFrame(perhit).to_csv(f"{OUT}/covariate_models_perhit.csv", index=False)
SH = pd.DataFrame(shrunk_rows)
SH.to_csv(f"{OUT}/lfc_shrinkage_frozen60.csv", index=False)

print("\n" + S.to_string(index=False))
ok = S[S.estimable == True]
print("\n=== frozen-60 retention by design ===")
print(ok.groupby("design").agg(evaluable=("n_evaluable", "sum"), sign=("sign_retained", "sum"),
                               padj01=("padj01_retained", "sum")).to_string())
print("\nG-8 REPRODUCTION GATE: ~grp must give 60/60, ~pmi + grp 30/60, ~arm + grp 42/60")

if len(SH):
    print(f"\n=== LFC shrinkage on the frozen {len(SH)} ===")
    print(f"MLE vs frozen CSV: max |diff| = {(SH.lfc_MLE - SH.lfc_frozen).abs().max():.3e}")
    print(f"median |MLE| {SH.lfc_MLE.abs().median():.3f} -> "
          f"median |shrunk| {SH.lfc_shrunk.abs().median():.3f}")
    print(f"|log2FC| > 1: MLE {(SH.lfc_MLE.abs() > 1).sum()}/{len(SH)} -> "
          f"shrunk {(SH.lfc_shrunk.abs() > 1).sum()}/{len(SH)}")
    print(f"sign preserved: {(np.sign(SH.lfc_MLE) == np.sign(SH.lfc_shrunk)).sum()}/{len(SH)}")
    nconv = int((SH.shrink_converged == False).sum())
    print(f"apeGLM L-BFGS-B did NOT converge for {nconv}/{len(SH)} frozen hits "
          f"(and {int((SH.shrink_converged_altparam == False).sum())}/{len(SH)} in the "
          f"alternative parameterisation)")
    print(f"the two parameterisations disagree by > 0.05 log2FC for "
          f"{int((SH.altparam_absdiff > 0.05).sum())}/{len(SH)} hits")
    bad = SH[(SH.shrink_converged == False) | (SH.altparam_absdiff > 0.05)]
    if len(bad):
        print("\n=== hits whose shrunken estimate is NOT trustworthy ===")
        print(bad[["celltype", "gene", "baseMean", "lfc_MLE", "lfcSE_MLE", "lfc_shrunk",
                   "lfc_shrunk_altparam", "shrink_converged", "shrink_converged_altparam",
                   "padj"]].round(4).to_string(index=False))
    named = ["SERPINE1", "IL15", "CCN2", "ANXA1", "TNC", "COL21A1", "APOD", "WNT5B", "KDR",
             "PLAT", "LAMA4", "APOC1", "LGALS9", "NPC2", "APOE", "ECM2", "IGFBP5",
             "COL12A1", "KITLG"]
    print("\n=== genes named in the Results prose ===")
    print(SH[SH.gene.isin(named)][["celltype", "gene", "baseMean", "lfc_MLE", "lfcSE_MLE",
                                   "lfc_shrunk", "shrink_ratio", "padj"]]
          .round(4).to_string(index=False))
print(f"\nwrote {OUT}/covariate_models_{{summary,perhit}}.csv, {OUT}/lfc_shrinkage_frozen60.csv")
