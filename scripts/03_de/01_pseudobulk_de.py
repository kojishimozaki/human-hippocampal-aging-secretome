#!/usr/bin/env python
"""Per-cohort donor pseudobulk DE. Pseudoreplication-safe (collapse cells -> donor).

Two modes:
  --mode age   : continuous age slope. design ~[sex+]age_years; numeric contrast on age_years
  --mode group : 2-level categorical. design ~[sex+]grp; contrast [grp, test, ref]
Per-cohort power is limited (few donors) -> real inference is the cross-cohort
meta-analysis / direction concordance (scripts/03_de/02_meta_analysis.R).

Output: results/de/de_<COHORT>_per_celltype.csv
  cols: gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, celltype, cohort, n_donors
"""
import argparse, os
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def get_counts(adata):
    return adata.layers["counts"] if "counts" in adata.layers else adata.X


def pseudobulk(adata, donor_col):
    X = get_counts(adata)
    if sp.issparse(X):
        X = X.tocsr()
    donors = adata.obs[donor_col].astype(str).values
    uniq = pd.unique(donors)
    mat = np.vstack([np.asarray(X[donors == d].sum(0)).ravel() for d in uniq])
    return pd.DataFrame(mat, index=uniq, columns=adata.var_names.astype(str))


def run(adata, cohort, celltype_col, donor_col, celltypes, min_cells, min_donors,
        mode, group_col, test_level, ref_level):
    obs = adata.obs
    cts = celltypes or sorted(obs[celltype_col].dropna().astype(str).unique())
    results = []
    for ct in cts:
        sub = adata[obs[celltype_col].astype(str) == ct]
        if sub.n_obs == 0:
            continue
        nc = sub.obs.groupby(donor_col, observed=True).size()
        sub = sub[sub.obs[donor_col].isin(nc[nc >= min_cells].index)]
        agg = {"sex": ("sex", "first")}
        if mode == "age":
            agg["age_years"] = ("age_years", "first")
        else:
            agg["grp"] = (group_col, "first")
        meta = sub.obs.groupby(donor_col, observed=True).agg(**agg)
        if mode == "age":
            meta["age_years"] = pd.to_numeric(meta["age_years"], errors="coerce")
            meta = meta.dropna(subset=["age_years"])
        else:
            meta["grp"] = meta["grp"].astype(str)
            meta = meta[meta["grp"].isin([test_level, ref_level])]
        if mode != "age":
            # Methods state "≥ min_donors donors/GROUP" — enforce per group, not just on the
            # total. The previous total-only guard could admit an unbalanced 2-level contrast
            # (e.g. 3 vs 1 still has len(meta)==4). This does NOT change the primary GSE268609
            # contrast: every reported niche type has ≥7 donors/group (Astro 9/8, Micro 9/7,
            # Oligo 9/8, OPC 9/8, Endo 9/7). (Audit fix 2026-06-05.)
            per_grp = meta["grp"].value_counts()
            if min(per_grp.get(test_level, 0), per_grp.get(ref_level, 0)) < min_donors:
                print(f"  skip {ct}: per-group donors {dict(per_grp)} (<{min_donors}/group)")
                continue
        if len(meta) < min_donors:
            print(f"  skip {ct}: n_donors={len(meta)} (<{min_donors})")
            continue
        pb = pseudobulk(sub, donor_col).loc[meta.index]
        keepg = (pb.sum(0) >= 10) & (pb.astype(bool).sum(0) >= max(3, int(0.5 * len(pb))))
        pb = pb.loc[:, keepg]
        if pb.shape[1] < 50:
            print(f"  skip {ct}: too few genes ({pb.shape[1]})")
            continue
        has_sex = meta["sex"].astype(str).nunique() >= 2
        try:
            if mode == "age":
                design = "~sex + age_years" if has_sex else "~age_years"
                dds = DeseqDataSet(counts=pb.astype(int), metadata=meta, design=design, quiet=True)
                dds.deseq2()
                cols = list(dds.varm["LFC"].columns)
                contrast = np.array([1.0 if c == "age_years" else 0.0 for c in cols])
            else:
                design = "~sex + grp" if has_sex else "~grp"
                dds = DeseqDataSet(counts=pb.astype(int), metadata=meta, design=design, quiet=True)
                dds.deseq2()
                contrast = ["grp", test_level, ref_level]
            ds = DeseqStats(dds, contrast=contrast, quiet=True)
            ds.summary()
            r = ds.results_df.copy()
            r["gene"] = r.index; r["celltype"] = ct; r["cohort"] = cohort; r["n_donors"] = len(meta)
            results.append(r.reset_index(drop=True))
            print(f"  {ct}: n_donors={len(meta)}, genes={pb.shape[1]}, design='{design}', "
                  f"sig(padj<0.1)={(r.padj < 0.1).sum()}")
        except Exception as e:
            print(f"  FAIL {ct}: {type(e).__name__} {str(e)[:160]}")
    return pd.concat(results, ignore_index=True) if results else pd.DataFrame()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--celltype-col", default="celltype_l1")
    ap.add_argument("--donor-col", default="donor_id")
    ap.add_argument("--mode", choices=["age", "group"], default="age")
    ap.add_argument("--age-col", default="age_years")
    ap.add_argument("--sex-col", default="sex")
    ap.add_argument("--group-col", default="group")
    ap.add_argument("--test-level", default=None)
    ap.add_argument("--ref-level", default=None)
    ap.add_argument("--celltypes", nargs="*", default=None)
    ap.add_argument("--min-cells", type=int, default=10)
    ap.add_argument("--min-donors", type=int, default=4)
    ap.add_argument("--min-age", type=float, default=None)
    ap.add_argument("--max-age", type=float, default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    PROJ = os.environ.get("PROJ", os.getcwd())
    adata = sc.read_h5ad(a.h5ad)
    adata.obs["sex"] = adata.obs[a.sex_col].astype(str) if a.sex_col in adata.obs else "U"
    if a.mode == "age":
        adata.obs["age_years"] = pd.to_numeric(adata.obs[a.age_col], errors="coerce")
        adata = adata[adata.obs["age_years"].notna()].copy()
        if a.min_age is not None:
            adata = adata[adata.obs["age_years"] >= a.min_age].copy()
        if a.max_age is not None:
            adata = adata[adata.obs["age_years"] <= a.max_age].copy()
    de = run(adata, a.cohort, a.celltype_col, a.donor_col, a.celltypes, a.min_cells,
             a.min_donors, a.mode, a.group_col, a.test_level, a.ref_level)
    # Enforce one row per (gene, celltype). Duplicate gene SYMBOLS can survive when two
    # Ensembl IDs map to the same HGNC symbol (e.g. read-through genes) and var_names were
    # not made unique upstream; pseudobulk() then carries both as same-named columns. Keep
    # the higher-baseMean (expressed) copy and drop the near-empty duplicate so downstream
    # gene-keyed joins are unambiguous. (Audit fix 2026-06-03; affected only ambiguous
    # non-secretome symbols, no effect on figure/validation results.)
    if len(de):
        before = len(de)
        de = (de.sort_values("baseMean", ascending=False)
                .drop_duplicates(subset=["gene", "celltype"], keep="first")
                .reset_index(drop=True))
        if len(de) < before:
            print(f"  collapsed {before - len(de)} duplicate (gene,celltype) rows (kept higher baseMean)")
    out = a.out or os.path.join(PROJ, "results", "de", f"de_{a.cohort}_per_celltype.csv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    de.to_csv(out, index=False)
    print(f"\nWrote {out}: {len(de)} rows, {de['celltype'].nunique() if len(de) else 0} celltypes")
