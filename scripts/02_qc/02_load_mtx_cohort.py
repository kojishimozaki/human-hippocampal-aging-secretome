#!/usr/bin/env python
"""Generic loader: MatrixMarket counts + cell metadata + gene list -> standardized .h5ad.

Raw counts go to layers['counts'] and X. Matrix orientation is auto-detected by
matching dims to (n_genes, n_cells). obs is standardized to:
  celltype_l2, donor_id, [region]  (+ whatever else is in the meta file)
Age/sex/diagnosis are mapped separately (per-cohort sample table) at DE time.

Usage:
  PROJ=<dir> python 02_load_mtx_cohort.py \
    --mtx raw/.../X.mtx.gz --genes raw/.../genes.txt.gz --meta raw/.../cell_meta.txt.gz \
    --celltype-col cluster --donor-col samplename --region-col region \
    --out processed/per_dataset/<cohort>.h5ad
"""
import argparse, gzip, os
import pandas as pd, anndata as ad
import scipy.io


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_lines(p):
    with opener(p) as f:
        return [l.rstrip("\n") for l in f]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mtx", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--cellname-col", default=None, help="default: first column of meta")
    ap.add_argument("--celltype-col", required=True)
    ap.add_argument("--donor-col", required=True)
    ap.add_argument("--region-col", default=None)
    ap.add_argument("--gene-col", type=int, default=0, help="column index if genes file is multi-col")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    print("reading mtx:", a.mtx, flush=True)
    with opener(a.mtx) as fh:
        X = scipy.io.mmread(fh).tocsr()
    print("  mtx shape:", X.shape, flush=True)

    genes = read_lines(a.genes)
    if genes and "\t" in genes[0]:
        genes = [g.split("\t")[a.gene_col] for g in genes]
    meta = pd.read_csv(a.meta, sep="\t")
    cellname = a.cellname_col or meta.columns[0]
    meta = meta.set_index(cellname)
    ng, nc = len(genes), len(meta)
    print(f"  n_genes={ng}, n_cells={nc}", flush=True)

    if X.shape == (ng, nc):
        X = X.T.tocsr()
    elif X.shape == (nc, ng):
        pass
    else:
        raise ValueError(f"mtx {X.shape} matches neither (genes,cells)=({ng},{nc}) nor its transpose")

    adata = ad.AnnData(X=X)
    adata.var_names = pd.Index(genes, dtype=str)
    adata.obs_names = meta.index.astype(str)
    adata.var_names_make_unique()
    adata.obs["celltype_l2"] = meta[a.celltype_col].astype(str).values
    adata.obs["donor_id"] = meta[a.donor_col].astype(str).values
    if a.region_col and a.region_col in meta.columns:
        adata.obs["region"] = meta[a.region_col].astype(str).values
    adata.layers["counts"] = adata.X.copy()

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    adata.write(a.out)
    print("wrote:", a.out, "->", adata, flush=True)
    print("celltype_l2 (top20):", adata.obs["celltype_l2"].value_counts().head(20).to_dict())
    if "region" in adata.obs:
        print("region:", adata.obs["region"].value_counts().to_dict())
    print("donor_id:", adata.obs["donor_id"].value_counts().to_dict())


if __name__ == "__main__":
    main()
