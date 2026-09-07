#!/usr/bin/env python
"""Compute and store the GSE268609 anchor UMAP embedding (presentation prep).

Owns the embedding that Figure 1 renders. Previously the UMAP was computed and
written back to the anchor INSIDE scripts/07_figures/01_figure1.py, so a figure
script mutated analysed data (reproducibility hazard). Moved here 2026-06-03
(audit fix) as an explicit, idempotent prep step that runs once, after
05_build_268609_anchor.py and before the figures.

Embedding params are unchanged from the original Fig 1 code (HVG=2000 batch=donor_id,
scale max=10, PCA=30, neighbors k=15/30 PCs, UMAP). PCA/neighbors/UMAP use scanpy's
default random_state=0, so the embedding is deterministic. Idempotent: if X_umap is
already present it does nothing unless --force is given.
"""
import os
import argparse
import scanpy as sc, numpy as np

P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
ANCHOR = f"{P}/processed/per_dataset/GSE268609_anchor.h5ad"
SEED = 42

ap = argparse.ArgumentParser()
ap.add_argument("--force", action="store_true", help="recompute even if X_umap exists")
a = ap.parse_args()

sc.settings.seed = SEED
np.random.seed(SEED)
A = sc.read_h5ad(ANCHOR)

if str(A.var_names[0]).startswith("ENSG"):
    raise SystemExit("anchor var_names are Ensembl IDs — run 05_build_268609_anchor.py "
                     "(symbol mapping) first")
if "X_umap" in A.obsm and not a.force:
    print("X_umap already present; nothing to do (use --force to recompute)")
    raise SystemExit(0)

# normalize a working copy for HVG/PCA (counts are preserved in layers['counts'])
A.X = A.layers["counts"].copy()
sc.pp.normalize_total(A, target_sum=1e4)
sc.pp.log1p(A)
sc.pp.highly_variable_genes(A, n_top_genes=2000, batch_key="donor_id")
A2 = A[:, A.var.highly_variable].copy()
sc.pp.scale(A2, max_value=10)
sc.tl.pca(A2, n_comps=30)
sc.pp.neighbors(A2, n_neighbors=15, n_pcs=30)
sc.tl.umap(A2)

# write ONLY the embedding back; keep X as raw counts so the stored anchor is unchanged
# except for the new obsm['X_umap'] (do not persist normalized X or .raw).
A.obsm["X_umap"] = A2.obsm["X_umap"]
A.X = A.layers["counts"].copy()
if A.raw is not None:
    del A.raw
A.write(ANCHOR)
print(f"wrote X_umap to {ANCHOR}  ({A.n_obs} cells)")
