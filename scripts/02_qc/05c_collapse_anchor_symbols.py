#!/usr/bin/env python
"""Bring an existing (pre-fix) GSE268609 anchor in line with the fixed build, in place.

The current anchor predates the 2026-06-03 audit fix: it carries 10 duplicate gene
symbols as GENE/GENE-1 (var_names_make_unique on read-through genes) and a normalized
X. A full from-scratch rebuild (05 + 05b) would instead collapse those 10 symbols to
their highest-count copy and store raw counts in X. This script applies exactly that
collapse to the existing file WITHOUT the heavy rebuild — and, crucially, WITHOUT
recomputing the UMAP, so the verified Fig 1 A/B embeddings are preserved.

Idempotent. Affects only the 10 colliding, non-secretome symbols (near-empty -1 copies
dropped); no figure/DE/validation number changes. Run once; figures then regenerate
identically. Backs nothing up itself — back up the anchor before running.
"""
import os
import gzip, re
import scanpy as sc, pandas as pd, numpy as np

P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
ANCHOR = f"{P}/processed/per_dataset/GSE268609_anchor.h5ad"
FEATURES = f"{P}/raw/GSE268609_lazarov2024/GSE268609_features.tsv.gz"

# The symbols that collided (>1 Ensembl 'Gene Expression' ID) -> these are the only ones
# var_names_make_unique() suffixed as SYMBOL / SYMBOL-1 in the pre-fix anchor.
e2s = {}
with gzip.open(FEATURES, "rt") as f:
    for line in f:
        p = line.rstrip().split("\t")
        if len(p) >= 3 and p[2] == "Gene Expression":
            e2s[p[0]] = p[1]
collided = {s for s, c in pd.Series(list(e2s.values())).value_counts().items() if c > 1}

A = sc.read_h5ad(ANCHOR)
n_before = A.n_vars
vn = pd.Index(A.var_names)
tot = np.asarray(A.layers["counts"].sum(0)).ravel()

drop_cols = []
for sym in sorted(collided):
    variants = [v for v in vn if v == sym or re.fullmatch(re.escape(sym) + r"-\d+", v)]
    if len(variants) <= 1:
        continue
    pos = [vn.get_loc(v) for v in variants]
    best = pos[int(np.argmax(tot[pos]))]                    # highest-count copy -> keep as bare symbol
    drop_cols += [p for p in pos if p != best]
    if vn[best] != sym:                                     # rename the survivor to the bare symbol
        names = list(A.var_names); names[best] = sym; A.var_names = pd.Index(names, dtype=str)

if drop_cols:
    keep = np.array([i for i in range(A.n_vars) if i not in set(drop_cols)])
    A = A[:, keep].copy()
A.X = A.layers["counts"].copy()                             # match from-scratch (raw counts in X)
if A.raw is not None:
    del A.raw
A.write(ANCHOR)
print(f"collapsed {len(collided)} colliding symbols: {n_before} -> {A.n_vars} vars "
      f"(dropped {n_before - A.n_vars}); X reset to raw counts; "
      f"X_umap preserved={'X_umap' in A.obsm}; dup symbols now={int(pd.Index(A.var_names).duplicated().sum())}")
