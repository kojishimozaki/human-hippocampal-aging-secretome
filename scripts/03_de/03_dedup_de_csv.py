#!/usr/bin/env python
"""Collapse duplicate (gene, celltype) rows in an existing per-celltype DE CSV.

Duplicate gene SYMBOLS arise when two Ensembl IDs map to one HGNC symbol (e.g.
read-through genes) and var_names were not made unique before pseudobulk, so the
same symbol appears as two columns and two DE rows. This applies the same
invariant now enforced by scripts/03_de/01_pseudobulk_de.py at write time
(one row per gene,celltype; keep the higher-baseMean / expressed copy) to a CSV
that was generated before that fix.

Audit fix 2026-06-03: affected only de_GSE268609_per_celltype.csv (26 ambiguous,
non-secretome symbol pairs). No figure or validation number depends on the
dropped low-baseMean copies. Usage: 03_dedup_de_csv.py results/de/<file>.csv
"""
import sys
import pandas as pd

path = sys.argv[1] if len(sys.argv) > 1 else \
    "/home/neurofuture/Bioanalysis/01.SGZ_aging_project/results/de/de_GSE268609_per_celltype.csv"
df = pd.read_csv(path)
before = len(df)
dups = df.duplicated(subset=["gene", "celltype"], keep=False)
n_pairs = df[dups].groupby(["gene", "celltype"]).ngroups
keep_idx = (df.sort_values("baseMean", ascending=False)
              .drop_duplicates(subset=["gene", "celltype"], keep="first")).index
out = df.loc[df.index.isin(keep_idx)].reset_index(drop=True)   # survivors, original order
out.to_csv(path, index=False)
print(f"{path}")
print(f"  {n_pairs} duplicate (gene,celltype) pairs; rows {before} -> {len(out)} "
      f"(dropped {before - len(out)} lower-baseMean copies)")
