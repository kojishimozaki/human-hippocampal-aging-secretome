#!/usr/bin/env python3
"""Multiple-testing scope: five per-cell-type BH families versus one genome-wide family.

The primary analysis applies Benjamini-Hochberg correction within each niche cell
type, so the 99,544 GSE268609 tests form five families rather than one. This script
reports what happens under the alternative scope: a single BH family across all five
cell types, intersected with the secretome universe.

Backs the Results sentence "Treating all five as a single genome-wide family retains
49 of the 60 reported pairs and admits 4 new ones (AZGP1, LTBP3, NRG1 and NXPE3, all
in oligodendrocytes)" and the corresponding Methods sentence.

Bare-clone runnable: reads only committed results/ and refs/ tables.

Usage:
    PROJ=$(pwd) python scripts/03_de/03_global_bh_scope.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests

PROJ = Path(os.environ.get("PROJ", Path(__file__).resolve().parents[2]))
DE = PROJ / "results/de/de_GSE268609_per_celltype.csv"
FROZEN = PROJ / "results/validation/frozen_primary_signature.csv"
SECRETOME = PROJ / "refs/secretome_union.csv"
OUT_DIR = PROJ / "results/de"
ALPHA = 0.1

de = pd.read_csv(DE)
frozen = pd.read_csv(FROZEN)
secretome = set(pd.read_csv(SECRETOME)["gene"])

# The deposited table is the frozen primary DE run; every row carries a p value.
assert de["pvalue"].notna().all(), "unexpected NaN p values in the primary DE table"
n_tests = len(de)

de["padj_global"] = multipletests(de["pvalue"].to_numpy(), method="fdr_bh")[1]
de["in_secretome"] = de["gene"].isin(secretome)
de["pair"] = list(zip(de["gene"], de["celltype"]))

frozen_pairs = set(zip(frozen["gene"], frozen["celltype"]))
assert len(frozen_pairs) == 60, f"expected 60 frozen pairs, got {len(frozen_pairs)}"

# Sanity check: the committed DESeq2 padj column must reproduce exactly the frozen 60.
per_ct = set(de.loc[(de["padj"] < ALPHA) & de["in_secretome"], "pair"])
assert per_ct == frozen_pairs, "per-cell-type BH does not reproduce the frozen 60"

global_hits = de[(de["padj_global"] < ALPHA) & de["in_secretome"]].copy()
global_pairs = set(global_hits["pair"])

retained = frozen_pairs & global_pairs
lost = frozen_pairs - global_pairs
entered = global_pairs - frozen_pairs

global_hits["scope_status"] = [
    "retained_from_frozen60" if p in frozen_pairs else "new_under_global_bh"
    for p in global_hits["pair"]
]
lost_rows = de[de["pair"].isin(lost)].copy()
lost_rows["scope_status"] = "lost_under_global_bh"

cols = [
    "gene", "celltype", "baseMean", "log2FoldChange", "lfcSE",
    "pvalue", "padj", "padj_global", "scope_status",
]
out = pd.concat([global_hits[cols], lost_rows[cols]], ignore_index=True)
out = out.sort_values(["scope_status", "celltype", "gene"]).reset_index(drop=True)

OUT_DIR.mkdir(parents=True, exist_ok=True)
out_path = OUT_DIR / "global_bh_scope_GSE268609.csv"
out.to_csv(out_path, index=False)

print(f"tests in the primary DE table          : {n_tests:,}")
print(f"per-cell-type BH < {ALPHA} (secretome)      : {len(per_ct)}")
print(f"single genome-wide BH < {ALPHA} (secretome) : {len(global_pairs)}")
print(f"  of the frozen 60, retained           : {len(retained)}")
print(f"  of the frozen 60, lost               : {len(lost)}")
print(f"  new pairs entering                   : {len(entered)}")
for gene, ct in sorted(entered):
    print(f"    + {gene} ({ct})")
print(f"wrote {out_path.relative_to(PROJ)}")
