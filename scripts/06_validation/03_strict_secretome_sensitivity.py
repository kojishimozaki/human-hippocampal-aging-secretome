#!/usr/bin/env python3
"""Sensitivity of the main results to a stricter definition of "secreted".

The secretome universe is the union of the Human Protein Atlas "predicted secreted"
annotation and reviewed human UniProtKB entries annotated as secreted (2,224 genes). A
union admits genes whose secreted status rests on one source only: FGF13, for instance,
is annotated as secreted by UniProt but is reported as intracellular in experimental work,
and KDR is a receptor. Editorial review raised this, and it matters because the gene
universe defines which pairs can be selected at all.

`refs/secretome_union.csv` already carries `is_core` — membership in the intersection of
both sources (1,813 genes) — so the stricter analysis needs no new annotation call and no
re-fitting: the differential expression was computed genome-wide per cell type and then
intersected with the universe, so tightening the universe only filters the hit list.

What is reported, fixed before running:
  1. how many of the 60 pairs survive the intersection universe, and which do not;
  2. whether external reproduction still holds on the surviving pairs, at the same
     direction-concordance metric used for the full set;
  3. whether the astrocyte gene-level RNA-chromatin concordance still holds;
  4. what the OPC accessibility result looks like with FGF13 removed, since the
     manuscript already reports that FGF13 contributes 13 of its 19 peaks.

This is a sensitivity analysis reported alongside the primary result, not a replacement
for it: the 60-pair set was frozen before any downstream analysis and seeded all of them.

Usage:
    PROJ=$(pwd) python scripts/06_validation/03_strict_secretome_sensitivity.py
"""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
from scipy import stats

PROJ = Path(os.environ.get("PROJ", Path(__file__).resolve().parents[2]))
OUT = PROJ / "results/validation/strict_secretome_sensitivity.csv"

frozen = pd.read_csv(PROJ / "results/validation/frozen_primary_signature.csv")
universe = pd.read_csv(PROJ / "refs/secretome_union.csv").set_index("gene")
frozen["is_core"] = universe.reindex(frozen.gene)["is_core"].to_numpy()
assert len(frozen) == 60, f"expected 60 frozen pairs, got {len(frozen)}"

rows: list[dict] = []


def rec(section: str, item: str, full: str, strict: str, note: str = "") -> None:
    rows.append({"section": section, "item": item, "full_union": full,
                 "core_intersection": strict, "note": note})


# --- 1. composition -------------------------------------------------------------------
n_core = int(frozen.is_core.sum())
dropped = frozen.loc[~frozen.is_core, ["gene", "celltype", "direction"]]
rec("composition", "pairs in the frozen set", "60", str(n_core),
    "core = HPA and UniProt both annotate the gene as secreted")
for ct, sub in frozen.groupby("celltype"):
    rec("composition", f"pairs in {ct}", str(len(sub)), str(int(sub.is_core.sum())))
rec("composition", "pairs dropped by the stricter universe", "-", str(len(dropped)),
    "; ".join(f"{g} ({c})" for g, c in zip(dropped.gene, dropped.celltype)))

# --- 2. external reproduction ---------------------------------------------------------
ext = pd.read_csv(PROJ / "results/de/de_GSE278576_per_celltype.csv")
ext_lfc = ext.set_index(["gene", "celltype"])["log2FoldChange"]
frozen["external_lfc"] = [ext_lfc.get((g, c)) for g, c in zip(frozen.gene, frozen.celltype)]
evaluable = frozen.dropna(subset=["external_lfc"]).copy()
evaluable["concordant"] = evaluable.log2FoldChange * evaluable.external_lfc > 0


def concord(df: pd.DataFrame) -> str:
    if not len(df):
        return "n/a"
    return f"{int(df.concordant.sum())}/{len(df)} = {df.concordant.mean() * 100:.0f}%"


strict_ev = evaluable[evaluable.is_core]
rec("external reproduction (GSE278576)", "evaluable pairs",
    str(len(evaluable)), str(len(strict_ev)))
rec("external reproduction (GSE278576)", "direction concordance",
    concord(evaluable), concord(strict_ev),
    "background concordance for the external secretome set is 58%")
b_full = stats.binomtest(int(evaluable.concordant.sum()), len(evaluable), 0.578,
                         alternative="greater").pvalue
b_core = stats.binomtest(int(strict_ev.concordant.sum()), len(strict_ev), 0.578,
                         alternative="greater").pvalue
rec("external reproduction (GSE278576)", "binomial P versus the 58% background",
    f"{b_full:.4f}", f"{b_core:.4f}",
    "one-sided; the matched-null analysis in the manuscript is the stricter comparison")

# --- 3. gene-level RNA-chromatin concordance -----------------------------------------
atac = pd.read_csv(PROJ / "results/de/atac_gene_vs_peak_level_per_gene.csv")
core_genes = set(universe.index[universe["is_core"].fillna(False).astype(bool)])
for ct in ("Astro", "Micro", "OPC"):
    sub = atac[(atac.celltype == ct) & (atac["set"] == "frozen_hit_linked")]
    sub_core = sub[sub.gene.isin(core_genes)]
    rec("gene-level RNA-chromatin concordance", f"{ct}: concordant genes",
        f"{int((sub.oriented_mean_lfc > 0).sum())}/{len(sub)}",
        f"{int((sub_core.oriented_mean_lfc > 0).sum())}/{len(sub_core)}")

# --- 4. the OPC result without FGF13 --------------------------------------------------
opc = atac[(atac.celltype == "OPC") & (atac["set"] == "frozen_hit_linked")]
opc_no = opc[opc.gene != "FGF13"]
rec("OPC accessibility", "linked peaks", str(int(opc.n_peaks.sum())),
    str(int(opc_no.n_peaks.sum())), "FGF13 removed")
rec("OPC accessibility", "concordant genes",
    f"{int((opc.oriented_mean_lfc > 0).sum())}/{len(opc)}",
    f"{int((opc_no.oriented_mean_lfc > 0).sum())}/{len(opc_no)}", "FGF13 removed")

out = pd.DataFrame(rows)
OUT.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT, index=False)
width = max(len(r["item"]) for r in rows)
print(f"{'':<26}{'item':<{width}}  {'full union':>16}  {'core intersection':>18}")
last = None
for r in rows:
    head = r["section"] if r["section"] != last else ""
    last = r["section"]
    print(f"{head:<26}{r['item']:<{width}}  {r['full_union']:>16}  {r['core_intersection']:>18}")
    if r["note"]:
        print(f"{'':<26}{'':<{width}}  -> {r['note']}")
print(f"\nwrote {OUT.relative_to(PROJ)}")
