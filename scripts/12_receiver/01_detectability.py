#!/usr/bin/env python
"""receiver-map STEP 01 — direction-blind per-cell-type detectability gate on the frozen 76.

Pre-registration: docs/RECEIVER_MAP_PLAN.md (sec 4), tag prereg/receiver-map-v1.  env: bio.
Locked input: results/de/de_GSE268609_neuron_aging_per_celltype.csv (~arm + grp, single contrast).

Detectable/testable in compartment C  <=>  the receptor has a row in C with non-NA padj
(the pipeline's expressed-gene / independent-filtering criterion). Gate is applied per compartment
and is BLIND to the sign/magnitude of log2FC and to padj/p magnitude. Reported: n_detectable/76 per
compartment ("can this receiver receive?"). No inference here.
"""
import os, json, pandas as pd

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
os.chdir(PROJ)
OUT = "results/receiver"
DE = "results/de/de_GSE268609_neuron_aging_per_celltype.csv"

CONFIRMATORY = ["Astro", "Micro", "Oligo", "OPC", "Endo", "DG_GC", "CA_ExN", "InN"]
NEUROGENIC = ["NSC", "Neuroblast", "Immature"]
COMPARTMENTS = CONFIRMATORY + NEUROGENIC

rec76 = pd.read_csv(f"{OUT}/_universe76_from_manifest.csv")["receptor"].tolist()
de = pd.read_csv(DE, usecols=["gene", "baseMean", "log2FoldChange", "padj", "celltype"])

rows = []
summ = []
for ct in COMPARTMENTS:
    d = de[de.celltype == ct]
    n_genes_in_table = len(d)
    n_detect_global = d["padj"].notna().sum()
    sub = d[d.gene.isin(rec76)].set_index("gene")
    n_in_table = 0; n_detect = 0; bm_detect = []
    for r in rec76:
        if r in sub.index:
            row = sub.loc[r]
            in_tab = True
            det = pd.notna(row["padj"])
            bm = float(row["baseMean"]); lfc = float(row["log2FoldChange"]) if pd.notna(row["log2FoldChange"]) else None
            padj = float(row["padj"]) if det else None
            n_in_table += 1
            if det: n_detect += 1; bm_detect.append(bm)
        else:
            in_tab = False; det = False; bm = 0.0; lfc = None; padj = None
        rows.append(dict(receptor=r, compartment=ct, in_table=in_tab, detectable=det,
                         baseMean=bm, log2FoldChange=lfc, padj=padj))
    med_bm = float(pd.Series(bm_detect).median()) if bm_detect else None
    summ.append(dict(compartment=ct,
                     family=("confirmatory" if ct in CONFIRMATORY else "neurogenic"),
                     n_in_table=n_in_table, n_detectable=n_detect, n_universe=76,
                     frac_detectable=round(n_detect / 76, 4),
                     median_baseMean_detectable=med_bm,
                     n_genes_in_table=int(n_genes_in_table),
                     global_detectable_frac=round(float(n_detect_global) / n_genes_in_table, 4)))

matrix = pd.DataFrame(rows)
matrix.to_csv(f"{OUT}/detectability_matrix.csv", index=False)
S = pd.DataFrame(summ)
S.to_csv(f"{OUT}/detectability_summary.csv", index=False)

print("STEP 01 detectability (frozen 76, direction-blind, per compartment)")
print(S.to_string(index=False))
print(f"\nwrote: {OUT}/detectability_matrix.csv ({len(matrix)} rows), {OUT}/detectability_summary.csv")
