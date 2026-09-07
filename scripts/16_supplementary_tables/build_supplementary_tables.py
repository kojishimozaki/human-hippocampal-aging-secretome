#!/usr/bin/env python
"""Build the three supplementary tables named in manuscript/manuscript.md.

The manuscript's "Supplementary tables" section describes S1, S2 and S(R1) and names their
source CSVs, but the deliverable tables themselves were never assembled. This script assembles
them from those committed sources -- it computes nothing new, so every number in the output is
traceable to a file already in the repository.

  Table S1    donor-label permutation calibration of the DE hit counts
  Table S2    peak-level versus gene-level accessibility-matched null
  Table S(R1) per-(TF x cell type) motif-accessibility versus TF-RNA concordance,
              HA-vs-YA primary AND AD-vs-HA projection stacked in one table with a `contrast`
              column, because the manuscript describes one table covering both

OUTPUT  manuscript/supplementary_tables/Table_S1.csv
        manuscript/supplementary_tables/Table_S2.csv
        manuscript/supplementary_tables/Table_SR1.csv      (no parentheses: Makefile-safe)
        manuscript/supplementary_tables/Supplementary_Tables.xlsx   (README + 3 sheets)

RUNTIME seconds.  Regenerate with:  make supplementary-tables
"""
import os
import pandas as pd

PROJ = os.environ.get("PROJ", os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
OUT = f"{PROJ}/manuscript/supplementary_tables"
os.makedirs(OUT, exist_ok=True)

# ---- Table S1 --------------------------------------------------------------------------
S1_SRC = f"{PROJ}/results/audit_reruns/de_calibration/donor_label_permutation_summary.csv"
s1 = pd.read_csv(S1_SRC)
S1_ORDER = ["cohort", "celltype", "status", "n_donors_YA", "n_donors_HA", "n_genes_tested",
            "n_perm", "n_distinct_label_splits",
            "obs_secretome_fdr01", "null_median_secretome_fdr01", "emp_p_secretome_fdr01",
            "obs_genomewide_fdr005", "null_median_genomewide_fdr005", "null_p95_genomewide_fdr005",
            "null_max_genomewide_fdr005", "emp_p_genomewide_fdr005",
            "null_frac_ge1_genomewide_fdr005", "null_frac_ge10_genomewide_fdr005",
            "obs_genomewide_fdr01", "emp_p_genomewide_fdr01", "verdict"]
assert set(S1_ORDER) == set(s1.columns), set(S1_ORDER) ^ set(s1.columns)
s1 = s1[S1_ORDER].sort_values(["cohort", "celltype"], kind="stable")

# ---- Table S2 --------------------------------------------------------------------------
S2_SRC = f"{PROJ}/results/de/atac_gene_vs_peak_level_matched_null.csv"
ENR_SRC = f"{PROJ}/results/de/atac_linkedpeak_enrichment.csv"
s2 = pd.read_csv(S2_SRC)
enr = pd.read_csv(ENR_SRC)[["celltype", "set", "n_fdr_concordant", "pass"]]
enr = enr.rename(columns={"n_fdr_concordant": "n_peaks_fdr_significant",
                          "pass": "prespecified_gonogo_pass"})
s2 = s2.merge(enr, on=["celltype", "set"], how="left", validate="one_to_one")
S2_ORDER = ["celltype", "set", "n_genes", "n_peaks", "peaks_per_gene",
            "top_gene", "top_gene_peaks", "top_gene_share",
            "obs_peak_level", "p_peak_level",
            "obs_gene_level", "p_gene_level", "p_gene_signperm",
            "genes_concordant", "obs_gene_LOO", "p_gene_LOO",
            "null_sd_peak", "null_sd_gene", "null_sd_ratio_gene_over_peak",
            "n_peaks_fdr_significant", "prespecified_gonogo_pass"]
assert set(S2_ORDER) == set(s2.columns), set(S2_ORDER) ^ set(s2.columns)
s2 = s2[S2_ORDER].sort_values(["set", "celltype"], kind="stable")

# ---- Table S(R1) -----------------------------------------------------------------------
R = f"{PROJ}/results/regulatory"
prim = pd.read_csv(f"{R}/motif_tf_rna_concordance_HAvYA_primary_pairs.csv")
proj = pd.read_csv(f"{R}/motif_tf_rna_concordance_ADvHA_projection_pairs.csv")
assert list(prim.columns) == list(proj.columns)
prim.insert(0, "contrast", "HAvYA_primary")
proj.insert(0, "contrast", "ADvHA_projection")
sr1 = pd.concat([prim, proj], ignore_index=True)
sr1 = sr1.sort_values(["contrast", "celltype", "gene"], kind="stable")

# ---- write ------------------------------------------------------------------------------
TABLES = [("Table_S1", s1, "S1 | donor-label permutation calibration of the DE hit counts",
           "results/audit_reruns/de_calibration/donor_label_permutation_summary.csv"),
          ("Table_S2", s2, "S2 | peak-level versus gene-level accessibility-matched null",
           "results/de/atac_gene_vs_peak_level_matched_null.csv + atac_linkedpeak_enrichment.csv"),
          ("Table_SR1", sr1, "S(R1) | motif-accessibility versus TF-RNA concordance",
           "results/regulatory/motif_tf_rna_concordance_{HAvYA_primary,ADvHA_projection}_pairs.csv")]
for name, df, _t, _s in TABLES:
    df.to_csv(f"{OUT}/{name}.csv", index=False)
    print(f"wrote {name}.csv  ({df.shape[0]} rows x {df.shape[1]} cols)")

readme = pd.DataFrame(
    [{"sheet": n, "table": t, "rows": len(d), "columns": d.shape[1], "source_data": s}
     for n, d, t, s in TABLES])
with pd.ExcelWriter(f"{OUT}/Supplementary_Tables.xlsx", engine="openpyxl") as xw:
    readme.to_excel(xw, sheet_name="README", index=False)
    for name, df, _t, _s in TABLES:
        df.to_excel(xw, sheet_name=name, index=False)
print("wrote Supplementary_Tables.xlsx (README + 3 sheets)")
