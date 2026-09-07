# Provenance — DE calibration tables (Supplementary Table S1)

## What these are
`donor_label_permutation_summary.csv` is the source table for **Supplementary Table S1** (donor-label
permutation nulls). `donor_label_permutation_null_draws.tsv` is the per-permutation draw table it
summarises (900 rows = 2 cohorts × cell types × 100 permutations).

Every number in the manuscript's DE-calibration paragraph traces to these files:

| Manuscript claim | Column / value here |
|---|---|
| secretome hit count above its null in 4 of 5 cell types, empirical p 0.020–0.049 | `emp_p_secretome_fdr01`: Micro 0.0198, Astro 0.0396, Endo 0.0396, OPC 0.0495 |
| but not oligodendrocytes (5 hits, p = 0.109) | `obs_secretome_fdr01` 5, `emp_p_secretome_fdr01` 0.1089 |
| the null median is zero in every cell type | `null_median_secretome_fdr01` = 0 throughout |
| genome-wide FDR<0.05: oligodendrocytes p = 0.139, vascular p = 0.069, inside the null's 95th percentile | `emp_p_genomewide_fdr005` 0.1386 / 0.0693 |
| GSE268609 (8 vs 9): ≥1 BH-significant gene 29–42%, ≥10 genes up to 16% | `null_frac_ge1_genomewide_fdr005` 0.29–0.42, `null_frac_ge10_genomewide_fdr005` max 0.16 |
| GSE278576 (10 vs 20): ≥1 gene 29–51%, ≥10 genes up to 17% | same columns, GSE278576 rows |
| ⚠ the two cohorts must not be pooled into one range — the manuscript quoted 29–51%/17% against the 8-vs-9 split until 2026-08-26 | |
| null draws reaching 1,331 genes here and 1,928 in GSE278576 | `null_max_genomewide_fdr005` |
| genome-wide astrocyte count 233 (→ 651 dropping one aged donor) | `obs_genomewide_fdr01` = 233 (the 651 is the leave-one-donor-out table, `results/de/aging_loo_robustness.csv` / Limitation 1) |
| empirical p floored at 1/101 = 0.0099 | `n_perm` = 100 |

## How they were made — and the gap
Produced by the **third-party audit run of 2026-08-21…25**, whose working tree lives OUTSIDE this
repository at `~/audit-runs/20260821-225834/audit/pass6/permutation/` (`summary.json` +
`null_distribution.tsv`). The files here are that output, imported verbatim on 2026-08-26 and
reshaped into one tidy row per cohort × cell type; no value was recomputed.

⚠ **The driver script for this permutation is not committed on any branch** — only its outputs are.
This is the same residual class as the consequence/receiver Phase-A engine (see
`results/resilience/_GENERATION_NOTE.md`), and it is disclosed rather than papered over. To re-derive:
100 donor-label permutations per cohort × cell type, refitting the full donor-pseudobulk pyDESeq2
model each time, counting BH-significant genes genome-wide (α = 0.05 and 0.1) and within the
secretome universe (α = 0.1). `n_distinct_label_splits` records how many distinct label assignments
the design admits (24,310 for the 8-versus-9 GSE268609 split), so 100 permutations do not exhaust it.

`GSE278576/Endo` is `GATE_FAIL`, not a null result: that cohort's endothelial compartment has 6 YA
and 2 HA donors and did not meet the run's minimum-donor gate. It is carried in the table as a row so
the absence is visible rather than silent.

## Reproduction check
The run re-fit each cohort × cell type before permuting and compared against the committed DE tables:
for GSE268609/Astro, `lfc_pearson` = 0.9999999996, `max_abs_lfc_diff` = 4.4×10⁻⁴, and the refit
recovered 233 genes at padj < 0.1 against the published 233. The per-cell-type checks are in the
source `summary.json`.
