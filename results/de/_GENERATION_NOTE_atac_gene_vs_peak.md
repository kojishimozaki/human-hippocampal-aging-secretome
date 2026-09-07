# Provenance — gene-level vs peak-level accessibility null (Supplementary Table S2)

`atac_gene_vs_peak_level_matched_null.csv` is the source table for **Supplementary Table S2**. One row
per (cell type × gene set); it carries the **peak-level** and the **gene-level** version of the same
accessibility-matched null side by side, which is what the Results mean by "peak-level p-values are
reported for comparison". `atac_gene_vs_peak_level_per_gene.csv` is the per-gene detail behind the
gene-level columns (5,018 rows).

Traces for the manuscript's regulatory paragraph:

| Manuscript claim | Row / column |
|---|---|
| astrocytes, 26 genes from 160 peaks | Astro/frozen_hit_linked `n_genes` 26, `n_peaks` 160 |
| matched-null p = 5.0×10⁻⁴, sign-permutation p = 1.0×10⁻³, leave-one-gene-out p = 5.0×10⁻⁴ | `p_gene_level` 0.0005, `p_gene_signperm` 0.001, `p_gene_LOO` 0.0005 |
| 20 of 26 genes concordant | `genes_concordant` 20 |
| all-secretome concordant in all five types (5.0×10⁻⁴ in Astro/Micro/Oligo/OPC; 0.027 in Endo) | `p_gene_level` on the `all_secretome_linked` rows (Endo 0.0265) |
| microglia and OPCs fail the gene-level sign permutation (p = 0.077, 0.072) | `p_gene_signperm` |
| the OPC foreground is effectively one gene (FGF13, 13 of 19 peaks) | OPC/frozen_hit_linked `top_gene` FGF13, `top_gene_peaks` 13, `n_peaks` 19 |
| the pre-specified go/no-go fails in all ten rows | `results/de/atac_linkedpeak_enrichment.csv` `pass` = False on all 10; `n_fdr_concordant` = 0 |

Why the gene is the unit: `peaks_per_gene` is 3.2–6.2 and one gene can supply a third to two thirds of
a set's peaks (`top_gene_share`; TAFA1 0.331 in Astro, FGF13 0.684 in OPC), so peak-level rows are not
independent. `null_sd_ratio_gene_over_peak` (1.23–1.85) quantifies how much the peak-level null
under-disperses as a result — that is the inflation the gene-level collapse removes.

## How it was made — and the gap
Produced by the **third-party audit run of 2026-08-21…25** (`~/audit-runs/20260821-225834/audit/
pass3b_gene_level_matched_null.csv` and `pass3b_gene_level_per_gene.csv`), imported verbatim on
2026-08-26; no value recomputed.

⚠ **The driver script is not committed on any branch** — only its outputs are (same residual class as
`results/audit_reruns/de_calibration/`). The peak-level half of the comparison IS reproducible from
committed code: `make atac-concordance` → `scripts/04_atac/01_atac_rna_concordance.py`, whose outputs
`results/de/atac_linkedpeak_enrichment{,_controls}.csv` agree with the `obs_peak_level` /
`p_peak_level` columns here.
