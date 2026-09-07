# scripts/10_resilience — CONSEQUENCE (Phase A) receiver-competence resilience driver

The generating driver for the **CONSEQUENCE / Phase A** arm (GSE325391 granule-lineage receiver-competence
resilience test). It produces `results/resilience/*` (Fig S(RT1) panel B + the resilience Methods/Results;
verdict **NOT SUPPORTED** — stage-averaged SAD−RES = −0.07 z, 95% CI [−0.335, 0.195], p = 0.54).

## ⚠ Provenance (read this)
This arm was originally run **inline** (Bash heredocs) in Claude Code session `027c5fdb-8f9a-4e9c-972a-a9a75b8704a6`
(2026-06-19→21) and was **never saved as a script** — only its outputs were committed (tags
`result/consequence-phaseA-v1-final` 33165fc, `freeze/consequence-receptor-set-v1` be864a3). These scripts were
**recovered VERBATIM from that session transcript** (2026-06-23) and split into per-step `.py`/`.R` files. Each
step was verified by an extraction agent to match the committed `results/resilience/*` byte-for-byte (and the cited
transcript commands were spot-checked against the agent's extraction). They are the authoritative source of the
committed outputs; they were **not re-run** during recovery (see "Re-running" below).
One honest caveat is flagged in `04c_nichenet_supporting_targets.R`.

## Environments (per step; inline, two envs)
- `.py` steps → **bio** (`/home/neurofuture/miniforge3/envs/bio/bin/python`; scanpy/anndata/scipy/statsmodels)
- `.R` steps → **sgz_r** (`/home/neurofuture/miniforge3/envs/sgz_r/bin/Rscript`; nlme, fgsea). NOT trajR (that's the trajectory arm).

## Run order (run from `$PROJ`; DATA-DEPENDENT — see Inputs)
1. `00_ligand_receptor_map.R`            → results/step0/ (UP secreted 34 → NicheNet ligands 22 → receptors 76)
2. `01_freeze_receptor_set.py`           → frozen_receiver_receptor_set.csv (76 → **56**) + manifest + excluded + edges
3. `02_receiver_score.py`                → donor_stage_receptor_scores.tsv (26 rows; RES 6 / SAD 7)
4. `03_primary_model.R`                  → primary_model_results.tsv (SAD−RES −0.07, p 0.54) + secondary_contrasts.tsv
5. `03b_robustness_build.py` → `03b_robustness_fit.R`              → robustness_summary.tsv
6. `03c_robustness_extended_build.py` → `03c_robustness_extended_fit.R` → robustness_{run_fixed,no_run,leave_one_run,leave_one_ligand_family,subtype_specific}.tsv
7. `04a_nichenet_targets.R` → `04b_nichenet_enrichment_build.py` → `04b_nichenet_enrichment_fgsea.R` → nichenet_top{20,50,100}_results.tsv + nichenet_method_sensitivity.tsv (0/22 ligands)
8. `04c_nichenet_supporting_targets.R` → `04c_nichenet_supporting_wilcox.py` → nichenet_supporting_results.tsv
9. `05_matrix_provenance_audit.py`       → matrix_provenance_audit.txt + matrix_value_diagnostics.tsv (VERDICT PASS)
10. `03d_robustness_index.sh`            → robustness_index.tsv (presence table; run last)

## Inputs (gitignored / not in a bare clone)
- `processed/per_dataset/GSE325391_resilience.h5ad` + `GSE325391_adultgc_metadata.csv` — built by the committed
  `scripts/02_qc/06_build_325391_resilience.py` from `raw/GSE325391_immature_resilience/*.RDS` (3.5 GB; on disk, not committed).
- `refs/nichenet/{lr_network_human_21122021.rds, ligand_target_matrix_nsga2r_final.rds}`.
- `results/validation/frozen_primary_signature.csv`, `results/de/de_GSE325391_RESvSAD.csv`.

## Re-running (heavy; feasibility-gated)
`make resilience` documents the order. It is data-dependent and requires first rebuilding the h5ad from the 3.5 GB
RDS (`scripts/02_qc/06_build_325391_resilience.py`) — a multi-GB computation. It was **not** run during recovery;
re-run to fully re-verify (esp. the `nichenet_supporting_results.tsv` per-ligand rows noted in 04c).
