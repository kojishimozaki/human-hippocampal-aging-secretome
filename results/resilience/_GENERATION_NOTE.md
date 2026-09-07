# results/resilience — generation provenance & a known reproducibility gap

Outputs of the **CONSEQUENCE (Phase A) receiver-competence resilience arm** (GSE325391), committed for the
manuscript (Fig S(RT1) panel B + the resilience Methods/Results) and for audit.

## Committed here (design / frozen inputs / outputs)
- **Pre-registered design:** `preregistration/manifests/consequence_design_v1.json` (tag `prereg/secretome-resilience-trajectory-v1`, 5b29c74).
- **Frozen receiver receptor set (pre-outcome freeze):** `frozen_receiver_receptor_set.csv`, `frozen_receiver_set_manifest.json` (tag `freeze/consequence-receptor-set-v1`, be864a3).
- **Primary result + full robustness / downstream tables:** `primary_model_results.tsv` (stage-averaged SAD−RES mixed model), per-donor scores `donor_stage_receptor_scores.tsv`; robustness `robustness_summary.tsv` + `robustness_{run_fixed,no_run,leave_one_run,leave_one_ligand_family,subtype_specific,index}.tsv` + `secondary_contrasts.tsv`; NicheNet downstream `nichenet_supporting_results.tsv` + `nichenet_{top20,top50,top100}_results.tsv` + `nichenet_method_sensitivity.tsv`; claim gate `claim_eligibility.md`. (All the variants/thresholds cited in `claim_eligibility.md` are now committed here — no longer tag-only.)
- **Audit provenance:** `matrix_provenance_audit.txt` + `matrix_value_diagnostics.tsv` (matrix-provenance PASS: X = raw int64 counts; CPM/log2 valid) and `deviation_log.md` (D1 nlme for the pre-registered REML mixed model; D2 NicheNet Wilcoxon + fgsea both reported).

## ✅ Driver recovered & committed (2026-06-23)
The **generating driver** for the receiver-competence model (the nlme `~group*stage + (1|Run/donor)` fit, the
receptor-set freeze, and the NicheNet downstream) was originally run **inline** (Bash heredocs) in origin session
`027c5fdb…` and was never saved as a file. It has now been **recovered VERBATIM from that session transcript** and
committed at **`scripts/10_resilience/`** (15 step scripts + README; bio `.py` + sgz_r `.R`). An extraction agent
verified each step matches the committed `results/resilience/*` byte-for-byte, and the cited transcript commands were
spot-checked. It was **not** re-run during recovery (data-dependent: needs the GSE325391 h5ad rebuilt from the 3.5 GB
raw RDS via `scripts/02_qc/06_build_325391_resilience.py`). Run order / envs / inputs and one honest caveat (the
`nichenet_supporting_results.tsv` per-ligand rows, STEP 04c) are in `scripts/10_resilience/README.md`. Re-verify
end-to-end with `make resilience`. (The TRAJECTORY drivers `scripts/11_trajectory/01–05` and the TF-RNA engine
`scripts/05_regulatory/06_motif_tf_rna_concordance.py` are likewise committed.)
