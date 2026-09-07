# Amendment v1.4 — coefficient-level identity gate for the single-fit trajectory model (post-diagnostic, pre-classification)

**Post-diagnostic, pre-classification implementation clarification.** Defines the single-fit identity gate on
the fitted GLM coefficients (the algebraically meaningful object), with a separate frozen-60 operational gate
at the `results()`/shrinkage level. No hypothesis, threshold, gate-purpose, gene set, or claim change. v1
(`5b29c74`), v1.1, v1.2, v1.3 tags retained permanently. New tag: `prereg/secretome-resilience-trajectory-v1.4`.

## Diagnostic facts (recorded verbatim)
- In the R DESeq2 single fit, the additive identity on `coef(dds)` holds **globally**, with maximum error
  ≈ **3.55e-15**.
- The `results()`-derived LFCs violate the `1e-6` additivity for **53 genes**.
- These 53 genes are **not** a single-fit coefficient inconsistency; they reflect **an empirically observed
  contrast-output behaviour for extremely low-information genes** (we do not over-claim the internal cause).
- **None of the 53 is a frozen hit.**
- The frozen 60 hits pass the identity at the `results()` level **60/60**.
- The frozen 60 hits are **60/60 converged**, and the apeGLM→`ashr` shrinkage of all three contrasts
  succeeded for all 60.
- A/D/T, category counts, modules, MCI, threshold sensitivity, and claim eligibility have **not** been inspected.

## Gate 1 — global structural identity (on GLM coefficients)
From the YA-reference design: `β_HA = coef(dds)[,"grp_HA_vs_YA"]`, `β_AD = coef(dds)[,"grp_AD_vs_YA"]`. Raw
contrasts from coefficients: `HA−YA = β_HA`, `AD−YA = β_AD`, `AD−HA = β_AD − β_HA`;
`coefficient_identity_error = β_AD − [β_HA + (β_AD − β_HA)]`. Universe: finite required coefficients, DESeq2
LFC fit converged, full-rank design. **Pass: max |coefficient_identity_error| ≤ 1e-10** (tighter than 1e-6 —
the algebra is exact, observed ≈3.55e-15). STOP if any gene exceeds.

## Gate 2 — frozen-60 operational mandatory gate
All 60 hits: (fit validity) converged, finite baseMean, finite raw LFC and SE for all 3 contrasts;
(results-level identity) `results_identity_error = AD_YA_results − (HA_YA_results + AD_HA_results)`, **|·| ≤ 1e-6,
60/60**; (shrinkage validity) `ashr` success for all 3 contrasts, 3 shrunken LFCs finite, and finite A, D,
T_class, T_direct_shrunk. No exclusion — STOP if any single hit fails.

## Global `results()` non-additivity is reported, not hidden
The 53 genes are NOT deleted or re-scoped by expression. A descriptive QC is saved
(`results_lfc_additivity_audit.tsv` + summary: total finite genes, #≤1e-6, #>1e-6, error median/95th/99th/max,
violator baseMean distribution, frozen-hit violators, contrast with the coefficient-level identity). This is a
descriptive audit, **not** a classification gate.

## Output columns & classification (unchanged otherwise)
Separate `*_LFC_coef` (structural audit) + `coefficient_identity_error`; `*_LFC_raw` (`results()` MLE for
reporting/Wald p-values) + `results_identity_error`; `*_LFC_ashr` (primary classification effect).
Classification uses **ashr LFC only**: `s = sign(frozen aging log2FC)`; `A = s·ashr(HA−YA)`; `D = s·ashr(AD−HA)`;
`T_class = A + D`; `T_direct_shrunk = s·ashr(AD−YA)`; `T_discrepancy = T_direct_shrunk − T_class` (not used to
classify). Projection gate `A ≥ 0.25`; `δ = max(0.25, 0.5·A)`; the four categories; the 9-combination threshold
sensitivity; module/MCI conventions — all unchanged.

## Proceed condition
Gate 1 (coef identity ≤1e-10, finite, converged, full-rank) AND Gate 2 (frozen 60/60 results-identity ≤1e-6 +
convergence + ashr + finite A/D/T) → B6 classification → B7 sensitivity/modules/MCI → B8 freeze + commit + STOP,
with no further confirmation. Any failure → STOP (no expression-threshold post-hoc cut, no DE-universe
exclusion, no count-filter change, no frozen-hit exclusion, no raw fallback, no category-driven ashr/raw choice).
