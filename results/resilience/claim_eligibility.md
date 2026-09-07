# CONSEQUENCE arm — claim eligibility (Phase A, post implementation-audit)

**Verdict: NOT SUPPORTED** (formally established after the complete implementation audit).

Of-record design `prereg/secretome-resilience-trajectory-v1` (5b29c74); receptor set frozen pre-outcome
(`freeze/consequence-receptor-set-v1`, be864a3); initial result (98dac20) unchanged; implementation audit
adds the remaining robustness/supporting analyses + a matrix-provenance check (no primary re-selection).

## Implementation audit
- **Matrix provenance = PASS.** `GSE325391_resilience.h5ad` `X` is raw int64 counts (min 1, max 2807, 0%
  non-integer, 0% negative; `X == layers['counts']`; raw not present). CPM/log2(CPM+1) in A1/A3/A6 valid →
  receptor freeze and primary stand; no amendment required.
- **Deviations logged:** D1 nlme(`~1|Run/donor`) for the pre-registered REML mixed model (lme4 absent in env);
  D2 NicheNet enrichment statistic not fully pre-fixed → both competitive Wilcoxon and fgsea reported,
  neither cherry-picked.

## Result (primary unchanged)
- **Primary (receiver-competence), stage-averaged SAD − RES = −0.070** z, 95% CI [−0.335, 0.195], df 6,
  **p = 0.544 (two-sided)**. Run variance ≈0 → donor-only fallback identical (p = 0.534). Interaction p=0.866.
- **Robustness — ALL pre-specified variants non-significant** (same-direction null): 3-stage −0.059 (p=0.62);
  raw-pooled −0.073 (0.51); ligand-balanced −0.244 (0.19); drop-ECM/GF +0.050 (0.68); leave-one-receptor-out
  range [−0.096,−0.051] all p>0.05; **Run-fixed −0.112 (0.52, estimable from mixed Runs)**; **no-Run −0.070
  (0.53)**; **leave-one-Run-out p∈[0.31,0.63]**; **per-family leave-one-out** apolipoprotein/collagen/
  cytokine/matricellular all p≥0.59, other_growth −0.279 (0.095); **subtype-specific** DiffN/NTF3/CHRM3 all
  p≥0.16. No variant trends toward the predicted SAD>RES; the weak non-significant trends are opposite.
- **NicheNet downstream supporting (network-prior-dependent), top-20/50/100, two methods:** competitive
  Wilcoxon union p = 0.48/0.44/0.35, **0/22 ligands BH<0.05 at every threshold**; fgsea NES −1.23/−1.17/−0.86,
  p = 0.12/0.14/0.87. No RES/SAD enrichment at any threshold or method.

## Ladder placement (per pre-registered rule)
All confirmatory + completion analyses are same-direction null → **NOT SUPPORTED** (maintained, not reclassified;
no significant variant exists to promote). Both legs null ⇒ **GSE325391 does not support the secretome-receiver
hypothesis** at the receptor or network-downstream level.

## What may / may not be said (unchanged ceiling)
- **May say:** the GSE325391 granule-lineage receptor-availability score for the aging niche secretome did not
  differ between resilient (RES) and susceptible (SAD) high-pathology donors (SAD−RES = −0.07 z, p = 0.54),
  robust to every pre-specified variant; no downstream RES/SAD enrichment.
- **Must NOT say:** glial secretome amount differs in RES; receptor expression measures ligand exposure; the
  secretome caused cognition; hippocampal glial sender state was observed.
- **Power caveat (honest null):** n = 13 (RES6/SAD7); age/sex/pathology not covariable; pathology-match
  study-defined; group×Run confound. A small true effect is not excluded, but there is no evidence for the
  predicted resilience association and the estimate is near zero and direction-discordant.

Implication: the "aging niche secretome ↔ cognitive resilience" link is **not supported on the granule-lineage
receiver side in this cohort**; the gated glial-sender SOURCE arm and the aging-to-AD TRAJECTORY arm remain
live. Reported as a negative result, not omitted.
