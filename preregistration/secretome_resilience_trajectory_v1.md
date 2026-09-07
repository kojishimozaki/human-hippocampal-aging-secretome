# PRE-REGISTRATION — Secretome resilience & aging-to-AD trajectory (v1)

**Status: LOCKED, pre-outcome.** Written and committed BEFORE any outcome analysis (no RES/SAD receiver
scores, no RES/SAD receptor effects, no GSE268609 HA−YA/AD−HA/AD−YA effects, no trajectory category
counts, no module group-differences, no SOURCE-cohort 60-hit directionality were inspected at P0).
Design locked via `/grill-me` (PI-confirmed, branches Q1–Q4 + execution procedure).

**Base commit:** `58e2a5b` (= tip of `analysis/methylome-hardening`).
**Base rationale:** stable shared core predating the `analysis/neuron-niche` and `analysis/tf-rna-concordance`
branches, carrying all frozen inputs required for the resilience and trajectory arms (frozen 60-hit
signature, GSE325391/GSE268609 processed data, NicheNet priors, GSE325391 RES/SAD DE).
**Branch:** `analysis/secretome-resilience-trajectory`.
**Registration commit / tag:** `5b29c74b91ba68886a239b066da10919899dfbf7` (registration of record; 2026-06-20) — tag `prereg/secretome-resilience-trajectory-v1`.
**Pre-stash of prior WIP:** `git stash a902e75f…` ("WIP before secretome-resilience branch: TF-RNA prose +
neuron work"); not applied on this branch.

---

## §0 Scientific question and frozen claims

The frozen niche aging-secretome signature is a **glial OUTPUT** (Astro/Micro/Endo/OPC/Oligo). This
project asks whether that secretome connects to **cognitive resilience/vulnerability** and to **AD state**,
across three arms with a fixed claim ceiling.

- **CONSEQUENCE (confirmatory, primary arm).** GSE325391 measures only the **SGZ granule-cell lineage**
  (no glia), so the glial secretome itself cannot be measured there. We therefore test the **receiver**
  side: *the granule-cell-lineage receptor state capable of receiving the aging niche secretome is
  associated with cognitive resilience/vulnerability.* The primary is a **frozen receptor-availability /
  receiver-competence score**; NicheNet downstream targets are a **supporting** network-prior-dependent
  analysis. We do NOT claim "high receptor expression = strong secretome stimulation"; actual
  transcriptional response is supported only if the downstream-target analysis concords.
- **TRAJECTORY (GSE268609 aging-to-AD).** Which components of the frozen aging-aligned niche signature are
  **amplified / maintained / AD-attenuated / AD-divergent** across a **cross-sectional** YA→HA→AD state
  axis. No longitudinal language (progression/conversion/evolution).
- **SOURCE (cohort-gated exploration).** Whether a public cohort with **niche glia + cognition + pathology**
  exists to test the glial secretome's source in resilience directly; analysed only if a cohort meets the
  pre-fixed eligibility criteria.

**Central claim ceiling (locked).** Until SOURCE is established, the paper's claim is NOT "the secretome is
itself a cognitive-vulnerability marker," but: **"the aging niche secretome and its neuronal-side receiver
response connect to AD state and cognitive resilience."**

---

## §1 Shared frozen inputs and provenance

All checksums in `manifests/frozen_input_checksums_v1.tsv` (SHA256). Key inputs:

| role | path | sha256 (head) |
|---|---|---|
| frozen 60-hit signature | `results/validation/frozen_primary_signature.csv` | `956e8743…` |
| NicheNet lr_network (2021-12-21) | `refs/nichenet/lr_network_human_21122021.rds` | `47c971d2…` |
| NicheNet ligand_target_matrix | `refs/nichenet/ligand_target_matrix_nsga2r_final.rds` | `699fce17…` |
| GSE325391 metadata | `processed/per_dataset/GSE325391_adultgc_metadata.csv` | `f215d10e…` |
| GSE325391 resilience h5ad | `processed/per_dataset/GSE325391_resilience.h5ad` | `ea814d17…` |
| GSE268609 metadata | `processed/per_dataset/GSE268609_metadata.csv` | `e13a5b3f…` |
| GSE268609 anchor h5ad | `processed/per_dataset/GSE268609_anchor.h5ad` | `1e467d76…` |
| GSE325391 RES-vs-SAD DE (supporting GSEA) | `results/de/de_GSE325391_RESvSAD.csv` | `2990367a…` |
| GENCODE v44 symbol map | `refs/gencode_v44_genes.tsv.gz` | `9f735983…` |

**Frozen signature (immutable):** 60 hits = Astro 26 / Micro 12 / Endo 11 / OPC 6 / Oligo 5; 36 UP / 24 DOWN.
Gene symbols, `direction`, `celltype`, `log2FoldChange` are NOT modified by any arm.

**GSE325391 known structure (pre-outcome):** granule-cell lineage only (`broad_cell_type = ExN_GC`); donors
CTRL 6 / RES 6 / MAD 5 / SAD 7; subtypes DiffN, DiffN_OTOF, MatN_OTOF, MatN_SGCZ_CHRM3, MatN_SGCZ_NTF3;
batch = Run (1–6). **No pathology / age / sex columns** (sex set "U"). RES/SAD pathology-matching depends on
the original study's group definitions and is not verifiable in-data; age recoverable, if at all, only from
external SOFT (not done at P0). The 3 main lineage stages (DiffN, MatN_SGCZ_NTF3, MatN_SGCZ_CHRM3) each have
≥50 cells in all 24 donors (feasibility-confirmed pre-outcome).

**GSE268609 known structure (pre-outcome):** YA 8 / HA 9 / AD 10 (+ MCI 6, SA 6); two tissue-prep arms
(whole-hippo orig.ident 1–14 / DG-microdissect 15–39) imbalanced across groups → model `~arm + grp`; sex "U".

**NicheNet feasibility (pre-outcome grounding, not an outcome):** of the 34 unique UP-secreted genes, 22 are
recognised NicheNet ligands (ANXA1, APOC1, APOD, APOE, AZGP1, CCN2, COL12A1, COL21A1, COL4A2, COL8A1, CXCL2,
CYTL1, DKK2, EFEMP1, FAP, IL15, KITLG, LGALS9, SEMA3D, SERPINE1, SPON1, TNC) → 76 distinct receptors before
expression filtering. The **final** receptor set is constructed and frozen at Phase A1 (post-tag) by the §A
rule; only the construction RULE is pre-registered here.

---

## §2 Analysis hierarchy and stopping rules

**Priority order (locked):** A) CONSEQUENCE confirmatory → STOP; B) GSE268609 trajectory → STOP; C) SOURCE =
eligibility exploration only, may run in parallel after the tag but never reshapes A or B.

**Phase order to the tag:** P0 write pre-reg (no outcomes) → P1 freeze artifacts (this doc + manifests +
checksums + planned scripts/outputs + env + full primary/secondary/exploratory hierarchy) → P2 commit +
annotated tag, verify clean tree. **Step-0 feasibility/data-integrity checks run AFTER the tag** and are
logged to `results/step0/`. Each arm ends at a STOP that commits outputs + a decision log regardless of
positive/negative/unstable outcome.

**Pre-tag freeze guarantees:** after seeing any arm's outcomes we do NOT change the frozen receptor set,
maturation grouping, Run handling, primary estimand, NicheNet target selection, the trajectory δ, or the
projection gate. Threshold sensitivities are pre-listed; the most-significant configuration is never chosen
post hoc.

---

## §A CONSEQUENCE confirmatory arm (GSE325391 receiver-competence)

**Claim:** *the granule-cell-lineage receptor state capable of receiving the aging niche secretome is
associated with cognitive resilience/vulnerability.* Primary = **frozen receptor-availability /
receiver-competence score**; NicheNet downstream = **supporting** (network-prior-dependent).

### A1. Frozen receptor set construction (executed post-tag; rule locked here)
- **Ligand origin:** frozen 60-hit, `direction == UP`, secreted categories = the 36 UP genes. aging-DOWN
  ligands are NOT mixed in (separate exploratory if ever).
- **ligand→receptor:** intersect the 36 UP genes with the NicheNet lr_network `from` column; take direct
  receptors (`to`). Record the gene-symbol convention. Select edges WITHOUT seeing RES/SAD. A receptor
  connected to several ligands is de-duplicated to one gene in the primary set. lr_network has coarse
  `source` (omnipath / nichenet_verschueren) and no fine confidence tier → use the whole 2021-12-21 network
  as-is, recording version + checksum. **Forbidden:** keeping only receptors that differ in RES/SAD.
- **Expression filter (diagnosis-blind):** computed on **CTRL donors only** (NOT RES/SAD). Per donor ×
  neuronal subtype pseudobulk; `CPM ≥ 1` = detectable; a receptor is adopted if detectable in **≥3 of 6 CTRL
  donors in ≥1 SGCZ maturation subtype**. A label-blind all-group pooled filter is a sensitivity only and
  never changes the primary set.
- **Freeze artifacts (before any RES/SAD test):** 36 UP origin genes; NicheNet-recognised ligands; L–R edge
  list; receptor list pre/post CTRL filter; excluded receptors + reasons; final frozen receptor set;
  resource version/checksum; expression threshold; score formula →
  `results/resilience/frozen_receiver_receptor_set.csv`, `results/resilience/frozen_receiver_set_manifest.json`.

### A2. Primary score
Unit = **donor × maturation subtype pseudobulk** (cells are NOT independent units). Equal-weight composite:
per-receptor pseudobulk logCPM → per-receptor z-standardisation across the primary (RES+SAD) donor×subtype
samples → mean over frozen receptor genes → one score per donor×subtype. NO weighting by network degree,
NicheNet weight, or RES/SAD effect size in the primary.

**Maturation (2-level primary):** immature = `DiffN`; mature = aggregate of `MatN_SGCZ_NTF3` and
`MatN_SGCZ_CHRM3`. The **mature score is the equal-weight mean of the donor-level NTF3 score and CHRM3 score**
(NOT a raw-count pool), so that group differences in mature-subtype composition are not mistaken for
receiver-expression differences. OTOF subtypes excluded from primary. Min-cell **≥20 per original subtype**
(DiffN, NTF3, CHRM3 each checked before pooling); donor×subtype below threshold excluded and logged
(`results/resilience/pseudobulk_cell_counts.csv`, `…_inclusion_log.csv`); no imputation; threshold never
lowered post hoc; group-wise exclusion counts reported.

### A3. Primary model, estimand and test
Primary fit on **RES + SAD only** (CTRL/MAD are a separate contextual model + figure, ordered CTRL→RES→MAD→SAD).
```
receiver_score ~ group * maturation_stage + (1 | Run) + (1 | donor)
```
REML; Kenward–Roger or Satterthwaite df/CI. **Primary estimand (single):** the **stage-averaged (equal-weight
across the 2 stages) marginal SAD − RES group contrast.** Biological pre-direction `SAD − RES > 0` (RES < SAD),
but the **test is two-sided**. Interpretation: RES<SAD ⇒ consistent with higher receiver-competence in SAD;
RES≈SAD ⇒ no receptor-level association; RES>SAD ⇒ inconsistent with a simple vulnerability model (treat as
compensatory expression / feedback). **Key secondary:** group × maturation_stage interaction; stage-specific
SAD−RES; within-group mature−immature. n RES 6 / SAD 7 → report effect size + CI; cell counts are NOT used as
df. Age/sex/pathology not adjustable; pathology-matching is study-defined; **Run is partially modelled, not
"batch-adjusted"** (Run nests donor; RES/SAD Run distribution is imbalanced; group–Run confound not fully
removable).

### A4. Robustness (pre-specified; do not pick the cleanest)
ligand-balanced score (mean within ligand, then across ligands); leave-one-receptor-out; leave-one-ligand-
family-out; ECM/growth-factor (many-receptor) drop; 3-level maturation; raw-count pooled mature; Run-fixed /
no-Run / leave-one-Run-out; subtype-stratified `score ~ group + Run`; OTOF-added; diagnosis-blind slingshot
pseudotime (exploratory). **Singular-fit policy:** report singular fit + Run variance; fall back to the
pre-defined donor-random-only model; compare with Run-fixed; never select the most significant.

### A5. NicheNet downstream supporting analysis (network-prior-dependent)
From the frozen UP ligands, NicheNet `ligand_target_matrix`; target rule fixed pre-outcome: **top-50 regulatory-
potential targets per ligand**, restricted to CTRL-granule-detectable targets; record per-ligand and union
sets; threshold sensitivity top-20/50/100 (never pick the most positive). NicheNet potential ≠ sign of
expression change, so the mean target expression is NOT called a "response score." Supporting test =
**competitive gene-set enrichment** of the frozen target set on the **donor-level RES-vs-SAD ranking**
(statistic = SAD − RES model statistic; positive ⇒ targets higher in SAD); per-ligand enrichment multiplicity-
corrected. Joint reading: receptor score AND downstream enrichment both SAD-ward ⇒ "receiver state + network-
predicted downstream program for aging-UP niche ligands relatively stronger in SAD"; receptor-only ⇒ receiver-
availability association; target-only ⇒ network-predicted-program association; both null ⇒ GSE325391 does not
support the secretome-response hypothesis. Exploratory only (never primary): per-ligand NicheNet ligand-
activity (previously non-discriminative; rerunning as primary = post-hoc rescue).

---

## §B GSE268609 aging-to-AD trajectory arm

Cross-sectional **state axis** across young aging, healthy aging, and AD (no progression/conversion language).

### B1. Model & unit
Donor × cell-type pseudobulk; **per-celltype separate fits** (Astro/Micro/Endo/Oligo/OPC); each frozen hit uses
only its assigned celltype's fit. One fit per celltype: `~ arm + grp`, `grp ∈ {YA, HA, AD}` (YA reference;
MCI/SA excluded from the primary fit). Extract contrasts **HA−YA, AD−HA, AD−YA** from the one fit. Pre-flight
(Step-0 B0): arm×grp crosstab, full-rank design, no group fully confined to one arm, per-celltype donor counts,
extreme library-size, donor exclusions+reasons. `arm` adjustment does not fully remove a strong arm–diagnosis
confound — stated. **Effect estimate:** shrunken log2FC for classification (raw + shrunken saved; shrinkage
method + contrast spec + any non-shrinkable contrast logged); padj reported as evidence strength, not used to
classify.

### B2. Directional definitions
`s = sign(frozen aging log2FC)` (UP→+1, DOWN→−1). Align GSE268609 contrasts to the frozen aging direction:
`A = s·LFC(HA−YA)` (aging recapitulation), `D = s·LFC(AD−HA)` (further AD movement), `T = s·LFC(AD−YA) = A + D`.

### B3. Aging-projection gate (Step 1; padj NOT used)
Classify only hits where GSE268609 recapitulates aging: **eligible if `A ≥ 0.25` log2FC.** Else:
`A < 0` → **aging direction not recapitulated**; `0 ≤ A < 0.25` → **weak / indeterminate aging projection**.
These are NOT forced into the four categories. (No significance filter — avoids power-based selection.)

### B4. Per-hit threshold & four mutually-exclusive categories (effect-size based)
`δᵢ = max(0.25, 0.5 × A)`. For gate-passing hits:
- **Amplified:** `D ≥ δᵢ` — AD-associated amplification of the aging-aligned change.
- **Maintained:** `−δᵢ < D < δᵢ` — healthy-aging state retained in AD (not "no change": the HA change persists).
- **AD-attenuated:** `D ≤ −δᵢ` AND `T > 0` — aging change attenuates in AD but AD still on the aging side of YA
  (redefinition of the ambiguous original "age-associated only"; "age-associated only" is descriptive prose
  only, never a class name).
- **AD-divergent:** `T ≤ 0` (usually also `D ≤ −δᵢ`) — aging-direction change lost; subflags: *reverted-toward-YA*
  `|T| < 0.25`, *direction-reversed* `T ≤ −0.25`.

### B5. Evidence, sensitivity, modules, MCI
Per hit save all three shrunken contrasts + SE/CI + nominal p + padj + category + margin distance (`D − δᵢ`).
Robustness (not used to re-tune rules): leave-one-donor-out category stability; bootstrap category proportions;
raw-vs-shrunken agreement; arm-adjusted vs not. **Threshold sensitivity:** relative {0.33, 0.50, 0.67} ×
absolute {0.20, 0.25, 0.30} (9 combos) reported for category counts + main-conclusion stability; **main result
fixed at relative 0.50 + absolute 0.25.** **Modules per celltype, UP and DOWN separately** (no cross-direction
averaging, no mixing all 5 celltypes); optional secondary signed direction-aligned score; `module_score ~ arm +
grp` per celltype with EMMs+CI, displayed YA→HA→MCI→AD but **primary inference YA/HA/AD only**. **MCI**: used only
in a secondary model/visualisation after the primary is frozen; "MCI provided intermediate clinical context but
was not used to define trajectory classes"; no transition/progression claim.

### B6. Outputs
`results/trajectory/frozen60_GSE268609_hit_trajectory.tsv` (columns: gene, celltype, frozen_direction,
frozen_log2FC, HA_vs_YA_LFC_raw/shrunken, AD_vs_HA_LFC_raw/shrunken, AD_vs_YA_LFC_raw/shrunken, A, D, T,
delta_relative, delta_absolute, delta_final, aging_projection_status, trajectory_category, trajectory_subflag,
pvalue/padj for the three contrasts, classification_stability), `…_module_scores.tsv`,
`…_threshold_sensitivity.tsv`, `…_classification_manifest.json`.

---

## §C SOURCE cohort-gated exploration

Gated Phase-0 search; may run in parallel after the tag; **never reshapes A or B.** Search for a public cohort
with niche glia + cognition + pathology that distinguishes a high-pathology-cognition-preserved group from
dementia. Candidates: ROSMAP, SEA-AD, Mathys, Green. **Eligibility (locked, `manifests/source_eligibility_v1.json`):**
(1) astrocyte/microglia/endothelial/OPC/oligodendrocyte analysable per donor; (2) pathology burden present (not
only cognition/diagnosis); (3) a high-pathology-cognition-preserved group definable; (4) donor count + covariates
adequate; (5) the frozen 60-hit signature projectable unchanged. Proceed to real analysis ONLY if a cohort meets
all criteria; freeze a separate analysis plan first (C3). If cortical, frame as **cross-brain-region glial
vulnerability association**, not an SGZ-niche reproduction; avoid hippocampal-niche extrapolation. Finding no
eligible cohort is itself a recorded result; no approximate adoption. Steps logged to
`results/step0/source_cohort_eligibility_audit.tsv`, `results/step0/step0_decisions.md`.

---

## §D Multiplicity and interpretation hierarchy

CONSEQUENCE: one pre-fixed global primary test (stage-averaged SAD−RES). Interaction + stage-specific +
within-group are key-secondary; NicheNet per-ligand enrichment is multiplicity-corrected within its family;
arms/families are corrected separately and never pooled. TRAJECTORY classification is descriptive point-estimate
based; padj/stability are confidence information, not class determinants. Claim ladder applied per result:
**supported / directionally consistent but uncertain / not supported / model-sensitive / not estimable /
feasibility failure** — evaluated against this ladder, not by a bare p-value.

---

## §E Negative, null, failed-feasibility and ambiguous outcomes

All reported, not dropped. CONSEQUENCE RES≈SAD or RES>SAD are reported with the locked interpretations.
TRAJECTORY hits failing the projection gate are reported as not-recapitulated / indeterminate (not silently
omitted). SOURCE "no eligible cohort" is a result. Feasibility failures (receptor set largely unmappable, rank-
deficient `arm+grp`, non-estimable primary contrast, missing donor×stage data, ineligible SOURCE) trigger §F
amendment, not a silent design change.

---

## §F Deviations and amendment policy

The v1 tag is never edited or deleted. **Minor implementation fixes** (symbol case → unique alias, path fixes,
syntax to implement the stated model) → logged in a deviation log; scientific design unchanged. **Design-affecting
fixes** (receptor set largely unmappable, rank-deficient design, non-estimable primary, missing required data,
SOURCE ineligibility) → written as `preregistration/amendments/amendment_v1.1.md` (reason, impact, alternative)
WITHOUT looking at outcomes, then a new commit + annotated tag `prereg/secretome-resilience-trajectory-v1.1`;
v1 retained permanently. Stash identity is tracked by commit hash + message (not `stash@{0}`).

---

## §G Reproducibility, software and output manifests

Env: Python 3.11.15, pandas 2.3.3, numpy 2.4.6, scipy 1.17.1, scanpy 1.11.5 (`bio`); R (`sgz_r`) for NicheNet
RDS, pyDESeq2/lme4/lmerTest for models. `seed = 42` (stochastic steps: bootstrap, permutation, any sampling
only). `PROJ` env-var portability; mpl/numba caches → /tmp. Machine-readable design: `manifests/consequence_design_v1.json`,
`trajectory_design_v1.json`, `source_eligibility_v1.json`, `frozen_input_checksums_v1.tsv`. Planned scripts:
`scripts/10_resilience/0X_*.py|R` (receiver set, pseudobulk, mixed model, NicheNet GSEA),
`scripts/11_trajectory/0X_*.py` (3-group fit, projection-gate classification, modules). Each STOP saves: input
checksums; inclusion/exclusion log; model matrix; exact formula; frozen set manifest; primary result; all
pre-specified sensitivities; warnings/singular-fit logs; session info; deviation log; claim-eligibility decision.

**Anti-rename guarantee:** the trajectory arm's `A ≥ 0.25` projection gate (applied before classification) makes
it a genuine freeze→test projection, not a renamed 3-group DE; the CONSEQUENCE CTRL-only expression filter makes
the receiver set diagnosis-blind, not a RES/SAD-selected set.
