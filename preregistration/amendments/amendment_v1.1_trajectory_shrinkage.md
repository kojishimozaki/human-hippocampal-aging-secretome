# Amendment v1.1 — trajectory contrast shrinkage & T=A+D classification (pre-outcome)

**Pre-outcome clarification.** This amends only the *implementation* of the TRAJECTORY arm (§B) shrinkage and
the additive classification identity. It does **not** change any hypothesis, threshold, gate, gene set, model,
or analysis hierarchy. The v1 registration (`prereg/secretome-resilience-trajectory-v1`, commit `5b29c74`) and
all v1 artifacts are retained permanently. New tag: `prereg/secretome-resilience-trajectory-v1.1`.

## 1. Premature B1 handling
- The earlier B1 attempt was **uncommitted and has been deleted**; nothing reached the record.
- Its `lfc_shrink()` silently fell back to **raw LFC** (wrong coefficient-name format). Any result containing
  that raw-LFC fallback is **not used** for the record, interpretation, or classification.
- **No** category counts, module trajectories, or individual-hit directions were used to alter the design.
- This amendment is a **technical fix of an implementation failure**, not an outcome-driven change of
  thresholds or hypotheses.

## 2. Fit structure (same counts, donors, filters, design for both fits)
Per celltype:
- **Primary YA-reference fit** `~ arm + grp`, grp reference = YA. **Of-record** for: model diagnostics,
  dispersions, raw LFC HA−YA and AD−YA, Wald p/padj for HA−YA, AD−YA, **and the AD−HA *direct contrast***.
- **Auxiliary HA-reference fit** `~ arm + grp`, grp reference = HA. Used **only** to obtain AD−HA as a named
  coefficient for apeGLM shrinkage. This is a re-parameterisation of the *same* design — not a different
  biological model — and does **not** replace the primary fit's inference or p-values.

## 3. Coefficient-name assertion (fail-closed)
The actual design-matrix LFC coefficient columns are saved per celltype
(`results/trajectory/design_coefficients_<celltype>.txt`) and the target coefficient is identified by
exact match in those columns (no hard-coded guessing). Expected (pydeseq2 patsy-style):
- YA-reference: `grp[T.HA]`, `grp[T.AD]`
- HA-reference: `grp[T.YA]`, `grp[T.AD]`

Mandatory assertions (any failure → **STOP, no classification**):
1. exactly one matching coefficient column is found;
2. `lfc_shrink()` completes without exception;
3. `DeseqStats.shrunk_LFCs == True`;
4. the shrunken LFC exists (non-NaN) for every frozen hit needing it;
5. **no silent fallback to raw LFC** anywhere.

## 4. AD−HA re-parameterisation equivalence audit
Compare, per gene, the **YA-reference primary direct contrast** AD−HA raw LFC against the **HA-reference
auxiliary named coefficient** (`grp[T.AD]`) raw LFC. Save
`results/trajectory/adha_reparameterization_audit.tsv` with columns: gene, celltype,
AD_HA_raw_direct_contrast, AD_HA_raw_HAref_coefficient, absolute_difference, relative_difference,
equivalence_pass. **Primary implementation check:** max absolute difference ≤ 1e-4. If exceeded → **STOP**,
do not run shrunken classification.

## 5. Final A/D/T definitions for classification
`s = sign(frozen aging log2FC)`.
- `A = s × shrunk_LFC(HA − YA)`
- `D = s × shrunk_LFC(AD − HA)`  (shrunken from the HA-reference auxiliary fit)
- **`T_class = A + D`** — used for the four-category classification.
- `T_direct_shrunk = s × shrunk_LFC(AD − YA)` — saved but **NOT** used to classify (per-contrast shrinkage is
  non-linear, so the additive identity is not guaranteed across independently shrunken contrasts; the
  pre-registered logical condition `T = A + D` is preserved by using `T_class`).
- `T_discrepancy = T_direct_shrunk − T_class` — shrinkage sensitivity / descriptive only; never changes a
  category.
- Raw identity audit: confirm `T_raw_direct ≈ A_raw + D_raw` and save the difference (`raw_identity_error`).

## 6. Roles of p-values vs shrunken LFC
- p-values / padj come from the **Wald contrasts** (YA-reference primary fit).
- trajectory classification is decided from **shrunken effect sizes**.
- significance alone never changes a category; p-values are **not** recomputed from shrunken LFCs.

## 7. Execution order (from scratch after this tag)
B0.1 amendment+freeze → B1 donor×celltype pseudobulk → B2 YA-ref primary fit → B3 HA-ref auxiliary fit →
B4 coefficient/shrinkage/equivalence assertions → B5 frozen-hit contrast table → B6 projection gate +
category assignment → B7 threshold + shrinkage sensitivity + modules + MCI context → B8 result freeze + STOP.
Fail-closed at B2–B4 (STOP, never continue on raw LFC) if: coefficient not uniquely resolved; shrinkage fails;
shrunken LFC missing for a frozen hit; AD−HA equivalence exceeds tolerance; design rank changed from Step-0;
donor inclusion mismatches the Step-0 audit; matrix not raw counts; or NaN enters a category computation.
