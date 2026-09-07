# Amendment v1.3 — single-fit DESeq2 + ashr direct-contrast shrinkage (post-diagnostic, pre-classification)

**Post-diagnostic, pre-classification methodological amendment.** Changes only the TRAJECTORY (§B) DE/shrinkage
*implementation*; no hypothesis, gate, threshold, gene set, or claim changes. v1 (`5b29c74`), v1.1, v1.2 tags
retained permanently. New tag: `prereg/secretome-resilience-trajectory-v1.3`.

## Why (diagnostic record)
- **v1.1 PyDESeq2 dual-full-fit** (YA-ref + HA-ref) failed the global AD−HA reparameterisation-equivalence gate:
  two independent `deseq2()` runs estimate dispersions independently, so 45/115,519 LFC-converged genes
  differed >1e-4 (max ≈0.37); all 60 frozen hits were exact (≤3.25e-5).
- **v1.2** scoped the gate to converged genes; the 45 both-converged failures persisted (real, not pure
  non-convergence — dispersions differed ~8% on extreme low-count genes).
- A **shared-nuisance coefficient-only transplant** was empirically tested and **rejected**: naive deepcopy +
  design-swap + `fit_LFC` gave max|Δ|≈3.36; a fresh object with transplanted (verified-identical) size-factors,
  dispersions and `non_zero_idx` + `fit_LFC` still gave max|Δ|≈1.40 even for non-Cook's-replaced converged
  genes. PyDESeq2 0.5.4 `fit_LFC()` depends on pipeline warm-start state that cannot be safely transplanted.
- The **v1.3 shared-nuisance proposal was never committed**; no raw-LFC fallback was used; and **A/D/T, category
  counts, modules, MCI, threshold sensitivity, and claim eligibility have not been inspected**.

## Resolution
Switch Phase B to **one R DESeq2 fit per cell type** and **`ashr`** shrinkage of *direct contrasts* from that
single fit. Unlike `apeglm` (named coefficients only), `ashr` shrinks a `results` object built from any
contrast — so all three contrasts come from one fit with one dispersion set, eliminating the second fit,
reference re-coding, dispersion re-estimation, and the equivalence problem entirely. Option A (two-full-fit
bound only by the frozen-60 gate) is **rejected** (the independent-dispersion divergence and reference-dependent
shrinkage remain undesirable for the primary classification).

## Unchanged
raw counts; donor inclusion; 5 cell types; `~ arm + grp`; YA/HA/AD primary population; MCI context-only; frozen
60 hits; projection gate `A ≥ 0.25`; `δ = max(0.25, 0.5×A)`; the four categories; `T_class = A + D`; separation
of p-values (Wald) from effect-size classification (shrunken); the 9-combination threshold sensitivity; the
claim ceiling.

## Changed
- DE implementation: **PyDESeq2 → R DESeq2** (env `trajR`: DESeq2 1.50.2, ashr 2.2.63, R 4.5.3; pinned,
  `envs/_trajR_pins.txt`; `sgz_r` untouched).
- Shrinker: **PyDESeq2 apeGLM → DESeq2 `ashr`**.
- AD−HA: **HA-reference second fit → single-fit direct contrast** `results(dds, contrast=c("grp","AD","HA"))`.
- Gate: **dual-fit equivalence → single-fit raw contrast-additivity identity**
  `raw(AD−YA) − [raw(HA−YA) + raw(AD−HA)] ≤ 1e-6` (exact within one fit), plus the frozen-60 mandatory gate.

## Execution
B1 fresh donor×celltype pseudobulk → B2 single `DESeq(dds)` per celltype → B3 three direct contrasts → B4
`ashr` shrinkage of each (no raw fallback; STOP if any frozen-hit shrinkage fails) → B5 single-fit identity gate
(Gate1 global ≤1e-6; Gate2 frozen-60 mandatory) → B6 projection gate + classification (`T_class`) → B7 threshold
sensitivity + modules + MCI context + R-vs-PyDESeq2 raw-LFC supporting audit (frozen-60, descriptive) → B8 freeze
+ commit + STOP. Fail-closed: STOP (no exclusion, no tolerance change, no raw fallback) on any gate failure.
