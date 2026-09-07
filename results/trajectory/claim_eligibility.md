# TRAJECTORY arm — claim eligibility (Phase B, B8)

> 📝 HISTORY (manuscript-integration branch): this file's body was **corrected** (codex re-audit, 2026-06-23) to within-cohort / effect-size language — the aging→AD projection is a **within-cohort** re-fit of the PRIMARY cohort GSE268609 (NOT independent replication), and the "maintained" calls are an **effect-size classification** (not a statistical maintenance / equivalence test; only 5/56 eligible hits reach AD-vs-HA padj < 0.1). The frozen-tag copy (`result/trajectory-phaseB-v1-final`) retains the original "independent / recapitulated / maintained" wording. See INTEGRATION_MANIFEST.md.

Of-record: `prereg/secretome-resilience-trajectory-v1` (5b29c74) + amendments v1.1–v1.4. Implementation:
R DESeq2 1.50.2 single fit per niche celltype (`~arm+grp`, YA/HA/AD) + `ashr` direct-contrast shrinkage
(env `trajR`). Gates PASS: coefficient additive identity 3.55e-15 (≤1e-10); frozen-60 mandatory (60/60
converged, finite, results-identity ≤1e-6, ashr valid).

## Result (cross-sectional aging-to-AD state axis; no longitudinal claim)
- **Aging-direction recovery (within-cohort self-consistency):** re-fitting the **primary** cohort GSE268609 (the
  cohort in which the signature was defined) with a 3-group model recovers the frozen aging direction in **56/60 hits
  eligible (A ≥ 0.25), 0 reversed**, 4 weak/indeterminate (0 ≤ A < 0.25). This is a self-consistency check under the
  model change, **NOT independent replication** (the genuinely independent transcriptomic replication is GSE278576, Fig A1).
- **Aging→AD effect-size category (ashr `T_class`, primary thresholds rel 0.5 / abs 0.25):** **category-maintained 54,
  amplified 1, AD-attenuated 1, AD-divergent 0** (of the 56 eligible). By this **effect-size rule** (NOT a statistical test
  of maintenance / equivalence) the aging-associated niche secretome is **not broadly amplified in the AD group**;
  individually only **5/56 eligible hits reach AD-vs-HA padj < 0.1**.
- **Threshold sensitivity (9 combos):** the maintained-dominant conclusion is stable (maintained 50–55;
  amplified 0–2; AD-attenuated 1–5; AD-divergent 0–2) across relative {0.33,0.50,0.67} × absolute {0.20,0.25,0.30}.
- **Shrinkage/implementation robustness:** raw-vs-ashr category concordance 0.83; `T_discrepancy` |median| 0.32;
  R-DESeq2-vs-PyDESeq2 raw LFC essentially identical (Pearson/Spearman 1.0000, max|Δ| ≤0.008, 100% direction).
- **Leave-one-donor-out stability (raw-LFC proxy):** median frac-stable **1.00**; **34/60 hits 100%-stable**.
  The single **amplified** hit is 100%-stable and the **AD-attenuated** hit 96%-stable; **maintained** hits are
  largely stable (median 1.0; 33/54 fully stable); the **4 weak/indeterminate (A near 0.25)** hits flip readily
  (expected boundary lability at the projection gate). The maintained-dominant headline is robust; the
  labile hits are exactly the gate-boundary ones, reported honestly.
- **ashr-primary vs raw-no-arm category (proxy; label corrected):** **50/60 (83%)** agreement (disagreements Astro 5 /
  Endo 2 / Oligo 2 / OPC 1). NOTE: this compares the ashr arm-adjusted **primary** category against a **RAW no-arm** refit, so
  it **conflates the `~arm` covariate AND ashr-vs-raw shrinkage** — it is NOT a clean arm-vs-no-arm comparison
  (`scripts/11_trajectory/04_stability.R`). A covariate/shrinkage-sensitivity note, not a change to the primary classes.
- **Caveat (raw-LFC proxy):** the LOO and ashr-primary-vs-raw-no-arm stability checks classify on the RAW `results()` LFC, not
  ashr — per-refit ashr over ~130 refits was computationally prohibitive. Raw-LFC classification is a validated
  proxy (raw-vs-ashr full-data category concordance 0.83). The **primary classification is unchanged and
  ashr-based**; these are robustness proxies only.

## Ladder placement
**SUPPORTED (descriptive, within-cohort):** in the **primary** cohort GSE268609 the frozen aging niche-secretome is, by an
**effect-size rule**, **not broadly amplified** across the aging→AD state axis (predominantly category-maintained; NOT a
statistical maintenance test; only 5/56 reach AD-vs-HA padj < 0.1). Reported as the component-level breakdown (amplify /
maintained / attenuate / diverge).

## What may / may not be said
- **May say:** within the primary cohort GSE268609, the frozen aging niche-secretome is, by effect size, **not broadly
  amplified in the AD group** (category-maintained 54/56 eligible; 5/56 reach AD-vs-HA padj < 0.1), with single hits
  amplified (1) or AD-attenuated (1) and none divergent, robust across thresholds — a within-cohort projection.
- **Must NOT say:** that GSE268609 is an independent / external replication cohort (it is the **primary** cohort); that
  "maintenance" is statistically proven (it is an **effect-size classification**); that individuals progressed from aging
  to AD; that HA is an AD-prodrome; that the secretome causes AD; that MCI is a temporal midpoint. No "progression /
  conversion / disease evolution" language.

Implication for the paper: complements the resilience arm — while the granule-lineage *receiver* side does not
show a resilience-associated effect (CONSEQUENCE = NOT SUPPORTED), the niche aging-secretome *output* signature, by a
within-cohort effect-size classification, is **not broadly amplified** in the AD group (predominantly category-maintained;
NOT a statistical maintenance/equivalence test; only 5/56 eligible hits reach AD-vs-HA padj < 0.1).
