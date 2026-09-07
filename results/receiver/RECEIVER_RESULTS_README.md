# results/receiver/ — receiver-map compute outputs (prereg/receiver-map-v1)

Computed 2026-06-23 under the frozen pre-registration `docs/RECEIVER_MAP_PLAN.md`
(tag `prereg/receiver-map-v1`). env: bio (matched-null), sgz_r (step 00b rds edges). seed=42, 2000×.
**No manuscript edits** — this is the STOP-gate-2 compute only.

## Headline
- **PRIMARY (matched-null-on-lfc, `03`): NO-GO.** 0/8 confirmatory compartments pass
  (mean signed lfc>0 AND BH-q<0.1). No coordinated age-UP receptor remodeling. Two-sided BH min-q=0.12
  (no coordinated shift in **either** direction survives FDR). Pre-authorized conclusion holds:
  *"The frozen age-UP ligand set is not accompanied by coordinated age-UP receptor regulation in
  candidate receiver compartments."* Claim tier: **NOT SUPPORTED** for coordinated up-regulation;
  the unified map = SUPPORTING/descriptive negative.
- **NSC repertoire (depth null, `02`): DEPTH-EXPLAINED.** v1 full-76 statistic: obs 22 of 76 testable
  (non-NA padj, i.e. testable after independent filtering) vs depth-matched null 21.99 (effect +0.01,
  p_low=0.69). The v1.1 decomposition (`02b`, below) resolves this into coverage (40 of 76 in the DE table)
  and in-table testability (22 of those 40), with the in-table matched-null at 21.99 (p_low=0.68). NSC's
  small receptor repertoire is **not separable from low-depth profiling of a rare population** — bounded
  caveat, NOT biology (as pre-specified §6).
- **04/05 not run** — no passer triggered LOO or donor-composite confirmation (correct per §8/§10).

## Files
| file | content |
|---|---|
| `_universe76_from_manifest.csv`, `ligands_22.csv` | frozen inputs (manifest-authoritative) |
| `frozen_receptor_universe_76.csv` | 76 receptors + n_ligands_connecting (rds-derived, cross-checked == manifest) |
| `ligand_receptor_edges.csv` | 109 ligand→receptor edges (22 ligands) |
| `detectability_matrix.csv` | 76 × 11 compartments: in_table / detectable / baseMean / lfc / padj |
| `detectability_summary.csv` | per-compartment n_detectable/76, global_detectable_frac |
| `depth_matched_null.csv` | depth-matched detectability null (`02`, v1 — original record) |
| `depth_testability_v1_1.csv` | **v1.1** decomposed depth/testability (`02b`: A coverage · B in-table detect · C in-table matched-null) |
| `matched_null_per_compartment.csv` | PRIMARY matched-null-on-lfc (`03`) + BH-q + calls |
| `_step00_provenance.json` | SHA pin verification |

## ⚠ Interpretation caveat — step 02 for high-global-detectability compartments
The depth-matched-null statistic is "# of the 76 detectable (in dds AND non-NA padj) vs baseMean-matched
genes drawn from the in-dds pool." It is **valid for low-global-detectability compartments** — notably
**NSC** (global detectable frac 0.34), the compartment §6 was designed to vet, where it correctly returns
depth-explained.

For **DG_GC, InN, Neuroblast, Immature** the null reports extreme "depletion" (p_low≈0.0005), but this is a
**known matching artifact, NOT a receiver-intrinsic finding:** these compartments have global detectable
frac ≈ 1.0 (0.997–0.9992), so the in-dds background pool is ~fully detectable at every baseMean decile and
the null predicts ~76. The observed shortfall is entirely the receptors that are **absent from the dds**
(below the pipeline's gene-inclusion / expression threshold; see `detectability_summary.csv` n_in_table:
Neuroblast 43, Immature 48, DG_GC 63, InN 64 of 76). The null cannot represent sub-inclusion-threshold
genes, so it cannot fairly test those receptors. **No claim is made from these four p-values.** They
reflect an expression-level descriptive (some receptors fall below inclusion) already visible in step 01.
A cleaner decomposition (inclusion vs independent-filtering, or restricting the null to in-dds receptors)
is available as an optional `prereg/receiver-map-v1.1` refinement if the neurogenic depletion question is
pursued — but the headline conclusions (NO-GO; NSC depth-explained) do not depend on it.

### v1.1 RESOLUTION (`02b_depth_testability_v1_1.py` → `depth_testability_v1_1.csv`)
The artifact above is resolved by decomposing testability into A) coverage `n_in_table/76`, B) in-table
detectability `n_detectable/n_in_table`, and C) a `baseMean`-matched null run **only within the in-table
universe** (absent receptors excluded), with saturated compartments flagged rather than p-valued.
- **NSC: interpretable, depth-explained** — 40 in-table receptors (of 76); 22 testable after independent
  filtering, matched-null 21.99 (≈ 22.0), effect +0.01, `p_low=0.68`. (The 22 of 76 testable overall is the
  coverage-inclusive descriptive value; the v1.1 null denominator is the in-table 40, not 76.) Same
  conclusion as v1, now isolated to the independent-filtering layer. Wording retained:
  *"NSC had limited receptor testability, but this was explained by low-depth profiling of a rare
  exploratory compartment rather than evidence for a receiver-intrinsic loss of receptor repertoire."*
- **DG_GC, InN, Neuroblast, Immature: `not_interpretable_saturated`** (in-table background detectable
  frac ≥ 0.98) — the spurious v1 `p_low≈0.0005` is replaced by an explicit "no resolving power" flag.
  Reported descriptively only (coverage `n_in_table/76`); **no biological depletion claim.**
- Other confirmatory compartments: interpretable, all depth/expression-explained (`p_low` 0.43–1.0).
- **STEP 03 (primary) is unchanged** (not re-run; no file-integrity dependency). NO-GO stands.
