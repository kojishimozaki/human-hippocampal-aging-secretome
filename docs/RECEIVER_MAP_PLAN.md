# RECEIVER_MAP_PLAN.md — Pre-registration: unified cross-cell-type receiver receptor-regulation map

**Arm:** `receiver-map` (Phase ②: secretome OUTPUT → candidate-receiver receptor regulation)
**Branch:** `analysis/receiver-map` (off `analysis/manuscript-integration`)
**Pre-registration date:** 2026-06-23
**Tag (applied after commit, before inferential compute):** `prereg/receiver-map-v1`
**Status when written:** descriptive scan + Step-0 feasibility completed; **all inferential thresholds and interpretation rules frozen here, before the formal matched-null run.**

---

## 0. Framing (binding — do not relax downstream)

- **Negative-expected, descriptive, report-regardless.** Three prior receptor analyses in this repo are
  already null/negative: CONSEQUENCE composite receiver score (GSE325391, SAD−RES = −0.07 z, p=0.544, all
  14 robustness variants same-direction null); the neuron→receptor map (`results/regulatory/neuron_nichenet_receptors.csv`,
  `age_regulated=FALSE`); and `scripts/05_regulatory/06_autocrine_validation.R` (niche receptors not age-regulated).
  The **prior expectation is that receiver receptors are NOT coordinately age-up-regulated.** This expectation
  is stated up front and preserved throughout. This arm is **not** an attempt to manufacture or rescue a
  positive OUTPUT→receiver axis.
- **Purpose:** complete the triangulation honestly by building a *unified, pre-specified* receiver
  receptor-regulation map across all available compartments in the primary aging contrast — including the
  **previously unmapped neurogenic compartment (NSC, Neuroblast, Immature)** — and report the result
  regardless of direction.
- **Claim tier:** **SUPPORTING / descriptive** if informative; **NOT SUPPORTED** for coordinated receiver
  receptor up-regulation if the null pattern persists. **Not PRIMARY.** A coordinated positive receiver axis
  may be claimed **only** if it clears the pre-specified GO bar in §9.
- **Authorized null conclusion (verbatim):** *"The frozen age-UP ligand set is not accompanied by
  coordinated age-UP receptor regulation in candidate receiver compartments."* This is a valid
  SUPPORTING/NOT-SUPPORTED outcome, **not** a failed analysis.

---

## 1. Frozen inputs (verbatim reuse — NOT re-derived, NOT re-tested)

The **ligand side is a fixed input**, not a discovery step. The receptor universe is reused verbatim.

### 1.1 Ligand anchor — 22 frozen age-UP secreted ligands
Source of truth: `results/resilience/frozen_receiver_set_manifest.json` (`n_up_ligands: 22`), itself frozen
under `prereg/secretome-resilience-trajectory-v1` (prereg_commit `5b29c74…`), freeze commit `d74bdfd…`.
Provenance: `frozen UP secreted (34 unique symbols) → NicheNet ligands (22)` via
`scripts/10_resilience/00_ligand_receptor_map.R`.

```
ANXA1, APOC1, APOD, APOE, AZGP1, CCN2, COL12A1, COL21A1, COL4A2, COL8A1,
CXCL2, CYTL1, DKK2, EFEMP1, FAP, IL15, KITLG, LGALS9, SEMA3D, SERPINE1, SPON1, TNC
```

### 1.2 Receptor universe — frozen 76 (pre-detectability)
Ligand→receptor pairing from the SHA-pinned NicheNet network
`refs/nichenet/lr_network_human_21122021.rds`
(`sha256 = 47c971d2fbba4ecd0ba7485d1846a74054432a0d1a97ea3e6a79ae27d0da8094`).

**The 76 is the pre-detectability universe** = manifest `frozen_receptors` (56) ∪ `excluded_receptors` (20).
(Verified 2026-06-23: 56 ∪ 20 = 76 distinct symbols.) The regenerable intermediate
`results/step0/_receptors_prefilter.txt` is data-dependent and may be absent from a bare tree; **the committed
manifest is the authoritative source of truth.**

**We use the full frozen 76 — NOT the GSE325391-detectability-filtered 56**, because the 56 bakes in another
dataset's detectability and must not be transferred to GSE268609. A **fresh GSE268609 per-cell-type
detectability gate** (§4) is applied to the 76 independently within each receiver compartment.

### 1.3 Prior-based pairing disclosure (honesty — binding)
Ligand→receptor *pairing* is **prior-based** (NicheNet curated network). Only the **prioritization/weighting
is measurement-based** (observed receptor expression + observed receptor age-regulation). Methods must state
this plainly; we never claim a "prior-free" receiver map.

---

## 2. Receiver-side DE input (locked; no recompute)

- **File (locked input):** `results/de/de_GSE268609_neuron_aging_per_celltype.csv`
- **Design:** `~arm + grp` (uniform across **all** receiver compartments). `arm` controls the
  whole-hippocampus-vs-DG dissection/batch structure identified in the neuron arm; `grp` = aged (HA) vs
  young-adult (YA).
- **Why this table:** (i) the only table containing the neurogenic compartment (NSC/Neuroblast/Immature);
  (ii) one uniform design across niche, neuronal, and neurogenic receivers → an internally comparable
  cross-cell-type matrix; (iii) more defensible receiver-side design (controls dissection batch); (iv) the
  ligand set is frozen, so the original `~group` ligand discovery does **not** require the receiver-side test
  to share that design; (v) the niche signature is documented as surviving under `~arm + grp`.
- **No DE recompute for the primary analysis.** The DE table is the locked input. Re-running receptor DE is
  permitted **only** on a specific file-integrity problem (not for tuning).
- **Niche-on-both-sides disclosure (binding):** niche cells (Astro/Micro/Oligo/OPC/Endo) appear **both** as
  ligand-producing senders (the 22 ligands) **and** as candidate receiver compartments. This is acceptable as
  **generalized paracrine/autocrine mapping**; `autocrine_validation.R` already tested a narrower version
  (4 ligands). Here we extend to all 22 ligands × 76 receptors × all compartments. Methods must disclose it.

---

## 3. Receiver compartments and multiple-testing families

| Family | Compartments | n_donor (table) | Role |
|---|---|---|---|
| **Confirmatory (well-powered)** | Astro, Micro, Oligo, OPC, Endo, DG_GC, CA_ExN, InN | 16–17 | Primary inference; BH across these 8 |
| **Exploratory (neurogenic)** | NSC, Neuroblast, Immature | 12, 13, 12 (`exploratory=True`) | Depth-caveated; stricter exception bar (§9) |

Neurogenic compartments are **never** primary evidence on a nominal signal alone.

---

## 4. Detectability gate (direction-blind; per compartment)

A frozen-76 receptor is **detectable/testable** in compartment C iff it has a row in the locked `~arm + grp`
table for C with **non-NA `padj`** (equivalently non-NA `log2FoldChange` & `pvalue`) — the pipeline's own
expressed-gene / independent-filtering criterion (the `baseMean` proxy already used by
`autocrine_validation.R`).

- Applied **independently per compartment**.
- **Direction-blind:** never gate on the magnitude or sign of `log2FoldChange`, nominal p, `padj`, or any
  apparent age-regulation. Detectability is defined independently of direction and effect size.
- A receptor undetectable in C is **absent from C's row**, not counted as a null result (no dropout-as-evidence).
- **Report `n_detectable / 76` per compartment as a first-class descriptive** ("can this receiver receive?").

> **Clarification (added 2026-06-24, post-hoc; original prereg text above unchanged).** The parenthetical
> "(equivalently non-NA `log2FoldChange` & `pvalue`)" is **not exact** in DESeq2/pyDESeq2: `padj` is set to
> `NA` by independent filtering (low-`baseMean`) and by Cook's-distance outlier removal **even when**
> `log2FoldChange` and `pvalue` are non-NA. The implemented gate is **non-NA `padj`**, which is therefore the
> **stricter "testable after independent filtering"** criterion — not literal detectability/expression. Where
> this document and the downstream surfaces say "detectable", read **"testable after independent filtering"**.
> This is a terminology/precision clarification only; it does **not** change the gate, the frozen design, the
> compute, or any result (the implementation always used non-NA `padj`). See
> `audit_log/2026-06-23_receiver_map/RESOLUTION.md` §10.

---

## 5. Primary statistic — matched-null on `log2FoldChange` (ATAC/methyl analogue)

Faithful analogue of `scripts/04_atac/05_differential_accessibility.py` and
`scripts/08_methyl/03_ageslope_matched_null.py`, operating on the locked DE effect sizes. **Per compartment C:**

- **Foreground:** detectable frozen-76 receptors in C.
- **Background:** non-receptor genes expressed in C (non-NA `padj`), excluding the foreground, **selected
  direction-blind** and **matched to the foreground on `baseMean`** (deciles / nearest-neighbour bins;
  exact matching procedure recorded in code). Rationale: low-expressed genes have inflated |lfc|, so
  expression-matching is required.
- **Primary statistic:** mean **signed** receptor `log2FoldChange`, **oriented age-UP** (ligands are age-UP →
  the coordination hypothesis is receptors age-UP). **One-sided** empirical matched-null p from **2000
  permutations, seed=42**. **Two-sided** empirical p **always reported** so a coordinated down-shift / other
  structure cannot be hidden.
- **Effect size reported:** observed mean receptor `lfc`, null mean, observed − null (in `lfc` units),
  one- and two-sided empirical p, and `n_detectable`.

### 5.1 Secondary statistic — age dynamism
Mean **|`log2FoldChange`|** vs the same matched-null (analogous to the methyl module): are receptors more
age-variable than matched genes, regardless of direction? **Descriptive only — must NOT be used to claim
age-UP receiver activation.**

### 5.2 Pooling and not-testable rule
- **Pooled** at the receptor-set level per compartment. **No per-ligand formal tests** as primary (most
  ligands map too few receptors → underpowered, multiplies tests). The ligand×receptor matrix is shown as
  **annotation/descriptive structure only**, not a significance-mining layer.
- **Not-testable-not-forced:** a compartment with **< 3 detectable receptors** is marked **not testable**,
  not forced to a result (same minimum-feature logic as the existing matched-null modules).

---

## 6. Depth-matched detectability null (neurogenic repertoire control)

Motivated by the feasibility scan (§11): NSC detectability of the receptor probe set was low (~39% of the
56-probe), but **confounded with low `baseMean`, low nucleus count, and rare/exploratory status**. We will
**not** interpret "NSC detects fewer receptors" as receiver-intrinsic biology unless it is significantly lower
than a depth/expression-matched expectation.

- **Per compartment:** compare the observed fraction of the frozen-76 that is detectable against
  **size-matched random gene sets drawn from that same compartment**, matched as closely as possible to the
  receptor set's `baseMean`/expression distribution. **Selected without using receptor age-regulation
  direction or significance.**
- **Report:** empirical null distribution, observed statistic, empirical p, effect size.
- **Interpretation rule (binding):**
  - If receptor detectability ≈ the depth/expression-matched expectation → report as a **bounded caveat**,
    e.g. *"The neurogenic stem-cell compartment had limited receptor testability, but this was not separable
    from low-depth profiling of a rare exploratory population."* — **not** a biological finding.
  - If receptors remain **significantly depleted** after the matched control → report as a
    **hypothesis-generating neurogenic receptor-repertoire observation**, still **not** a signaling-activation
    claim unless receptor age-regulation also clears §9.
- This control is a **methodological safeguard and reportable audit feature**, not a stand-alone biological
  main finding.

---

## 7. Multiple-testing handling

- **Confirmatory family (8 compartments):** one-sided matched-null empirical p → **BH across the 8**; GO bar
  **q < 0.1** (§9). Matches the primary OUTPUT's `padj < 0.1` threshold and gives proper FDR control across
  compartments.
- **Neurogenic family (3 compartments):** report **nominal** empirical p and, as a conservative reference,
  the **all-compartment BH-q**. Nominal neurogenic signals are **never** primary evidence.
- **Per-receptor `padj < 0.1`:** matrix **annotation only**. Set-level positivity does **not** license naming
  individual receptors as significant unless a separately pre-specified receptor-level bar is met (not the
  goal here).

---

## 8. Leave-one-ligand-out (LOO) robustness (passers only)

For any compartment that passes §9: remove each ligand's mapped receptor block in turn and re-run the
matched-null. **Operational rule:** after removing each ligand's receptor block, the compartment must
**retain the same positive direction** and **must not collapse into an obvious single-ligand-driven effect**
(e.g. an integrin/collagen-heavy block). Report LOO empirical p-values as robustness diagnostics. **LOO is a
robustness check, not a new fishing layer.**

---

## 9. GO / NO-GO bar (frozen)

### 9.1 Confirmatory compartment — "coordinated age-UP receptor remodeling" requires ALL of:
1. **≥ 3 detectable receptors** (testability);
2. observed **mean oriented (age-UP) receptor `lfc` > 0** (direction);
3. one-sided matched-null empirical p passes **BH q < 0.1** across the 8 confirmatory compartments;
4. **LOO-robust** (§8);
5. **two-sided** empirical p reported alongside.

### 9.2 Neurogenic exception — requires 9.1(1–5) **AND**:
6. **survives the depth-matched detectability control** (§6);
7. **confirmed by the donor-level composite-score model** (§10), if feasible;
8. **rare-cell / low-depth caveat carried in every interpretation.**

### 9.3 NO-GO (expected)
If no compartment passes → the §0 authorized null conclusion; SUPPORTING / NOT-SUPPORTED; descriptive;
lightweight integration (§13).

---

## 10. Donor-level composite-score model — CONDITIONAL confirmation only

Run **only if** a compartment passes §9. CONSEQUENCE-style: per-donor pseudobulk logCPM of the detectable
receptor set → per-receptor z → mean composite per donor → `~ arm + grp` (group effect = age). This is the
**only** step permitted to recompute from counts, and **only** when there is a passer to confirm. Preserves
the locked-DE / no-recompute constraint and avoids threshold-shopping.

---

## 11. Feasibility scan disclosure (transparency; non-binding on thresholds)

A Step-0 feasibility/power scan was run on **2026-06-23** (frozen-56 probe subset, `~arm + grp` table) to
size the analysis and motivate the depth-matched-null control. It showed: (a) per-receptor age-regulation
essentially null — 4 receptor×compartment cells at `padj < 0.1` across 11 compartments (Astro 1, Micro 1,
Endo 1, NSC 1), i.e. ≈ the per-transcriptome baseline rate, receptors behaving like random genes; (b) a
detectability gradient with **NSC lowest (22/56 ≈ 39% detectable, mean `baseMean` ≈ 26)**, confounded with
depth. **This scan motivated §6 but did not set any threshold:** all §4–§10 thresholds are chosen on
principle (match primary `padj < 0.1`; house matched-null procedure) and frozen here before the formal run.
The formal analysis uses the **full frozen 76** (not the 56 probe) with a fresh per-cell-type gate.

---

## 12. Staged compute (`scripts/12_receiver/`; outputs → `results/receiver/`)

| Step | Script | env | Role |
|---|---|---|---|
| 00 | `00_freeze_universe.py` | bio | Reconstruct 76 = manifest(56 ∪ 20); assert count==76; verify rds `sha256` pin; emit `results/receiver/frozen_receptor_universe_76.csv` + per-receptor n_ligands |
| 01 | `01_detectability.py` | bio | Per-compartment detectability gate (§4); emit `detectability_matrix.csv` (76 × compartments) + `n_detectable` summary |
| 02 | `02_depth_matched_null.py` | bio | Depth-matched detectability null (§6), seed=42; emit `depth_matched_null.csv` |
| 03 | `03_matched_null_on_lfc.py` | bio | Primary matched-null-on-lfc (§5), seed=42, 2000×; emit `matched_null_per_compartment.csv` (+ dynamism) |
| 04 | `04_leave_one_ligand_out.py` | bio | LOO (§8) — **passers only** |
| 05 | `05_donor_composite.py` | bio/sgz_r | Donor composite model (§10) — **conditional, passers only** |

Seeds: `seed=42` everywhere (permutations + any jitter). No magic slices; iterate all compartments + a QC
predicate, log "N of M used".

---

## 13. Integration (default lightweight; STOP gate 2 first)

Default (expected null): **one SUPPORTING/descriptive paragraph in `manuscript.md` + matching paragraph in
`manuscript_ja.md`**; **one supplementary figure in `figures/figureS_receiver/`** (detectability matrix +
matched-null forest, both showing the null); legend + Source-Data CSV; **`make receiver`** committed-data
target (redraws from committed `results/receiver/*.csv`; the matched-null *compute* is data-dependent on the
gitignored DE table); `audit_log/2026-06-23_receiver_map/RESOLUTION.md`. **Expand only** if a neurogenic
exception clears §9.2. At integration: **鉄則9 grep-all-surfaces** across EN manuscript, JA manuscript,
figures, legends, Makefile, README, stdout, Source-Data, audit logs.

---

## 14. STOP gates (鉄則3 — no auto-chaining)

- **STOP gate 1** — after this plan is committed + tagged `prereg/receiver-map-v1`, **before** any formal
  inferential compute. Report: files committed, git hash, tag, ready-to-run status.
- **STOP gate 2** — after the formal compute (steps 00–04/05), **before** manuscript integration. Report:
  GO/NO-GO outcome, pass/fail table by compartment, neurogenic caveats, proposed claim tier. **Do not
  auto-chain into manuscript edits.**

---

## 15. Reproducibility

`seed=42`; envs `bio` (pandas/numpy matched-null) and `sgz_r` (only if re-deriving the rds, optional — the
committed manifest is authoritative); `PROJ` from environment; inputs SHA-referenced above. This plan is the
binding contract; deviations are logged in the eventual `audit_log/2026-06-23_receiver_map/RESOLUTION.md` with
rationale (鉄則9: a correction is a new analysis).

---

## AMENDMENT v1.1 (2026-06-23) — STEP 02 depth/testability decomposition (narrow, technical)

**Scope (binding): STEP 02 only.** This amendment does **NOT** change the primary matched-null-on-lfc
(STEP 03) result, the PASS bar (§9), the multiple-testing families (§7), the receptor universe (76), the
ligand universe (22), the locked DE table (§2), `seed=42`, the permutation count (2000), or any
interpretation rule. **STEP 03 remains exactly as computed under v1** (compute commit `f9d3af3`,
`results/receiver/matched_null_per_compartment.csv`) and is **not re-run** — there is no file-integrity
dependency requiring it. This is a clarification of the neurogenic-repertoire safeguard, **not** an attempt
to rescue a positive result. The accepted v1 conclusions stand: confirmatory 8/8 FAIL; no coordinated
age-UP receiver receptor remodeling; NSC depth-explained; claim tier NOT SUPPORTED / SUPPORTING-descriptive.

**v1 STEP 02 result (preserved as the original record):** `results/receiver/depth_matched_null.csv`.
Statistic = number of the full frozen 76 detectable (in-dds AND non-NA padj) vs `baseMean`-matched genes
drawn from the in-dds pool, absent receptors mapped to the lowest decile. NSC: obs 22 vs null 21.99,
effect +0.01, `p_low=0.69` (depth-explained). DG_GC/InN/Neuroblast/Immature: `p_low≈0.0005`.

**Observed artifact:** in compartments with near-complete global detectability
(`global_detectable_frac ≥ 0.98`: DG_GC 0.997, InN 0.998, Neuroblast 0.9992, Immature 0.9957), the in-dds
background pool is ~fully detectable at every `baseMean` decile, so the null predicts ~76; receptors that are
**absent from the DE table** (below the pipeline's gene-inclusion threshold) cannot be fairly matched and
generate an **uninterpretable depletion p-value**. **No headline claim depends on these four p-values.**

**v1.1 refinement (`scripts/12_receiver/02b_depth_testability_v1_1.py` → `depth_testability_v1_1.csv`):**
decompose receptor testability into explicit layers, per compartment —
- **A. Coverage / table inclusion:** `n_in_table / 76` (descriptive; contextualized by the compartment's
  total `n_genes_in_table`, since rare low-depth compartments include fewer genes overall).
- **B. In-table detectability:** among in-table (eligible) receptors, `n_detectable / n_in_table`
  (direction-blind non-NA padj gate; descriptive).
- **C. Expression/depth-matched null WITHIN the in-table universe only:** in-table receptors vs
  `baseMean`-matched in-table **non-receptor** genes (`seed=42`, 2000×, decile-matched). Absent receptors are
  **excluded** from the null (handled descriptively in layer A). For compartments with **saturated
  detectability** (background-pool detectable fraction `≥ 0.98`), layer C has no resolving power and is
  reported **"not interpretable / saturated detectability"**, not a p-value.

**Interpretation preserved:** for **NSC**, if layer C again shows detectability ≈ depth/expression-matched
expectation, retain: *"NSC had limited receptor testability, but this was explained by low-depth profiling of
a rare exploratory compartment rather than evidence for a receiver-intrinsic loss of receptor repertoire."*
For **Neuroblast/Immature**, make **no** biological depletion claim unless layer C is interpretable and
supports it; if saturated, report only the descriptive `n_detectable/76` + caveat.

Tag: `prereg/receiver-map-v1.1` (applied to this amendment commit, before the refined STEP 02 recompute).
