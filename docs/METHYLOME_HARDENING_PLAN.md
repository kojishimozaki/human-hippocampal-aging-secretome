# Methylome hardening plan (M1) — instructions for a future analysis session

**Status:** pre-specified analysis plan. **Lock the GO/NO-GO criteria BEFORE running**, report nulls as
results (this project's identity is rigor, not positive findings).
**Author context:** drafted 2026-06-05 in response to the PI's proposal to extend the regulatory layer
from ATAC to **DNA methylation / enhancer methylation** — testing *"secretome-linked regulatory regions
are age-dynamic in the methylome beyond a feature-matched genomic background"* (set-level, not
single-factor). The proposal was scrutinised and judged **sound and high-value, feasibility-gated on the
methylome arm's scale** (see §0.1). This plan is the natural successor to A2 (`docs/IN_SILICO_HARDENING_PLAN.md`),
which the cohort search already flagged as a "possible future independent epigenetic check"
(`audit_log/2026-06-04_in_silico_hardening/A1_COHORT_SEARCH.md:65`).

> Read this with the rigor-skill mindset (`single-cell-reanalysis-rigor` / `public-data-de-statistical-audit`):
> donor-level units only, matched nulls, FDR, within-celltype, universe-aware, pre-specified criteria,
> honest negatives. **This is publication-bound public-data reanalysis — follow the rigor skill.**

---

## 0. Orientation — read these first (the future session has no chat history)

- **Current claims (authoritative):** `manuscript/manuscript.md` (+`_ja.md`). (This plan originally pointed at `MANUSCRIPT_DRAFT.md`, which was superseded and moved to `docs/archive/` on 2026-08-26; corrected here so the pointer does not send a reader to withdrawn claims.) Fig 2 (secretome
  OUTPUT signature + external references) = PRIMARY; Fig 3 (chromatin) = SUPPORTING; mechanism =
  **honest dual-null** (cis-TF and autocrine both "not supported"). Read `results/regulatory/PHASE_REPORT.md`.
- **What's already done on the regulatory layer:** A2 (`scripts/04_atac/05_differential_accessibility.py`)
  — genome-wide donor-pseudobulk DA over 369,393 peaks = **0/369k at FDR** (underpowered at n=8/9), but
  secretome-linked peaks are enriched for **RNA-direction-concordant accessibility above a
  GC+width+accessibility+TSS-distance matched null** (Astro/Micro/OPC p≤8e-3; direction-specific).
  Record: `audit_log/2026-06-04_in_silico_hardening/A2_RESULTS.md`. **This methylome analysis is the
  orthogonal-modality version of that matched-null test.**
- **The external cohort:** GSE278576 (Zemke et al. 2024, PMID 39463924) = the snRNA+ATAC arm used for A1
  replication (40 donors, YA<40 n=10 / HA≥60 n=20; Astro/Micro/OPC/Oligo well-powered, **Endo
  NOT_EVALUABLE**). See `audit_log/2026-06-04_in_silico_hardening/A1_COHORT_SEARCH.md`.
- **Memory:** `~/.claude/.../memory/gse268609-covariate-and-fig3.md`.
- **Envs:** `conda activate bio` (scanpy/pydeseq2/numpy/scipy) and `sgz_r` (Seurat/Signac/GenomicRanges).
  Methylome processing likely needs **ALLCools** (`pip`/conda; methylpy ecosystem) — install into a
  dedicated env (`envs/methyl.yml`) and pin it. Set `PROJ=$(pwd)`. Seed 42. New scripts in `scripts/`
  (numbered, e.g. `scripts/08_methyl/`); write `audit_log/<date>_methylome_hardening/RESOLUTION.md`.

### 0.1 DATA-PROVENANCE CORRECTION (read before downloading anything)

The PI's note attributed the snm3C methylome to **GSE278576**. That is the snRNA+ATAC arm. **The
methylome/3D-genome arms are deposited separately** under the same study (PMID 39463924):

- **GSE299139 / GSE299899** — the epigenomic arms (snm3C / CUT&Tag), per this repo's own cohort search
  (`A1_COHORT_SEARCH.md:26`). **Which accession holds the base-level snm3C methylation (ALLC) must be
  confirmed in step 0** (I could not web-verify this session: the gstack `/browse` binary is the wrong
  architecture on this host and won't run).
- **4DN `4DNESATHQR75`** (data.4dnucleome.org) — the 3D-genome / contact arm of snm3C (snm3C yields BOTH
  per-nucleus methylation AND chromatin contacts; 4DN may host ALLC + `.cool`/contacts).

→ Confirm the **primary publication via `!Series_pubmed_id` = 39463924** for whichever accession you use,
and record the exact accession + file types in the RESOLUTION note (citation-provenance rule). **Do not
cite GSE278576 as the methylome source.**

### 0.2 Hard constraints (apply everywhere)

- **Donor is the unit.** donor×celltype **pseudobulk methylation**; no cell-level tests for inference.
- **mCG is primary; mCH is a glial caveat, not the lead.** The niche senders are **glia**
  (Astro/Micro/Oligo/OPC/Endo). Non-CpG methylation (mCH) is a **neuronal** feature (abundant in
  neurons, low in glia) — so for these cell types **lead with mCG**; report mCH only as a secondary
  check (and expect it to be low/uninformative in glia). Do not headline glial mCH.
- **"Modified cytosine", not pure 5mC.** Bisulfite/snm3C cannot separate 5mC from 5hmC. 5hmC is
  non-trivial in brain (higher in neurons; lower in glia). Read all signal as **modified cytosine** and
  keep 5mC/5hmC conflation as a stated caveat (GSE159671-type 5hmC data is bulk / weak cell-type
  resolution — at most a side note, not a dependency).
- **Matched null + FDR**, tested **above the empirical background**, never vs 0.5. **Match within
  region-class** on the methylation-specific confounders (§M1 method).
- **Within-celltype** for any methylation↔RNA / methylation↔ATAC covariation; never pool across celltypes.
- **Honest-null:** a clean "secretome regions are NOT more age-dynamic than matched background" is a
  publishable result and consistent with the dual-null identity. Pre-register that outcome as acceptable.

---

## M1 — Secretome regulatory-region methylation: age-dynamics vs a matched genomic background

**The question (set-level, pre-specified):** Are the **regulatory regions of secretome genes**
(promoters + linked enhancers) **more age-dynamic in DNA methylation** than feature-matched non-secretome
regions — and, where they change, do they change in the direction expected from the RNA signature
(aging-UP secretome → enhancer hypomethylation; aging-DOWN → hypermethylation)?

**Claim it would license (if PASS) — and the calibration that matters:**
- PASS licenses: *"the secretome's regulatory regions are **coordinately age-dynamic in DNA methylation**,
  beyond a feature-matched background, and directionally coupled to the RNA (and ATAC) changes"* — i.e. a
  **stronger, multi-modal SUPPORTING regulatory-coordination layer** (could upgrade Fig 3 from
  "ATAC direction-concordance" to "ATAC + methylation matched-null coordination across modalities").
- PASS does **NOT** license a causal/upstream-regulator claim. **Methylation change can be a
  consequence of expression/accessibility change, not a cause.** Mechanism (which TF/signal drives it)
  stays the open question. Also: mCG and ATAC are biologically anti-correlated, so a concordant
  methylation result is **partly corroborative-of-ATAC, not fully independent** — state that. So this
  *advances "the output is epigenetically coordinated"*, it does **not** "resolve mechanism".
- FAIL: report honestly ("secretome regulatory regions are not more age-dynamic in methylation than
  matched background"); the regulatory layer stays as A2 leaves it.

### Step 0 — FEASIBILITY GATE (do this first; stop and report before any heavy compute)

The make-or-break is the **methylome arm's scale** (snm3C cohorts are usually a *subset* of the RNA
cohort). Verify and **lock GO/NO-GO before processing**:

1. **Locate the data:** confirm which accession (GSE299139/GSE299899/4DN 4DNESATHQR75) holds base-level
   snm3C methylation (**ALLC** files or an equivalent per-cell mCG/CH count matrix) + cell metadata.
   Confirm PMID 39463924 via `!Series_pubmed_id`.
2. **Count donors with methylome**, split young vs old (need a young-vs-old contrast OR a usable
   continuous-age span). Record per-donor age/sex.
3. **Per-celltype nucleus coverage:** for Astro/Micro/Oligo/OPC/(Endo), how many nuclei per donor have
   adequate genome coverage? (snmC covers ~5–10% of CpGs per cell → per-region mCG is only stable after
   pseudobulk over **many** nuclei.) Expect **Endo to be unevaluable** (it was even in the 40-donor RNA
   arm: 641 cells).
4. **Cell-type annotations** for methylome nuclei (author labels, or transferable from the RNA arm).

**GO** if: ≥ ~8–10 donors spanning young **and** old (or a usable continuous-age range), with ≥3–4 niche
celltypes having enough nuclei/donor for stable region-level pseudobulk mCG. **NO-GO / downgrade** if the
methylome arm is tiny, all-old, or lacks glial annotation → say so and stop (do not force an underpowered
matched-null). **Report the step-0 findings and wait for PI confirmation before heavy compute.**

### Region definitions (build once, per secretome gene; same for the matched background)

Use the frozen secretome universe (`refs/secretome_union.csv`) and the frozen 60-hit signature
(`results/validation/frozen_primary_signature.csv`). Per gene, define and label by **region class**:
- **promoter:** TSS ± 1 kb (and ±2 kb as a sensitivity);
- **gene body;**
- **proximal regulatory:** ATAC peaks within gene body ± 2 kb (reuse the existing proximity links) and a
  ±10 kb sensitivity;
- **distal enhancer:** cCRE / ATAC peaks **linked** to the gene — prefer **contact-based linking from the
  snm3C 3D-genome** (ABC-style: accessibility × contact) over proximity, since this cohort *has* the
  contacts (a genuine upgrade over the project's current proximity-only links). If contacts are
  unavailable/too sparse, fall back to proximity and **say so**.

### Method (pre-specify; donor-level, within-celltype)

1. **Pseudobulk methylation:** per (donor × celltype), aggregate per-nucleus mC/coverage calls into
   region-level **mCG = methylated/total CpG basecalls** (and mCH separately, glia-caveat). Apply a
   **per-region coverage threshold** (e.g. ≥ N CpG basecalls per donor×celltype) and drop regions below
   it; report how many regions survive per celltype.
2. **Age model (donor = unit):** per region, per celltype, `mCG ~ age + sex + global_mCG/coverage_QC`
   (young-vs-old if groups exist, else continuous age slope). Extract the age coefficient/slope.
   *No cell-level tests.* BH-FDR within celltype.
3. **Matched background (the core test):** for the secretome-linked region set, build a matched
   non-secretome region set **within each region-class**, matching on **GC content, CpG density (the
   dominant methylation-specific confounder — add this to A2's covariates), region length, baseline mCG,
   mean accessibility, and distance-to-TSS** (KDTree nearest-neighbour, as in A2). 2,000× permutation.
4. **Set-level enrichment tests (pre-specified):**
   - (i) **Age-dynamism:** is |age-slope of mCG| in secretome regulatory regions **greater** than the
     matched background? (rank/AUC or 2-sample test, per celltype + pooled).
   - (ii) **Direction-concordance:** among regions linked to the 60 frozen hits, is mCG change
     **RNA-direction-concordant** (aging-UP gene → enhancer **hypo**methylation; aging-DOWN →
     **hyper**methylation) above the matched-background concordance rate? (binomial/Fisher **vs the
     empirical background rate, not 0.5** — same logic as Fig 3 / A2).
   - (iii) **Sign-permutation negative control:** the direction effect must vanish under label/sign
     permutation (as in A2).
5. **Cross-modal concordance (secondary, descriptive):** at gene level, does mCG age-change correlate
   with (a) RNA log2FC and (b) ATAC age-change, within celltype? Report ρ; **frame as coordination, not
   causation**, and note mCG↔ATAC are intrinsically anti-correlated.

### Pre-specified GO/NO-GO

- **PASS** ("methylation supports the regulatory coordination") if: (i) secretome regulatory regions are
  significantly **more age-dynamic** than the matched null (p<0.05) in ≥2 niche celltypes **AND** (ii)
  the change is **RNA-direction-concordant above background** **AND** (iii) the effect is
  **direction-specific** (dies under sign-permutation). Per-region FDR survivors are a bonus, not
  required (expect sparsity).
- **NO-GO / honest null** otherwise: report "secretome regulatory regions are not more age-dynamic in
  methylation than matched background"; Fig 3 stays as A2 leaves it. **Either outcome is reportable.**

**Output:** `results/methyl/methyl_ageslope_<celltype>.csv`, `results/methyl/secretome_region_methyl_enrichment.csv`
(+ matched-null + permutation p), a matched-null + cross-modal scatter panel. **Scripts:**
`scripts/08_methyl/{01_build_allc_pseudobulk,02_region_methylation,03_matched_null_enrichment}.py`.

---

## Caveats to embed in any write-up (pre-stated, so they're not "found" by a reviewer)

1. **Consequence vs cause:** methylation change may follow expression change; this is a coordination
   signal, not a causal mechanism.
2. **mCG↔ATAC not independent:** a concordant mCG result partly recapitulates the A2 ATAC result.
3. **5mC/5hmC conflation:** "modified cytosine"; 5hmC non-trivial in brain.
4. **mCH uninformative in glia:** lead with mCG; mCH is a neuronal mark.
5. **snmC sparsity:** per-region mCG needs many pooled nuclei; rare celltypes (Micro/Endo) may be
   unevaluable — say which, don't silently drop.
6. **Methylome cohort scale:** likely smaller than the 40-donor RNA arm; the contrast's power is the
   gate (§step 0).
7. **Linking:** prefer contact-based (snm3C) enhancer→gene links; if proximity-fallback, state it.

## Global rules (inherited from `docs/IN_SILICO_HARDENING_PLAN.md` §Global rules)

Donor-level units; pre-specify metric + GO/NO-GO before peeking; matched nulls + FDR above empirical
background (never vs 0.5); within-celltype only; universe-aware (secretome vs genome-wide, lead with
within-universe/directional); report nulls as results; provenance via `!Series_pubmed_id`; seed 42; new
code in `scripts/`; write `audit_log/<date>_methylome_hardening/RESOLUTION.md`; **commit/push only when
the PI says so.**

---

## Kick-off prompt (paste as the first message of a fresh session in this repo)

```text
SGZ 海馬加齢 niche-secretome プロジェクトの追加 in-silico 解析（M1: メチローム）を実行する。
作業ディレクトリは本リポジトリ。これは論文化を目指す公開単一細胞データの再解析なので、
single-cell-reanalysis-rigor スキルに従って進めること。

■ まず必ずこの順で読む（新セッションは履歴を持たない）:
  1. docs/METHYLOME_HARDENING_PLAN.md            ← 本タスクの指示書。GO/NO-GO 付き。これが正典
  2. docs/IN_SILICO_HARDENING_PLAN.md §Global rules + A2  ← 同型の matched-null を ATAC で既に実施済み
  3. manuscript/MANUSCRIPT_DRAFT.md（Fig2=PRIMARY / Fig3=SUPPORTING / 機構=honest dual-null）
  4. results/regulatory/PHASE_REPORT.md（機構は既に null）
  5. audit_log/2026-06-04_in_silico_hardening/{A1_COHORT_SEARCH,A2_RESULTS}.md

■ スコープ = M1 のみ:「secretome の調節領域（promoter + linked enhancer）が、matched background
  を超えて加齢で methylation 変化を受けるか／その方向は RNA signature と一致するか」を donor 単位・
  matched-null で検定する。set-level（特定因子ではない）。

■ 厳守事項（指示書の §0.2 と Global rules）:
  - データ出所の訂正: snm3C メチロームは GSE278576 ではなく GSE299139/GSE299899（+4DN 4DNESATHQR75,
    同一study PMID 39463924）。!Series_pubmed_id で確認し、ALLC（塩基レベル）の在処を特定する。
  - donor×celltype pseudobulk methylation のみ（cell-level 禁止）。mCG を主、mCH は glia では副次/注意。
  - matched-null は region-class 内で GC + CpG密度 + length + baseline mCG + accessibility + TSS距離 を
    マッチ（A2 に CpG密度を追加）。enrichment は empirical background に対して（0.5 帰無は使わない）。
    sign-permutation 対照。within-celltype。seed 42。
  - "modified cytosine"（5mC/5hmC 不分離）として慎重に。PASS でも因果機構ではなく「協調的に age-dynamic」
    までしか主張しない（methylation は結果の可能性／ATAC と非独立）。null も結果として報告。

■ 最初のチェックポイント（重い処理の前に必ず停止して報告）:
  指示書 §Step 0 のフィージビリティ・ゲート — メチローム arm の (a) donor 数と若齢/高齢構成、
  (b) niche celltype 別の nuclei coverage、(c) cell-type annotation の有無 を確認し、GO/NO-GO を
  audit_log/<今日の日付>_methylome_hardening/RESOLUTION.md にロックしてから私の確認を待つこと。

■ 運用: env は bio/sgz_r ＋ methylome 用に ALLCools を envs/methyl.yml に pin。新規コードは
  scripts/08_methyl/ に番号付きで。新規ブランチを切ってから着手。コミット/push は私が指示するまで保留。
```
