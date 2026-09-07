# Neuron–niche secretome comparison plan (N1–N3) — instructions for a future analysis session

**Status:** pre-specified analysis plan. **Lock the GO/NO-GO criteria and freeze the neuron analysis
BEFORE running any neuron DE.** Report nulls as results (this project's identity is rigor, not positive
findings). The niche side is already frozen — see `results/validation/frozen_primary_signature.csv`
(60 hits) and `results/validation/frozen_primary_secretome_logfc.csv` (per-celltype log2FC vector).

**Author context:** drafted 2026-06-18 from a `/grill-me` design session with the PI. The published
manuscript (`manuscript/manuscript.md`, GeroScience) is deliberately **niche-only, neurogenesis-agnostic,
AD-separate**: it analyses the secreted-protein output of the five niche *sender* cell types
(Astro/Micro/Oligo/OPC/Endo) and explicitly **excludes neurons**. The PI's proposal **inverts exactly that
exclusion**: characterise how the *neuron* (its secretome **and** its internal state) changes in normal
aging and in AD, and **position it against the niche signature** — is the niche aging-secretome
niche-specific, and do neurons respond to it? The matched material already exists: GSE268609 carries niche
and neurons in the **same donors**.

> Read this with the rigor-skill mindset (**`public-omics-reanalysis-protocol`**, which supersedes
> `single-cell-reanalysis-rigor` / `public-data-de-statistical-audit` and auto-fires for this work):
> donor-level units only, matched/permutation nulls, FDR within family, within-celltype, universe-aware,
> pre-specified criteria, honest negatives. **This is publication-bound public-data reanalysis — follow the
> rigor skill, and the project CLAUDE.md (no discovery verbs; exact p + effect size + n; trace-before-claim).**

---

## 0. Orientation — read these first (the future session has no chat history)

- **Current claims (authoritative):** `manuscript/manuscript.md` (+`_ja.md`). The niche **secretome-OUTPUT**
  signature (Fig 2) is PRIMARY; chromatin/methylation are SUPPORTING; **mechanism is an honest dual-null**
  (cis-TF and autocrine both "not supported"). `manuscript/implementation_notes.md` = the decision ledger.
  `manuscript/MANUSCRIPT_DRAFT*.md` are **superseded** narrative drafts — do not cite their numbers.
- **The matched cohort:** GSE268609 (Disouky A, …, Lazarov O. *Nature* 2026; PMID 41741649;
  doi:10.1038/s41586-026-10169-4). snRNA+ATAC multiome, 39 donors. Anchor h5ad with author cell-type labels
  and raw counts: `processed/per_dataset/GSE268609_anchor.h5ad` (`celltype_l1`, `donor_id`, `Group`).
  Build script: `scripts/02_qc/05_build_268609_anchor.py` (cluster→label map there).
- **The DE engine:** `scripts/03_de/01_pseudobulk_de.py` — donor-pseudobulk pyDESeq2, generic over
  `--celltypes`. The niche DE used only the 5 sender types; **neurons just need their labels added**.
  Existing outputs: `results/de/de_GSE268609_per_celltype.csv` (YA-vs-HA, niche only),
  `results/de/de_GSE268609_ADvHA.csv` (AD-vs-HA, niche only).
- **The replication engine (reuse for the niche↔neuron test):** `scripts/03_de/05_external_replication.py`
  — the niche signature was *tested* in GSE278576 with this (ρ≈0.24). **The same machinery, pointed at
  neuron cell types instead of an external cohort, is the niche↔neuron specificity/correspondence test.**
- **Secretome lens (frozen, reuse as-is for symmetry):** `refs/secretome_union.csv` (2,224 genes;
  columns gene/n_sources/sources/is_core/category) and `refs/secretome_core.csv` (1,813, sensitivity).
- **External references:** SenMayo `refs/senmayo.txt`, SASP Atlas `refs/sasp_atlas_core.txt`,
  Wingo cognition β `refs/wingo_cogtraj_beta.csv`, CSF aging `refs/csf_aging_estimate.csv`,
  NicheNet priors `refs/nichenet/`.
- **Power probe (read-only, already run):** `scripts/03_de/_power_probe.py` → the donors/group×celltype
  table embedded in §0.2.
- **Memory:** `~/.claude/.../memory/neuron-niche-secretome-plan.md` and `…/gse268609-covariate-and-fig3.md`.
- **Envs:** `conda activate bio` (scanpy 1.11.5 / pydeseq2 0.5.4 / numpy / scipy) for DE & comparison;
  `sgz_r` (Seurat/Signac) for NicheNet (N3-B). Pathway enrichment (N2) may need **gseapy** — install into
  `bio` (or a pinned env) and record the version. Set `PROJ=$(pwd)`. **Seed 42.** New code under
  `scripts/` (numbered, e.g. `scripts/09_neuron/`). Write `audit_log/<date>_neuron_niche/RESOLUTION.md`.

### 0.1 Hard constraints (apply everywhere)

- **Donor is the unit.** donor×celltype pseudobulk; no cell-level tests for inference. ≥10 cells/donor,
  ≥4 donors/group (per-group, not total — the engine already enforces this).
- **Even-handed, whole-hippocampus, no a-priori DG privilege.** The published niche signature is itself
  **pan-hippocampal** (Visium: median |ρ|<0.04 vs the DG marker PROX1), so symmetry *forbids* centring the
  neuron analysis on the dentate gyrus. Analyse DG_GC / CA_ExN / InN **equally**. DG only becomes the focus
  if the pre-registered **DG-trigger** fires (see that section) — otherwise report a mature-neuron panel.
- **"Control" is dual.** CA_ExN/InN are **not SGZ-lineage** but **are paracrine-exposed** to the niche
  secretome (which bathes the whole tissue). They control the *developmental-lineage* axis, **not** the
  *exposure* axis. So "correspondence in all three neuron types" can mean **paracrine niche exposure**
  (still niche-coupled) OR generic neuronal aging — disambiguate with N3 (ligand→receptor map), never read
  "all three move" as automatically "global, niche-irrelevant".
- **Disease axis kept clean (mirror the niche discipline).** Normal aging = **YA-vs-HA**; AD =
  **AD-vs-HA** (age-matched aged reference, isolates the AD increment). Keep them separate, then correlate
  the two log2FC vectors (aging↔AD). AD-vs-YA is descriptive only.
- **SA group excluded from primary.** SA is 100% DG-microdissect arm + extreme age (86–100) + n=6 ⇒
  exclude from primary; exploratory only.
- **Tissue-prep / two-arm batch structure — CHARACTERIZED (N0 gate run 2026-06-18; covariate required, NOT
  blocking; supersedes the earlier "blocking" note).** GSE268609 has two arms: a **whole-hippocampus
  multiome** (GEO titles `1_GEX`–`14_GEX`, `orig.ident` 1–14, 14 donors) and a **dentate-gyrus-microdissected**
  arm (titles `15`–`39 RNA/ATAC`, `orig.ident` 15–39, 25 donors); the anchor pools both. **The feared
  DG_GC-vs-CA_ExN confound is REFUTED data-intrinsically** (`_tissueprep_probe.py`): every one of the 39
  donors contains both DG_GC and CA_ExN (≥10 each; zero DG-only donors; CA fraction 0.01–0.30) — the
  microdissection retained CA pyramidal neurons, so both cell types share the donor pool and the
  even-handed design is valid. **But the arm is a real batch variable, imbalanced across contrast groups**
  (DG-arm fraction: YA 75%, HA 56%, AD 40%, MCI 67%, SA 100%). → N0 MUST recover the arm label (anchor lacks
  it; map cell→`orig.ident` via `processed/per_dataset/GSE268609_metadata.csv`, `orig.ident`≤14 = whole-hippo
  / ≥15 = DG), add it to the anchor, and **model it as a batch covariate** (`~arm + group`); whole-hippo-only
  is too thin to restrict to (YA n=2), so covariate-on-full-cohort + a thin whole-hippo sensitivity. **The
  same unmodeled imbalance sits in the published niche YA-vs-HA → run a retrospective `~arm+group`
  robustness check on the current submission (intellectual-honesty item).** Pre-register all of this.
  (Earlier `metadata_master` GSM-level counts that looked SA-only / then all-mixed described the full GEO
  deposit, not the analysed anchor — corrected by the data-intrinsic probe.)
- **Do not re-run Disouky's neuron work as a new claim.** The primary paper (GEO `!Series_summary`) already
  reported: PV+ inhibitory-neuron loss in aging; immature-neuron transcriptional shutdown
  (ribosomal/mitochondrial down) in AD; neuroblast/immature/GABAergic inverse-correlation with pathology;
  and **CellChat/NeuronChat ligand–receptor inference for neurogenic cells**. Our lane is **orthogonal**:
  secretome-OUTPUT lens, **mature** neurons head-to-head with the **niche in the same donors**, the
  normal-aging axis, and external SASP/proteome anchoring. Engage Disouky's findings only as
  "consistent-with" cross-checks, never as our discovery.

### 0.2 Locked design (the `/grill-me` output, 2026-06-18 — all PI-confirmed)

| # | Decision |
|---|---|
| 1 Frame | **C**: de-novo neuron DE foundation + niche-specificity contrast. Narrative weight = niche-specificity / positioning. |
| 2 Structure | **3 layers**: **L1** neuron secretome (sender-parallel) / **L2** internal factors (genome-wide intrinsic) / **L3** correspondence (**A** statistical correspondence = core; **B** NicheNet niche→neuron receiver = supplemental). |
| 3 Cohort | discovery = **GSE268609** (niche+neuron, same donors). Replication is **finding-driven**: GSE325391 (granule-AD; `de_GSE325391_*` exist), GSE278576 (aging-neuron). |
| 4 Cell types | **DG_GC / CA_ExN / InN even-handed.** immature (NSC/Neuroblast/Immature) = normal-aging exploratory only. **SA excluded.** |
| 5 Disease axis | **HA-anchored two-axis** (YA-vs-HA aging / AD-vs-HA disease) + **aging↔AD log2FC correlation** (accelerated-aging vs divergent, compared niche-vs-neuron). **Extend the niche side to AD-vs-HA secretome too** (the paper only published YA-vs-HA). MCI-vs-HA exploratory; AD-vs-YA descriptive. |
| 6 Comparison metric | **even-handed global = primary** (de-novo neuron DE + niche↔neuron log2FC correlation via `05_external_replication.py`) + **anchored 60-hit = confirmatory** + category-level (esp. **neuropeptide / growth-factor**). *L1 specificity ≡ L3-A correspondence = one metric, two framings.* |
| 7 L2 content | genome-wide DE → **unbiased pathway + targeted "niche-response" readout + donor-level coupling** (same-donor strength). Disouky's turf = cross-check only. |
| 8 Lens | same **2,224 union** (symmetry is mandatory for the comparison), core 1,813 as sensitivity. |
| 9 Pre-reg defaults | freeze before running; BH per (celltype × axis); **frozen padj<0.1 + extended padj<0.2**; 2,000× permutation nulls; **PMI sensitivity `~pmi+group` both axes**; leave-one-donor-out; seed 42; exact p + effect size + n. |
| 10 DG-trigger | **conjunctive strict** (DG significant **AND** significantly > CA_ExN *and* InN). Else → mature-neuron panel. Spin-out to a separate claim additionally requires **GSE325391 replication**. |
| 11 External anchors | **SASP/SenMayo always** (as a specificity test: "niche aligns, neuron does not" = external evidence of niche-specificity). Wingo / CSF / neuron-specific neurotrophic refs = **finding-driven** (verify new refs via literature-review first). |
| 12 Output | **data-gated**: (a) supplement that armors the current submission, or (b) a separate neuron paper. Decide *after* seeing N1–N3 + the DG-trigger + replication. |

---

## N0 — Build the DE and pass the verification gate (do first; stop and report)

1. **Cell-type labels** (already in `celltype_l1`): neurons = `DG_GC` (mGC, mature dentate granule),
   `CA_ExN` (CA pyramidal), `InN` (GABAergic); immature = `NSC`, `Neuroblast`, `Immature`. Niche =
   `Astro`/`Micro`/`Oligo`/`OPC`/`Endo`.
2. **Verify the tissue-prep field** per donor for YA/HA/MCI/AD (see §0.1). Lock how it is handled
   (covariate / restrict / pooled-with-caveat) and record in the RESOLUTION note **before** DE.
3. **Run the DE** (`scripts/03_de/01_pseudobulk_de.py`, mode `group`, donor-col `donor_id`,
   celltype-col `celltype_l1`) for the matched niche+neuron cell set, three contrasts:
   - **YA-vs-HA** (normal aging) — test=HA ref=YA (match the existing niche convention/sign).
   - **AD-vs-HA** (AD) — test=AD ref=HA.
   - **MCI-vs-HA** (exploratory).
   Output `results/de/de_GSE268609_neuron_{aging,ADvHA,MCIvHA}_per_celltype.csv` (or extend the existing
   per-celltype files). Run the **PMI-adjusted** (`~pmi+group`) and **LOO** robustness variants per the
   existing `scripts/03_de/04_aging_robustness.py` pattern.

**Power (from `_power_probe.py`; donors/group with ≥10 cells; DE needs ≥4/group):**

| celltype | n_cells | YA | HA | MCI | AD | SA |
|---|---|---|---|---|---|---|
| DG_GC | 16,882 | 7 | 9 | 6 | 10 | 6 |
| CA_ExN | 12,679 | 8 | 9 | 6 | 10 | 6 |
| InN | 7,650 | 8 | 9 | 6 | 9 | 6 |
| NSC | 877 | 8 | 4 | 4 | 8 | 3 |
| Neuroblast | 1,309 | 7 | 6 | 4 | **1** | 5 |
| Immature | 1,366 | 5 | 7 | 5 | **0** | 6 |

→ **DG_GC / CA_ExN / InN are fully powered for YA-vs-HA and AD-vs-HA.** Immature lineage: normal-aging
exploratory only; **Neuroblast/Immature AD is infeasible** (1/0 donors); NSC is borderline (HA=4) and rare
— exploratory with the thinness stated. Donor groups: YA n=8 (30 y), HA n=9 (81 y), MCI n=6 (86 y),
AD n=10 (84 y); PMI YA 5.84 / HA 7.67 / MCI 7.76 / AD 5.65 (so AD-vs-HA carries a PMI gap → the PMI
sensitivity is required, same as the niche YA-vs-HA).

---

## N1 — Neuron secretome (Layer 1): is the niche aging-secretome niche-specific?

**Question (pre-specified):** Does the neuron reprogram its **secreted-protein output** with aging/AD, and
does that program **match or differ from** the niche's? Read low niche↔neuron alignment as *niche-specific*
(and low correspondence); high alignment as *shared program* (and high correspondence) — same number, both
framings.

**Method.**
1. **De-novo neuron secretome DE** — intersect the N0 neuron DE (both axes) with `secretome_union.csv`;
   call hits at **padj<0.1 (frozen) and <0.2 (extended)**, per (celltype × axis). This is the neuron's own
   program, with **no niche anchoring**.
2. **Global niche↔neuron log2FC correlation** *(primary, even-handed)* — Spearman of per-celltype secretome
   log2FC, every niche-type × neuron-type pair, over the shared secretome universe, with a 2,000×
   permutation null. **Reuse `05_external_replication.py`** with neuron cell types in place of the external
   cohort. (Contrast: GSE278576 gave ρ≈0.24 as a genuine *replication*; a low neuron ρ means "neuron is not
   a replica of the niche" = niche-specific, shown with validated code.)
3. **Anchored 60-hit test** *(confirmatory)* — do the frozen 60 niche hits move concordantly in neurons?
   direction-concordance vs the neuron's **own** secretome background (binomial/Fisher vs background, not
   0.5 — the project's standard).
4. **Category-level** — cytokine / ECM / growth-factor / **neuropeptide** × celltype heatmap. Neuropeptide
   and growth-factor are the candidate **neuron-specific** modules (neurons secrete neuropeptides,
   BDNF/NTF3, Wnts) — a biologically meaningful asymmetry, not noise.

**PASS/report framing (every outcome reportable):**
- **niche-specific** if niche↔neuron secretome ρ is low/non-significant AND the 60 hits do not move
  concordantly in neurons (above their own background) → external-validity for the niche claim; armors the
  current paper.
- **shared/global** if ρ is high AND hits move concordantly → the signature is a tissue-wide aging
  secretome; reframe (and route to N3 to test paracrine exposure vs generic aging).
- Either way report the **neuron's own** secretome hits (the de-novo program) as a result.

---

## N2 — Neuron internal factors (Layer 2): what changes inside the neuron, and does it track the niche?

**Question (pre-specified):** Beyond secretion, what **intrinsic** programs shift in neurons with
aging/AD, and is the shift **consistent with a response to the niche's SASP-like output**?

**Method (donor-level, within-celltype).**
1. **Genome-wide neuron DE** — the unfiltered N0 neuron DE (all genes), both axes, BH per (celltype × axis).
2. **Unbiased pathway enrichment** — Hallmark / Reactome / GO on the ranked log2FC (gseapy GSEA or similar);
   let the data name the programs (may surface neuron-specific aging biology Disouky did not frame as
   secretome-response).
3. **Targeted "niche-response" readout** *(pre-specified gene sets)* — test whether the neuron intrinsic DE
   is enriched, in the niche-predicted direction, for: inflammatory/NF-κB & interferon response, complement
   and cytokine **receptors**, integrated stress response (ATF4/ISR), and senescence-like markers
   (e.g. CDKN1A). Descriptive/correlative — **no causal language**.
4. **Donor-level coupling** *(the matched-design strength)* — niche and neuron are the **same donors**:
   per donor compute a niche-SASP-output score and a neuron-response score (and a neuron-secretome score),
   then **correlate across donors** (Spearman, with the PMI/age caveats). This is a within-subject coupling
   test the niche paper could not do for methylation (it had to settle for cross-cohort effect-size
   concordance).

**Report framing:** "neuron intrinsic state is / is not enriched for a niche-SASP-response program, and
does / does not co-vary with niche output across donors." Disouky's PV+ loss / immature ribo-mito shutdown
appear only as "consistent-with" notes.

---

## N3 — Correspondence (Layer 3): how the neuron relates to the aging niche secretome

**A — statistical correspondence *(core; implement first)*.** This is the N1 global correlation + the N2
donor-level coupling, read as the niche↔neuron *coupling* (not just specificity). **No ligand–receptor
mechanism here** — deliberately disjoint from Disouky's CellChat/NeuronChat.

**B — NicheNet niche→neuron receiver *(supplemental data; after A)*.** Take the niche aging-secretome
ligands (60 hits + broader UP/DOWN), map ligand→receptor with `refs/nichenet/lr_network_human_*.rds`, and
test on neurons: (i) are the receptors expressed; (ii) do receptor genes move with aging/AD on the neuron
side; (iii) NicheNet ligand-activity predicting neuron aging DE, with the existing 2,000× permutation null.
**Caveats (state them):** the niche paper's *autocrine* NicheNet was non-discriminative (null passed
~82–91% of ligands) — the paracrine niche→neuron version inherits that weak calibration; snRNA receptors
are sparse; this is correlative, not causal; and it overlaps Disouky's method-space, so frame strictly as
honest exploration.

---

## Cross-cutting: the aging↔AD contrast (the PI's "通常加齢と AD")

Per cell type (niche and neuron), correlate the **YA-vs-HA** and **AD-vs-HA** secretome (N1) and
genome-wide (N2) log2FC vectors:
- positive, same-direction, amplified ⇒ **AD continues/accelerates aging**;
- orthogonal/anti-correlated ⇒ **AD is a distinct program**.
The novel comparison is **whether this aging↔AD relationship differs between the niche and the neuron** —
which Disouky did not do. MCI sits as the exploratory intermediate (monotonicity check, n=6).

## DG-trigger (pre-registered — lock before peeking at neuron results)

Deep-dive DG_GC **only if** it is BOTH (1) significant on its own (de-novo secretome hit at padj<0.1, or
N3 correspondence / donor-coupling significant vs null) **AND** (2) significantly stronger than **both**
CA_ExN and InN (Δ tested by permutation). Otherwise report the **even mature-neuron panel** — no DG
privilege. A **separate paper/claim** about DG additionally requires **GSE325391 granule replication**.

## Output decision gate (data-gated)

After N1–N3 + DG-trigger: (a) **supplement** armoring the current GeroScience submission (most likely if
the result is "niche-specific, neuron does not mirror it"), or (b) **separate neuron paper** (only with a
standalone neuron/DG finding **and** GSE325391/GSE278576 replication). Decide with the PI; do not force.

## Caveats to embed in any write-up (pre-stated, so a reviewer does not "find" them)

1. **Disouky overlap** — mature-neuron + secretome-output + niche-symmetry is our lane; PV+ loss, immature
   shutdown, and CellChat/NeuronChat are theirs (cross-check only).
2. **Post-mitotic SASP** — "neuronal senescence/SASP" is contested; neuron non-alignment to SASP is the
   expected, specificity-supporting outcome, not a failure.
3. **Paracrine vs lineage** — correspondence in CA_ExN/InN can be paracrine exposure, not "global/aging".
4. **Correlative** — donor-coupling and ligand-activity are coordination, not causation (no discovery verbs).
5. **Nuclear RNA** — snRNA is an imperfect proxy for secreted protein (mitigated by the protein refs).
6. **Tissue-prep** — state how whole-hippo vs DG-microdissect was handled; SA excluded.
7. **n and PMI** — small per-group n; PMI differs (report the `~pmi+group` sensitivity and LOO).

## Global rules (inherited from `docs/IN_SILICO_HARDENING_PLAN.md` + project CLAUDE.md)

Donor-level units; freeze + pre-specify metric and GO/NO-GO before peeking; matched/permutation nulls + FDR
above the empirical background (never vs 0.5); within-celltype; universe-aware (secretome vs genome-wide;
lead with within-universe/directional); report nulls as results; provenance via `!Series_pubmed_id`;
exact p + effect size + n; seed 42; new code in `scripts/`; write
`audit_log/<date>_neuron_niche/RESOLUTION.md`; **commit/push (and branch) only when the PI says so.**

---

## Kick-off prompt (paste as the first message of a fresh session in this repo)

```text
SGZ 海馬加齢プロジェクトの新規 in-silico 解析（N: ニューロン vs ニッチ secretome/内的因子）を実行する。
作業ディレクトリは本リポジトリ。論文化を目指す公開単一細胞データの再解析なので、
public-omics-reanalysis-protocol スキルに従って進めること（解析開始前に発動）。

■ まず必ずこの順で読む（新セッションは履歴を持たない）:
  1. docs/NEURON_NICHE_PLAN.md                         ← 本タスクの指示書（正典）。§0.2 にロック済み設計
  2. manuscript/manuscript.md（niche secretome-OUTPUT=PRIMARY / 機構=honest dual-null / AD は分離）
  3. scripts/03_de/01_pseudobulk_de.py（DE エンジン）, 05_external_replication.py（niche↔neuron 比較に流用）
  4. scripts/02_qc/05_build_268609_anchor.py（cluster→celltype マップ）

■ スコープ = 同一ドナー(GSE268609)で成熟ニューロン(DG_GC/CA_ExN/InN)の secretome(L1) と
  内的因子(L2, genome-wide) の加齢(YA-vs-HA)・AD(AD-vs-HA)変化を「無偏に」DE し、既凍結の niche signature と
  head-to-head 対比(L3-A) する。L3-B(NicheNet niche→neuron) は supplemental。DG は特権化しない
  （DG-trigger が発火した時のみ深掘り）。SA 除外。

■ 厳守事項（§0.1 と Global rules）:
  - donor 単位 pseudobulk pyDESeq2、≥10 cells/donor・≥4 donors/group、padj<0.1+0.2、perm 2000×、
    PMI 感度 ~pmi+group 両軸、LOO、seed 42、正確 p+効果量+n。
  - even-handed（DG 偏重禁止＝niche signature は pan-hippocampal, |ρ|<0.04 vs PROX1）。
  - Disouky 本体の縄張り（PV+喪失/immature shutdown/CellChat-NeuronChat）は再走しない＝cross-check のみ。
  - null も結果として報告（"niche 揃う/neuron 揃わない" は niche 特異性の積極的証拠）。因果語は使わない。

■ 最初のチェックポイント（DE を回す前に必ず停止して報告）:
  N0 の検証ゲート — (a) tissue-prep（whole-hippo vs DG-microdissect）の混在を donor 別 tissue フィールドで
  確認し処理方針をロック、(b) SA 除外の確認、を audit_log/<今日の日付>_neuron_niche/RESOLUTION.md に
  記録してから私の確認を待つこと。

■ 運用: env は bio（scanpy/pydeseq2）＋ sgz_r（NicheNet）。pathway 用に gseapy を必要なら pin。
  新規コードは scripts/09_neuron/ に番号付きで。新規ブランチを切ってから着手。コミット/push は私が指示するまで保留。
```
