# In-silico hardening plan — instructions for a future analysis session

**Status:** pre-specified analysis plan. **Lock the GO/NO-GO criteria in each section BEFORE running**,
and report nulls as results (this project's identity is rigor, not positive findings).
**Author context:** drafted 2026-06-04 after the codex v2.4/v2.5 publication-readiness audits, in
response to a codex proposal of 8 in-silico analyses to "make stronger claims". The 8 were triaged
against this dataset's real constraints; only the high-value/feasible ones are specced below.

> Read this with the rigor skill mindset (`single-cell-reanalysis-rigor` /
> `public-data-de-statistical-audit`): donor-level units only, matched nulls, FDR, within-celltype,
> universe-aware, pre-specified criteria, honest negatives.

---

## 0. Orientation — read these first (the future session has no chat history)

- **Current claims (authoritative):** `manuscript/manuscript.md` (+`_ja.md`). (This plan originally pointed at `MANUSCRIPT_DRAFT.md`, which was superseded and moved to `docs/archive/` on 2026-08-26; corrected here so the pointer does not send a reader to withdrawn claims.) Fig 2 (secretome
  OUTPUT signature + 4 protein-level external references) is PRIMARY; Fig 3 (chromatin) is SUPPORTING
  (68% direction-concordance, borderline p≈0.06, 0/4,971 genome-wide ATAC FDR); mechanism = **honest
  dual-null**.
- **Mechanism is already exhausted:** `results/regulatory/PHASE_REPORT.md` — Layer A (motif
  enrichment) FDR-null; B (NicheNet) non-discriminative; C (donor footprinting) 0/21, nominal
  *against*; D (L-R validation) fails. **Do not expect a positive by re-running these.**
- **Audit trail:** `audit_log/2026-06-03_codex_audit_v2_statistical/RESOLUTION.md` and
  `audit_log/2026-06-04_codex_audit_v2.4_publication_readiness/RESOLUTION.md`.
- **Memory:** `~/.claude/.../memory/gse268609-covariate-and-fig3.md`.
- **Hard constraints to respect everywhere:**
  - **n = 8 YA / 9 HA donors** (GSE268609). Donor is the unit. Cell-level p-values are inflated and
    forbidden for inference. Most per-feature FDR tests are underpowered at this n.
  - Secretome universe vs genome-wide background changes every enrichment result — report both,
    lead with the within-secretome-universe / directional statistic.
  - ATAC accessibility must be depth-normalised by donor total fragments-in-peaks (`nCount_Peaks`),
    not a peak subset (audit v2.4).
  - Say "externally validated" / "candidate replication"; never "replicated" unless A1 PASSES.
- **Envs:** `conda activate bio` (scanpy/pydeseq2/numpy/scipy) and `sgz_r` (Seurat/Signac/ArchR/
  chromVAR/nichenetr). Set `PROJ=$(pwd)`. Seed 42. New scripts go in `scripts/` (numbered);
  write a `RESOLUTION`-style note per analysis in `audit_log/<date>_in_silico_hardening/`.
- **On-disk data already downloaded (raw/):** GSE268609 (multiome incl. 8.9 G peak MTX + 77 G
  fragments), **GSE186538 Franjic** (human/rhesus/pig HPC, subregion dissections, UNTAPPED),
  GSE199243 (glia lifespan — already the signed-null replication attempt), GSE185277/185553 (Zhou
  lifespan/adult), GSE163737 (aged human+macaque, **all-old**), GSE160189 (HPC A/P, epilepsy),
  GSE264624/264692 (Tran snRNA + Visium), GSE325391 (resilience).

### Triage summary of the 8 codex proposals

| codex # | proposal | verdict | why |
|---|---|---|---|
| 1 | independent sc/snRNA replication | **TIER A — do** | only 1 attempt so far (GSE199243, signed-null); Franjic on disk is untapped; strengthens the core claim |
| 2 | genome-wide DA peak + secretome-linked enrichment | **TIER A — do** | full peak MTX on disk; completes the *deferred* accessibility-matched null; can have power at set level even if per-peak FDR is sparse |
| 3 | peak-gene link re-build (ArchR/Cicero/scGLUE) + matched null + FDR | **TIER B — do if useful** | current link test lacks matched-null/FDR; this credibly upgrades "proximity"→"computationally linked" (or confirms null with a standard tool) |
| 5 | SGZ specificity (region×age) | **TIER B — REDIRECT** | Visium can't resolve DG → do it on **Franjic subregion dissections** (DG vs CA/EC/SUB) instead; honest ceiling |
| 4 | TF mechanism triangulation | **TIER C — already done** | Layers A–D already triangulate TF-RNA/chromVAR/footprint/motif/target/link → dual-null. No new positive without new data |
| 7 | L-R sender/receiver fixed | **TIER C — already done** | Layer D already used permutation-null + L-R expression filter (exactly codex's spec) → failed. Null stands |
| 6 | Visium deconvolution correction | **TIER C — defer** | contingent on the resolution-limited Visium (#5); low yield vs effort |
| 8 | mediation-style analysis | **TIER C — against (ranking only)** | n=8/9 makes mediation statistically hopeless; use at most for candidate ranking, never for any claim |

---

## TIER A — highest value, feasible, strengthens the existing core claim

### A1 — Independent human-cohort replication of the niche-secretome aging signature  (codex #1)

**Claim it would license (if PASS):** upgrade "externally validated at the protein level" to "+ an
independent *transcriptomic* replication in human hippocampus" — i.e. closer to multi-cohort, without
saying "replicated" unless criteria are met.

**Data (in priority order):**
1. **GEO search first (step 0):** look for *post-2023* human hippocampus/DG **aging** snRNA/scRNA with
   **both young-adult and aged donors**, donor-resolved (not pooled), niche glia present. Search terms:
   "human hippocampus snRNA aging", "dentate gyrus single-nucleus aging", check GEO/CELLxGENE/Synapse.
   Record primary publication via `!Series_pubmed_id` (citation-provenance rule).
2. **GSE186538 Franjic (on disk, UNTAPPED):** human HPC with subregion dissections. **Caveat: age
   range is compressed (~44–79 y, no true young-adult)** → this is a *mid-to-old age-slope* test, not
   YA-vs-HA. Still a legitimate independent human substrate, and it doubles for B2.
3. **GSE185277 Zhou lifespan (secondary):** broad age incl. aged; check donor/age composition and
   whether ≥4 aged + ≥4 younger adult donors exist per niche celltype.

**Method (pre-specify, no peeking):**
- Re-annotate Astro/Micro/Endo/OPC/Oligo with the **same markers** used for the GSE199243 glia
  re-annotation (`scripts/02_qc/07_reannotate_glia.py`) — or use author labels if present (Franjic
  has them). Document marker provenance.
- **Donor pseudobulk** DE: YA-vs-HA if a young group exists, else **continuous-age** slope (design
  `~age`), pyDESeq2, ≥10 cells/donor, ≥4 donors/group. Donor = unit (no cell-level).
- **Freeze the primary signature** (the 60 RNA-significant secretome gene×celltype hits and the
  per-celltype log2FC vector from `results/de/de_GSE268609_per_celltype.csv`) and *test it* in the
  external cohort — do not re-discover.

**Pre-specified metrics + GO/NO-GO (all donor-level):**
- (i) **Signature-level directional test:** one-sided Mann-Whitney / linear test that the external
  log2FC of the primary aging-UP set > aging-DOWN set (or >0 for UP). PASS if p<0.05.
- (ii) **log2FC correlation:** Spearman(primary log2FC, external log2FC) over shared secretome genes,
  per niche celltype and pooled. PASS if ρ>0 with p<0.05 in ≥2 niche celltypes.
- (iii) **Sign-concordance of the 60 hits**, tested **above the external cohort's own background
  concordance rate** (same logic as Fig 3 — NOT vs 0.5), binomial + Fisher. PASS if significantly
  above background.
- (iv) **Null control:** sign-flip / label-permutation null (≥1,000×) for (i)–(iii).
- **OVERALL GO** (call it "independent transcriptomic replication") only if ≥2 of (i)/(ii)/(iii) pass
  with the permutation null respected, in ≥1 genuine cohort. **Otherwise** report honestly:
  "transcriptomic replication not established (GSE199243 signed-null; [new cohort] [result])" and keep
  the protein-level external validation as the load-bearing external evidence.

**Output:** `results/de/de_<COHORT>_per_celltype.csv`, `results/validation/replication_<COHORT>.csv`
(metrics + null p-values), a scatter panel. **Script:** `scripts/03_de/05_external_replication.py`.

---

### A2 — Genome-wide differential accessibility + secretome-linked-peak matched-null enrichment  (codex #2)

**Claim it would license (if PASS):** "ATAC supports the RNA signature beyond proximity direction-
concordance — secretome-linked peaks are enriched among donor-level DA peaks in the RNA-concordant
direction." (If FAIL, ATAC stays SUPPORTING, as now.) Also completes the **accessibility-matched null
that PHASE_REPORT explicitly deferred** (line 41).

**Data:** full peak matrix `processed/per_dataset/GSE268609_peaks_counts.mtx` (8.9 G, on disk) +
`nCount_Peaks` depth. Heavy — budget memory/time; consider per-celltype chunking.

**Method (pre-specify):**
- Per niche celltype, **donor pseudobulk over ALL peaks** (not just linked). Design `~group` (YA/HA);
  **depth-normalise by donor total fragments-in-peaks** (`nCount_Peaks`). Optional sensitivity:
  add a depth/QC covariate (`~depth + group`) and report stability.
- Test every peak (Mann-Whitney or pyDESeq2 on the peak pseudobulk), BH-FDR genome-wide → `DA peaks`.
- **Build a matched background** for the secretome-linked peak set: match on **GC content, peak width,
  mean accessibility, and distance-to-TSS** (now feasible — the full matrix gives mean accessibility,
  the field-standard 4th covariate the prior matched-null lacked).
- **Set-enrichment test:** are secretome-linked peaks shifted toward DA (and toward *opening* near
  RNA-UP genes / *closing* near RNA-DOWN genes) vs the matched background? Use a rank/AUC or
  2-sample test on the DA statistic, with the matched-null.

**Pre-specified GO/NO-GO:**
- PASS ("strong-ish ATAC") if: secretome-linked peaks are significantly shifted vs matched null
  (p<0.05) **AND** the shift is RNA-direction-concordant (UP→opening, DOWN→closing) **AND** ≥3 peaks
  reach DA FDR<0.1 with RNA-concordant direction.
- Otherwise: report the set-level result honestly; ATAC remains a SUPPORTING direction-concordance
  (the genome-wide per-peak family is expected to be underpowered at n=8/9 — say so).

**Output:** `results/de/atac_DA_peaks_<celltype>.csv`, `results/de/atac_linkedpeak_enrichment.csv`.
**Script:** `scripts/04_atac/05_differential_accessibility.py`. **Env:** bio (or sgz_r/ArchR if you
prefer its DA machinery — then document the tool).

---

## TIER B — do if Tier A is encouraging / moderate value

### B1 — Peak-gene link hardening: published tool + matched null + FDR  (codex #3)

**Why:** `scripts/04_atac/04_coaccess_links.py` already does within-celltype metacell peak↔RNA
Spearman, but with a fixed ρ>0.1/p<0.01 threshold and **no matched random peak-gene background and no
FDR**. Redo with a standard method to make the "proximity-only" caveat either go away or be credibly
confirmed.

**Method:** ArchR `addPeak2GeneLinks` (or Signac `LinkPeaks`, or Cicero) **within each niche
celltype** (NEVER pool across celltypes — this was the audit-fixed artifact). Donor/metacell-aware.
Compare validated-link rate to a **matched random peak-gene background**; control FDR.

**GO/NO-GO:** PASS ("computationally linked") if the secretome peak-gene links validate at FDR
significantly above the matched-random rate, within-celltype, for ≥k genes incl. ≥1 Fig-3 hit.
Otherwise links stay **proximity-based** (current wording). **Script:**
`scripts/04_atac/06_peak2gene_links.R`.

### B2 — Subregion (DG vs non-DG) specificity via Franjic dissections  (codex #5, redirected)

**Why redirect:** the GSE264692 Visium is anterior-HPC and **does not crisply resolve DG/SGZ** (stated
in the manuscript). A spot-level `region×age` mixed model would rest on unreliable region labels.
**Franjic has explicit subregion dissections** (DG, CA1, CA2-4, EC, SUB) → a cleaner substrate.

**Method:** score the frozen niche aging-UP signature per Franjic cell/pseudobulk; test
`signature ~ subregion (+ age where range allows) + celltype_composition + (1|donor)`. Question: is the
signature **higher in DG than other subregions**? **Honest ceiling:** with Franjic's compressed age
range and modest n, the defensible outcomes are "DG-enriched" vs "not DG-restricted / pan-HPC" — same
honesty bar as the current Visium conclusion. Keep the Visium analysis as-is for spatial context.
**Script:** `scripts/06_validation/03_franjic_subregion.py`. (If #6 deconvolution is ever pursued, it
attaches here, not to the Visium.)

---

## TIER C — do NOT (re-)run for a positive (rationale, so a future session doesn't waste effort)

- **#4 TF triangulation — ALREADY COMPLETE & NULL.** Layers A–D in `PHASE_REPORT.md` already align
  TF-RNA / chromVAR / footprint / motif-enrichment / target-RNA / peak-link. Result: FDR-null, nominal
  *against*. Re-triangulating null signals cannot manufacture a positive. Only revisit with **new data**
  (e.g. matched bulk-ATAC depth, or CUT&RUN/ChIP for NF-κB/IRF) — which is wet-lab, not in-silico.
- **#7 L-R fixed sender/receiver — ALREADY COMPLETE & NULL.** Layer D already required ligand
  aging-UP+expressed, receptor expressed+age-regulated, target coherence, cross-cohort robustness, with
  a permutation null — exactly codex's spec. CCN2/SPON1/IL15/APOD all fail. The autocrine null stands.
- **#6 Visium deconvolution (RCTD/cell2location)** — only meaningful if the Visium could resolve DG,
  which it can't; B2 (Franjic) is the better route. Defer unless a higher-resolution spatial dataset
  is obtained.
- **#8 mediation-style analysis** — with n=8/9 donors this is statistically hopeless and will mislead.
  At most: an *exploratory candidate-ranking* heuristic, clearly labelled, **never** used for any
  claim of "mechanism-like" causation.

---

## Global rules (inherited; apply to every analysis above)

1. **Donor is the unit.** Pseudobulk; no cell-level p-values for inference.
2. **Pre-specify** the hypothesis, metric, and GO/NO-GO **before** looking at the result; log it.
3. **Matched nulls + FDR**, not raw p or vs-0.5 nulls. Test enrichment **above the empirical
   background**, not against an implausible 50%.
4. **Within-celltype** for any peak↔RNA covariation; never pool across celltypes.
5. **Universe-aware:** report secretome-universe AND genome-wide; lead with the within-universe /
   directional statistic.
6. **Wording:** "externally validated", "candidate/attempted replication", "supporting" — reserve
   "replicated" / "computationally linked" / "strong ATAC evidence" for analyses that PASS their
   pre-set criteria.
7. **Report nulls as results.** A clean external-replication null or a clean DA null is publishable
   and consistent with the paper's honest-dual-null identity.
8. **Provenance:** primary publication via `!Series_pubmed_id` for any new cohort; seed 42; new code in
   `scripts/`; write `audit_log/<date>_in_silico_hardening/RESOLUTION.md`.

## Minimum bars to actually upgrade a claim (adopted from codex, sharpened)

- **Multi-cohort replication (A1):** signature-level significant in an independent cohort, logFC
  correlation positive, primary-hit sign-concordance above the cohort's own background, permutation
  null respected.
- **Strong ATAC (A2):** secretome-linked peaks enriched vs a GC+width+accessibility+distance-matched
  null, RNA-direction-concordant, with ≥3 peaks at FDR.
- **Computationally linked (B1):** within-celltype peak-gene links validate at FDR above a
  matched-random background.
- **SGZ/DG specificity (B2):** subregion×age effect with the signature higher in DG than other
  subregions (Franjic); else state "not DG-restricted".
- **Mechanism:** would require TF-activity + motif + footprint + peak-gene link + target-RNA all
  concordant **in the same celltype** — already shown NOT to hold; needs orthogonal/wet-lab data.

## Suggested order
A1 (external replication) and A2 (DA + matched-null enrichment) in parallel → if A2 encouraging, B1
(peak-gene) → B2 (Franjic subregion) alongside A1 since it reuses the same Franjic load. Skip Tier C.
