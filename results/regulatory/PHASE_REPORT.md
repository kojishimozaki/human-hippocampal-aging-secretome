# Phase 2 (regulatory) — results vs the locked pre-specified plan

Plan: `PHASE1_5_PRE_REGISTRATION.md` (locked 2026-06-04, before running). Honest verdict first.

## Layer A — cis TF-motif enrichment in aging-secretome peaks: **NULL at FDR**
- Evaluable arms (≥4 sig secretome genes): Astro-UP (18 genes), Astro-DOWN (8), Micro-UP (10),
  Endo-DOWN (8), OPC-DOWN (4). Others NOT_EVALUABLE.
- **0 motif×arm passes enrichment FDR<0.10. 0 CONVERGENT.** The pre-specified hypothesis
  (UP~SASP-TFs, DOWN~homeostatic-TFs) is **not statistically established** in this cohort.
- Exploratory (nominal raw p<0.05, NOT FDR-significant): a **direction-coherent SASP-TF signal in
  the astro+micro UP arm** — IRF1, CEBPA, STAT1::STAT2 (Astro), STAT3, NFKB1, IRF1 (Micro) — 6/6
  prediction-matched AND chromVAR-dz-sign-concordant. Best FDR=0.16 (Astro-DOWN STAT3, raw p=0.003,
  but direction-discordant). Underpowered: foreground peak sets are only 28–88 peaks.
- chromVAR donor-level (n=8/9) remains FDR-null (unchanged) — used here only as a direction axis.

## Layer B — NicheNet ligand activity (prior-based, geneset-permutation null): **signal in Endo+Astro**
- Controlling hub-ligand bias by permuting the response geneset (2000×): in **Endothelial** and
  **Astrocyte** receivers the permutation null is **non-discriminative**: **~82–91% of ALL 1,226
  ligands** pass perm-FDR<0.1 (Astro 1006/1226, Endo 1117/1226; Micro/Oligo/OPC 0/1226). NicheNet
  therefore does **not** select specific drivers — the aging response is a large co-regulated set most
  ligands' target profiles correlate with. CCN2/CTGF and SPON1 are merely **top-ranked** within that
  mass, not specific hits (the `z`-vs-all-ligand column was additionally degenerate/NA).
- This is a **prior/co-regulation signal, not evidence of signalling**; the real adjudication is the
  ligand–receptor validation (Layer D). Hypothesis-generating only.

## Synthesis — what is and isn't established
- **Established (FDR): nothing new at the cis-regulatory level.** The chromatin/TF mechanism of the
  aging niche secretome is **not resolvable in this single multiome cohort at this power.**
- **Coherent, hypothesis-generating, underpowered:** a two-arm structure —
  (i) an **astro/micro SASP-cytokine arm** with a nominal NF-κB/C-EBP/STAT/IRF cis-signature
  (motif raw p<0.05, chromVAR-concordant) and IL15/CXCL2 ligands; and
  (ii) an **astro-ECM / endothelial matricellular arm** where CCN2/CTGF (Astro) and SPON1 (Endo)
  rank top in a **non-discriminative** NicheNet permutation (most ligands pass; see Layer B/D).
- This is consistent with — but does not prove — the SASP-program hypothesis. It is the honest
  ceiling of the no-new-data analysis.

## Deviations from the pre-specified plan (logged)
- NicheNet receiver (ii) NSC/Neuroblast **dropped**: the niche DE has no NSC/Neuroblast contrast
  (rare cells, underpowered, not load-bearing by design). Receivers restricted to niche cells.
- Matched-null used **GC+width** (not GC+accessibility): accessibility-matching needs the 8.9 G peak
  matrix; deferred. GC is the dominant motif-frequency confound (field-standard control).

## Layer C — donor-level TF footprinting (Stage 2): **no evidence for increased SASP-TF footprints (null; nominal trend opposite)**
Donor-level (n=8 YA / 9 HA) differential footprint at all genomic motif sites within accessible peaks,
from 277 M / 81 M / 25 M Tn5 insertions (Astro/Micro/Endo). **Caveats:** cut sites are *unshifted* (no
+4/-5 Tn5 offset), there is *no Tn5-bias model*, and the score is a simple *flank/center ratio* — but
the YA-vs-HA *same-site* differential cancels both the constant offset and the sequence bias, and the
ratio is depth-normalised within donor, so the differential is interpretable (donor = unit, no pooling
pseudoreplication). It is higher-powered than the 28–88-peak motif test (millions of insertions
aggregated) but still **donor-limited (n=8/9)**.
- **0 / 21 TF×celltype tests pass footprint FDR<0.1.** 5 reach nominal p<0.05 and **all five go the
  *opposite* way to the SASP hypothesis** (footprints SHALLOWER with age): Endo IRF1 (dz=−0.066,
  p=0.011), Astro RELA (−0.025, p=0.015), Astro TEAD4 (−0.028, p=0.021), Astro TEAD1 (−0.026, p=0.036),
  Astro IRF1 (−0.052, p=0.046). Accessibility around these motifs is flat (n.s.), so not a depth artifact.
- **Read as a clean negative: NO evidence for increased SASP-TF (NF-κB/IRF/STAT/C-EBP) or TEAD chromatin
  binding with age, and the nominal trend is the opposite. We do NOT claim a positive "refutation" — a
  null + nominal-reverse at n=8/9 is "not supported", not proof of absence.**

## Layer D — autocrine L-R validation of the NicheNet leads: **does NOT hold up**
Requiring more than prior-based activity — ligand actually aging-UP + expressed, receptor expressed AND
age-regulated in the receiver, coherent targets, ligand robust in an independent cohort (GSE199243):
- **CCN2/CTGF (Astro):** ligand aging-UP (lfc=+2.3, padj=0.013) but very low expression (baseMean=8),
  not testable cross-cohort; its 8 receptors are expressed but NONE age-regulated (all padj>0.4); 8/250
  top targets aging-UP (incl. CCN2 itself = circular). WEAK.
- **SPON1 (Endo):** ligand is NOT aging-UP in Endo (lfc=−0.2, padj=1.0); receptor APP nominally DOWN
  (padj=0.08). FAILS — the NicheNet activity was target co-regulation, not a real L-R axis.
- **IL15 (Micro):** ligand strongly aging-UP (lfc=+3.1, padj<1e-3; a real microglial SASP cytokine) but
  NOT replicated cross-cohort (GSE199243 lfc=+0.04), receptor IL15RA barely expressed (baseMean=9,
  static), 1 target → the hit is real but the autocrine LOOP is unsupported.
- **APOD (Astro):** nominal up (padj=0.06), flat cross-cohort, receptor LEPR static.
=> The NicheNet permutation signal was **prior/co-regulation-inflated** (the pre-specified caveat).
   No robust autocrine program survives L-R validation.

## FINAL SYNTHESIS (the truth, across all layers)
**Both** mechanism hypotheses fail rigorous validation: cis (chromVAR/motif/footprint — FDR-null,
footprint nominally anti-SASP-direction) AND trans (NicheNet L-R — prior-inflated, fails L-R checks).
**The aging niche secretome is a real, externally-validated RNA-level OUTPUT signature
(SASP/proteome/CSF-concordant) whose mechanism — cis-TF activation OR autocrine signalling — is NOT
resolved in this single multiome cohort.** Individual hits are real and strong (microglial IL15
lfc=+3.1; astrocyte ECM TNC/CCN2/COL), but as a *program* the regulator is unestablished: the SASP-TF
cis-activation route is footprint-**unsupported** (FDR-null, nominal trend against — not a positive
refutation), and the autocrine route fails L-R validation.

**Headline implication:** reframe around the **solid core** — the externally-validated, pan-hippocampal,
SASP-aligned niche-secretome OUTPUT signature — and present BOTH mechanism hypotheses (cis-TF, autocrine)
as **honest negatives** ("no evidence for", not "disproven"), with mechanism stated as the open question
for future perturbation work. The footprint null and L-R null are themselves clean results.
This is the rigorous, defensible, genuinely-informative position. (NB: the analysis plan was
*pre-specified / locked before running* this session but is not independently git-timestamped ahead of
the results — all committed together — so "pre-specified", not "pre-registered", is the claim made.)
