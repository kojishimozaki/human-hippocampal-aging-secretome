#!/usr/bin/env python3
"""Single source of truth for the bioRxiv v2 revision.

The submitted manuscript (`submission/shimozaki-aging-secretome-final_draft.docx`,
identical in body text to the submitted PDF) is the base. This module declares every
change as an operation on that document's ordered paragraph list, so the DOCX, the
English markdown of record and the Japanese mirror are all generated from one
definition rather than edited independently (CLAUDE.md iron rule 9).

Provenance for every number introduced here is recorded in
`audit_log/2026-09-04_biorxiv_v2_reframe/RESOLUTION.md`.
"""

from __future__ import annotations

# --------------------------------------------------------------------------------------
# Anchors: the first 60 characters of each base paragraph an operation targets. Verified
# against the base document at build time so the plan fails loudly if the base changes.
# --------------------------------------------------------------------------------------

TITLE = (
    "A support-cell secretome program in human hippocampal aging is directionally "
    "reproducible but not gene-specific"
)

ABSTRACT = (
    "Cells that build the hippocampal extracellular environment are a plausible but "
    "poorly specified locus of human brain aging, and the single-nucleus designs used to "
    "study them are rarely calibrated before they are interpreted. We applied a "
    "pre-specified analysis of 2,224 secretome-related genes to a human hippocampal "
    "multiome (8 young and 9 older, clinically defined control donors), calibrating the "
    "contrast by donor-label permutation before reading any biological result and pairing "
    "every claim with a matched null designed to break it. Three findings survived. Older "
    "age was associated with a directionally stable secretome-related program across "
    "astrocytes, microglia, a mixed vascular compartment, oligodendrocyte precursor cells, "
    "and oligodendrocytes: all 60 gene–cell-type pairs kept their sign in every "
    "leave-one-donor-out refit, all 57 testable directions were preserved after "
    "ambient-RNA correction, and the direction reproduced in an independent hippocampal "
    "cohort (Spearman ρ = 0.24; 81% sign concordance versus a 58% background, "
    "P = 5×10−4). The same program was not reproduced above background in mature "
    "dentate granule, CA excitatory, or inhibitory neurons from the same donors, "
    "identifying it as support-cell-weighted rather than tissue-wide. RNA and chromatin "
    "accessibility were directionally coupled at the gene level in astrocytes, surviving "
    "an end-to-end donor-label permutation that repeated gene selection itself "
    "(P = 0.025). Four claims did not survive their controls. Expression- and "
    "effect-matched control genes reproduced as well as the selected genes (P = 0.146 and "
    "P = 0.187), so gene-level specificity is not supported; an apparent RNA–methylation "
    "coordination failed all four specificity controls; no upstream regulator was "
    "identifiable by motif, footprint, or transcription-factor RNA analysis; and no "
    "candidate receiver compartment showed coordinated receptor up-regulation. Donor-label "
    "permutation also showed that under a complete null this 8-versus-9 design returns at "
    "least one gene significant at FDR < 0.05 in 29–42% of draws, against a bound of "
    "0.05, so the per-gene P values are anti-conservative here and hit counts cannot be "
    "read as false-discovery estimates. Three "
    "depositor-assigned rare cell labels were not supported by marker or classifier "
    "re-analysis and were excluded from every conclusion. The result is a directionally "
    "reproducible, mechanistically unresolved support-cell aging program, reported with an "
    "explicit account of which inferences human multiome data at this scale can and cannot "
    "sustain."
)

KEYWORDS = (
    "Keywords: hippocampal aging; niche-supporting cells; secretome-related transcription; "
    "extracellular environment; single-nucleus multiome; cognitive aging; matched-null "
    "controls; analytical reproducibility"
)

# --------------------------------------------------------------------------------------
# Methods: new "Reproducibility" subsection, inserted after "Statistical procedures".
# --------------------------------------------------------------------------------------

METHODS_REPRO_HEADING = "Reproducibility"

METHODS_REPRO_1 = (
    "Analysis code, the pre-specified analysis documents, the software environments, and "
    "the result tables underlying every number reported here are available in the project "
    "repository (Code availability). Environments are pinned in envs/, and analyses use "
    "random seed 42 unless stated otherwise. Workflow targets are separated into those "
    "that run from the committed result tables alone — the figures, the validation tests, "
    "the supplementary tables, and the multiple-testing-scope check — and those that "
    "additionally require the deposited matrices, so that a reader without the primary data "
    "can still regenerate every figure and validation table from a clone of the repository."
)

METHODS_REPRO_2 = (
    "Three limits on that reproducibility are stated rather than implied. First, the script "
    "that exported count matrices from the deposited Seurat object was not retained, so the "
    "ingestion step — cell filtering, quality thresholds, and RNA–ATAC barcode matching — "
    "cannot be re-executed from the repository; the workflows can be repeated from the "
    "processed-matrix stage onward. Second, re-running the primary differential-expression "
    "target in the pinned environment returns 99,543 rows rather than the 99,544 of the "
    "table used here, because the low-count oligodendrocyte transcript LINC01238 (two "
    "counts in total, detected in 2 of 17 donors) falls below the current gene filter. All "
    "60 predefined pairs are present in the regenerated table and all 60 retain adjusted "
    "P < 0.1; the largest absolute difference in log2 fold change across those 60 is "
    "4.3×10−5. Table S1 reports the regenerated per-cell-type test counts and its "
    "oligodendrocyte entry is accordingly 22,844 rather than 22,845, so the five entries sum "
    "to 99,543. Third, the scripts that produced the donor-label permutation calibration of "
    "Table S1, the effect-size-matched reproduction background, and the external donor-label "
    "null were run in the pre-submission audit environment and are not part of the "
    "repository; their output tables and a provenance record for each are included, and "
    "every manuscript number drawn from them was recomputed from those tables."
)

# --------------------------------------------------------------------------------------
# Results: new opening subsection. Replaces the calibration paragraph that currently sits
# third in "Older age reshapes secretome-related transcription...", which is deleted.
# --------------------------------------------------------------------------------------

RESULTS_CALIB_HEADING = "Calibrating what this design can support"

RESULTS_CALIB_1 = (
    "Deposited single-nucleus data now supply much of what is claimed about human brain "
    "aging, yet the inferential capacity of the resulting designs is rarely measured. We "
    "calibrated the primary contrast — 8 young adults versus 9 older, clinically "
    "defined control donors — by donor-label permutation, and report that calibration here, "
    "after the results it governs, because it sets the terms on which every count above "
    "should be read."
)

RESULTS_CALIB_2 = (
    "Donor labels were permuted 100 times within each cohort×cell-type combination and the "
    "complete donor-pseudobulk model was refitted (Table S1). Two properties of the design "
    "emerged. The first is favourable. The observed number of secretome-related associations "
    "exceeded the permutation null in microglia (empirical P = 0.020), astrocytes "
    "(P = 0.040), the vascular compartment (P = 0.040), and OPCs (P = 0.050), and the null "
    "median was zero in every compartment; oligodendrocytes did not exceed their null (five "
    "pairs, P = 0.109). The second is not. Under a complete null, Benjamini–Hochberg "
    "control at level α bounds the probability that any gene is called significant by α, "
    "so at the genome-wide FDR < 0.05 scope that probability should not exceed 0.05. The "
    "permutations returned at least one significant gene in 29–42% of draws and at least "
    "ten in 8–16%; GSE278576’s 10-versus-20 split gave 29–51% and 3–17%. Exceeding the "
    "bound by that margin is not a property of the procedure but evidence that the "
    "per-gene P values are anti-conservative at this donor count, so the "
    "Benjamini–Hochberg guarantee does not hold here. Counts of significant genes are "
    "accordingly not an interpretable currency, and no expected number of false "
    "associations can be read off them."
)

RESULTS_CALIB_3 = (
    "Applying the same calibration to the genome-wide count at FDR < 0.05 separates the "
    "strong arms from the weak. Astrocytes (P = 0.030), microglia (P = 0.030), and OPCs "
    "(P = 0.050) exceeded their nulls, whereas oligodendrocytes (P = 0.139) and the vascular "
    "compartment (P = 0.069) fell inside the null’s 95th percentile. We therefore treat "
    "astrocytes, microglia, and OPCs as the compartments that carry demonstrable signal at "
    "both scopes, and report the oligodendrocyte and vascular arms as directionally "
    "consistent but not separable from their permutation nulls genome-wide. Empirical P "
    "values cannot fall below 1/101 = 0.0099, and 100 draws do not exhaust the 24,310 "
    "distinct label assignments available (11,440 for the microglial and vascular "
    "comparisons, which are 7 versus 9 after the 10-nucleus threshold), so these values "
    "bound the evidence rather than estimate it precisely."
)

RESULTS_CALIB_4 = (
    "Multiple-testing scope was fixed in advance and is reported in both forms. "
    "Benjamini–Hochberg correction was applied within each cell type, so the 99,544 tests "
    "across the five niche classes form five families rather than one. Treating all five as "
    "a single genome-wide family retains 49 of the 60 reported pairs and admits 4 new ones, "
    "all in oligodendrocytes (AZGP1, LTBP3, NRG1, and NXPE3). Because the predefined "
    "60-pair set seeded every downstream accessibility, methylation, reproduction, receptor, "
    "and neuronal analysis, we retain the original set and report the alternative scope "
    "rather than silently replacing it."
)

RESULTS_CALIB_4B = (
    "The gene universe was defined as a union, and its stricter intersection is reported "
    "beside it. The 2,224-gene universe unions the Human Protein Atlas predicted-secreted "
    "annotation with reviewed UniProtKB secreted entries; 1,813 genes carry both. Ten of "
    "the 60 pairs rest on one annotation only — USH2A, CHST9, CRB2, KITLG and ANTXR2 in "
    "astrocytes, KDR and KITLG in the vascular compartment, ASAH1 in microglia, FGF13 in "
    "OPCs and EPHA3 in oligodendrocytes — and FGF13 is reported as intracellular in "
    "experimental work while KDR is a receptor. Restricting to the 50 pairs annotated by "
    "both sources leaves external direction concordance at 32 of 40 evaluable pairs (80%, "
    "against 39 of 48, or 81%, for the full set) and astrocyte gene-level RNA–accessibility "
    "concordance at 16 of 21 genes (against 20 of 26). Removing FGF13 alone cuts the OPC "
    "linked-peak count from 19 to 6 while leaving its remaining genes concordant. The "
    "reported results therefore do not depend on the more permissive annotations, although "
    "the OPC accessibility evidence is thin once its dominant gene is set aside."
)

RESULTS_CALIB_5 = (
    "Two consequences follow, and we hold to them throughout. Direction, not per-gene "
    "significance, is the unit of evidence this design can support; and no claim in this "
    "paper rests on a count of significant genes."
)

# --------------------------------------------------------------------------------------
# Discussion: softened opening heading, and the convergence subsection moved to the front.
# --------------------------------------------------------------------------------------

DISCUSSION_HEADING_SOFTENED = (
    "A directionally stable but mechanistically unresolved support-cell program"
)

CONVERGENCE_BRIDGE = (
    "We consider first what the external comparisons do and do not license, because that "
    "boundary governs how the rest of the program should be read. "
)

# --------------------------------------------------------------------------------------
# Conclusions: full replacement of the two-paragraph "Conclusions and future directions".
# --------------------------------------------------------------------------------------

CONCLUSIONS_HEADING = "Conclusions: what this design can and cannot sustain"

CONCLUSIONS_1 = (
    "Three findings survived every control we applied. Human hippocampal aging is "
    "accompanied by a directionally stable secretome-related program in support-cell "
    "compartments: all 60 gene–cell-type pairs retained their sign in every "
    "leave-one-donor-out refit, 57 of 57 testable directions survived ambient-RNA "
    "correction, and the direction reproduced in an independent hippocampal cohort "
    "(Spearman ρ = 0.24; 81% sign concordance against a 58% background, P = 5×10−4). The "
    "program is support-cell-weighted rather than tissue-wide: mature dentate granule, CA "
    "excitatory, and inhibitory neurons from the same donors did not reproduce it above "
    "their own secretome backgrounds and yielded 4, 3, and 0 de novo associations against "
    "17, 11, 7, 6, and 6 in the five support-cell compartments. And in astrocytes, RNA and "
    "chromatin accessibility were directionally coupled at the gene level — 20 of 26 genes "
    "across 160 linked peaks — surviving a gene-matched null, a sign-permutation null, "
    "leave-one-gene-out analysis, and an end-to-end donor-label permutation that repeated "
    "gene selection itself (P = 0.025)."
)

CONCLUSIONS_2 = (
    "Four claims did not survive. Genes matched for expression and primary-cohort effect "
    "size reproduced as well as the selected pairs (P = 0.146 and P = 0.187 at tight "
    "matching), so the external agreement supports a broader secretome-related pattern and "
    "not the specificity of these 60 pairs. An apparent RNA–methylation coordination failed "
    "all four specificity controls, including orientation by a non-matching cell type and a "
    "matched diagonal indistinguishable from 10,000 scrambled grids. No upstream regulator "
    "was identifiable: 0 of 300 motif tests and 0 of 21 footprint tests reached FDR "
    "significance, and cognate transcription-factor RNA followed its motif in 11 of 22 "
    "testable combinations against a 52.4% matched background. No candidate receiver "
    "compartment showed coordinated receptor up-regulation (8 of 8 confirmatory tests "
    "failed; directional Q = 0.986). Joint adjustment for sex, preparation arm, post-mortem "
    "interval, and group retained 0 of 60 pairs at adjusted P < 0.1, and the three-term "
    "model without sex retained 12 of 60, while all 60 directions remained unchanged — a "
    "result that neither establishes nor excludes confounding, but that marks precisely "
    "where per-gene inference stops."
)

CONCLUSIONS_3 = (
    "Two further results concern the data rather than the biology. Donor-label permutation "
    "showed that under a complete null this design returns at least one gene significant at "
    "FDR < 0.05 in 29–42% of draws, against a bound of 0.05: the per-gene P values are "
    "anti-conservative at this donor count, so hit counts at this scale cannot be read as "
    "false-discovery estimates. And three depositor-assigned rare cell labels were not "
    "supported by marker-panel or classifier re-analysis: 0.9% of NSC-labelled nuclei "
    "retained their label, against 83–99.8% retention for established cell types. We "
    "excluded those populations from every conclusion, and note that any reanalysis "
    "inheriting such labels uncritically will inherit the error."
)

CONCLUSIONS_4 = (
    "What remains is a directionally reproducible, mechanistically unresolved support-cell "
    "aging program, delivered together with an explicit account of its evidentiary "
    "boundary. That boundary defines the experiments that would move it: direct protein "
    "measurement in the compartments nominated here, cell-type-resolved profiling at depth "
    "sufficient to model enhancer–gene relationships within a compartment, a formal "
    "cell-type-by-age interaction test in a cohort powered for it, and perturbation with "
    "single-cell readouts in human multicellular models. We report the negative results at "
    "the same length as the positive ones because, at this sample size, knowing which "
    "inferences fail is the more transferable result."
)

# --------------------------------------------------------------------------------------
# Availability statements.
# --------------------------------------------------------------------------------------

DATA_AVAIL_ADDENDUM = (
    " Processed intermediates are not redistributed: the analysis matrices exceed 50 GB and "
    "are regenerable from the deposits above by the documented workflow, subject to the "
    "ingestion limit stated in Methods. The code archive contains the pre-specified analysis "
    "documents and the result tables underlying the reported numbers, including the "
    "donor-label permutation calibration, the effect-size-matched reproduction background, "
    "and the external donor-label permutation null, each with a provenance record."
)

CODE_AVAIL_HEADING = "Code availability"

CODE_AVAIL_BODY = (
    "The analysis code, the pre-specified analysis documents, and the result tables that "
    "support the reported numbers are available at {REPO_URL} and archived at {ZENODO_DOI}; "
    "the version used for this manuscript is the release tagged {RELEASE_TAG}. Software "
    "environments are "
    "pinned in envs/ (pyDESeq2 0.5.4, scanpy 1.11.5, Python 3.11.15; R with DESeq2, Signac, "
    "chromVAR, and ashr), and analyses use random seed 42 unless stated otherwise. The "
    "figure, validation, supplementary-table, and multiple-testing-scope targets run from a "
    "clone of the repository using only the committed result tables; the targets that "
    "recompute chromatin-accessibility concordance, differential-expression robustness, the "
    "regulatory analyses, and the primary differential-expression table additionally require "
    "the deposited matrices, which are not redistributed here. The limits of that "
    "reproducibility are stated in Methods."
)

# Filled in by the author at submission; the build refuses to emit a final document while
# any of these is still the placeholder unless --allow-placeholders is passed.
# A commit hash cannot be quoted inside the commit that contains it, so the manuscript
# names a git tag instead: create the tag on the finished commit and the reference resolves
# to exactly one revision with no circularity. Only the two bracketed values need filling.
PLACEHOLDERS = {
    "REPO_URL": "[repository URL to be inserted at submission]",
    "ZENODO_DOI": "[Zenodo DOI to be inserted at submission]",
    "RELEASE_TAG": "biorxiv-v2",
}

# --------------------------------------------------------------------------------------
# Rule-13 corrections. The Word rewrite replaced four measured values with approximations
# ("approximately zero", "approximately 25,000 ...", "approximately 82-91%",
# "approximately 22.0") and softened a fifth with "almost exactly". The exact values are
# restored here. Sources: signed reproducibility 0.02 = replication_GSE199243v2 (also
# Fig. S7 legend, rho = +0.02); nucleus counts recounted from the deposited metadata
# (Cluster column of processed/per_dataset/GSE268609_metadata.csv: Astrocytes 25,218,
# Microglia 12,165, mOli 59,100, OPCs 12,597 -- note the draft's 12,600 for OPCs was
# already one step removed from the data); 82-91% and 21.99 are unchanged values that
# merely lost their exactness in the wording.
# --------------------------------------------------------------------------------------

EXACT_VALUE_FIXES = [
    ("because its signed reproducibility was approximately zero [20].",
     "because its signed reproducibility was 0.02 [20]."),
    ("at depth\u2014approximately 25,000 astrocyte, 12,000 microglial, 59,000 oligodendrocyte, "
     "and 12,600 oligodendrocyte precursor cell (OPC) nuclei",
     "at depth\u201425,218 astrocyte, 12,165 microglial, 59,100 oligodendrocyte, "
     "and 12,597 oligodendrocyte precursor cell (OPC) nuclei"),
    ("null, approximately 82\u201391% of all ligands passed",
     "null, 82\u201391% of all ligands passed"),
    ("This was almost exactly the depth- and expression-matched expectation of 21.99",
     "This matched the depth- and expression-matched expectation of 21.99"),
    ("expectation of 21.99 (approximately 22.0; P = 0.68)",
     "expectation of 21.99 (P = 0.68)"),
    ("showed a signed-magnitude reproducibility near zero",
     "showed a signed-magnitude reproducibility of 0.02"),
]

# --------------------------------------------------------------------------------------
# Corrections arising from external editorial review (2026-09-07), each verified against
# source data before being applied. The review's remaining points are recorded in
# audit_log/2026-09-07_editorial_review/RESOLUTION.md, including the ones not accepted.
# --------------------------------------------------------------------------------------

REVIEW_FIXES = [
    # The CSF comparison is made against the whole measured proteome, not against a
    # secretome-matched background. The manuscript criticises unmatched backgrounds
    # elsewhere, so both comparisons are now reported.
    ("Agreement with the CSF aging proteome [16] was weaker and depended on the threshold. "
     "TNC, ANXA1, SERPINE1, and IGFBP5 all increased with age in CSF, but the predefined "
     "signature showed only a non-significant trend (60% up, n = 25, P = 0.21). The extended "
     "set was concordant (69% up, n = 42, P = 0.001; Fig. S4).",
     "Agreement with the CSF aging proteome [16] was weaker and depended on both the "
     "threshold and the comparison background. TNC, ANXA1, SERPINE1, and IGFBP5 all "
     "increased with age in CSF. Against the whole measured CSF proteome (6,174 features) "
     "the predefined signature showed only a non-significant trend (60% up, n = 25, "
     "P = 0.21) while the extended set was concordant (69% up, n = 42, P = 0.001; Fig. S4). "
     "Restricting the background to the 1,325 CSF features that lie inside the secretome "
     "universe — the matched comparison — moves those values to P = 0.44 and P = 0.014."),

    # Nuclear RNA is not secretion; the figure title claimed the latter.
    ("Fig. 2 | Aging reprograms the niche\u2019s secreted-protein output.",
     "Fig. 2 | Age-associated secretome-related transcription across niche cell classes."),

    ("CSF age estimates [16] are compared among aging-UP, aging-DOWN, and background genes "
     "within the niche-secretome universe. The predefined adjusted-P < 0.1 set shows a "
     "non-significant upward trend (60% up, n = 25, P = 0.21). The extended adjusted-P < 0.2 "
     "set is concordant (69% up, n = 42, P = 0.001). The result supports the extended UP "
     "component, not the core 60-pair set.",
     "The aging-UP and aging-DOWN sets are drawn from the niche-secretome universe and their "
     "CSF age estimates [16] are compared with the whole measured CSF proteome as background "
     "(6,174 features). The predefined adjusted-P < 0.1 set shows a non-significant upward "
     "trend (60% up, n = 25, P = 0.21) and the extended adjusted-P < 0.2 set is concordant "
     "(69% up, n = 42, P = 0.001). Against the 1,325 background features that lie inside the "
     "secretome universe the same comparisons give P = 0.44 and P = 0.014. The result is "
     "therefore background-dependent, and supports the extended UP component rather than the "
     "core 60-pair set."),

    # Fig. 4A's own panel carries the concordance numbers and the scope of the one
    # significant gene; the legend said only what the panel showed. Verified against
    # results/de/atac_rna_concordance.csv: 41/60 concordant, background 2,269 of the
    # 3,942 pairs with an RNA adjusted P (57.6%), COL21A1 significant only at
    # padj_atac_rnasig = 0.0099 while its genome-wide padj_atac is 0.164.
    ("A, Directional agreement between RNA and linked-peak ATAC effects. B, Age-associated "
     "chromVAR motif activity.",
     "A, Directional agreement between RNA and linked-peak ATAC effects: 41 of the 60 pairs "
     "(68%) are concordant, against a background of 57.6% among the 3,942 RNA-non-significant "
     "pairs. COL21A1 is the only gene to reach ATAC FDR significance, and it does so within "
     "the RNA-significant family alone: none of the 369,393 genome-wide peak tests survives "
     "correction. B, Age-associated chromVAR motif activity."),

    # "independent" invited the reading that a separate investigator had audited the work.
    ("MAD is included only in C to show age matching and is not part of the resilience "
     "contrast.",
     "MAD is included only in C to show the age distribution and is not part of the "
     "resilience contrast."),

    ("Prior to submission, Claude (Anthropic) was also used for an independent AI-assisted "
     "robustness audit and for manuscript proofreading.",
     "Prior to submission, separate AI sessions were used for adversarial robustness "
     "re-checks (Claude, Anthropic) and for editorial review and proofreading (Claude, "
     "Anthropic; GPT, OpenAI). These were additional automated passes over the author\u2019s "
     "own work, not review by an independent investigator, and their findings were verified "
     "against source data by the author before any were acted on."),
]

TABLE_S1_ADDENDUM = (
    " These calibration results are summarised in the first Results subsection, which "
    "governs how the counts reported in later sections should be read."
)
