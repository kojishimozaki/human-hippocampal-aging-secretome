# bioRxiv resubmission — note for the screening team

*Draft. Paste into the "comments to the editor" box, or send as the appeal message.
Fill the three bracketed items first.*

---

Dear bioRxiv screening team,

Our earlier submission (manuscript [SUBMISSION ID]) was declined as out of scope, on the
ground that simple automated or computational analyses of public data are generally not
sufficient. We have revised the manuscript and would like to resubmit.

We think the original submission gave a misleading impression of what the study contains,
in two concrete ways that we have now fixed.

**1. The submitted version had lost its code and reproducibility statements.** A
Code availability section and a Methods subsection on reproducibility were present in our
working draft and were dropped when the manuscript was converted to its final format. That
omission is corrected. The revised manuscript states where the analysis code, the
pre-specified analysis documents, the pinned software environments and the result tables
underlying every reported number are deposited ([REPOSITORY URL], archived at
[ZENODO DOI]), which of the workflow targets run from a clone of the repository alone,
which additionally require the deposited matrices, and — explicitly — the three points at
which our own reproducibility is limited, including a one-row difference between the
deposited differential-expression table and its regeneration in the pinned environment.

**2. The abstract described the results but not the design.** The work is not a pipeline
run over public data. Every claim in it is paired with a control built to break it: the
8-versus-9 donor contrast is first calibrated by 100 donor-label permutations with the
complete model refitted each time, which is what tells us that this design yields at least
one Benjamini–Hochberg-significant gene in 29–42% of permutations under the null and that
hit counts are therefore not an interpretable currency here; reproduction in a second
cohort is tested against expression- and effect-size-matched control genes, which removes
our own gene-level specificity claim; the chromatin result is tested against a matched
null, a sign-permutation null, leave-one-gene-out analysis and an end-to-end permutation
that repeats gene selection itself; an apparent RNA–methylation coordination is reported
as failing all four of its specificity controls; and three cell-type labels carried by the
public deposit are re-analysed and excluded. Four of the study's substantive claims are
negative, and they are reported at the same length as the positive ones. The revised title
and abstract now say this, and a new first Results subsection reports the calibration of
the design before any biological result is interpreted.

The revision also corrects four supplementary figure numbers that had drifted out of step
with their legends, and restores five measured values that had been replaced by
approximations during formatting.

We would be glad to answer any further questions.

Yours sincerely,
Koji Shimozaki
Brain Science Research Unit, Graduate School of Biomedical Sciences, Nagasaki University
