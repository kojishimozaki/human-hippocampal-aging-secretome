# Deviation log — CONSEQUENCE arm

## D1 (minor implementation fix; §F): mixed-model package nlme instead of lme4
- Pre-registration A3/A4 specified `receiver_score ~ group*maturation_stage + (1|Run) + (1|donor)`, REML,
  Kenward–Roger or Satterthwaite df/CI (lme4/lmerTest).
- **lme4/lmerTest/emmeans/pbkrtest are not installed in any conda env** (sgz_r is the only R; lme4 deps
  minqa/nloptr absent). Installing was not attempted (pinned env; offline-uncertain).
- **Same model, different package:** each donor belongs to exactly one Run, so donor is nested in Run and
  `(1|Run)+(1|donor)` ≡ `nlme::lme(random = ~1|Run/donor)`. The REML linear mixed model and the random
  structure are UNCHANGED. Sum-coding `maturation_stage` makes the `groupSAD` fixed-effect coefficient equal
  the pre-specified stage-averaged SAD−RES estimand directly.
- **Only difference:** df/CI use nlme's containment (within-stratum) approximation rather than KR/Satterthwaite.
  For the group (between-donor) contrast this yields the donor-stratum df. This is an implementation detail,
  logged here; scientific design (model, estimand, two-sided test, hierarchy) is unchanged. Revisit if lme4
  becomes installable.

## D2 (implementation judgment; §F): NicheNet enrichment test not fully fixed pre-outcome
- The pre-registration/manifest fixed "competitive gene-set enrichment of the frozen target set on the
  donor-level SAD−RES ranking" but did NOT fix the exact statistic. To avoid presenting a post-hoc-chosen
  test as registered, BOTH are reported and neither is cherry-picked:
  (a) competitive Wilcoxon rank-sum of the SAD−RES statistic, target vs non-target (registered "competitive"
      reading);
  (b) ranked-list GSEA via fgsea (rank-based enrichment) as a method sensitivity.
- Fixed details (documented): ranking = SAD−RES = −(RES-vs-SAD DESeq2 Wald stat) from de_GSE325391_RESvSAD.csv
  (donor-pseudobulk DESeq2; verified orientation: stat>0 = up in RES, e.g. EFHB/TNS3); duplicate genes
  collapsed (drop_duplicates, keep first); gene universe = DE genes with non-NA stat; targets restricted to
  CTRL-granule-detectable (CPM≥1 in ≥3/6 CTRL donors in ≥1 primary subtype); top-N by regulatory potential
  with 0-weight entries excluded; ties at the N-th rank broken by matrix row order (logged).
