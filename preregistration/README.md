# Pre-specified analysis documents

These are the documents the manuscript refers to when it says an analysis was
pre-specified, pre-registered or separately registered. Each was written and committed
**before** the analysis it governs was run, so its git timestamp in the working repository
is the evidence that the decision rules were fixed in advance.

| Document | Governs | Manuscript |
|---|---|---|
| `secretome_resilience_trajectory_v1.md` + `amendments/v1.1`–`v1.4` | Within-cohort AD projection and the independent resilience receiver analysis | Fig. S9; "Two pre-registered projections" |
| `tf_rna_concordance_v1.md` | Whether transcription of each cognate factor follows its chromVAR motif change | Fig. S6; "A separately registered analysis" |
| `manifests/consequence_design_v1.json`, `frozen_input_checksums_v1.tsv`, `source_eligibility_v1.json`, `trajectory_design_v1.json` | Frozen inputs, eligibility rules and design for the projection analyses | Methods, "Receiver map, Alzheimer's-disease projection, and resilience analysis" |
| `prespec/*.md` | Decision criteria fixed before each round of corrective re-analysis, including the gate the regenerated primary differential-expression table had to pass | Methods, "Reproducibility" |
| `../docs/RECEIVER_MAP_PLAN.md`, `NEURON_NICHE_PLAN.md`, `METHYLOME_HARDENING_PLAN.md`, `IN_SILICO_HARDENING_PLAN.md` | The receiver-map, neuron-specificity, methylation and in-silico hardening arms | Results and Methods for those arms |
| `../results/regulatory/PHASE1_5_PRE_REGISTRATION.md` | The candidate regulatory-mechanism search | "Before analysis, we specified the hypothesis that…" |

## A note on paths

The files under `prespec/` and `tf_rna_concordance_v1.md` sit inside `audit_log/` in the
author's working repository, and a few script docstrings still refer to them by those
original paths (`audit_log/<date>_<topic>/`). The audit narrative itself is not part of
this release — no script reads it and it is not needed to reproduce any result — so those
documents were moved here rather than dropped. The content is unchanged.
