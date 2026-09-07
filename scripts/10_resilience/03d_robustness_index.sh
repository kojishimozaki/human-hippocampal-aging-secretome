#!/usr/bin/env bash
# CONSEQUENCE (Phase A) — STEP 03d: consolidated robustness index (presence table).
# Provenance: VERBATIM (non-git fragment) from origin session 027c5fdb (transcript line 905). In the original this
#   was embedded in the commit command (which also wrote claim_eligibility.md via a heredoc + git add/commit/tag);
#   only the index-building loop is the reproducible artifact and is preserved here. Run from $PROJ. Run last.
# Out: results/resilience/robustness_index.tsv
set -euo pipefail
# consolidated robustness index
{ echo -e "file\tcontent"
  for f in robustness_summary robustness_run_fixed robustness_no_run robustness_leave_one_run robustness_leave_one_ligand_family robustness_subtype_specific; do
    echo -e "results/resilience/${f}.tsv\t$( [ -f results/resilience/${f}.tsv ] && echo present || echo MISSING )"; done
} > results/resilience/robustness_index.tsv
