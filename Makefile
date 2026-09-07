# Batch regeneration for the hippocampal aging niche-secretome reanalysis.
# Run from the repo root with the `bio` conda env active (and `sgz_r` for the R/* targets).
#
#   conda activate bio
#   make help
#
# PROJ defaults to the current checkout; every script reads it from the environment.
# Override interpreters if your envs are named differently, e.g.  make figures PY=python

PROJ ?= $(CURDIR)
export PROJ
# Keep scanpy/numba + matplotlib quiet and deterministic on a cacheless/read-only HOME (e.g. a bare clone).
MPLCONFIGDIR    ?= /tmp/mpl_cache
NUMBA_CACHE_DIR ?= /tmp/numba_cache
export MPLCONFIGDIR NUMBA_CACHE_DIR
PY      ?= python      # expects the `bio` env on PATH (conda activate bio)
RSCRIPT ?= Rscript     # expects the `sgz_r` env on PATH (conda activate sgz_r)
# `make trajectory` spans TWO conda envs (bio for the .py steps, trajR for the .R steps); pinned absolute
# interpreters so it runs regardless of the active env (override if your envs live elsewhere). Data-dependent.
TRAJ_PY      ?= $(HOME)/miniforge3/envs/bio/bin/python
TRAJ_RSCRIPT ?= $(HOME)/miniforge3/envs/trajR/bin/Rscript

.PHONY: help figures replication validation all-light atac-concordance de-robustness supplementary-tables \
        figures-data regulatory trajectory tf-rna-concordance resilience receiver \
        resilience-de primary-de bh-scope biorxiv-v2 verify-manuscript strict-secretome

help:
	@echo "Runnable from a BARE CLONE (only committed results/ + refs/ -- no raw/ or processed/ data needed):"
	@echo "  make figures       Fig 2,4,6 (+ composites), Fig 4D footprints, CSF S3,"
	@echo "                     Fig 3 replication, Fig S(E1) external-cohort replication spectrum, Fig S(N) neuron specificity, Fig S(R1) TF-RNA concordance, Fig S(RT1) AD-trajectory/resilience, Fig S(Rcv) receiver-map (regenerated from committed results/*.csv)"
	@echo "  make validation    SASP-overlap + proteome direction tests -> results/validation/"
	@echo "  make all-light     figures + validation"
	@echo "  NOTE: figures/figureN/ map 1:1 to manuscript Fig N (physically renumbered 2026-06-29)."
	@echo "  supplementary-tables  Assemble Supplementary Tables S1/S2/S(R1) -> manuscript/supplementary_tables/ (committed data)"
	@echo "  make bh-scope      Five per-cell-type BH families vs one genome-wide family -> results/de/global_bh_scope_GSE268609.csv"
	@echo "  make strict-secretome  Sensitivity of the main results to the HPA-and-UniProt intersection universe"
	@echo "  make verify-manuscript  Check the manuscript against its source tables (56 checks; runs from a clone)"
	@echo "  make biorxiv-v2    Rebuild the manuscript (DOCX + EN/JA markdown), then verify. Needs the"
	@echo "                     superseded submission DOCX as its edit base, which is not redistributed."
	@echo "        Build SCRIPT filenames are historical: 03_figure3.py builds Fig 4; 04b_visium_improved.py builds Fig 5 (Visium, data-dependent);"
	@echo "        15_replication_spectrum_supp.py builds Fig S(E1); 05_figure5.py builds Fig 6; 03_de/05_external_replication.py builds Fig 3 (was figureA1)."
	@echo ""
	@echo "DATA-DEPENDENT -- need the large, GITIGNORED raw/ + processed/ data (NOT in this repo;"
	@echo "                  these RECOMPUTE the result tables that 'make figures' merely re-plots):"
	@echo "  make atac-concordance recompute RNA<->ATAC concordance   (needs raw/ features + processed/ peaks + linked-peak cache)"
	@echo "  make de-robustness    PMI-adjusted + leave-one-donor-out  (needs raw/ SOFT + processed/ anchor h5ad)"
	@echo "  make figures-data     Fig 1 (+composite) + Fig 5 Visium/composite (Visium-only scope) + methyl Fig S-M1 (needs processed/ data)"
	@echo "  make regulatory       motif / NicheNet / footprint (R)    (needs processed/ ; conda activate sgz_r)"
	@echo "  make tf-rna-concordance  motif<->TF-RNA concordance -> Fig S(R1) tables   (bio; needs processed/)"
	@echo "  make trajectory          within-cohort aging->AD projection -> Fig S(RT1) (DUAL env: bio .py + trajR .R, pinned abs paths)"
	@echo "  make resilience          CONSEQUENCE Phase-A receiver-competence -> results/resilience/ (DUAL env bio+sgz_r; needs GSE325391 h5ad)"
	@echo "  make resilience-de       GSE325391 pairwise group DE -> results/de/de_GSE325391_*.csv (bio; needs GSE325391 h5ad)"
	@echo "  make primary-de          PRIMARY DE (origin of the frozen 60) -> results/de/de_GSE268609_per_celltype.regenerated.csv (bio; needs anchor h5ad)"
	@echo "  make receiver            OUTPUT->receiver receptor-regulation map -> results/receiver/ + Fig S(Rcv) tables (bio+sgz_r; needs refs/nichenet/*.rds)"

# ---- runnable from a bare clone (committed results/ + refs/ only) --------------------------------
# Each script writes its legend-referenced single-panel PDFs AND, for Fig 2/4/6, an assembled
# composite figures/figureN/figureN.pdf. Shared publication style: scripts/07_figures/_pubstyle.py.
# Assemble Supplementary Tables S1/S2/S(R1) from committed result CSVs (no new computation).
supplementary-tables:
	$(PY) scripts/16_supplementary_tables/build_supplementary_tables.py

# Multiple-testing scope check behind the Results calibration subsection (committed data only).
bh-scope:
	$(PY) scripts/03_de/03_global_bh_scope.py

# Do the results depend on genes annotated as secreted by only one of the two sources?
strict-secretome:
	$(PY) scripts/06_validation/03_strict_secretome_sensitivity.py

# Rebuild the current manuscript from scripts/17_biorxiv_v2/revision_edits.py and check it.
# Figure renumbering is a one-off over submission/figures/ and is not re-run here.
biorxiv-v2:
	$(PY) scripts/17_biorxiv_v2/build_revision.py
	$(PY) scripts/17_biorxiv_v2/build_japanese_mirror.py
	$(MAKE) verify-manuscript

# Checks only: reads the committed manuscript, tables and figures, so it runs from a clone.
verify-manuscript:
	$(PY) scripts/17_biorxiv_v2/verify_revision.py

figures:
	$(PY) scripts/07_figures/02_figure2.py
	$(PY) scripts/07_figures/02b_csf_direction.py
	$(PY) scripts/07_figures/03_figure3.py
	$(PY) scripts/05_regulatory/05_footprint_figure.py
	$(PY) scripts/07_figures/15_replication_spectrum_supp.py
	$(PY) scripts/07_figures/05_figure5.py
	$(PY) scripts/09_neuron/06_supp_figures.py
	$(PY) scripts/07_figures/12_tf_rna_concordance_supp.py
	$(PY) scripts/07_figures/13_resilience_trajectory_supp.py
	$(PY) scripts/07_figures/14_receiver_supp.py
	$(MAKE) replication
	@echo "Composites figures/figureN/figureN.pdf for Fig 2,4,6 and figureSE1 (Fig S(E1)), figureS_neuron (Fig S(N)), figureSR1 (Fig S(R1)), figureSRT1 (Fig S(RT1)), figureS_receiver (Fig S(Rcv)). Fig 1 + Fig 5 (Visium, and their composites)"
	@echo "and methyl Fig S-M1 are data-dependent — see 'make figures-data' / 'make help'."

# Fig 3 (was Fig A1) — test the frozen signature in each external cohort (committed external DE; seed 42).
replication:
	$(PY) scripts/03_de/05_external_replication.py --cohort GSE278576        --external results/de/de_GSE278576_per_celltype.csv
	$(PY) scripts/03_de/05_external_replication.py --cohort GSE278576_age    --external results/de/de_GSE278576_age_per_celltype.csv
	$(PY) scripts/03_de/05_external_replication.py --cohort GSE278576_sexadj --external results/de/de_GSE278576_sexadj_per_celltype.csv
	$(PY) scripts/03_de/05_external_replication.py --cohort GSE186538        --external results/de/de_GSE186538_per_celltype.csv
	$(PY) scripts/03_de/05_external_replication.py --cohort GSE186538_drop79 --external results/de/de_GSE186538_drop79_per_celltype.csv

validation:
	$(PY) scripts/06_validation/01_sasp_overlap.py
	$(PY) scripts/06_validation/02_proteome_validation.py

all-light: validation figures

# ---- DATA-DEPENDENT: require the gitignored raw/ + processed/ data (will FAIL on a bare clone) ----
atac-concordance:
	$(PY) scripts/04_atac/01_atac_rna_concordance.py

de-robustness:
	$(PY) scripts/03_de/04_aging_robustness.py

figures-data:
	$(PY) scripts/07_figures/01_figure1.py
	$(PY) scripts/07_figures/04b_visium_improved.py
	$(PY) scripts/08_methyl/05_figure.py

regulatory:
	$(RSCRIPT) scripts/05_regulatory/01_tf_motif_enrichment.R
	$(RSCRIPT) scripts/05_regulatory/02_nichenet.R
	$(RSCRIPT) scripts/05_regulatory/03_footprint_sites.R
	$(RSCRIPT) scripts/05_regulatory/04_footprint_test.R
	$(PY)      scripts/05_regulatory/05_footprint_figure.py
	$(RSCRIPT) scripts/05_regulatory/06_autocrine_validation.R

# TF-RNA cross-modal concordance engine -> Fig S(R1) source tables (bio env; needs processed chromVAR + donor-pseudobulk DE)
tf-rna-concordance:
	$(PY) scripts/05_regulatory/06_motif_tf_rna_concordance.py

# Within-cohort aging->AD projection -> Fig S(RT1) trajectory tables. DUAL-ENV + data-dependent (needs processed/):
# .py steps run under bio (TRAJ_PY), .R steps under trajR (TRAJ_RSCRIPT: DESeq2 1.50.2 / ashr). Pinned absolute paths.
trajectory:
	$(TRAJ_PY)      scripts/11_trajectory/01_pseudobulk.py
	$(TRAJ_RSCRIPT) scripts/11_trajectory/02_deseq2_ashr.R
	$(TRAJ_PY)      scripts/11_trajectory/03_classify_modules.py
	$(TRAJ_RSCRIPT) scripts/11_trajectory/04_stability.R
	$(TRAJ_PY)      scripts/11_trajectory/05_cross_impl_audit.py

# CONSEQUENCE (Phase A) receiver-competence resilience driver -> results/resilience/* (Fig S(RT1) panel B).
# Recovered verbatim from the origin session transcript (see scripts/10_resilience/README.md). DUAL-ENV +
# data-dependent: .py under bio (RESIL_PY), .R under sgz_r (RESIL_RSCRIPT); needs the GSE325391 h5ad rebuilt
# from the raw RDS first (scripts/02_qc/06_build_325391_resilience.py). Pinned absolute interpreters.
RESIL_PY      ?= $(HOME)/miniforge3/envs/bio/bin/python
RESIL_RSCRIPT ?= $(HOME)/miniforge3/envs/sgz_r/bin/Rscript
resilience:
	mkdir -p results/step0 results/resilience
	$(RESIL_RSCRIPT) scripts/10_resilience/00_ligand_receptor_map.R
	$(RESIL_PY)      scripts/10_resilience/01_freeze_receptor_set.py
	$(RESIL_PY)      scripts/10_resilience/02_receiver_score.py
	$(RESIL_RSCRIPT) scripts/10_resilience/03_primary_model.R
	$(RESIL_PY)      scripts/10_resilience/03b_robustness_build.py
	$(RESIL_RSCRIPT) scripts/10_resilience/03b_robustness_fit.R
	$(RESIL_PY)      scripts/10_resilience/03c_robustness_extended_build.py
	$(RESIL_RSCRIPT) scripts/10_resilience/03c_robustness_extended_fit.R
	$(RESIL_RSCRIPT) scripts/10_resilience/04a_nichenet_targets.R
	$(RESIL_PY)      scripts/10_resilience/04b_nichenet_enrichment_build.py
	$(RESIL_RSCRIPT) scripts/10_resilience/04b_nichenet_enrichment_fgsea.R
	$(RESIL_RSCRIPT) scripts/10_resilience/04c_nichenet_supporting_targets.R
	$(RESIL_PY)      scripts/10_resilience/04c_nichenet_supporting_wilcox.py
	$(RESIL_PY)      scripts/10_resilience/05_matrix_provenance_audit.py
	bash             scripts/10_resilience/03d_robustness_index.sh

# Pairwise group DE for GSE325391 -> results/de/de_GSE325391_{RESvCTRL,RESvSAD,SADvCTRL}.csv,
# the source of the Fig 6 "RES vs CTRL is 2 genes" statement. The audit (B-8 / F-0-001) found these
# three CSVs had NO write site anywhere in the repository and only an initial-commit history, so the
# numbers could not be reproduced. The commands below were reconstructed from the CSVs' own `cohort`
# column and the h5ad's `group` levels, then verified on 2026-08-26: all three outputs reproduce the
# committed files with a maximum absolute difference of 0 across baseMean, log2FoldChange, lfcSE,
# stat, pvalue and padj (29,688 / 29,521 / 29,563 rows, identical gene sets). See
# audit_log/2026-08-26_remaining_still_open/RESOLUTION.md section 2-1.
# ROW ORDER DIFFERS and that is expected: this target emits baseMean-descending order because
# 01_pseudobulk_de.py:L139-L144 sorts and de-duplicates (the 2026-06-03 audit fix), whereas the
# committed CSVs were generated before that block was applied to them and retain the h5ad
# var_names order. Content is unaffected -- no duplicate (gene,celltype) rows exist in this
# cohort, so the block only reorders. Running this target therefore REWRITES the three committed
# files with the same content in a different order; commit that only if you intend to.
# Data-dependent: needs processed/per_dataset/GSE325391_resilience.h5ad (gitignored; rebuild with
# scripts/02_qc/06_build_325391_resilience.py). NOTE this closes only the first half of B-8 -- the
# Seurat RDS -> MTX export that produced the count matrices is still not in the repository.
resilience-de:
	bash scripts/14_audit_reruns/12_recover_gse325391_de.sh

# [B-12] The PRIMARY differential expression table -- the origin of the frozen 60 and of every
# downstream leg. The audit found no target invoked its engine. This writes BESIDE the committed
# table (…_per_celltype.regenerated.csv), not over it, because the two are not identical: the
# committed file has one extra row (LINC01238/Oligo) that the current gene filter rejects, so it
# was produced with a looser filter. The frozen 60 are unaffected -- 60/60 keep padj < 0.1 with a
# maximum log2FC difference of 4.3e-05. See audit_log/2026-08-27_provenance_b12_b13_b14/.
# Data-dependent: needs processed/per_dataset/GSE268609_anchor.h5ad (gitignored).
primary-de:
	bash scripts/14_audit_reruns/13_primary_de.sh

# OUTPUT->receiver receptor-regulation map -> results/receiver/* (Fig S(Rcv) source tables). Pre-reg tags
# prereg/receiver-map-v1 (+ v1.1); seed 42. DUAL-ENV: 00b reads the SHA-pinned NicheNet rds under sgz_r
# (RECV_RSCRIPT); 00/01/02/02b/03 under bio (RECV_PY). DATA-DEPENDENT on refs/nichenet/*.rds (gitignored) for
# 00b only; 01-03 read the committed donor-pseudobulk DE table + committed manifest. The FIGURE re-plots from
# committed results/receiver/*.csv on a bare clone via `make figures`.
RECV_PY      ?= $(HOME)/miniforge3/envs/bio/bin/python
RECV_RSCRIPT ?= $(HOME)/miniforge3/envs/sgz_r/bin/Rscript
receiver:
	mkdir -p results/receiver
	$(RECV_PY)      scripts/12_receiver/00_freeze_universe.py
	$(RECV_RSCRIPT) scripts/12_receiver/00b_receptor_edges.R
	$(RECV_PY)      scripts/12_receiver/01_detectability.py
	$(RECV_PY)      scripts/12_receiver/02_depth_matched_null.py
	$(RECV_PY)      scripts/12_receiver/02b_depth_testability_v1_1.py
	$(RECV_PY)      scripts/12_receiver/03_matched_null_on_lfc.py
