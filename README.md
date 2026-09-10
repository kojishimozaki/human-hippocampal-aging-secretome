# A support-cell secretome program in human hippocampal aging is directionally reproducible but not gene-specific

Analysis code, pre-specified analysis documents, audit trail and result tables for the
manuscript of that title.

Koji Shimozaki, Brain Science Research Unit, Graduate School of Biomedical Sciences,
Nagasaki University.

The manuscript of record is [`manuscript/biorxiv_v2/manuscript_biorxiv_v2.md`](manuscript/biorxiv_v2/manuscript_biorxiv_v2.md).
This repository is where it is published: the typeset manuscript with its figures is
[`manuscript/biorxiv_v2/shimozaki-aging-secretome-v2.pdf`](manuscript/biorxiv_v2/shimozaki-aging-secretome-v2.pdf)
(54 pages: 39 of text, then Figures 1–6 and S1–S9, which are also here one per file under
[`manuscript/biorxiv_v2/figures/`](manuscript/biorxiv_v2/figures)). The tables every number
in it comes from are under `results/`. It has not been through peer review.

## What this is

A pre-specified reanalysis of public human hippocampal single-nucleus multiome data
through a 2,224-gene secretome lens, with **donors — not nuclei — as the unit of
replication**. The design is calibrated by donor-label permutation before any biological
result is interpreted, and every claim is paired with a control built to break it.

Four of the study's substantive claims are negative, and the repository contains the
tables behind those as well as the positive ones.

## Reproducing the results

```bash
conda env create -f envs/bio.yml       # Python 3.11: scanpy 1.11.5, pyDESeq2 0.5.4, ...
conda activate bio
export PROJ=$(pwd)

make figures              # figure panels and composites, from the committed result tables
make validation           # senescence / proteome direction tests
make supplementary-tables # Tables S1, S2, S3
make bh-scope             # five per-cell-type BH families vs one genome-wide family
make verify-manuscript    # 76 checks of the manuscript against its source tables
```

These five targets run from a clone of this repository using only what is committed here.
Random seed is 42 throughout; environments are pinned in `envs/`.

`make verify-manuscript` is the one to run first if you want to know whether the numbers
in the manuscript are the numbers in the tables. It re-derives each calibration value from
`Table_S1`, recounts the 60-pair composition from the differential-expression table,
checks that every supplementary figure's footer matches its legend, confirms the Word file
and the markdown carry identical text, and confirms that the typeset PDF sends readers to
the same three addresses the text of record does. It does not take the manuscript's word
for anything.

`make help` lists the remaining targets, which recompute the upstream tables and need the
deposited matrices (see below).

## What is deliberately not here

| Not included | Why |
|---|---|
| `raw/` (113 GB), `processed/` (59 GB) | The primary data are public — see Data sources. These are regenerable, and too large to redistribute. |
| Publisher supplementary files for the SenMayo, SASP Atlas, brain-proteome and CSF references | Third-party material. The small derived tables the pipeline actually reads are included; `refs/README.md` says where to obtain the originals. |
| A separate project arm (donor-dispersion pilot) | Not cited in this manuscript. |
| Superseded manuscript drafts | The current text is under `manuscript/biorxiv_v2/`. |
| The audit narrative (`audit_log/`) | Records what each audit round found and changed. Read by no script and not needed to reproduce any result; the pre-specification documents it contained are in `preregistration/`. |
| The Japanese mirror of the manuscript | Maintained alongside the English in the working repository; it is a translation, not a source of results. |
| Generated figure panels under `figures/` | `make figures` produces them from the committed tables. |

## Limits of reproducibility, stated rather than implied

1. The script that exported count matrices from the deposited Seurat object was not
   retained, so the ingestion step — cell filtering, quality thresholds, RNA–ATAC barcode
   matching — cannot be re-executed here. The workflows repeat from the processed-matrix
   stage onward.
2. Re-running the primary differential-expression target in the pinned environment returns
   99,543 rows rather than the 99,544 of the table used in the manuscript, because the
   low-count oligodendrocyte transcript `LINC01238` falls below the current gene filter.
   All 60 predefined pairs are present and retain adjusted *P* < 0.1; the largest absolute
   difference in log2 fold change across those 60 is 4.3×10⁻⁵.
3. The scripts that produced the donor-label permutation calibration (Table S1), the
   effect-size-matched reproduction background and the external donor-label null were run
   in the pre-submission audit environment and are not part of this repository. Their
   output tables and a provenance record for each are in
   `results/audit_reruns/calibration/`, and every manuscript number drawn from them was
   recomputed from those tables.

## Data sources

All primary and secondary data are public.

| Accession | Role | Primary paper |
|---|---|---|
| **GSE268609** | Primary hippocampal multiome (RNA + ATAC); the young-versus-older control contrast | Disouky et al., *Nature* 2026 |
| **GSE278576** | Independent reproduction cohort (40 donors); the frozen 60-pair signature is *tested*, not redefined | Zemke et al., *Science* 2026 |
| GSE299139 | Matched snm3C methylation arm of the same study | Zemke et al., *Science* 2026 |
| GSE186538 | Supporting continuous-age cohort covering the endothelial arm | Franjic et al., *Neuron* 2022 |
| GSE199243 | Boundary comparison; signed reproducibility 0.02 | Su et al., *Cell Stem Cell* 2022 |
| GSE264692 | Visium spatial | Thompson et al., *Nat Neurosci* 2025 |
| GSE325391 | Cognitive-resilience granule lineage | Tosoni et al., *Cell Stem Cell* 2026 |

## Layout

```
scripts/      numbered pipeline: 01_download -> 02_qc -> 03_de -> 04_atac -> 05_regulatory
              -> 06_validation -> 07_figures -> ... -> 17_biorxiv_v2 (manuscript build)
results/      committed analysis outputs; the source of every number in the manuscript
refs/         the secretome definition and the derived external reference tables
manuscript/   biorxiv_v2/ (text of record in English, the typeset PDF, the Word file,
              and figures/ one PDF per figure) and supplementary_tables/
preregistration/   the pre-specified analysis documents, with a README mapping each one
              to the analysis it governs and to where the manuscript relies on it
docs/*_PLAN.md     the arm-level analysis plans
envs/         pinned conda environments
```

The manuscript reached its final form through six adversarial audit rounds. The
pre-specification each round was held to is in `preregistration/`; the audit narrative
itself is kept in the author's working repository, since no script reads it and it is not
needed to reproduce any result here.

## Licence and citation

Code and result tables are released under the MIT Licence (`LICENSE`). The manuscript text
and figures are © the author, all rights reserved, and are not covered by that licence; if
they are later accepted by a journal, that journal's terms will apply to them instead.
Third-party reference data retain their own licences.

Citation metadata is in `CITATION.cff`; the Zenodo deposit is described by `.zenodo.json`.
The manuscript and the code behind it are one Zenodo record, registered as a preprint, so
there is one thing to cite: concept DOI **10.5281/zenodo.22640702**, which always resolves
to the newest version.
