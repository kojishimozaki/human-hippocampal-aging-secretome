# External reference data

The pipeline reads the small derived tables in this directory. The publisher supplementary
files they were built from are **not redistributed here** — they belong to their
publishers, and no script in this repository reads them. Obtain them from the sources below
if you want to rebuild a derived table from scratch.

| Derived table (included) | Built from | Source |
|---|---|---|
| `csf_aging_estimate.csv` | Supplementary tables of the CSF aging proteome | Seo et al., *Nat Aging* 2025;5(10):2125–2141. doi:10.1038/s43587-025-00971-6 |
| `wingo_cogtraj_beta.csv` | Supplementary tables of the brain proteome cognitive-trajectory study | Wingo et al., *Nat Commun* 2019;10(1):1619. doi:10.1038/s41467-019-09613-z (CC BY 4.0) |
| `sasp_atlas_core.txt` | SASP Atlas supplementary tables | Basisty et al., *PLoS Biol* 2020;18(1):e3000599. doi:10.1371/journal.pbio.3000599 (CC0) |
| `senmayo.txt` | MSigDB gene set `SAUL_SEN_MAYO` | Saul et al., *Nat Commun* 2022;13(1):4827. doi:10.1038/s41467-022-32552-1 |
| `secretome_union.csv`, `secretome_core.csv` | Union of the Human Protein Atlas "predicted secreted" annotation and reviewed human UniProtKB entries annotated as secreted | `hpa_secretome.tsv`, `uniprot_secreted_human.tsv` (both included); built by `scripts/01_download/02_build_secretome.py` |
| `metadata_master.csv`, `metadata_raw.csv` | GEO SOFT sample metadata for the cohorts used | Public GEO records; de-identified |

`enrichr/` holds gene-set libraries used by the neuronal gene-set analysis.
`gencode_v44_genes.tsv.gz` is the GENCODE v44 gene annotation.

Nothing in this directory contains donor-identifying information: the GEO metadata here is
the de-identified sample annotation as deposited.
