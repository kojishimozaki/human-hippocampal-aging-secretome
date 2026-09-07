# NCBI citation audit: `manuscript_proofread_tracked_v3.docx`

Audit date: 2026-08-29 (JST)

## Overall conclusion

- The tracked DOCX contains **36 numbered reference entries**. Two entries contain an additional DOI-bearing work, giving **38 DOI-bearing records** in total.
- **All 38 DOI records resolve uniquely in NCBI PubMed.** No fabricated or non-existent reference was detected.
- A PubMed search restricted to `Retracted Publication` and `Retraction of Publication` returned **0 records** among the 38 records.
- The title, first author, journal, year, volume/issue, and DOI are essentially correct throughout. Reference 3 has one minor pagination omission: PubMed gives `589–599.e5`, whereas the manuscript gives `589–599`.
- The in-text uses of references 1–36 are contextually appropriate. No clear case was found in which a cited paper is unrelated to the sentence it supports.
- The bibliography is nevertheless **not yet submission-ready**: two analysed GEO datasets have published primary papers that are not cited, and references 18 and 22 contain avoidable version/record mixing.

## Changes recommended before submission

### 1. Add the primary paper for GSE264692

The manuscript uses GSE264692 for the Visium analysis but does not cite its published primary paper. NCBI GEO currently has no PubMed link in the series record, but a PubMed title search identifies the exact published article:

> Thompson JR, Nelson ED, Tippani M, et al. An integrated single-nucleus and spatial transcriptomics atlas reveals the molecular landscape of the human hippocampus. *Nat Neurosci.* 2025;28(9):1990–2004. doi:10.1038/s41593-025-02022-0. PMID: [40739059](https://pubmed.ncbi.nlm.nih.gov/40739059/).

Recommended citation locations include the Methods paragraph introducing GSE264692, the spatial-localisation Methods paragraph, Data availability, and the Fig. S4/associated figure legend if the journal expects citations in legends.

### 2. Add the primary paper for GSE325391

NCBI GEO links GSE325391 to PMID 42034060, but that paper is absent from the reference list:

> Tosoni G, Ayyildiz D, Snoeck S, et al. Transcriptional profiles of immature neurons in aged human hippocampus track Alzheimer's pathology and cognitive resilience. *Cell Stem Cell.* 2026;33(5):763–783.e9. doi:10.1016/j.stem.2026.04.002. PMID: [42034060](https://pubmed.ncbi.nlm.nih.gov/42034060/).

Recommended citation locations include the Methods paragraph introducing GSE325391, the Results sections using the resilience cohort, and Data availability.

### 3. Complete reference 3 pagination

Change `589–599` to `589–599.e5` to match PubMed PMID [29625071](https://pubmed.ncbi.nlm.nih.gov/29625071/).

### 4. Simplify reference 18

The final Science paper and the bioRxiv preprint both exist and PubMed explicitly links them as update/preprint versions of the same work:

- Final article: PMID [42490474](https://pubmed.ncbi.nlm.nih.gov/42490474/)
- Preprint: PMID [39463924](https://pubmed.ncbi.nlm.nih.gov/39463924/)

Because the final article is now published and indexed, retain the Science citation and normally remove the parenthetical bioRxiv citation. Keeping both is not factually wrong, but it is redundant in a final bibliography.

### 5. Do not place two independent papers under reference 22

Reference 22 currently combines:

- *The human secretome*, PMID [31772123](https://pubmed.ncbi.nlm.nih.gov/31772123/), which directly supports the “predicted secreted” annotation; and
- *Proteomics. Tissue-based map of the human proteome*, PMID [25613900](https://pubmed.ncbi.nlm.nih.gov/25613900/).

Both papers are real and the metadata are correct, but placing two independent papers under one reference number makes the in-text citation ambiguous. The 2019 secretome paper is the direct source for the stated annotation. Delete the 2015 paper if it is not needed, or give it a separate reference number and full title.

## GEO-to-publication checks

| GEO series | Use in manuscript | NCBI-linked/corresponding PubMed record | Judgment |
|---|---|---|---|
| [GSE268609](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268609) | Primary RNA/ATAC multiome | Disouky et al., PMID [41741649](https://pubmed.ncbi.nlm.nih.gov/41741649/) | Correct mapping |
| [GSE278576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278576) | Reproduction RNA/ATAC multiome | Zemke et al.; GEO still links the preprint PMID [39463924](https://pubmed.ncbi.nlm.nih.gov/39463924/), which PubMed links to the final Science PMID [42490474](https://pubmed.ncbi.nlm.nih.gov/42490474/) | Correct mapping; cite final article |
| [GSE299139](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299139) | Matched snm3C arm | Same Zemke study as above | Correct mapping |
| [GSE186538](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186538) | Supporting continuous-age cohort | Franjic et al., PMID [34798047](https://pubmed.ncbi.nlm.nih.gov/34798047/) | Correct mapping |
| [GSE199243](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199243) | Boundary comparison | Su et al., PMID [36332572](https://pubmed.ncbi.nlm.nih.gov/36332572/) | Correct mapping |
| [GSE264692](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264692) | Visium spatial analysis | Thompson et al., PMID [40739059](https://pubmed.ncbi.nlm.nih.gov/40739059/) | Primary-paper citation missing |
| [GSE325391](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325391) | Cognitive-resilience analysis | Tosoni et al., PMID [42034060](https://pubmed.ncbi.nlm.nih.gov/42034060/) | Primary-paper citation missing |

## Item-by-item assessment

“Accurate” means the DOI uniquely resolves in PubMed and the manuscript's core bibliographic fields match the PubMed record. “Appropriate” means that the paper's subject and study design support the claim for which it is cited; it does not mean that the cited study proves the present manuscript's own findings.

| Ref. | PubMed | Bibliographic accuracy | Contextual appropriateness | Note |
|---:|---:|---|---|---|
| 1 | [27814520](https://pubmed.ncbi.nlm.nih.gov/27814520/) | Accurate | Appropriate | Review supports the adult-hippocampal-neurogenesis framework. |
| 2 | [29513649](https://pubmed.ncbi.nlm.nih.gov/29513649/) | Accurate | Appropriate | Directly supports the report of sharply declining/undetectable adult human hippocampal neurogenesis. |
| 3 | [29625071](https://pubmed.ncbi.nlm.nih.gov/29625071/) | Minor correction | Appropriate | Add the electronic supplement suffix: `589–599.e5`. |
| 4 | [30911133](https://pubmed.ncbi.nlm.nih.gov/30911133/) | Accurate | Appropriate | Directly supports persistence in healthy adults and a marked reduction in Alzheimer's disease. |
| 5 | [25611508](https://pubmed.ncbi.nlm.nih.gov/25611508/) | Accurate | Appropriate | Direct human-hippocampus evidence for age-associated blood–brain barrier breakdown. |
| 6 | [30944478](https://pubmed.ncbi.nlm.nih.gov/30944478/) | Accurate | Appropriate | Experimental CD22 blockade/microglial phagocytosis precedent. |
| 7 | [33473210](https://pubmed.ncbi.nlm.nih.gov/33473210/) | Accurate | Appropriate | Experimental myeloid-metabolism/cognitive-outcome precedent. |
| 8 | [31086348](https://pubmed.ncbi.nlm.nih.gov/31086348/) | Accurate | Appropriate | Directly supports aged-blood, endothelial VCAM1, microglial, and neural-precursor effects. |
| 9 | [31413369](https://pubmed.ncbi.nlm.nih.gov/31413369/) | Accurate | Appropriate | Directly supports niche stiffness as an aging mechanism in CNS progenitor cells. |
| 10 | [34880498](https://pubmed.ncbi.nlm.nih.gov/34880498/) | Accurate | Appropriate | Experimental exercise-plasma/clusterin effects on memory and neuroinflammation. |
| 11 | [35545674](https://pubmed.ncbi.nlm.nih.gov/35545674/) | Accurate | Appropriate | Directly supports young-CSF effects on oligodendrogenesis and memory in aged mice. |
| 12 | [20078217](https://pubmed.ncbi.nlm.nih.gov/20078217/) | Accurate | Appropriate | Foundational SASP review; suitable for SASP definition and biological heterogeneity. |
| 13 | [35974106](https://pubmed.ncbi.nlm.nih.gov/35974106/) | Accurate | Appropriate | Primary SenMayo gene-set source. The manuscript appropriately states its non-brain validation context. |
| 14 | [31945054](https://pubmed.ncbi.nlm.nih.gov/31945054/) | Accurate | Appropriate | Primary SASP Atlas/proteomic secretome source. The manuscript appropriately limits this to cultured-cell conditioned media. |
| 15 | [30962425](https://pubmed.ncbi.nlm.nih.gov/30962425/) | Accurate | Appropriate | Direct source for the human-brain proteomic cognitive-trajectory comparison. |
| 16 | [41087722](https://pubmed.ncbi.nlm.nih.gov/41087722/) | Accurate | Appropriate with stated limits | Supports age-related CSF proteomic comparison. The manuscript correctly treats it as suggestive and not as measured hippocampal secretion. |
| 17 | [41741649](https://pubmed.ncbi.nlm.nih.gov/41741649/) | Accurate | Appropriate | Correct primary paper for GSE268609; GEO links the same PMID. |
| 18 | [42490474](https://pubmed.ncbi.nlm.nih.gov/42490474/) / [39463924](https://pubmed.ncbi.nlm.nih.gov/39463924/) | Accurate but redundant | Appropriate | Final Science article and its preprint are the same study. Prefer the final article alone. |
| 19 | [34798047](https://pubmed.ncbi.nlm.nih.gov/34798047/) | Accurate | Appropriate | Correct primary paper for GSE186538; GEO links the same PMID. |
| 20 | [36332572](https://pubmed.ncbi.nlm.nih.gov/36332572/) | Accurate | Appropriate | Correct primary paper for GSE199243; GEO links the same PMID. |
| 21 | [30914743](https://pubmed.ncbi.nlm.nih.gov/30914743/) | Accurate | Appropriate | Primary Leiden-method paper for the clustering method. |
| 22 | [31772123](https://pubmed.ncbi.nlm.nih.gov/31772123/) / [25613900](https://pubmed.ncbi.nlm.nih.gov/25613900/) | Accurate but structurally ambiguous | Appropriate | Both works exist; the 2019 secretome paper is the direct annotation source. Split or delete the secondary 2015 citation. |
| 23 | [36408920](https://pubmed.ncbi.nlm.nih.gov/36408920/) | Accurate | Appropriate | Suitable UniProtKB database citation for the stated reviewed-entry annotation. |
| 24 | [37669147](https://pubmed.ncbi.nlm.nih.gov/37669147/) | Accurate | Appropriate | Primary pyDESeq2 paper; version number in the manuscript is a software-version statement, not a bibliographic claim. |
| 25 | [29409532](https://pubmed.ncbi.nlm.nih.gov/29409532/) | Accurate | Appropriate | Primary Scanpy paper. |
| 26 | [28825706](https://pubmed.ncbi.nlm.nih.gov/28825706/) | Accurate | Appropriate | Primary chromVAR paper. |
| 27 | [34725479](https://pubmed.ncbi.nlm.nih.gov/34725479/) | Accurate | Appropriate | Primary Signac paper. |
| 28 | [31701148](https://pubmed.ncbi.nlm.nih.gov/31701148/) | Accurate | Appropriate | Primary JASPAR 2020 database paper. |
| 29 | [31819264](https://pubmed.ncbi.nlm.nih.gov/31819264/) | Accurate | Appropriate | Primary NicheNet paper. |
| 30 | [36426870](https://pubmed.ncbi.nlm.nih.gov/36426870/) | Accurate | Appropriate | Primary GSEApy paper. |
| 31 | [28099414](https://pubmed.ncbi.nlm.nih.gov/28099414/) | Accurate | Appropriate | Directly supports induction of reactive astrocyte states by activated microglia. |
| 32 | [30232451](https://pubmed.ncbi.nlm.nih.gov/30232451/) | Accurate | Appropriate | Direct tauopathy-model precedent for senescent-glial-cell clearance and cognitive/pathology outcomes; the manuscript correctly avoids extrapolating this to control aging. |
| 33 | [27984732](https://pubmed.ncbi.nlm.nih.gov/27984732/) | Accurate | Appropriate | Suitable precedent for pooled perturbation with single-cell readouts. |
| 34 | [30301888](https://pubmed.ncbi.nlm.nih.gov/30301888/) | Accurate | Appropriate | Supports cerebral organoid systems incorporating microglia. |
| 35 | [31591580](https://pubmed.ncbi.nlm.nih.gov/31591580/) | Accurate | Appropriate | Supports brain organoids with vascular-like components. |
| 36 | [30911168](https://pubmed.ncbi.nlm.nih.gov/30911168/) | Accurate | Appropriate | Supports transcriptome-scale, spatially resolved RNA imaging as a technical route. |

## Scope and interpretation

This audit used the final tracked-text view of the DOCX: deleted text was excluded and inserted text was included. Existence and metadata were checked through NCBI Entrez/PubMed by DOI, and dataset–publication relationships were checked through NCBI GEO and PubMed. Contextual appropriateness was assessed against the exact in-text sentence, PubMed metadata/abstract information, and GEO study descriptions. It was not a full-text reproduction of every cited experiment.
