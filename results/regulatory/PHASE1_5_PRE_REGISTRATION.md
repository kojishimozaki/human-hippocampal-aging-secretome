# Phase 1.5 pre-specified analysis plan — TF-regulatory layer of the aging niche secretome
**Locked 2026-06-04, BEFORE running the enrichment/convergence pipeline** (pre-*specified*: this plan was written before the analysis in-session, but it is NOT git-timestamped ahead of the results — plan, scripts and results were committed together — so external auditors should read this as a *locked analysis plan*, not an independently-verifiable pre-registration). Post-hoc edits to
this plan after seeing results = confirmation bias; any deviation is logged in a DEVIATIONS section.

## Discovery hypothesis (directional, pre-specified from prior biology — NOT post-hoc)
The aging-secretome reprogramming (Fig 2) is driven by a defined, cell-type-resolved TF program:
- **UP arm** (ECM↑ + inflammation↑) is regulated by the **canonical SASP transcriptional program**:
  NF-κB, C/EBP, AP-1, STAT, IRF, GATA4, plus ECM/mechano TEAD & TGF-β/SMAD.
- **DOWN arm** (endothelial/vascular↓) reflects **loss of homeostatic/vascular TFs**:
  ETS family, Wnt-TCF, endothelial SOX/KLF, FOXO.
H0 (must be reported if it holds): the aging niche secretome shows **no** detectable enrichment of
the pre-specified TF panel above a GC/width-matched background at this power.

## Candidate TF motif panel (JASPAR2020 CORE; all verified present in the 633-motif set)
UP-predicted: NFKB1,NFKB2,RELA,RELB,REL | CEBPA,CEBPB,CEBPD,CEBPG | JUN,JUNB,JUND,FOS,FOSB,FOSL1,FOSL2,BATF,BATF3 |
STAT1,STAT2,STAT3,STAT1::STAT2 | IRF1,IRF2,IRF3,IRF8,SPI1 | GATA2,GATA4,GATA6 | TEAD1,TEAD3,TEAD4 |
SMAD2,SMAD3,SMAD4,SMAD2::SMAD3::SMAD4 | RUNX2,RUNX3
DOWN-predicted: ETS1,ETV1,ETV2,ETV4,ETV6,ERG,FLI1,ELK1,ELK3,ELK4,ELF1,EHF,FEV | TCF7,TCF7L1,TCF7L2,LEF1 |
SOX18 | KLF2,KLF4,KLF5,SP1 | FOXO3,FOXO4
(The UP/DOWN label is the *prediction*; the test reports the actual arm, so a prediction-discordant
hit is flagged honestly, not hidden.)

## Inputs (all on hand; no new data)
- RNA DE (YA→HA, per niche celltype): `results/de/de_GSE268609_per_celltype.csv` (secretome-annotated).
- Motif×peak matrix (binary, JASPAR2020 CORE 9606): cached `processed/per_dataset/GSE268609_peaks_motifs.rds`.
- chromVAR motif deviations (donor-level dz, nominal p): `results/de/chromvar_motifs_YAvsHA.csv`.
- Genome / GC: BSgenome.Hsapiens.UCSC.hg38. Secretome gene coords: GSE268609 features.tsv.gz.

## Gene/peak sets (locked)
- Per niche celltype C ∈ {Astro,Micro,Oligo,OPC,Endo}:
  - UP set = secretome genes with padj<0.1 & log2FC>0 in C; DOWN set = padj<0.1 & log2FC<0.
  - A celltype-arm is **evaluable** only if its set has ≥4 genes (else reported NOT_EVALUABLE).
- Foreground peaks = peaks within gene body ±2 kb of the set's genes (same proximity rule as Fig 3).
- Background (matched null) = peaks sampled from the full motif-annotated peak universe, matched to
  the foreground on **GC content (20 bins) × log10 width (5 bins)**, 10× foreground size. GC is the
  dominant confound for motif frequency (field-standard AME/HOMER control). Accessibility-matching is
  a noted refinement (would need the 8.9 G peak matrix; deferred).

## Three-layer evidence + locked judgment
For each (celltype C, arm D∈{UP,DOWN}, motif M in panel):
- **L1 enrichment**: Fisher one-sided (motif+ in foreground vs matched background); BH-FDR over the
  panel *within* each (C,D). **PASS = FDR < 0.10.**
- **L2 chromVAR concordance** (orthogonal, donor-level): the motif's chromVAR dz in C has the arm's
  sign (UP→dz>0, DOWN→dz<0). **PASS = sign-concordant** (nominal; dz & p reported, not thresholded —
  chromVAR is known underpowered at FDR, used here only as a direction axis).
- **L3 targets**: ≥3 secretome genes in the (C,D) set carry M in ≥1 linked peak. **PASS = ≥3.**
Verdicts (locked): **CONVERGENT = L1 & L2 & L3 all PASS**; **SUGGESTIVE = exactly 2 PASS**; else not called.
"Prediction-match" = arm D equals the panel's predicted direction for M (reported as a separate column).

## Anti-bias rules (locked)
- No motif added/removed after seeing results. No threshold changes post-hoc.
- Direction-discordant hits (e.g., a "DOWN-predicted" ETS enriched in the UP arm) are reported, not dropped.
- If 0 CONVERGENT across all evaluable arms → report H0 (null) as the result.
- Matched-null seed fixed (42); enrichment is a peak-set property (no donor pseudoreplication);
  L2 is the only donor-level layer and is treated as nominal/exploratory.

## Parallel: NicheNet ligand→target (orthogonal, prior-based — disclosed as [WEAK] external source)
- Priors: `refs/nichenet/{lr_network,ligand_target_matrix,weighted_networks}` (v2, Zenodo 7074291).
- Senders = niche cells; candidate ligands = expressed aging-UP secretome ligands in lr_network.
- Receivers (two pre-specified): (i) niche cells (auto/paracrine), (ii) NSC+Neuroblast (neurogenic readout).
- Geneset of interest = receiver aging DE (padj<0.1); background = expressed genes.
- Metric = NicheNet ligand activity (Pearson corr of ligand-target potential vs the response), ranked
  vs the full ligand distribution; report top 15 + predicted targets among the aging DE.
- Framing (locked): prior-based, **hypothesis-generating, not causal**; inherits NicheNet prior
  limitations (disclose in Limitations).

## Outputs
`results/regulatory/tf_motif_enrichment.csv`, `tf_convergence_summary.csv`, `nichenet_ligand_activity.csv`,
and a `PHASE_REPORT.md` with the CONVERGENT/SUGGESTIVE table or the honest null.
