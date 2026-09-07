#!/usr/bin/env python
"""[F-3b-012] Apply the RNA leg's ">= 10 nuclei/donor" gate to the ATAC leg and re-measure.

WHY
---
scripts/03_de/01_pseudobulk_de.py:L43 drops donors with < `min_cells` (=10) nuclei before
pseudobulk. The ATAC leg does not: scripts/04_atac/05_differential_accessibility.py:L156-L163
gates only on donors-per-group (>=4), and scripts/04_atac/01_atac_rna_concordance.py:L110 only
on >=3. Two donor x celltype units therefore enter the ATAC rank tests with almost no data --
Micro/donor 39/YA (5 nuclei) and Endo/donor 47/YA (8 nuclei) -- each carrying 1/8 of the weight
of a Mann-Whitney over 8 YA donors. Manuscript L74 states ">= 10 cells/donor" once, in the RNA
Methods, where a reader will read it as applying to the whole pipeline.

The third-party audit (run 20260821-225834) raised this as F-3b-012 and left it [UNVERIFIED]:
`compute/pass3b_03_atac_multiplet_upperbound.py` was supposed to carry a cells-per-donor
variant and never did. This script supplies it.

WHAT IS RECOMPUTED, AND HOW IT MATCHES THE PUBLISHED PIPELINE
-------------------------------------------------------------
Two arms per niche cell type, over YA+HA only:
  all_donors    every donor with >=1 nucleus of that cell type   (the published behaviour)
  min10_cells   only donors with >= MIN_CELLS (=10) nuclei       (the RNA leg's rule)

and two statistics, each built exactly as the committed script that publishes it:

  (1) Fig 4A gene-level concordance  -- 04_atac/01_atac_rna_concordance.py:L105-L125
        acc  = (summed linked-peak counts) / (donor's total nCount_Peaks) * 1e6
        lfc  = log2((mean acc in HA + 1) / (mean acc in YA + 1))
        concordant iff sign(lfc_atac) == sign(lfc_rna); scored over the 60 RNA-significant
        gene x celltype pairs.  Published value: 41/60 (Astro 21, Micro 8, Oligo 2, OPC 6, Endo 4).

  (2) A2 matched-null foreground obs -- 04_atac/05_differential_accessibility.py:L183-L204
        foreground = peaks linked to a secretome gene (gene body +/- PROM) AND *tested* in that
        cell type, oriented by that gene's RNA log2FC sign; obs = mean(lfc_atac * sign).
        `tested` = detected in >=3 donors and non-zero CPM variance, recomputed inside each arm.
        Published obs: Astro +0.1834, Micro +0.2391, Oligo -0.0151, OPC +0.2863, Endo +0.1118.

        NOTE the audit kit's pass3b_02 / pass3b_03 dropped the "AND tested" half and so did not
        reproduce the published foreground (OPC: 55 peaks vs the published 19, obs +0.033 vs
        +0.286). Restoring it is what makes the two arms comparable to the paper.

        Only the observed statistic is recomputed, not its matched-null p: the null resamples
        KNN-matched NON-linked peaks genome-wide, whose lfc would also have to be refitted under
        the gate, and that needs the 8.9 GB peak MTX. The pre-registered criterion for F-3b-012
        (audit_log/2026-08-25_p0_cheap_reruns/RESOLUTION.md, section 0-3) is stated on obs sign
        and magnitude for exactly this reason.

Inputs are the linked-peak cache (32 MB) plus the committed DA tables, so the peak MTX is never
loaded.  seed 42.  Runtime: a few minutes.

Out: results/audit_reruns/atac_mincells_gate_summary.csv     per celltype x arm
     results/audit_reruns/atac_mincells_gate_per_gene.csv    per gene x celltype x arm
     results/audit_reruns/atac_mincells_gate_donors.csv      which donors the gate removes
"""
import os, gzip
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_cache")
import numpy as np, pandas as pd, scipy.sparse as sp

SEED = 42
rng = np.random.default_rng(SEED)
PROJ = os.environ.get("PROJ", os.getcwd())
PD = f"{PROJ}/processed/per_dataset"
OUT = f"{PROJ}/results/audit_reruns"
os.makedirs(OUT, exist_ok=True)

CL = {"Astrocytes": "Astro", "Microglia": "Micro", "mOli": "Oligo", "OPCs": "OPC",
      "Endothelial": "Endo"}
NICHE = ["Astro", "Micro", "Oligo", "OPC", "Endo"]
PROM = 2000           # gene body +/- 2 kb, as 05_differential_accessibility.py:L115
MIN_CELLS = 10        # the RNA leg's rule, 01_pseudobulk_de.py:L43
MIN_DONORS_A2 = 4     # 05_differential_accessibility.py:L163
MIN_DONORS_CONC = 3   # 01_atac_rna_concordance.py:L110

# ---------------------------------------------------------------- peaks and gene links
peaks = np.array([l.strip() for l in open(f"{PD}/GSE268609_peaks_names.txt")])
parts = [p.rsplit("-", 2) for p in peaks]
pkdf = pd.DataFrame({"chr": [x[0] for x in parts], "start": [int(x[1]) for x in parts],
                     "end": [int(x[2]) for x in parts], "idx": np.arange(len(peaks))})
pk_by_chr = {c: g for c, g in pkdf.groupby("chr")}

sec = set(pd.read_csv(f"{PROJ}/refs/secretome_union.csv").gene)
sec_coord = {}
with gzip.open(f"{PROJ}/raw/GSE268609_lazarov2024/GSE268609_features.tsv.gz", "rt") as f:
    for line in f:
        q = line.rstrip().split("\t")
        if len(q) >= 6 and q[2] == "Gene Expression" and q[1] in sec:
            try:
                sec_coord[q[1]] = (q[3], int(q[4]), int(q[5]))
            except ValueError:
                pass
g2peaks = {}
for gene, (ch, gs, ge) in sec_coord.items():
    if ch not in pk_by_chr:
        continue
    s = pk_by_chr[ch]
    hit = s[(s.end >= gs - PROM) & (s.start <= ge + PROM)].idx.values
    if len(hit):
        g2peaks[gene] = hit
print(f"secretome genes with linked peaks: {len(g2peaks)}", flush=True)

# ---------------------------------------------------------------- linked-peak cache
Xl = sp.load_npz(f"{PD}/GSE268609_linkedpeaks.npz").tocsr()      # linked_peaks x cells
linked_idx = np.load(f"{PD}/GSE268609_linkedpeaks_idx.npy")
peak_local = pd.Series(np.arange(len(linked_idx)), index=linked_idx)
print(f"linked-peak cache: {Xl.shape[0]} peaks x {Xl.shape[1]} cells", flush=True)

barcodes = [l.strip() for l in open(f"{PD}/GSE268609_peaks_barcodes.txt")]
m = pd.read_csv(f"{PD}/GSE268609_metadata.csv", index_col=0, low_memory=False)
m.index = m.index.astype(str)
mr = m.reindex(pd.Index(barcodes, dtype=str))
donor_all = mr["SampleNumber"].values.astype(object)
group_all = mr["Group"].values.astype(object)
ct_all = np.array([CL.get(x, None) for x in mr["Cluster"].values], dtype=object)
ncp_all = mr["nCount_Peaks"].values.astype(float)
assert Xl.shape[1] == len(barcodes), "cache columns must align to peaks_barcodes.txt"

rna = pd.read_csv(f"{PROJ}/results/de/de_GSE268609_per_celltype.csv")
rna = rna[rna.gene.isin(sec)]

rows, gene_rows, donor_rows = [], [], []
for celltype in NICHE:
    base = (ct_all == celltype) & np.isin(group_all, ["YA", "HA"])
    if base.sum() == 0:
        continue
    ncell = pd.Series(donor_all[base]).value_counts()
    small = sorted(ncell[ncell < MIN_CELLS].index.tolist(), key=str)
    for d in ncell.index:
        donor_rows.append(dict(celltype=celltype, donor=str(d), n_cells=int(ncell[d]),
                               group=str(group_all[base][donor_all[base] == d][0]),
                               dropped_by_min10=bool(ncell[d] < MIN_CELLS)))
    print(f"\n{celltype}: {base.sum()} nuclei, {len(ncell)} donors; "
          f"< {MIN_CELLS} nuclei: {[(str(d), int(ncell[d])) for d in small] or 'none'}", flush=True)

    # the published DA table gives the *tested* peak set on all donors (G-8 reference)
    da = pd.read_csv(f"{PROJ}/results/de/atac_DA_peaks_{celltype}.csv",
                     usecols=["peak_idx", "lfc_atac", "is_sec_linked"])
    pub_lfc = pd.Series(da.lfc_atac.values, index=da.peak_idx.values)

    rna_ct = rna[rna.celltype == celltype].set_index("gene")

    for arm in ["all_donors", "min10_cells"]:
        cells = base.copy()
        if arm == "min10_cells" and small:
            cells &= ~np.isin(donor_all, small)
        dc = donor_all[cells]
        gc_grp = group_all[cells]
        ncp_c = ncp_all[cells]
        donors = pd.unique(dc)
        code = pd.Categorical(dc, categories=donors).codes
        OH = sp.csr_matrix((np.ones(len(dc)), (code, np.arange(len(dc)))),
                           shape=(len(donors), len(dc)))
        pb = np.asarray((Xl[:, cells] @ OH.T).todense())          # linked_peaks x donors
        depth = np.asarray(OH @ ncp_c).ravel() + 1.0
        dgrp = np.array([gc_grp[dc == d][0] for d in donors])
        ha, ya = dgrp == "HA", dgrp == "YA"
        cpm = pb / depth * 1e6

        # ---- (2) A2 foreground: linked AND tested, as 05:L166-L197 -------------------
        detected = (pb > 0).sum(1)
        testable = (detected >= 3) & (cpm.std(1) > 0)
        meanHA, meanYA = cpm[:, ha].mean(1), cpm[:, ya].mean(1)
        lfc_local = np.log2((meanHA + 1.0) / (meanYA + 1.0))       # 05:L71, identical formula

        a2_ok = (ha.sum() >= MIN_DONORS_A2) and (ya.sum() >= MIN_DONORS_A2)
        fg_vals, fg_pubs = [], []
        for gene, pks in g2peaks.items():
            if gene not in rna_ct.index:
                continue
            rp = rna_ct.loc[gene, "padj"]
            if not (pd.notna(rp) and rp < 0.1):
                continue
            s = float(np.sign(rna_ct.loc[gene, "log2FoldChange"]))
            for pk in pks:
                if pk not in peak_local.index:
                    continue
                loc = int(peak_local[pk])
                if not testable[loc]:
                    continue
                fg_vals.append(lfc_local[loc] * s)
                if pk in pub_lfc.index:
                    fg_pubs.append(pub_lfc[pk] * s)
        a2_obs = float(np.mean(fg_vals)) if (fg_vals and a2_ok) else np.nan
        a2_obs_published_lfc = float(np.mean(fg_pubs)) if (fg_pubs and a2_ok) else np.nan

        # ---- (1) Fig 4A gene-level concordance, as 01_atac_rna_concordance.py:L118-L126 --
        conc_ok = (ha.sum() >= MIN_DONORS_CONC) and (ya.sum() >= MIN_DONORS_CONC)
        n_hit = n_conc = 0
        if conc_ok:
            for gene, pks in g2peaks.items():
                if gene not in rna_ct.index:
                    continue
                lps = [int(peak_local[pk]) for pk in pks if pk in peak_local.index]
                if not lps:
                    continue
                acc = pb[lps].sum(0) / depth * 1e6
                if acc.sum() == 0:
                    continue
                lfc_g = np.log2((acc[ha].mean() + 1) / (acc[ya].mean() + 1))
                lfc_r = rna_ct.loc[gene, "log2FoldChange"]
                rp = rna_ct.loc[gene, "padj"]
                is_hit = bool(pd.notna(rp) and rp < 0.1)
                concordant = bool(np.sign(lfc_g) == np.sign(lfc_r))
                gene_rows.append(dict(celltype=celltype, arm=arm, gene=gene, n_peaks=len(lps),
                                      lfc_atac=round(float(lfc_g), 5),
                                      lfc_rna=round(float(lfc_r), 5), padj_rna=float(rp),
                                      rna_hit=is_hit, concord=concordant))
                if is_hit:
                    n_hit += 1
                    n_conc += int(concordant)

        rows.append(dict(
            celltype=celltype, arm=arm, n_donor=len(donors), n_YA=int(ya.sum()),
            n_HA=int(ha.sum()), n_cells=int(cells.sum()),
            donors_dropped=";".join(map(str, small)) if arm == "min10_cells" else "",
            a2_evaluable=a2_ok, n_fg_peaks=len(fg_vals), a2_obs_mean_oriented=a2_obs,
            a2_obs_using_published_lfc=a2_obs_published_lfc,
            conc_evaluable=conc_ok, n_rna_hits=n_hit, n_concordant=n_conc,
            concordance=round(n_conc / n_hit, 4) if n_hit else np.nan))
        print(f"  [{arm:11s}] donors {len(donors)} (YA {ya.sum()} / HA {ha.sum()}), "
              f"fg peaks {len(fg_vals)}, A2 obs {a2_obs:+.4f}, "
              f"concordance {n_conc}/{n_hit}", flush=True)

S = pd.DataFrame(rows)
S.to_csv(f"{OUT}/atac_mincells_gate_summary.csv", index=False)
pd.DataFrame(gene_rows).to_csv(f"{OUT}/atac_mincells_gate_per_gene.csv", index=False)
pd.DataFrame(donor_rows).to_csv(f"{OUT}/atac_mincells_gate_donors.csv", index=False)

print("\n" + S.to_string(index=False))
tot = S.groupby("arm").agg(hits=("n_rna_hits", "sum"), conc=("n_concordant", "sum"))
print("\n=== Fig 4A concordance, totalled over the five niche cell types ===")
print(tot.assign(frac=(tot.conc / tot.hits).round(4)).to_string())
print("\nG-8 REPRODUCTION GATE (all_donors arm must match the committed pipeline):")
print("  Fig 4A published 41/60 -> Astro 21, Micro 8, Oligo 2, OPC 6, Endo 4")
print("  A2 obs published      -> Astro +0.1834, Micro +0.2391, Oligo -0.0151, "
      "OPC +0.2863, Endo +0.1118")
print(f"\nwrote {OUT}/atac_mincells_gate_{{summary,per_gene,donors}}.csv")
