#!/usr/bin/env python
"""M1 methylome — build candidate regions (the donor x celltype mCG counting target).

Region classes (pre-specified, RESOLUTION.md / plan M1), each with secretome FOREGROUND + matched
non-secretome BACKGROUND (matching itself happens at the enrichment stage, within region-class):
  - promoter      : TSS +/- 1 kb (strand-aware; TSS = gene start[+]/end[-]).  (2 kb = sensitivity, sep.)
  - gene_body     : [gene_start, gene_end].
  - proximal_peak : GSE268609 ATAC peaks overlapping gene-body +/- 2 kb (the A2 link rule; reuse peaks).
Foreground = secretome-universe genes (refs/secretome_union.csv); background = non-secretome genes /
non-secretome-linked peaks. Direction labels (UP/DOWN) are celltype-specific and joined at test time
from the per-celltype DE, NOT baked here.

hg38, chr1-22,X only. Gene coords + strand from refs/gencode_v44_genes.tsv.gz (canonical entry per
symbol). seed 42 (no randomness here). Out: results/methyl/regions_meta.csv (+ gene_tss.csv for §dist).
"""
import os
import numpy as np
import pandas as pd

SEED = 42
PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
RES = f"{PROJ}/results/methyl"
PD = f"{PROJ}/processed/per_dataset"
PROM = 1000
PROX = 2000
MAIN = [f"chr{c}" for c in list(range(1, 23)) + ["X"]]
os.makedirs(RES, exist_ok=True)


def canonical_genes():
    """One GENCODE entry per symbol on main chroms: prefer protein_coding, then longest span."""
    g = pd.read_csv(f"{PROJ}/refs/gencode_v44_genes.tsv.gz", sep="\t")
    g = g[g.chrom.isin(MAIN)].copy()
    g["span"] = g.end - g.start
    g["pc"] = (g.gene_type == "protein_coding").astype(int)
    g = (g.sort_values(["symbol", "pc", "span"], ascending=[True, False, False])
           .drop_duplicates("symbol", keep="first").reset_index(drop=True))
    g["tss"] = np.where(g.strand == "+", g.start, g.end)
    return g


def load_peaks():
    names = np.array([l.strip() for l in open(f"{PD}/GSE268609_peaks_names.txt")])
    parts = [p.rsplit("-", 2) for p in names]
    df = pd.DataFrame({"peak": names,
                       "chrom": [x[0] for x in parts],
                       "start": [int(x[1]) for x in parts],
                       "end": [int(x[2]) for x in parts],
                       "pidx": np.arange(len(names))})
    return df[df.chrom.isin(MAIN)].reset_index(drop=True)


def link_peaks_to_genes(peaks, genes):
    """peak -> set of gene symbols whose [start-2k, end+2k] window overlaps the peak (per chrom)."""
    pk_gene = {i: [] for i in peaks.pidx}
    for ch, gsub in genes.groupby("chrom"):
        psub = peaks[peaks.chrom == ch]
        if psub.empty:
            continue
        pe = psub.end.values; ps = psub.start.values; pidx = psub.pidx.values
        order = np.argsort(ps); ps_s = ps[order]; pidx_s = pidx[order]; pe_s = pe[order]
        for _, gr in gsub.iterrows():
            lo, hi = gr.start - PROX, gr.end + PROX
            # peaks overlapping [lo,hi]: start<=hi and end>=lo
            j0 = np.searchsorted(ps_s, hi, side="right")
            cand = np.where(pe_s[:j0] >= lo)[0]
            for k in cand:
                pk_gene[pidx_s[k]].append(gr.symbol)
    return pk_gene


def main():
    genes = canonical_genes()
    sec = set(pd.read_csv(f"{PROJ}/refs/secretome_union.csv").gene) & set(genes.symbol)
    genes["is_sec"] = genes.symbol.isin(sec)
    genes[["symbol", "chrom", "tss", "strand", "start", "end", "is_sec"]].to_csv(f"{RES}/gene_tss.csv", index=False)
    print(f"canonical genes (main chrom): {len(genes)} | secretome matched: {genes.is_sec.sum()}")

    rows = []
    # promoter + gene_body, foreground (secretome) and background (non-secretome)
    for _, g in genes.iterrows():
        ps, pe = g.tss - PROM, g.tss + PROM
        rows.append((f"promoter:{g.symbol}", g.chrom, max(1, ps), pe, g.strand, "promoter", g.symbol, g.is_sec))
        rows.append((f"genebody:{g.symbol}", g.chrom, g.start, g.end, g.strand, "gene_body", g.symbol, g.is_sec))

    # proximal_peak: link peaks to genes; classify secretome-linked vs non-secretome-linked
    peaks = load_peaks()
    pk_gene = link_peaks_to_genes(peaks, genes)
    secset = sec
    sym2tss = dict(zip(genes.symbol, genes.tss))                    # precompute (avoid in-loop set_index)
    pcoord = {int(p): (c, int(s), int(e)) for p, c, s, e
              in zip(peaks.pidx, peaks.chrom, peaks.start, peaks.end)}
    for pidx, gl in pk_gene.items():
        if not gl:
            continue                                   # unlinked peak: not in the proximal_peak universe
        ch, s, e = pcoord[int(pidx)]
        sec_links = [x for x in gl if x in secset]
        is_sec = len(sec_links) > 0
        # representative gene = nearest-TSS secretome gene if foreground, else nearest gene overall
        cand = sec_links if is_sec else gl
        gtss = np.array([sym2tss[x] for x in cand])
        pmid = (s + e) // 2
        rep = cand[int(np.argmin(np.abs(gtss - pmid)))]
        rows.append((f"peak:{pidx}", ch, s, e, ".", "proximal_peak", rep, is_sec))

    reg = pd.DataFrame(rows, columns=["region_id", "chrom", "start", "end", "strand",
                                      "region_class", "gene", "is_secretome"])
    reg = reg[reg.chrom.isin(MAIN)].copy()
    reg["start"] = reg.start.clip(lower=1).astype(int)
    reg["end"] = reg.end.astype(int)
    reg["length"] = reg.end - reg.start + 1
    reg = reg.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
    reg.to_csv(f"{RES}/regions_meta.csv", index=False)

    print("\n== regions per class (foreground secretome / background) ==")
    tab = reg.groupby(["region_class", "is_secretome"]).size().unstack(fill_value=0)
    tab.columns = ["background", "secretome"]
    print(tab.to_string())
    print(f"\nTOTAL regions: {len(reg)}  (counting target)")
    print(f"wrote {RES}/regions_meta.csv  and  {RES}/gene_tss.csv")


if __name__ == "__main__":
    main()
