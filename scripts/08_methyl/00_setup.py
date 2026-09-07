#!/usr/bin/env python
"""M1 methylome — shared config + download manifest / donor table builder.

Dataset: GSE299139 (snm3C DNA-methylation arm of the aging human hippocampus study).
Primary publication: Zemke et al. 2024, "Epigenetic and 3D genome reprogramming during the aging
  of human hippocampus", PMID 39463924 (bioRxiv doi:10.1101/2024.10.14.618338).
Verified: !Series_pubmed_id = 39463924 (GEO SOFT), 2026-06-05. Assembly: hg38 / GRCh38.
Pipeline (authors'): YAP cemba-data 1.6.9, bismark 0.20, methylome profiles by ALLCools 1.0.23.
Data form: 503 per-(donor x celltype) ALLC files (.allc.tsv.gz, 7-col: chrom,pos,strand,context,
  mc_count,cov,methylated) + 40 per-donor 3C contacts + per-nucleus metadata CSV (22,239 nuclei).

This module: (a) holds the shared M1 config (paths, seed, celltype map, region params);
(b) parses the GEO RAW filelist into a per-(donor x celltype) ALLC download manifest with per-GSM
FTP URLs; (c) parses the per-nucleus metadata into a donor table (age/sex/group + per-class nucleus
counts). Donor = unit. mCG primary. Endo built for the record but flagged NOT_EVALUABLE for the
donor-level age contrast (~2-8 nuclei/donor). seed 42.

Run: `python scripts/08_methyl/00_setup.py`  ->  results/methyl/{donor_table,allc_manifest}.csv
"""
import os
import numpy as np
import pandas as pd

SEED = 42
PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
RAW = f"{PROJ}/raw/GSE299139_snm3C"
RES = f"{PROJ}/results/methyl"
PD = f"{PROJ}/processed/per_dataset"
os.makedirs(RES, exist_ok=True)

# GEO per-GSM supplementary FTP path (all samples are GSM9034317..GSM9034356 -> folder GSM9034nnn)
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9034nnn/{gsm}/suppl/{name}"

# class-level aggregation: target niche celltype -> set of deposited per-cluster ALLC suffixes.
# (the deposit splits Astro/Micro into states; we sum mc+cov per CpG to the RNA-arm 'class' level.)
CLASS_MAP = {
    "Astro": {"Astro1", "Astro2"},
    "Micro": {"Micro1", "Micro2"},
    "Oligo": {"Oligo"},
    "OPC": {"OPC"},
    "Endo": {"Endo"},          # built for record; NOT_EVALUABLE for donor-level age contrast
}
EVALUABLE = ["Astro", "Micro", "Oligo", "OPC"]   # locked in RESOLUTION.md (Endo too sparse)

# region parameters (pre-specified, RESOLUTION.md / plan M1)
PROM = 1000          # promoter = TSS +/- 1 kb (primary); 2 kb sensitivity built separately
PROM_SENS = 2000
PROX = 2000          # proximal-regulatory: ATAC peaks within gene-body +/- 2 kb (reuse A2 link rule)
PROX_SENS = 10000
MIN_CG_COV = 20      # per-region coverage floor: >= this many CpG basecalls per donor x celltype
NPERM = 2000


def _parse_allc_name(name):
    """'GSM9034317_hc11_Astro1.allc.tsv.gz' -> (GSM9034317, hc11, Astro1)."""
    base = name[:-len(".allc.tsv.gz")]
    parts = base.split("_")
    return parts[0], "_".join(parts[1:-1]), parts[-1]


def load_manifest():
    """Parse RAW filelist.txt into the ALLC download manifest (target glial+Endo classes only)."""
    fl = pd.read_csv(f"{RAW}/filelist.txt", sep="\t")
    fl = fl[fl["Name"].str.endswith(".allc.tsv.gz")].copy()
    recs = []
    for _, r in fl.iterrows():
        gsm, donor, src = _parse_allc_name(r["Name"])
        target = next((t for t, s in CLASS_MAP.items() if src in s), None)
        if target is None:
            continue
        recs.append(dict(donor=donor, gsm=gsm, src_celltype=src, target_class=target,
                         size_bytes=int(r["Size"]), name=r["Name"],
                         url=FTP.format(gsm=gsm, name=r["Name"])))
    m = pd.DataFrame(recs).sort_values(["target_class", "donor", "src_celltype"]).reset_index(drop=True)
    return m


def load_donor_table():
    """Per-donor age/sex/group + per-class nucleus counts from the per-nucleus metadata CSV."""
    cols = ["donor", "age", "group", "gender", "cluster", "class", "mCGFrac", "mCGCov", "mCHFrac"]
    md = pd.read_csv(f"{RAW}/GSE299139_human_aging_snm3C_metadata.csv.gz", usecols=cols)
    don = (md.groupby("donor")
             .agg(age=("age", "first"), sex=("gender", "first"), group=("group", "first"),
                  n_nuclei=("class", "size"))
             .reset_index())
    # per-class nucleus counts per donor (wide)
    cls = (md.groupby(["donor", "class"]).size().unstack(fill_value=0)
             .rename(columns=lambda c: f"n_{c}"))
    don = don.merge(cls, on="donor", how="left")
    don["age_grp"] = np.where(don.age < 40, "young", np.where(don.age >= 60, "old", "mid"))
    return don.sort_values("age").reset_index(drop=True), md


if __name__ == "__main__":
    man = load_manifest()
    don, md = load_donor_table()
    man.to_csv(f"{RES}/allc_manifest.csv", index=False)
    don.to_csv(f"{RES}/donor_table.csv", index=False)

    # ---- sanity checks vs the locked feasibility findings ----
    assert don.donor.nunique() == 40, f"expected 40 donors, got {don.donor.nunique()}"
    cls_tot = md["class"].value_counts()
    print("== per-nucleus class tally (expect Oligo 13188, Micro 2604, Astro 2497, OPC 1118, Endo_VLMC 240):")
    print(cls_tot.to_string())
    print(f"\n== donors: {don.donor.nunique()} | age {don.age.min()}-{don.age.max()} | "
          f"young<40 {sum(don.age_grp=='young')} / mid {sum(don.age_grp=='mid')} / old>=60 {sum(don.age_grp=='old')} | "
          f"sex {dict(don.sex.value_counts())}")

    print("\n== ALLC manifest: per target_class, donors covered + GB to download ==")
    g = man.groupby("target_class").agg(n_files=("name", "size"),
                                        n_donors=("donor", "nunique"),
                                        GB=("size_bytes", lambda s: round(s.sum()/1e9, 1)))
    print(g.to_string())
    print(f"\nTOTAL files {len(man)} | total {round(man.size_bytes.sum()/1e9,1)} GB | "
          f"EVALUABLE target download {round(man[man.target_class.isin(EVALUABLE)].size_bytes.sum()/1e9,1)} GB")

    # per (donor x evaluable class) nuclei -> coverage expectation
    print("\n== nuclei per donor for evaluable classes (median / min / max) ==")
    for t in EVALUABLE:
        col = f"n_{t}"
        if col in don:
            v = don[col][don[col] > 0]
            print(f"  {t:6s} donors_with_cells={ (don[col]>0).sum() }  nuclei/donor med={int(v.median())} "
                  f"min={int(v.min())} max={int(v.max())}")
    print(f"\nwrote {RES}/allc_manifest.csv  and  {RES}/donor_table.csv")
