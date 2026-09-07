#!/usr/bin/env python
"""Build a secretome gene set from HPA + UniProt (+ Matrisome if present).

- union  = gene appears in >=1 source (the FILTER applied to genome-wide DE results)
- core   = gene supported by >=2 sources (high-confidence flag)
Adds rough functional categories (pattern + curated anchors) for the
category-resolved analysis (cytokine / growth factor / neuropeptide / ECM).
NOTE: categories are regex/keyword-defined, NOT database memberships. The ECM
category is therefore labelled `ecm_matrisome_like` (a symbol-pattern proxy);
it is NOT the curated MatrisomeDB unless Hs_Matrisome_Masterlist.xlsx is present
(it adds a source for set MEMBERSHIP only, never the category label). Manuscript
text must say "regex-defined ECM/matrisome-like", not "matrisome".

Outputs: refs/secretome_union.csv, refs/secretome_core.csv
Run:  PROJ=<projdir> /home/neurofuture/miniforge3/envs/bio/bin/python 02_build_secretome.py
"""
import os, re
import pandas as pd
from collections import defaultdict

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
REFS = os.path.join(PROJ, "refs")
sources = {}

# --- HPA predicted secreted ---
hpa = pd.read_csv(os.path.join(REFS, "hpa_secretome.tsv"), sep="\t")
sources["HPA"] = {str(g).strip() for g in hpa["Gene"].dropna() if str(g).strip() not in ("", "nan")}

# --- UniProt reviewed secreted (primary symbol = first token of Gene Names) ---
uni = pd.read_csv(os.path.join(REFS, "uniprot_secreted_human.tsv"), sep="\t")
uni_primary = set()
for names in uni["Gene Names"].dropna().astype(str):
    toks = names.split()
    if toks:
        uni_primary.add(toks[0].strip())
sources["UniProt"] = {g for g in uni_primary if g}

# --- Matrisome (optional, if the masterlist xlsx is present) ---
matri_path = os.path.join(REFS, "Hs_Matrisome_Masterlist.xlsx")
if os.path.exists(matri_path):
    m = pd.read_excel(matri_path)
    col = next((c for c in m.columns if "gene symbol" in str(c).lower()), None)
    if col:
        sources["Matrisome"] = {str(g).strip() for g in m[col].dropna()}

# --- membership ---
counts, src_map = defaultdict(int), defaultdict(list)
for sname, gset in sources.items():
    for g in gset:
        counts[g] += 1
        src_map[g].append(sname)
all_genes = sorted(counts)
core = {g for g in all_genes if counts[g] >= 2}

# --- rough functional categories ---
GF = {"NGF","BDNF","NTF3","NTF4","GDNF","NRTN","ARTN","PSPN","CNTF","IGF1","IGF2",
      "VEGFA","VEGFB","VEGFC","PGF","HGF","EGF","TGFA","PDGFA","PDGFB","PDGFC","PDGFD"}
NEUROPEP = {"NPY","SST","VIP","CRH","CCK","PENK","PDYN","PNOC","TAC1","TAC3","GAL",
            "CARTPT","ADCYAP1","GRP","NMB","NTS","AGRP","POMC","OXT","AVP","NMU","NPW","PYY"}
def categorize(g):
    if re.match(r"^(IL\d|CXCL|CCL|CX3CL|XCL|TNF|IFNG|IFNB|IFNA|LTA|LTB|CSF1|CSF2|CSF3)", g): return "cytokine_chemokine"
    if re.match(r"^(FGF|IGF|VEGF|PDGF|TGFB|BMP|GDF|WNT|NGF|NTF|NRG|INHB|INHA|EGF)", g) or g in GF: return "growth_factor"
    if g in NEUROPEP: return "neuropeptide"
    if re.match(r"^(COL\d|LAMA|LAMB|LAMC|FN1|FBN|ELN|MMP|ADAMTS|TIMP|SPARC|THBS|TNC|TNR|POSTN|LOX|SERPIN|HSPG2|NID|VCAN|ACAN|BGN|DCN|LUM|AGRN|MGP)", g): return "ecm_matrisome_like"
    return "other_secreted"

df = pd.DataFrame([{"gene": g, "n_sources": counts[g], "sources": "|".join(sorted(src_map[g])),
                    "is_core": g in core, "category": categorize(g)} for g in all_genes])
df = df.sort_values(["is_core","n_sources","gene"], ascending=[False,False,True])
df.to_csv(os.path.join(REFS, "secretome_union.csv"), index=False)
df[df.is_core].to_csv(os.path.join(REFS, "secretome_core.csv"), index=False)

print("Sources:", {k: len(v) for k, v in sources.items()})
print(f"Union genes: {len(all_genes)}   Core (>=2 sources): {len(core)}")
print("\nUnion category counts:\n", df["category"].value_counts().to_string())
print("\nSanity - known aging/niche-secreted markers:")
for g in ["CHI3L1","CCL2","GDF15","CLU","APOE","SPP1","C3","C1QA","SERPINA3","IGF1",
          "BDNF","VEGFA","TGFB1","NPY","SST","B2M","CST3","LGALS3","TIMP1","IGFBP5"]:
    print(f"  {g:10s} in_union={g in counts!s:5s} core={g in core!s:5s} "
          f"cat={categorize(g):18s} src={'|'.join(sorted(src_map.get(g,[]))) or '-'}")
