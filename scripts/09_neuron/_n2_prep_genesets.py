#!/usr/bin/env python
"""N2 prep: cache pathway libraries locally (reproducibility) + build the pre-registered targeted
'niche-response' gene sets. Run once before the N2 analyses. Records download provenance.

Caches (Enrichr via gseapy):  refs/enrichr/{MSigDB_Hallmark_2020,Reactome_2022,GO_Biological_Process_2023}.gmt
Builds (pre-registered N2.3): refs/neuron_response_sets.gmt   (objective Hallmark subsets + SenMayo +
                                                               cited curated ISR/senescence/receptor sets)
Provenance:                   refs/enrichr/PROVENANCE.txt
"""
import os, datetime
import gseapy
P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
OUT = f"{P}/refs/enrichr"; os.makedirs(OUT, exist_ok=True)
LIBS = ["MSigDB_Hallmark_2020", "Reactome_2022", "GO_Biological_Process_2023"]

def write_gmt(path, d):
    with open(path, "w") as f:
        for name, genes in d.items():
            f.write("\t".join([name, "na"] + list(genes)) + "\n")

libdicts = {}
prov = [f"# N2 pathway-library cache  (download date: {datetime.date.today().isoformat()})",
        f"# gseapy {gseapy.__version__}; source Enrichr (maayanlab.cloud). Versions are encoded in the library names."]
for lib in LIBS:
    d = gseapy.get_library(name=lib, organism="human")
    libdicts[lib] = d
    write_gmt(f"{OUT}/{lib}.gmt", d)
    prov.append(f"{lib}: {len(d)} sets -> refs/enrichr/{lib}.gmt")
    print(f"cached {lib}: {len(d)} sets")

# ---- pre-registered targeted 'niche-response' sets (RESOLUTION 2026-06-18) ----
HALL = libdicts["MSigDB_Hallmark_2020"]
def hall(k):
    # Hallmark_2020 keys look like 'TNF-alpha Signaling via NF-kB'; match by normalised token
    norm = lambda s: s.upper().replace("-", "").replace(" ", "").replace("_", "")
    for key in HALL:
        if norm(key) == norm(k):
            return list(HALL[key])
    hits = [key for key in HALL if norm(k) in norm(key)]
    if hits:
        return list(HALL[hits[0]])
    raise KeyError(f"Hallmark set not found: {k} (avail e.g. {list(HALL)[:3]})")

def read_list(fn):
    return [l.strip() for l in open(f"{P}/refs/{fn}") if l.strip() and not l.startswith("#")]

targeted = {
    # objective (Hallmark)
    "NFkB_TNFA":      hall("TNF-alpha Signaling via NF-kB"),
    "IFN_gamma":      hall("Interferon Gamma Response"),
    "IFN_alpha":      hall("Interferon Alpha Response"),
    "Inflammatory":   hall("Inflammatory Response"),
    "Complement":     hall("Complement"),
    "UPR":            hall("Unfolded Protein Response"),
    # objective (local published)
    "SenMayo":        read_list("senmayo.txt"),
    # supporting curated (cited in RESOLUTION)
    "ISR_ATF4":       ["ATF4","ATF3","DDIT3","TRIB3","ASNS","CHAC1","PPP1R15A","DDIT4","SLC7A11",
                       "SESN2","WARS1","PSAT1","MTHFD2"],
    "Senescence_core":["CDKN1A","CDKN2A","TP53","GLB1","SERPINE1","LMNB1","MKI67","IGFBP3","GADD45A","CCND1"],
    "Cyt_Compl_Receptors":["IL6R","IL6ST","IL1R1","IL1RAP","TNFRSF1A","TNFRSF1B","IFNGR1","IFNGR2",
                       "IFNAR1","IFNAR2","IL10RA","IL10RB","CSF1R","OSMR","LIFR","C3AR1","C5AR1","C5AR2",
                       "CR1","CXCR4","CCR5","IL4R","IL13RA1"],
}
write_gmt(f"{P}/refs/neuron_response_sets.gmt", targeted)
prov.append(f"\n# targeted niche-response sets -> refs/neuron_response_sets.gmt")
for k, v in targeted.items():
    prov.append(f"  {k}: {len(v)} genes")
    print(f"targeted {k}: {len(v)} genes")

with open(f"{OUT}/PROVENANCE.txt", "w") as f:
    f.write("\n".join(prov) + "\n")
print(f"\nwrote {OUT}/PROVENANCE.txt and refs/neuron_response_sets.gmt")
