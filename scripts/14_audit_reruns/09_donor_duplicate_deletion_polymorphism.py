#!/usr/bin/env python
"""[F-8-025 / F-1-021] Can the three suspected duplicate-donor pairs be REFUTED from the data?

WHAT IS ALREADY ESTABLISHED (audit Pass 8 addendum -- not re-derived here)
-------------------------------------------------------------------------
Disouky et al. Supplementary Table 1 has 39 rows, and three pairs agree on all 11 columns
other than ID / Cohort / Hippocampal Sclerosis:

    A  orig.ident 5 / 23   Healthy Aging  78  male    PMI 10.0  BRAAK 0  LATE 0  APOE 3/3
    B  orig.ident 6 / 29   Healthy Aging  92  female  PMI  7.5  BRAAK 4  LATE 1  APOE 2/3
    C  orig.ident 14 / 25  Young Aging    33  male    PMI  5.0  BRAAK 0  LATE 0  APOE unknown

The paper's own text says 38 participants while the table has 39 rows, with exactly one
excess in HA. Each pair splits across arm (WH/DG) and sequencer (NovaSeq6000/NovaSeqX).

WHY A NEW TEST
--------------
The audit concluded that a decisive call needs genotype concordance, which needs FASTQ/BAM
from SRA: GSE268609's ATAC is a coordinate-only fragment TSV, so no alleles are available.
That is correct and is NOT attempted here.

This script runs a different, strictly LOCAL test, and it is deliberately ASYMMETRIC:

    Common gene-deletion polymorphisms are homozygous-null in a sizeable fraction of people.
    If two records are the same individual, the presence/absence of expression at such a
    locus MUST agree. A single disagreement refutes that pair. Full agreement proves nothing.

So: mismatch -> pair refuted (CLOSED-CHANGE). Full agreement -> STILL-OPEN, never CLOSED-OK.
And full agreement is only interpretable against a calibration: how often do two UNRELATED
donors also agree on all markers? That is computed here too, and reported alongside.

Markers verified present in this reference: GSTM1, UGT2B17, HLA-DRB5, RHD, FCGR3B, NAT2,
HLA-DRB1. (GSTT1 and LILRA3 are absent from the annotation and cannot be used.)
XIST / RPS4Y1 / UTY / DDX3Y are carried as a positive control: sex must agree within a pair.

Detection rule, fixed in PRESPEC.md before any number was computed:
    "expressed" = donor pseudobulk CP10K >= 1.0 AND detected in >= 3 nuclei.

INPUTS   processed/per_dataset/GSE268609_anchor.h5ad   (counts, donor_id, Group)
         processed/per_dataset/GSE268609_metadata.csv  (SampleNumber -> orig.ident)
OUTPUT   results/audit_reruns/donor_duplicate_marker_matrix.csv
         results/audit_reruns/donor_duplicate_pair_calls.csv
RUNTIME  ~5 min.
"""
import itertools
import os

import numpy as np
import pandas as pd
import scipy.sparse as sp

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
PD = f"{PROJ}/processed/per_dataset"
OUT = f"{PROJ}/results/audit_reruns"
os.makedirs(OUT, exist_ok=True)

# pre-registered detection rule -- do not change after seeing results (CLAUDE.md rule 7)
CP10K_MIN = 1.0
NUCLEI_MIN = 3

DELETION_MARKERS = ["GSTM1", "UGT2B17", "HLA-DRB5", "RHD", "FCGR3B", "NAT2", "HLA-DRB1"]
SEX_MARKERS = ["XIST", "RPS4Y1", "UTY", "DDX3Y"]

# candidate duplicate pairs, keyed by orig.ident (from the published donor table)
CANDIDATES = {"A": (5, 23), "B": (6, 29), "C": (14, 25)}

import scanpy as sc  # noqa: E402

meta = pd.read_csv(f"{PD}/GSE268609_metadata.csv", index_col=0)
oid_of_sample = (meta[["SampleNumber", "orig.ident"]].drop_duplicates()
                 .assign(SampleNumber=lambda d: d.SampleNumber.astype(str))
                 .set_index("SampleNumber")["orig.ident"].astype(int).to_dict())

A = sc.read_h5ad(f"{PD}/GSE268609_anchor.h5ad")
var = A.var_names.astype(str)
gpos = {g: i for i, g in enumerate(var)}
markers = [g for g in DELETION_MARKERS + SEX_MARKERS if g in gpos]
missing = [g for g in DELETION_MARKERS + SEX_MARKERS if g not in gpos]
if missing:
    print(f"[note] not in reference, excluded: {missing}")

X = A.layers["counts"].tocsc()[:, [gpos[g] for g in markers]].tocsr()
lib = np.asarray(A.layers["counts"].sum(1)).ravel()
don = A.obs.donor_id.astype(str).values
grp = A.obs.Group.astype(str).values

donors = pd.unique(don)
rows = []
for d in donors:
    m = don == d
    sub = X[m]
    tot = lib[m].sum()
    umi = np.asarray(sub.sum(0)).ravel()
    ncell = np.asarray((sub > 0).sum(0)).ravel()
    rec = dict(donor=d, oid=oid_of_sample.get(d, -1), Group=grp[m][0], n_nuclei=int(m.sum()),
               total_umi=float(tot))
    for j, g in enumerate(markers):
        cp = umi[j] / max(tot, 1) * 1e4
        rec[f"{g}_cp10k"] = round(float(cp), 4)
        rec[f"{g}_ncell"] = int(ncell[j])
        rec[f"{g}_expressed"] = bool(cp >= CP10K_MIN and ncell[j] >= NUCLEI_MIN)
    rows.append(rec)

M = pd.DataFrame(rows).sort_values("oid").reset_index(drop=True)
M.to_csv(f"{OUT}/donor_duplicate_marker_matrix.csv", index=False)
print(f"wrote {OUT}/donor_duplicate_marker_matrix.csv  ({len(M)} donors x {len(markers)} markers)")

del_markers = [g for g in DELETION_MARKERS if g in markers]
sex_markers = [g for g in SEX_MARKERS if g in markers]
print("\nper-marker null-rate across donors (how informative each marker is):")
for g in del_markers:
    frac = float(M[f"{g}_expressed"].mean())
    print(f"  {g:10s} expressed in {M[f'{g}_expressed'].sum():2d}/{len(M)} donors "
          f"({frac:.2f})   {'informative' if 0.05 < frac < 0.95 else 'UNINFORMATIVE (no variation)'}")


def agree(a, b, cols):
    va = M.loc[M.oid == a, [f"{g}_expressed" for g in cols]].iloc[0].values
    vb = M.loc[M.oid == b, [f"{g}_expressed" for g in cols]].iloc[0].values
    return int((va == vb).sum()), len(cols), [c for c, x, y in zip(cols, va, vb) if x != y]


# ---- calibration: unrelated pairs, matched on Group and on arm where possible -----------
oids = M.oid.tolist()
cal = []
cand_flat = {frozenset(p) for p in CANDIDATES.values()}
for a, b in itertools.combinations(oids, 2):
    if frozenset((a, b)) in cand_flat:
        continue
    ga = M.loc[M.oid == a, "Group"].iloc[0]
    gb = M.loc[M.oid == b, "Group"].iloc[0]
    n_ok, n_tot, diff = agree(a, b, del_markers)
    cal.append(dict(oid_a=a, oid_b=b, same_group=bool(ga == gb),
                    n_agree=n_ok, n_markers=n_tot, full_agreement=bool(n_ok == n_tot)))
C = pd.DataFrame(cal)
rate_all = float(C.full_agreement.mean())
rate_same = float(C[C.same_group].full_agreement.mean()) if C.same_group.any() else np.nan
print(f"\ncalibration on {len(C)} unrelated pairs: full agreement on all "
      f"{len(del_markers)} deletion markers in {C.full_agreement.sum()}/{len(C)} "
      f"({rate_all:.3f});  same-Group pairs {rate_same:.3f}")

# ---- the three candidate pairs ----------------------------------------------------------
out = []
print("\ncandidate pairs:")
for name, (a, b) in CANDIDATES.items():
    if a not in oids or b not in oids:
        print(f"  {name}: oid {a}/{b} not both in the analysed donors -> untestable")
        out.append(dict(pair=name, oid_a=a, oid_b=b, testable=False))
        continue
    n_ok, n_tot, diff = agree(a, b, del_markers)
    s_ok, s_tot, sdiff = agree(a, b, sex_markers)
    refuted = n_ok < n_tot
    out.append(dict(pair=name, oid_a=a, oid_b=b, testable=True,
                    n_agree=n_ok, n_markers=n_tot, disagreeing=";".join(diff) or "",
                    sex_agree=s_ok, sex_markers=s_tot, sex_disagreeing=";".join(sdiff) or "",
                    refuted_as_same_individual=bool(refuted),
                    calibration_full_agreement_rate=round(rate_all, 4)))
    verdict = ("REFUTED (different individuals)" if refuted
               else "not refuted -- STILL-OPEN (agreement proves nothing)")
    print(f"  {name}: oid {a} vs {b}  deletion markers {n_ok}/{n_tot}"
          f"{' [differ: ' + ','.join(diff) + ']' if diff else ''}"
          f"  sex controls {s_ok}/{s_tot}  -> {verdict}")

pd.DataFrame(out).to_csv(f"{OUT}/donor_duplicate_pair_calls.csv", index=False)
C.to_csv(f"{OUT}/donor_duplicate_calibration_pairs.csv", index=False)
print(f"\nwrote {OUT}/donor_duplicate_pair_calls.csv")
