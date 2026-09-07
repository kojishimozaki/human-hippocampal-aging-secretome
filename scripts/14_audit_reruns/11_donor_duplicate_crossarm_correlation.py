#!/usr/bin/env python
"""[F-8-025 / F-1-021] Second attempt: are the three suspected duplicate pairs outliers in
cross-arm donor-to-donor expression similarity?

WHY A SECOND ATTEMPT
--------------------
Script 09 ran the pre-registered gene-deletion-polymorphism test and it came back
UNINFORMATIVE: all seven markers (GSTM1, UGT2B17, HLA-DRB5, RHD, FCGR3B, NAT2, HLA-DRB1) are
expressed in 0/39 donors in these hippocampal nuclei, so all 738 unrelated pairs "agree" on all
seven. The pre-registered calibration is what exposed that -- the test has no power here, and
PRESPEC.md §2-3 says an uninformative calibration means STILL-OPEN, not agreement.

This script tries the one remaining local angle. Each candidate pair straddles the WH/DG arm
boundary AND the NovaSeq6000/NovaSeqX boundary:

    A  oid  5 (WH, NovaSeq6000) / 23 (DG, NovaSeqX)
    B  oid  6 (WH, NovaSeq6000) / 29 (DG, NovaSeqX)
    C  oid 14 (WH, NovaSeq6000) / 25 (DG, NovaSeqX)

That geometry is what makes a similarity test possible at all: batch pushes cross-arm pairs
APART, so if two records are the same person their similarity should stand out ABOVE the other
cross-arm pairs despite the batch difference. The comparison set is therefore restricted to
other cross-arm pairs, never to within-arm pairs.

STILL ASYMMETRIC. An outlier is suggestive, not proof: two donors can be similar for reasons
other than identity. A NON-outlier is likewise not proof of distinctness. Whatever comes out,
the verdict cannot be CLOSED-OK; the decisive test remains genotype concordance from SRA
FASTQ/BAM, which the coordinate-only fragment deposit cannot supply.

Pre-registered here, before running (PRESPEC.md amendment 2):
  * similarity  = Spearman rho between donor pseudobulk log-CP10K over the 2,000 most variable
                  genes, computed per cell type and then averaged over cell types with >=10
                  nuclei in BOTH donors
  * comparison  = all other cross-arm donor pairs in the same Group
  * statistic   = z of the candidate pair within that comparison distribution
  * SUGGESTIVE  = z >= 2 for a candidate pair
  * calibration = report how many non-candidate cross-arm pairs also reach z >= 2

INPUTS   processed/per_dataset/GSE268609_anchor.h5ad
         results/de/gse268609_donor_covariates.csv   (oid -> SampleNumber/arm/seqbatch/Group)
OUTPUT   results/audit_reruns/donor_duplicate_crossarm_similarity.csv
RUNTIME  ~5 min.
"""
import itertools
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
PD = f"{PROJ}/processed/per_dataset"
OUT = f"{PROJ}/results/audit_reruns"
os.makedirs(OUT, exist_ok=True)

N_HVG = 2000
MIN_NUCLEI = 10
CANDIDATES = {"A": (5, 23), "B": (6, 29), "C": (14, 25)}

import scanpy as sc  # noqa: E402

cov = pd.read_csv(f"{PROJ}/results/de/gse268609_donor_covariates.csv")
cov["SampleNumber"] = cov.SampleNumber.astype(str)
oid2sample = dict(zip(cov.oid.astype(int), cov.SampleNumber))
sample2oid = {v: k for k, v in oid2sample.items()}
arm_of = dict(zip(cov.SampleNumber, cov.arm))
grp_of = dict(zip(cov.SampleNumber, cov.Group))

A = sc.read_h5ad(f"{PD}/GSE268609_anchor.h5ad")
don = A.obs.donor_id.astype(str).values
ct = A.obs.celltype_l1.astype(str).values
lib = np.asarray(A.layers["counts"].sum(1)).ravel()
X = A.layers["counts"].tocsr()
donors = [d for d in pd.unique(don) if d in arm_of]
print(f"donors with covariates: {len(donors)}", flush=True)

celltypes = [c for c in pd.unique(ct) if (ct == c).sum() >= 500]
print(f"cell types used: {celltypes}", flush=True)

# per (celltype, donor) pseudobulk log-CP10K, restricted to that cell type's HVGs
prof = {}
for c in celltypes:
    m_ct = ct == c
    rows, keep = [], []
    for d in donors:
        m = m_ct & (don == d)
        if m.sum() < MIN_NUCLEI:
            continue
        v = np.asarray(X[m].sum(0)).ravel()
        tot = lib[m].sum()
        rows.append(v / max(tot, 1) * 1e4)
        keep.append(d)
    if len(keep) < 4:
        continue
    P = np.log1p(np.vstack(rows))
    hv = np.argsort(P.var(0))[::-1][:N_HVG]
    prof[c] = (keep, P[:, hv])
print(f"profiles built for {len(prof)} cell types", flush=True)


def similarity(da, db):
    vals = []
    for c, (keep, P) in prof.items():
        if da in keep and db in keep:
            ia, ib = keep.index(da), keep.index(db)
            vals.append(spearmanr(P[ia], P[ib]).statistic)
    return (float(np.mean(vals)), len(vals)) if vals else (np.nan, 0)


rows = []
for da, db in itertools.combinations(donors, 2):
    if arm_of[da] == arm_of[db] or grp_of[da] != grp_of[db]:
        continue                                  # cross-arm, same Group only
    s, n = similarity(da, db)
    if not np.isfinite(s):
        continue
    oa, ob = sample2oid.get(da, -1), sample2oid.get(db, -1)
    pair = next((k for k, v in CANDIDATES.items() if set(v) == {oa, ob}), "")
    rows.append(dict(donor_a=da, donor_b=db, oid_a=oa, oid_b=ob, Group=grp_of[da],
                     n_celltypes=n, mean_spearman=round(s, 5), candidate_pair=pair))

R = pd.DataFrame(rows)
if not len(R):
    raise SystemExit("no cross-arm same-Group pairs found")
bg = R[R.candidate_pair == ""]
mu, sd = bg.mean_spearman.mean(), bg.mean_spearman.std(ddof=1)
R["z_vs_crossarm_background"] = ((R.mean_spearman - mu) / sd).round(3)
R = R.sort_values("z_vs_crossarm_background", ascending=False)
R.to_csv(f"{OUT}/donor_duplicate_crossarm_similarity.csv", index=False)

print(f"\ncross-arm same-Group background: n={len(bg)} pairs, "
      f"mean rho={mu:.4f} sd={sd:.4f}")
print(f"non-candidate pairs reaching z>=2: {int((bg.mean_spearman - mu).div(sd).ge(2).sum())}/{len(bg)}")
print("\ncandidate pairs:")
for name, (a, b) in CANDIDATES.items():
    r = R[R.candidate_pair == name]
    if not len(r):
        print(f"  {name}: oid {a}/{b} not both present as cross-arm same-Group -> untestable")
        continue
    r = r.iloc[0]
    verdict = "SUGGESTIVE (z>=2)" if r.z_vs_crossarm_background >= 2 else "not an outlier"
    print(f"  {name}: oid {a} vs {b}  mean rho={r.mean_spearman:.4f} over {r.n_celltypes} cell types  "
          f"z={r.z_vs_crossarm_background:+.2f}  -> {verdict}")
print("\ntop 5 cross-arm pairs overall (candidates marked):")
print(R.head(5)[["oid_a", "oid_b", "Group", "mean_spearman",
                 "z_vs_crossarm_background", "candidate_pair"]].to_string(index=False))
print(f"\nwrote {OUT}/donor_duplicate_crossarm_similarity.csv")
