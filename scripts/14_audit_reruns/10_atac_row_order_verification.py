#!/usr/bin/env python
"""[F-3b-009 / B-8] Is peaks_counts.mtx really in peaks_names.txt order?

WHAT IS ALREADY SETTLED (read-only checks, no computation)
----------------------------------------------------------
* `processed/per_dataset/GSE268609_peaks_counts.mtx` header is 369,393 x 153,530, and
  `peaks_names.txt` has 369,393 lines while `peaks_barcodes.txt` has 153,530 -- dimensions agree.
* After normalising the separator, `peaks_names.txt` is row-for-row IDENTICAL to
  `raw/GSE268609_lazarov2024/GSE268609_recall_peaks.bed.gz` (369,393/369,393). The provenance of
  the peak NAMES is therefore settled.
* The deposit carries TWO disjoint peak sets: `features.tsv.gz` lists 141,330 peaks (the RNA
  MTX's 177,931 rows = 36,601 genes + those 141,330 peaks), while `recall_peaks.bed.gz` has
  369,393. Zero coordinates are shared. Confusing the two would silently break any ATAC analysis.

WHAT IS NOT SETTLED
-------------------
Whether the ROWS of `peaks_counts.mtx` were written in that same order. The Seurat RDS -> MTX
export script is absent from the repository (that is B-8), so this cannot be followed in code.
It can, however, be checked against the data.

THE TEST (pre-registered in audit_log/2026-08-26_remaining_still_open/PRESPEC.md §1-3)
-------------------------------------------------------------------------------------
Recount one chromosome's fragments per peak straight from the deposited fragment file and
compare with the matrix's row sums for the same rows.
    statistic : Spearman rho over that chromosome's peaks
    CLOSED-OK : rho >= 0.90 AND above the 95th percentile of a null in which that row
                block is randomly permuted (100 shuffles)

AMENDMENT 1 (declared before running, results unseen): PRESPEC.md named chr22. The target is
chr1 instead, for two reasons that do not depend on any result: chr1 is the FIRST block in the
coordinate-sorted fragment file, so awk can exit after it and read ~a tenth of the 82 GB rather
than all of it; and chr1 carries ~35k peaks against chr22's 6,331, which makes the test more
discriminating, not less. The statistic, the 0.90 threshold and the shuffled null are unchanged.
A stated limit, fixed before running: the fragment file covers all 366,175 deposited barcodes
while the matrix keeps the 153,530 that passed QC, so counts are NOT expected to match. This
tests ORDER, not magnitude.

INPUTS   raw/GSE268609_lazarov2024/GSE268609_Aggregated_atac_fragments.tsv.gz  (82 GB, streamed)
         processed/per_dataset/GSE268609_peaks_counts.mtx                      (9.5 GB, streamed)
         processed/per_dataset/GSE268609_peaks_names.txt
OUTPUT   results/audit_reruns/atac_row_order_<chrom>.csv
         results/audit_reruns/atac_row_order_summary.csv
RUNTIME  ~10-20 min (awk does both heavy scans). Memory stays under ~2 GB.
"""
import os
import subprocess

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
PD = f"{PROJ}/processed/per_dataset"
RAW = f"{PROJ}/raw/GSE268609_lazarov2024"
OUT = f"{PROJ}/results/audit_reruns"
os.makedirs(OUT, exist_ok=True)
CHROM = os.environ.get("ROWCHK_CHR", "chr1")   # chr1 is the first block in the fragment file
NSHUF = 100
rng = np.random.default_rng(42)

# ---- peaks on the target chromosome, keeping their global row index --------------------
names = [l.strip() for l in open(f"{PD}/GSE268609_peaks_names.txt")]
rowidx, starts, ends = [], [], []
for i, nm in enumerate(names):
    ch, s, e = nm.rsplit("-", 2)
    if ch == CHROM:
        rowidx.append(i); starts.append(int(s)); ends.append(int(e))
rowidx = np.array(rowidx); starts = np.array(starts); ends = np.array(ends)
order = np.argsort(starts)
rowidx, starts, ends = rowidx[order], starts[order], ends[order]
n_pk = len(rowidx)
print(f"{CHROM}: {n_pk} peaks (global rows {rowidx.min()}..{rowidx.max()})", flush=True)
assert n_pk > 100, "too few peaks to be discriminating"

# ---- 1. matrix row sums for those rows -------------------------------------------------
# awk does the 636M-line scan; python only reads the per-row totals back.
# NOTE the mtx is 1-based and its rows are in whatever order the export wrote them -- that is
# exactly the thing under test, so we address rows by NUMBER and never by name here.
lo, hi = int(rowidx.min()) + 1, int(rowidx.max()) + 1
awk_mtx = (r"NR>2 && $1>=%d && $1<=%d {s[$1]+=$3} END {for (r in s) print r, s[r]}" % (lo, hi))
res = subprocess.run(["awk", awk_mtx, f"{PD}/GSE268609_peaks_counts.mtx"],
                     capture_output=True, text=True, check=True)
rowsum = np.zeros(len(names), np.int64)
for ln in res.stdout.splitlines():
    r, v = ln.split()
    rowsum[int(r) - 1] = int(v)
mat = rowsum[rowidx].astype(float)
print(f"matrix row sums: total {mat.sum():,.0f}, nonzero {int((mat > 0).sum())}/{n_pk}", flush=True)

# ---- 2. recount the same intervals straight from the fragments -------------------------
# The fragment file is coordinate-sorted and chr1 is its FIRST block, so awk can exit as soon
# as the chromosome changes -- roughly a tenth of the 82 GB instead of all of it.
# Binary search over the peak starts is done inside awk; peaks are disjoint after re-calling.
frag = np.zeros(n_pk, np.int64)
pk_tmp = f"/tmp/_rowchk_{CHROM}_peaks.tsv"
with open(pk_tmp, "w") as fh:
    for s_, e_ in zip(starts, ends):
        fh.write(f"{s_}\t{e_}\n")
awk_frag = r"""
  BEGIN { while ((getline line < PK) > 0) { split(line, a, "	"); n++; S[n]=a[1]; E[n]=a[2] } }
  /^#/ { next }
  $1 != CH { if (seen) exit; next }
  { seen=1
    lo=1; hi=n; j=0
    while (lo <= hi) { m=int((lo+hi)/2); if (S[m] <= $3) { j=m; lo=m+1 } else { hi=m-1 } }
    if (j > 0 && E[j] >= $2) C[j]++ }
  END { for (k in C) print k, C[k] }
"""
p1 = subprocess.Popen(["zcat", f"{RAW}/GSE268609_Aggregated_atac_fragments.tsv.gz"],
                      stdout=subprocess.PIPE)
p2 = subprocess.run(["awk", "-v", f"PK={pk_tmp}", "-v", f"CH={CHROM}", awk_frag],
                    stdin=p1.stdout, capture_output=True, text=True)
p1.stdout.close(); p1.wait()
for ln in p2.stdout.splitlines():
    k, v = ln.split()
    frag[int(k) - 1] = int(v)
os.remove(pk_tmp)
print(f"fragment recount: total {frag.sum():,}, nonzero {int((frag > 0).sum())}/{n_pk}", flush=True)

# ---- 3. compare -------------------------------------------------------------------------
rho = float(spearmanr(frag, mat).statistic)
null = np.empty(NSHUF)
for i in range(NSHUF):
    null[i] = spearmanr(frag, rng.permutation(mat)).statistic
p95 = float(np.quantile(null, 0.95))
passed = bool(rho >= 0.90 and rho > p95)
print(f"\nSpearman(fragment recount, matrix row sums) = {rho:.4f}")
print(f"row-shuffled null: mean {null.mean():+.4f}  sd {null.std():.4f}  p95 {p95:+.4f}")
print(f"pre-registered gate (rho >= 0.90 AND rho > null p95): {'PASS' if passed else 'FAIL'}")

pd.DataFrame({"global_row": rowidx, "peak": [names[i] for i in rowidx],
              "start": starts, "end": ends,
              "fragment_recount": frag, "matrix_rowsum": mat}).to_csv(
    f"{OUT}/atac_row_order_{CHROM}.csv", index=False)
pd.DataFrame([dict(chrom=CHROM, n_peaks=n_pk, spearman_rho=round(rho, 6),
                   null_mean=round(float(null.mean()), 6), null_p95=round(p95, 6),
                   n_shuffles=NSHUF, gate_pass=passed,
                   fragment_total=int(frag.sum()), matrix_total=int(mat.sum()),
                   note="fragments cover all 366,175 barcodes; matrix keeps 153,530 QC-passing "
                        "cells, so magnitudes are not expected to match. This tests ORDER.")]
             ).to_csv(f"{OUT}/atac_row_order_summary.csv", index=False)
print(f"wrote {OUT}/atac_row_order_{CHROM}.csv and _summary.csv")
