#!/usr/bin/env python
# CONSEQUENCE (Phase A) — STEP 05: matrix provenance audit (confirms X = raw non-negative integer counts -> CPM valid).
# Provenance: VERBATIM from origin session 027c5fdb (transcript line 878). Produced the committed
#   matrix_provenance_audit.txt + matrix_value_diagnostics.tsv (VERDICT: PASS) byte-for-byte (agent-verified).
# env: bio.  DATA-DEPENDENT (GSE325391 h5ad) -> run from $PROJ. seed 42 (rule 10; sampling only).
# Out: results/resilience/{matrix_provenance_audit.txt, matrix_value_diagnostics.tsv}
import anndata as ad, numpy as np, scipy.sparse as sp, pandas as pd
A = ad.read_h5ad("processed/per_dataset/GSE325391_resilience.h5ad", backed='r')
lines=[]; w=lambda s:(lines.append(s),print(s))
w("=== GSE325391 matrix provenance audit (outcome-blind) ===")
w(f"adata.X dtype (backed): {A.X.dtype if hasattr(A.X,'dtype') else 'n/a'}")
w(f"adata.layers keys: {list(A.layers.keys())}")
w(f"adata.raw present: {A.raw is not None}")
# sample 8000 cells across all groups for value diagnostics
rng=np.random.default_rng(42); n=A.n_obs; idx=np.sort(rng.choice(n, size=min(8000,n), replace=False))
sub=A[idx].to_memory()
def diag(M,name):
    M=M.tocsr() if sp.issparse(M) else sp.csr_matrix(M)
    dat=M.data
    nonint=np.mean(np.abs(dat-np.round(dat))>1e-9) if dat.size else 0.0
    neg=np.mean(dat<0) if dat.size else 0.0
    libs=np.asarray(M.sum(1)).ravel()
    return dict(layer=name,min=float(dat.min()) if dat.size else 0,max=float(dat.max()) if dat.size else 0,
                nonint_frac=float(nonint),neg_frac=float(neg),
                lib_min=float(libs.min()),lib_med=float(np.median(libs)),lib_max=float(libs.max()),
                dtype=str(M.dtype))
rows=[diag(sub.X,"X")]
if "counts" in sub.layers: rows.append(diag(sub.layers["counts"],"counts"))
# is X == counts?
if "counts" in sub.layers:
    Xc=sub.X.tocsr() if sp.issparse(sub.X) else sp.csr_matrix(sub.X)
    Cc=sub.layers["counts"].tocsr() if sp.issparse(sub.layers["counts"]) else sp.csr_matrix(sub.layers["counts"])
    same=(Xc!=Cc).nnz==0
    w(f"X identical to layers['counts']: {same}")
df=pd.DataFrame(rows); df.to_csv("results/resilience/matrix_value_diagnostics.tsv",sep="\t",index=False)
w("\n"+df.to_string(index=False))
xr=rows[0]
verdict = "PASS (X is raw non-negative integer counts)" if (xr["nonint_frac"]<1e-6 and xr["neg_frac"]==0 and xr["max"]>30) else "FLAG (X not raw integer counts -> CPM invalid; use counts layer / rebuild)"
w(f"\nVERDICT: {verdict}")
w("build script (scripts/02_qc/06_build_325391_resilience.py): X <- counts.mtx (raw); layers['counts'] <- X.copy()")
w(f"count matrix to use for receptor set + score: {'X (== counts; raw)' if xr['nonint_frac']<1e-6 else 'layers[counts] or rebuild'}")
open("results/resilience/matrix_provenance_audit.txt","w").write("\n".join(map(str,lines)))
