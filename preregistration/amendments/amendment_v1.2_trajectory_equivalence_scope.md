# Amendment v1.2 — converged-universe scope for the trajectory equivalence audit (SUPERSEDED by v1.3)

**Status: recorded retroactively for a complete audit trail; SUPERSEDED by v1.3.** This adjudication occurred
between v1.1 and v1.3: the global AD−HA reparameterisation-equivalence gate was scoped, using the *actual*
PyDESeq2 gene-level LFC convergence flag (`dds.var["_LFC_converged"]`), to a **dual-fit converged + finite gene
universe** (Gate 1), with the **frozen-60 mandatory no-exclusion gate** (Gate 2). The gene universe, counts,
donors, `~arm+grp`, Wald tests, p-value universe, frozen set, projection gate, δ, and categories were NOT
changed — only the QC set for testing reparameterisation equivalence.

**Outcome (why it was superseded, not adopted):** the dual-gate was implemented and run. The 94 original
all-gene failures resolved against the real convergence flag as: 32 both-non-converged, 8 YA-only, 9 HA-only,
and **45 BOTH-converged genes still failing >1e-4** (dual-converged-universe max|Δ| ≈ 0.366). Per the
pre-specified rule ("両fit収束なのに >1e-4 → reparameterisation equivalence unresolved → STOP"), **Gate 1
FAILED** (Gate 2 passed, 60/60). Root cause: two independent `deseq2()` fits estimate dispersions
independently (~8% rel-diff on extreme low-count genes); a shared-nuisance coefficient-only transplant was
then tested and rejected as unsafe in PyDESeq2 0.5.4.

Therefore the converged-universe dual-fit gate is **abandoned**, and Phase B moves to a single R DESeq2 fit +
`ashr` direct-contrast shrinkage with a single-fit raw contrast-additivity identity gate (amendment v1.3). The
v1.2 commit/tag are recorded after the v1.3 commit only because the gate failed before v1.2 was formally
committed; the logical order is v1.1 → v1.2 → v1.3. v1, v1.1, v1.3 retained.
