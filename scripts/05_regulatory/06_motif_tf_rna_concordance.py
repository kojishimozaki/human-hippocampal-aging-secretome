#!/usr/bin/env python
"""06 — Motif <-> TF-RNA same-direction cross-modal concordance.

Pre-registration: audit_log/2026-06-19_tf_rna_concordance/PRE_REGISTRATION.md
(registration commit 26ea769, tag tf-rna-concordance-prereg — locked BEFORE this code).

Headline (HAvYA primary, GSE268609): do aging-associated motif chromatin-accessibility
changes come with SAME-DIRECTION TF RNA change at higher frequency than a matched empirical
background?  Framed strictly as same-direction cross-modal concordance ENRICHMENT — NOT TF
activity, binding, occupancy, coherence, or causation. Descriptive ceiling (branch A): the
cis-mechanism conclusion (footprint/chromVAR/motif FDR-null) is UNCHANGED regardless of result.

ADvHA = descriptive projection of the FROZEN HAvYA foreground (no 3-group refit, no formal
heterogeneity, no foreground reselection). GSE278576 = RNA-direction robustness ONLY (not
replication; no external chromVAR leg).

env: bio. seed=42 (only the permutation/sampling steps are stochastic). PROJ for portability.
"""
import os, re, json, hashlib, subprocess, datetime
import numpy as np
import pandas as pd
from scipy import stats

PROJ = os.environ.get("PROJ", os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
SEED = 42
N_PERM = 10000
NICHE = ["Astro", "Micro", "Oligo", "OPC", "Endo"]
OUT = f"{PROJ}/results/regulatory"
os.makedirs(OUT, exist_ok=True)
PFX = f"{OUT}/motif_tf_rna_concordance"

IN = {
    "chromvar": f"{PROJ}/results/de/chromvar_motifs_strengthened.csv",
    "rna_HAvYA": f"{PROJ}/results/de/de_GSE268609_per_celltype.csv",
    "rna_ADvHA": f"{PROJ}/results/de/de_GSE268609_ADvHA.csv",
    "rna_278576": f"{PROJ}/results/de/de_GSE278576_sexadj_per_celltype.csv",
}
# 1:1 celltype map GSE268609 -> GSE278576 (Endo has no 278576 counterpart). Frozen pre-analysis.
CT_MAP_278576 = {"Astro": "Astro", "Micro": "Micro", "OPC": "OPC", "Oligo": "Oligo"}  # Endo -> mismatch

LOG = []
def log(m): print(m); LOG.append(str(m))

# ----------------------------------------------------------------------------- mapping
def motif_to_genes(tf):
    """Direct mapping: split dimers on '::', strip '(var.N)'. Returns list of gene symbols."""
    return [re.sub(r"\(var\.\d+\)", "", p).strip() for p in re.split("::", str(tf))]

def classify_family(g):
    g = str(g).upper()
    rules = [
        ("ETS", r"^(ELK\d|ELF\d|ETV\d|ETS\d|ERG|FLI1|FEV|EHF|GABPA|SPI[1BC]|SPDEF| ELG)"),
        ("AP-1", r"^(JUN|JUNB|JUND|FOS|FOSB|FOSL\d|BATF\d?|JDP2)$"),
        ("bHLH-neural", r"^(MYF\d|MYOD1|MYOG|ATOH\d|NEUROD\d|NEUROG\d|NHLH\d|BHLHE\d+|BHLHA\d+|HAND\d|HES\d|HEY\d|ID\d|TCF[34]$|TCF12|OLIG\d|PTF1A|FERD3L|MSC|MYF6)"),
        ("bHLH-other", r"^(MYC|MAX|MXI1|MNT|MLX|USF\d|SREBF\d|TFE[3B]|MITF|AHR|ARNT|HIF1A|EPAS1|NPAS\d|CLOCK|BMAL|ARNTL)"),
        ("TWIST", r"^TWIST\d$"),
        ("STAT", r"^STAT\d"),
        ("IRF", r"^IRF\d"),
        ("NFKB-REL", r"^(NFKB\d|RELA|RELB|REL)$"),
        ("CEBP", r"^CEBP[ABDEG]$"),
        ("FOX", r"^FOX"),
        ("SOX", r"^SOX\d"),
        ("TCF-LEF", r"^(TCF7|TCF7L1|TCF7L2|LEF1)$"),
        ("GATA", r"^GATA\d"),
        ("KLF-SP", r"^(KLF\d+|SP\d)$"),
        ("TEAD", r"^TEAD\d"),
        ("SMAD", r"^SMAD\d"),
        ("RUNX", r"^RUNX\d"),
        ("TP53fam", r"^TP(53|63|73)$"),
        ("E2F-DP", r"^(E2F\d|TFDP\d)$"),
        ("EGR-WT1", r"^(EGR\d|WT1)$"),
        ("HOX", r"^HOX[A-D]\d"),
        ("homeobox-other", r"^(BARHL\d|LHX\d|DLX\d|NKX|PAX\d|POU\d|PROX1|MEIS\d|PBX\d|PKNOX\d|EMX\d|VAX\d|EN\d|GSX\d|VSX\d|RAX|SIX\d|OTX\d|GBX\d)"),
        ("NR", r"^(NR\d|RAR[ABG]|RXR[ABG]|ESR\d|ESRR[ABG]|AR$|PGR$|PPAR[ABDG]|THR[AB]|VDR$|NR2C1|NR2C2|HNF4[AG]|RORA|RORB|RORC)"),
        ("ZBTB", r"^ZBTB\d"),
        ("RFX", r"^RFX\d"),
        ("CTCF", r"^CTCF"),
        ("ZNF-KRAB", r"^(ZNF\d|ZSCAN\d|ZKSCAN\d)"),
    ]
    for fam, pat in rules:
        if re.search(pat, g):
            return fam
    return "other"

# ----------------------------------------------------------------------------- load
chromvar = pd.read_csv(IN["chromvar"])
rna = {k: pd.read_csv(IN[k]) for k in ("rna_HAvYA", "rna_ADvHA", "rna_278576")}

# Gene-ID harmonization: de_GSE268609_ADvHA.csv is Ensembl-keyed (ENSG…); HAvYA & 278576 are
# symbol-keyed (and chromVAR TF names are symbols). Map ADvHA ENSG -> symbol via GENCODE v44,
# collapse (symbol,celltype) collisions by max baseMean (dominant transcript). Logged, traceable.
_gc = pd.read_csv(f"{PROJ}/refs/gencode_v44_genes.tsv.gz", sep="\t")
ENSG2SYM = dict(zip(_gc.ensg.astype(str), _gc.symbol.astype(str)))
_ad = rna["rna_ADvHA"].copy()
_ad["ensg"] = _ad["gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)
_n_ensg = _ad["ensg"].str.startswith("ENSG").mean()
_ad["gene"] = _ad["ensg"].map(ENSG2SYM)
_unmapped = int(_ad["gene"].isna().sum())
_ad = _ad.dropna(subset=["gene"])
_collide = int(_ad.duplicated(["gene", "celltype"], keep=False).sum())
_ad = _ad.sort_values("baseMean", ascending=False).drop_duplicates(["gene", "celltype"], keep="first")
rna["rna_ADvHA"] = _ad
ADVHA_IDNOTE = (f"ADvHA Ensembl->symbol via GENCODE v44 (input {_n_ensg:.0%} ENSG): "
                f"unmapped ENSG rows dropped={_unmapped}, (symbol,celltype) collisions collapsed by max baseMean={_collide}")
log(f"[harmonize] {ADVHA_IDNOTE}")

# sign-convention asserts (pre-reg §2): +log2FC == aging-UP (HA>YA / old>young)
# (ADvHA convention = AD − HA, +=disease-up, per repo `test vs ref` naming; documented in manifest.)
def lfc(df, g, c):
    r = df[(df.gene == g) & (df.celltype == c)]
    return None if r.empty else float(r.iloc[0].log2FoldChange)
assert lfc(rna["rna_HAvYA"], "IL15", "Micro") > 0 and lfc(rna["rna_HAvYA"], "COL21A1", "Astro") > 0, "HAvYA sign convention!"
assert lfc(rna["rna_278576"], "IL15", "Micro") > 0, "GSE278576 sign convention!"
log("[sign] verified: +log2FC == aging-UP in 268609 HAvYA and 278576 (old>young).")

# ----------------------------------------------------------------------------- build (gene,celltype) units per contrast
def build_units(contrast):
    """Collapse chromVAR motifs to (gene,celltype) canonical units for a contrast.
    Returns DataFrame: gene, celltype, motif (representative), dz, cv_p, cv_padj, family, is_foreground."""
    cv = chromvar[chromvar.contrast == contrast].copy()
    rows = []
    for _, r in cv.iterrows():
        for g in motif_to_genes(r.tf):
            rows.append((g, r.celltype, r.tf, r.dz, r.p, r.padj))
    ex = pd.DataFrame(rows, columns=["gene", "celltype", "motif", "dz", "cv_p", "cv_padj"])
    # canonical representative per (gene,celltype): smallest cv_p; tie-break prefer non-dimer, non-variant, lexical
    ex["_dimer"] = ex.motif.str.contains("::").astype(int)
    ex["_var"] = ex.motif.str.contains(r"\(var\.").astype(int)
    ex = ex.sort_values(["gene", "celltype", "cv_p", "_dimer", "_var", "motif"])
    rep = ex.drop_duplicates(["gene", "celltype"], keep="first").drop(columns=["_dimer", "_var"]).reset_index(drop=True)
    # foreground = a (gene,celltype) reachable by ANY motif with nominal cv_p<0.05 in that contrast's chromVAR
    fgkeys = set(zip(*[ex[ex.cv_p < 0.05][c] for c in ("gene", "celltype")]))
    rep["is_foreground"] = [ (g, c) in fgkeys for g, c in zip(rep.gene, rep.celltype) ]
    rep["family"] = rep.gene.map(classify_family)
    return rep

# ----------------------------------------------------------------------------- attach RNA + concordance
def attach_rna(units, rna_df, ct_filter=None):
    """Join RNA DE (gene,celltype). testable = log2FC & pvalue estimable. Adds concordance vs dz sign."""
    r = rna_df[["gene", "celltype", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"]].rename(
        columns={"log2FoldChange": "rna_lfc", "pvalue": "rna_p", "padj": "rna_q", "lfcSE": "rna_se", "baseMean": "rna_baseMean"})
    m = units.merge(r, on=["gene", "celltype"], how="left")
    if ct_filter is not None:
        m = m[m.celltype.isin(ct_filter)]
    m["rna_ci_lo"] = m.rna_lfc - 1.96 * m.rna_se
    m["rna_ci_hi"] = m.rna_lfc + 1.96 * m.rna_se
    m["testable"] = m.rna_lfc.notna() & m.rna_p.notna()
    cv_sign = np.sign(m.dz)
    rna_sign = np.sign(m.rna_lfc)
    m["cv_sign"] = cv_sign
    m["rna_sign"] = rna_sign
    m["concordant"] = np.where(m.testable & (cv_sign != 0) & (rna_sign != 0), cv_sign == rna_sign, np.nan)
    def cat(row):
        if not row.testable or pd.isna(row.concordant): return "D_not_testable"
        if row.concordant and pd.notna(row.rna_q) and row.rna_q < 0.05: return "A_concordant_q<0.05"
        if row.concordant: return "B_concordant_q>=0.05"
        return "C_discordant"
    m["category"] = m.apply(cat, axis=1)
    return m

# ----------------------------------------------------------------------------- enrichment (matched/stratified permutation + MH-OR + sensitivities)
def expr_bins(df, n=3):
    """log(baseMean) quantile bin within celltype over testable units."""
    out = pd.Series(0.0, index=df.index, dtype="float64")
    for ct, g in df.groupby("celltype"):
        v = np.log1p(g.rna_baseMean)
        try:
            b = pd.qcut(v, n, labels=False, duplicates="drop")
        except (ValueError, IndexError):
            b = pd.Series(0.0, index=g.index)
        out.loc[g.index] = pd.Series(b, index=g.index).astype("float64").fillna(0.0)
    return out.astype(int)

def wilson_ci(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z*z/n
    c = (p + z*z/(2*n)) / d
    h = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
    return (max(0, c-h), min(1, c+h))

def mh_or(fg, bg, strata_col):
    """Mantel-Haenszel OR (concordant vs discordant, fg vs bg) + RBG 95% CI."""
    num = den = 0.0; R = S = 0.0; rbg_pr = rbg_ps_qr = rbg_qs = 0.0
    for s in set(fg[strata_col]) | set(bg[strata_col]):
        f = fg[fg[strata_col] == s]; b = bg[bg[strata_col] == s]
        a = int(f.concordant.sum()); bb = int((~f.concordant.astype(bool)).sum())
        c = int(b.concordant.sum()); d = int((~b.concordant.astype(bool)).sum())
        nT = a+bb+c+d
        if nT == 0: continue
        R += a*d/nT; S += bb*c/nT
        rbg_pr += ((a+d)*a*d)/(nT*nT)
        rbg_ps_qr += ((a+d)*bb*c + (bb+c)*a*d)/(nT*nT)
        rbg_qs += ((bb+c)*bb*c)/(nT*nT)
    if S == 0 or R == 0:
        return (np.inf if S == 0 and R > 0 else (0.0 if R == 0 and S > 0 else np.nan), np.nan, np.nan)
    orr = R/S
    var_ln = rbg_pr/(2*R*R) + rbg_ps_qr/(2*R*S) + rbg_qs/(2*S*S)
    se = np.sqrt(var_ln)
    return (orr, orr*np.exp(-1.96*se), orr*np.exp(1.96*se))

def strat_perm(fg_tab, bg_tab, strata_cols, label, seed=SEED, n=N_PERM):
    """Stratified permutation: fix fg joint composition over strata, draw equal n from bg per stratum."""
    rng = np.random.default_rng(seed)
    fg = fg_tab.dropna(subset=["concordant"]).copy()
    bg = bg_tab.dropna(subset=["concordant"]).copy()
    fg["_s"] = list(zip(*[fg[c].astype(str) for c in strata_cols]))
    bg["_s"] = list(zip(*[bg[c].astype(str) for c in strata_cols]))
    pools = {k: g.concordant.astype(float).values for k, g in bg.groupby("_s")}
    matchable = fg[fg["_s"].map(lambda s: s in pools)]
    unmatched = len(fg) - len(matchable)
    if len(matchable) == 0:
        log(f"[{label}] no matchable foreground (all strata absent in bg) — test skipped.")
        return None
    obs = float(matchable.concordant.mean())
    counts = matchable.groupby("_s").size()
    total = int(counts.sum())
    contrib = np.zeros(n); repl = []
    for s, k in counts.items():
        pool = pools[s]; L = len(pool)
        if L >= k:
            order = rng.random((n, L)).argsort(1)[:, :k]
        else:
            order = rng.integers(0, L, size=(n, k)); repl.append((s, L, int(k)))
        contrib += pool[order].sum(1)
    perm = contrib / total
    p = (1 + int(np.sum(perm >= obs))) / (n + 1)
    # matched background rate = stratum-weighted bg concordance over matchable strata
    bg_rate = float(np.average([pools[s].mean() for s in counts.index],
                               weights=[counts[s] for s in counts.index]))
    orr, lo, hi = mh_or(matchable.assign(_s=matchable["_s"]), bg[bg["_s"].isin(counts.index)], "_s")
    k_c = int(matchable.concordant.sum()); ci = wilson_ci(k_c, total)
    if repl: log(f"[{label}] with-replacement strata (bg<fg): {repl}")
    if unmatched: log(f"[{label}] unmatched foreground dropped (stratum absent in bg): {unmatched}")
    return dict(label=label, strata=";".join(strata_cols), fg_testable=total, fg_concordant=k_c,
                fg_rate=obs, fg_ci_lo=ci[0], fg_ci_hi=ci[1], matched_bg_rate=bg_rate,
                risk_diff=obs-bg_rate, mh_or=orr, mh_or_lo=lo, mh_or_hi=hi,
                perm_p_onesided_greater=p, n_perm=n, seed=seed,
                n_unmatched_fg=unmatched, n_repl_strata=len(repl))

def enrichment_suite(units_rna, contrast_label):
    """Run primary stratified test + sensitivities for one contrast."""
    df = units_rna.copy()
    df = df[df.testable & df.concordant.notna()].copy()
    df["expr_bin"] = expr_bins(df)
    df["absdz_bin"] = (df.groupby("celltype").dz.transform(lambda s: pd.qcut(s.abs(), 3, labels=False, duplicates="drop"))
                       .fillna(0).astype(int))
    fg = df[df.is_foreground].copy(); bg = df[~df.is_foreground].copy()
    res = []
    # PRIMARY: celltype x cv_sign x expr_bin
    r = strat_perm(fg, bg, ["celltype", "cv_sign", "expr_bin"], f"{contrast_label}|PRIMARY_matched(ct,sign,expr)")
    if r: res.append(r)
    # SENS 3: + |dz| bin
    r = strat_perm(fg, bg, ["celltype", "cv_sign", "expr_bin", "absdz_bin"], f"{contrast_label}|SENS_absdz_matched")
    if r: res.append(r)
    # SENS 4: family-matched (ct, sign, family)
    r = strat_perm(fg, bg, ["celltype", "cv_sign", "family"], f"{contrast_label}|SENS_family_matched")
    if r: res.append(r)
    # SENS 1: pooled Fisher (one-sided greater)
    a = int(fg.concordant.sum()); b = int((~fg.concordant.astype(bool)).sum())
    c = int(bg.concordant.sum()); d = int((~bg.concordant.astype(bool)).sum())
    orf, pf = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    ci = wilson_ci(a, a+b)
    res.append(dict(label=f"{contrast_label}|SENS_pooled_Fisher", strata="none(pooled)", fg_testable=a+b,
                    fg_concordant=a, fg_rate=a/(a+b) if a+b else np.nan, fg_ci_lo=ci[0], fg_ci_hi=ci[1],
                    matched_bg_rate=c/(c+d) if c+d else np.nan, risk_diff=(a/(a+b)-c/(c+d)) if (a+b and c+d) else np.nan,
                    mh_or=orf, mh_or_lo=np.nan, mh_or_hi=np.nan, perm_p_onesided_greater=pf,
                    n_perm=0, seed=SEED, n_unmatched_fg=0, n_repl_strata=0))
    # SENS 2: binomial vs 0.5 (REFERENCE ONLY, not headline)
    bt = stats.binomtest(a, a+b, 0.5, alternative="greater") if a+b else None
    res.append(dict(label=f"{contrast_label}|SENS_binomial_vs0.5_REFERENCE_ONLY", strata="none", fg_testable=a+b,
                    fg_concordant=a, fg_rate=a/(a+b) if a+b else np.nan, fg_ci_lo=ci[0], fg_ci_hi=ci[1],
                    matched_bg_rate=0.5, risk_diff=(a/(a+b)-0.5) if a+b else np.nan, mh_or=np.nan,
                    mh_or_lo=np.nan, mh_or_hi=np.nan, perm_p_onesided_greater=(bt.pvalue if bt else np.nan),
                    n_perm=0, seed=SEED, n_unmatched_fg=0, n_repl_strata=0))
    # SENS 5: signed continuous — Spearman(rna_lfc, dz) fg & bg
    for name, sub in (("fg", fg), ("bg", bg)):
        if len(sub) >= 3:
            rho, pp = stats.spearmanr(sub.rna_lfc, sub.dz)
            res.append(dict(label=f"{contrast_label}|SENS_spearman_lfc_vs_dz_{name}", strata="continuous",
                            fg_testable=len(sub), fg_concordant=np.nan, fg_rate=rho, fg_ci_lo=np.nan, fg_ci_hi=np.nan,
                            matched_bg_rate=np.nan, risk_diff=np.nan, mh_or=np.nan, mh_or_lo=np.nan, mh_or_hi=np.nan,
                            perm_p_onesided_greater=pp, n_perm=0, seed=SEED, n_unmatched_fg=0, n_repl_strata=0))
    return pd.DataFrame(res), fg, bg

# ----------------------------------------------------------------------------- family (paralog) sensitivity per branch C
def family_paralog_sensitivity(units_rna, contrast_label):
    """For each FOREGROUND (gene,celltype), does ANY expressed family paralog move in the motif dz direction?"""
    fg = units_rna[units_rna.is_foreground].copy()
    rows = []
    for _, r in fg.iterrows():
        fam, ct, want = r.family, r.celltype, np.sign(r.dz)
        pool = units_rna[(units_rna.celltype == ct) & (units_rna.family == fam) & units_rna.testable]
        conc = pool[np.sign(pool.rna_lfc) == want]
        rows.append(dict(contrast=contrast_label, celltype=ct, family=fam, motif_gene=r.gene,
                         named_gene_expressed=bool(r.testable), dz=r.dz, dz_sign=int(want),
                         n_expressed_paralogs=int(len(pool)), n_paralogs_motif_direction=int(len(conc)),
                         paralogs_in_direction=";".join(f"{g}({l:+.2f})" for g, l in zip(conc.gene, conc.rna_lfc))))
    return pd.DataFrame(rows)

# ============================================================================= RUN
log(f"PROJ={PROJ}  seed={SEED}  n_perm={N_PERM}")

# --- HAvYA primary + ADvHA projection share the SAME frozen foreground (defined on HAvYA chromVAR)
units_HA = build_units("HAvYA")
fg_keys = set(zip(units_HA[units_HA.is_foreground].gene, units_HA[units_HA.is_foreground].celltype))
log(f"[foreground] frozen on HAvYA chromVAR nominal p<0.05: {len(fg_keys)} (gene,celltype) units")
assert len(fg_keys) == 31, f"foreground != 31 ({len(fg_keys)}) — frozen set changed!"

# ADvHA units: take SAME frozen foreground keys, do NOT reselect on ADvHA chromVAR
units_AD = build_units("ADvHA")
units_AD["is_foreground"] = [ (g, c) in fg_keys for g, c in zip(units_AD.gene, units_AD.celltype) ]

# attach RNA per contrast
HA = attach_rna(units_HA, rna["rna_HAvYA"])
AD = attach_rna(units_AD, rna["rna_ADvHA"])

n_fg_test_HA = int(HA[HA.is_foreground].testable.sum())
log(f"[HAvYA] foreground testable in RNA: {n_fg_test_HA}/31")
assert n_fg_test_HA == 22, f"HAvYA foreground testable != 22 ({n_fg_test_HA})"

# frozen_foreground.csv (with 278576 testability annotation)
r278 = rna["rna_278576"]; g278 = set(zip(r278.gene, r278.celltype))
def t278(g, c):
    if c not in CT_MAP_278576: return "celltype_mismatch"
    return "YES" if (g, CT_MAP_278576[c]) in g278 else "low_expr"
fz = HA[HA.is_foreground][["gene", "celltype", "motif", "family", "dz", "cv_p", "cv_padj", "testable",
                           "rna_lfc", "rna_q", "rna_baseMean", "category"]].copy()
fz["test_268609_HAvYA"] = fz.testable
fz["test_278576"] = [t278(g, c) for g, c in zip(fz.gene, fz.celltype)]
fz.to_csv(f"{PFX}_frozen_foreground.csv", index=False)

# mapping audit (all motifs -> genes, HAvYA)
ex_rows = []
for _, r in chromvar[chromvar.contrast == "HAvYA"].iterrows():
    for g in motif_to_genes(r.tf):
        ex_rows.append((r.tf, g, classify_family(g), r.celltype, r.dz, r.p, r.padj, (g, r.celltype) in fg_keys))
pd.DataFrame(ex_rows, columns=["motif", "mapped_gene", "family", "celltype", "dz", "cv_p", "cv_padj", "is_foreground"]
             ).to_csv(f"{PFX}_mapping_audit.csv", index=False)

# per-pair tables
cols_pair = ["gene", "celltype", "motif", "family", "is_foreground", "dz", "cv_p", "cv_padj",
             "rna_lfc", "rna_se", "rna_ci_lo", "rna_ci_hi", "rna_p", "rna_q", "rna_baseMean",
             "testable", "cv_sign", "rna_sign", "concordant", "category"]
HA[HA.is_foreground][cols_pair].to_csv(f"{PFX}_HAvYA_primary_pairs.csv", index=False)
AD[AD.is_foreground][cols_pair].to_csv(f"{PFX}_ADvHA_projection_pairs.csv", index=False)

# enrichment suites
enr_HA, fgHA, bgHA = enrichment_suite(HA, "HAvYA_primary")
enr_AD, fgAD, bgAD = enrichment_suite(AD, "ADvHA_projection")
enr_HA.to_csv(f"{PFX}_HAvYA_primary_enrichment.csv", index=False)
enr_AD.to_csv(f"{PFX}_ADvHA_projection_enrichment.csv", index=False)

# family paralog sensitivity
family_paralog_sensitivity(HA, "HAvYA").to_csv(f"{PFX}_HAvYA_family_sensitivity.csv", index=False)
family_paralog_sensitivity(AD, "ADvHA").to_csv(f"{PFX}_ADvHA_family_sensitivity.csv", index=False)

# cross-contrast descriptive (6 categories)
xc = HA[HA.is_foreground][["gene", "celltype", "dz", "rna_lfc", "rna_q", "concordant", "testable"]].merge(
    AD[AD.is_foreground][["gene", "celltype", "dz", "rna_lfc", "rna_q", "concordant", "testable"]],
    on=["gene", "celltype"], suffixes=("_HAvYA", "_ADvHA"))
def xcat(r):
    if not (r.testable_HAvYA and r.testable_ADvHA) or pd.isna(r.concordant_HAvYA) or pd.isna(r.concordant_ADvHA):
        return "6_not_testable"
    ch, ca = bool(r.concordant_HAvYA), bool(r.concordant_ADvHA)
    same_dir = np.sign(r.rna_lfc_HAvYA) == np.sign(r.rna_lfc_ADvHA)
    if ch and not ca: return "5_age-concordant_disease-discordant"
    if ch and ca and same_dir: return "3_shared_same-direction"
    if ch and ca and not same_dir: return "4_shared_opposite-direction"
    if (not ch) and ca: return "2_disease-only_concordant"
    if ch and not ca: return "1_age-only_concordant"
    return "C_both_discordant"
xc["cross_category"] = xc.apply(xcat, axis=1)
xc.to_csv(f"{PFX}_crosscontrast_descriptive.csv", index=False)

# GSE278576 RNA-direction robustness (frozen foreground; Endo -> mismatch)
rob_rows = []
ha_rna = rna["rna_HAvYA"]
for g, c in sorted(fg_keys):
    rec = dict(gene=g, celltype_268609=c, family=classify_family(g))
    r1 = ha_rna[(ha_rna.gene == g) & (ha_rna.celltype == c)]
    rec["lfc_268609"] = float(r1.iloc[0].log2FoldChange) if len(r1) else np.nan
    rec["q_268609"] = float(r1.iloc[0].padj) if len(r1) and pd.notna(r1.iloc[0].padj) else np.nan
    if c not in CT_MAP_278576:
        rec["status"] = "NOT_TESTABLE_CELLTYPE_MISMATCH"; rob_rows.append(rec); continue
    c2 = CT_MAP_278576[c]; r2 = r278[(r278.gene == g) & (r278.celltype == c2)]
    if not len(r1) or not len(r2) or pd.isna(r1.iloc[0].log2FoldChange) or pd.isna(r2.iloc[0].log2FoldChange):
        rec["status"] = "NOT_TESTABLE_LOW_EXPRESSION"; rob_rows.append(rec); continue
    rec.update(celltype_278576=c2, lfc_278576=float(r2.iloc[0].log2FoldChange),
               se_278576=float(r2.iloc[0].lfcSE), p_278576=float(r2.iloc[0].pvalue),
               q_278576=float(r2.iloc[0].padj) if pd.notna(r2.iloc[0].padj) else np.nan,
               baseMean_278576=float(r2.iloc[0].baseMean),
               direction_agree=bool(np.sign(rec["lfc_268609"]) == np.sign(r2.iloc[0].log2FoldChange)),
               status="TESTABLE")
    rob_rows.append(rec)
rob = pd.DataFrame(rob_rows)
rob.to_csv(f"{PFX}_GSE278576_RNA_robustness.csv", index=False)
rt = rob[rob.status == "TESTABLE"]
agree_k = int(rt.direction_agree.sum()); agree_n = len(rt)
rho278 = stats.spearmanr(rt.lfc_268609, rt.lfc_278576) if agree_n >= 3 else (np.nan, np.nan)
ci278 = stats.binomtest(agree_k, agree_n).proportion_ci(method="exact") if agree_n else (np.nan, np.nan)

# ----------------------------------------------------------------------------- manifest
def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
import scipy, numpy
manifest = dict(
    git_commit=subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJ, capture_output=True, text=True).stdout.strip(),
    prereg_tag="tf-rna-concordance-prereg", prereg_commit="26ea769bf2f9bdcff0e1abeeca48301e1070ac65",
    run_datetime=datetime.datetime.now().isoformat(timespec="seconds"),
    python=os.sys.version.split()[0], pandas=pd.__version__, numpy=numpy.__version__, scipy=scipy.__version__,
    seed=SEED, n_perm=N_PERM, seed_note="seed governs ONLY stratified-permutation sampling; all joins/tests deterministic",
    inputs={k: dict(path=v, sha256=sha256(v)) for k, v in IN.items()},
    foreground_units=len(fg_keys), foreground_testable_HAvYA=n_fg_test_HA,
    mapping_version="direct-name v1 (split '::', strip '(var.N)')",
    celltype_map_278576=CT_MAP_278576, celltype_map_note="Endo -> NOT_TESTABLE_CELLTYPE_MISMATCH",
    expr_bins="log1p(baseMean) tertile within celltype",
    bh_families="per-contrast RNA q from input DE universe (NOT recomputed on 22); enrichment families separated: HAvYA primary, ADvHA projection, family-sensitivity, 278576 robustness; never pooled",
    outputs=sorted(os.path.basename(p) for p in os.listdir(OUT) if p.startswith("motif_tf_rna_concordance")),
)
with open(f"{PFX}_run_manifest.json", "w") as f:
    json.dump(manifest, f, indent=2)

# ----------------------------------------------------------------------------- STOP-report summary
def pick(df, lab):
    r = df[df.label == lab]
    return r.iloc[0] if len(r) else None
log("\n================= STOP REPORT SUMMARY =================")
for lab, df in (("HAvYA primary", enr_HA), ("ADvHA projection", enr_AD)):
    prim = df.label[df.label.str.contains("PRIMARY")]
    if len(prim) == 0:
        log(f"\n[{lab}] PRIMARY matched test unavailable (no matchable strata)."); continue
    p = pick(df, prim.iloc[0])
    log(f"\n[{lab}]  matched-permutation (ct×sign×expr):")
    log(f"  foreground concordant {int(p.fg_concordant)}/{int(p.fg_testable)} = {p.fg_rate:.3f} "
        f"(95%CI {p.fg_ci_lo:.2f}–{p.fg_ci_hi:.2f}) | matched bg {p.matched_bg_rate:.3f} | "
        f"RD {p.risk_diff:+.3f} | MH-OR {p.mh_or:.2f} | perm p(greater)={p.perm_p_onesided_greater:.4f}"
        f" | unmatched_fg={int(p.n_unmatched_fg)}")
    fis = pick(df, f"{('HAvYA_primary' if 'HAvYA' in lab else 'ADvHA_projection')}|SENS_pooled_Fisher")
    log(f"  [sens] pooled Fisher OR={fis.mh_or:.2f} p={fis.perm_p_onesided_greater:.4f}")
log(f"\n[category counts HAvYA foreground]\n{HA[HA.is_foreground].category.value_counts().to_string()}")
log(f"\n[GSE278576 RNA-direction robustness] testable={agree_n} (Endo mismatch={int((rob.status=='NOT_TESTABLE_CELLTYPE_MISMATCH').sum())}, "
    f"low-expr={int((rob.status=='NOT_TESTABLE_LOW_EXPRESSION').sum())}); "
    f"direction-agree {agree_k}/{agree_n}" + (f" (95%CI {ci278[0]:.2f}–{ci278[1]:.2f}); Spearman ρ={rho278[0]:.2f} p={rho278[1]:.3f}" if agree_n else ""))
log(f"\n[outputs] {len(manifest['outputs'])} files at {OUT}/")
with open(f"{PFX}_run_log.txt", "w") as f:
    f.write("\n".join(LOG))
print("\nDONE.")
