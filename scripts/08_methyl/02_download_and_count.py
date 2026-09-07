#!/usr/bin/env python
"""M1 methylome — per-(donor x celltype) region mCG counting from GSE299139 ALLC files.

For each ALLC in results/methyl/allc_manifest.csv:
  1. download to a temp .gz  (curl --retry -C-, resumable)
  2. VERIFY integrity: exact byte size == manifest AND `gzip -t` CRC ok  (guards silent truncation ->
     a partial stream would corrupt mCG; we never trust an unverified file)
  3. parse CG-context rows on chr1-22,X (mawk, LC_ALL=C) -> per-region sum(mc), sum(cov) via per-chrom
     cumulative-sum + searchsorted (regions may overlap; each summed independently)
  4. save results/methyl/counts/<gsm>_<donor>_<src>.npz  (mc_cg, cov_cg aligned to regions_meta order,
     + genome-wide CG totals for the global-methylation/coverage QC covariate)
  5. delete the temp .gz
Resumable (skips existing npz). Donor = unit; this is pure per-(donor x celltype) pseudobulk methylation
(cell-level is structurally impossible — the deposit is already aggregated). mCG primary. seed-free.

hg38. Reference pipeline for the same quantity = ALLCools `allc-to-region-count` (Zemke 2024); we
reimplement transparently (sum mc_count / cov over CG-context cytosines, both strands, per region).

Usage: python scripts/08_methyl/02_download_and_count.py [--workers K] [--limit N] [--only CLASS]
"""
import os, sys, subprocess, tempfile, argparse, time
import numpy as np
import pandas as pd

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
RES = f"{PROJ}/results/methyl"
CNT = f"{RES}/counts"
TMP = os.environ.get("METHYL_TMP", f"{PROJ}/raw/GSE299139_snm3C/_tmp")
MAIN = [f"chr{c}" for c in list(range(1, 23)) + ["X"]]
os.makedirs(CNT, exist_ok=True)
os.makedirs(TMP, exist_ok=True)

_AWK = ('BEGIN{LC_ALL="C"; n=split("%s",a," "); for(i=1;i<=n;i++) keep[a[i]]=1}'
        ' keep[$1] && substr($4,1,2)=="CG"{print $1"\\t"$2"\\t"$5"\\t"$6}') % " ".join(MAIN)

# region geometry, loaded once per worker
_R = {}


def _init():
    reg = pd.read_csv(f"{RES}/regions_meta.csv")
    _R["n"] = len(reg)
    _R["chrom"] = reg.chrom.values
    by = {}
    for ch in MAIN:
        m = np.where(reg.chrom.values == ch)[0]
        if len(m):
            by[ch] = (m.astype(np.int64),
                      reg.start.values[m].astype(np.int64),
                      reg.end.values[m].astype(np.int64))
    _R["by"] = by


def _download_verify(url, size_bytes, dst, tries=5):
    for t in range(tries):
        r = subprocess.run(["curl", "-fsS", "--retry", "8", "--retry-delay", "5", "-C", "-",
                            "--speed-limit", "50000", "--speed-time", "60",   # abort+retry stalled (<50KB/s/60s)
                            "-o", dst, url])
        ok_size = os.path.exists(dst) and os.path.getsize(dst) == int(size_bytes)
        ok_gz = ok_size and subprocess.run(["gzip", "-t", dst]).returncode == 0
        if r.returncode == 0 and ok_size and ok_gz:
            return True
        if os.path.exists(dst) and not ok_size:
            # partial/corrupt: restart clean next attempt
            try: os.remove(dst)
            except OSError: pass
        time.sleep(3 * (t + 1))
    return False


def _count_one(row):
    gsm, donor, src = row["gsm"], row["donor"], row["src_celltype"]
    out = f"{CNT}/{gsm}_{donor}_{src}.npz"
    if os.path.exists(out):
        return (out, "skip", 0.0)
    t0 = time.time()
    dst = os.path.join(TMP, f"{gsm}_{donor}_{src}.allc.tsv.gz")
    if not _download_verify(row["url"], row["size_bytes"], dst):
        return (out, "DOWNLOAD_FAIL", time.time() - t0)
    try:
        p = subprocess.Popen(f"gzip -dc {dst!r} | LC_ALL=C mawk '{_AWK}'",
                             shell=True, stdout=subprocess.PIPE)
        df = pd.read_csv(p.stdout, sep="\t", header=None, names=["chrom", "pos", "mc", "cov"],
                         dtype={"chrom": "category", "pos": np.int64, "mc": np.int32, "cov": np.int32})
        p.wait()
        if p.returncode != 0:
            return (out, "PARSE_FAIL", time.time() - t0)
        n = _R["n"]
        mc_cg = np.zeros(n, np.float64); cov_cg = np.zeros(n, np.float64)
        g_mc = g_cov = 0
        for ch, sub in df.groupby("chrom", observed=True):
            ch = str(ch)
            if ch not in _R["by"]:
                continue
            pos = sub["pos"].values; o = np.argsort(pos, kind="stable")
            pos = pos[o]; mc = sub["mc"].values[o].astype(np.int64); cov = sub["cov"].values[o].astype(np.int64)
            cmc = np.concatenate([[0], np.cumsum(mc)]); ccov = np.concatenate([[0], np.cumsum(cov)])
            g_mc += int(cmc[-1]); g_cov += int(ccov[-1])
            ridx, rs, re = _R["by"][ch]
            li = np.searchsorted(pos, rs, "left"); ri = np.searchsorted(pos, re + 1, "left")  # [start,end] incl
            mc_cg[ridx] = cmc[ri] - cmc[li]; cov_cg[ridx] = ccov[ri] - ccov[li]
        np.savez_compressed(out, mc_cg=mc_cg.astype(np.int32), cov_cg=cov_cg.astype(np.int32),
                            global_mc=g_mc, global_cov=g_cov,
                            gsm=gsm, donor=donor, src=src)
        return (out, "ok", time.time() - t0)
    finally:
        try: os.remove(dst)
        except OSError: pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--only", default="")          # comma target_class filter
    ap.add_argument("--row", type=int, default=-1)  # single manifest row (test)
    a = ap.parse_args()
    man = pd.read_csv(f"{RES}/allc_manifest.csv")
    # process small/most-useful classes first (OPC validates the stats path fast; Oligo's 181 GB last)
    ORDER = {"OPC": 0, "Endo": 1, "Micro": 2, "Astro": 3, "Oligo": 4}
    man = (man.assign(_o=man.target_class.map(lambda c: ORDER.get(c, 9)))
              .sort_values(["_o", "donor", "src_celltype"]).drop(columns="_o").reset_index(drop=True))
    if a.only:
        man = man[man.target_class.isin(a.only.split(","))]
    if a.row >= 0:
        man = man.iloc[[a.row]]
    if a.limit:
        man = man.iloc[:a.limit]
    rows = man.to_dict("records")
    print(f"counting {len(rows)} ALLC files, workers={a.workers}", flush=True)
    _init()
    if a.workers <= 1:
        res = [_count_one(r) for r in rows]
    else:
        from multiprocessing import Pool
        with Pool(a.workers, initializer=_init) as pool:
            res = []
            for i, r in enumerate(pool.imap_unordered(_count_one, rows), 1):
                res.append(r)
                if r[1] != "ok" or i % 10 == 0 or i == len(rows):
                    print(f"  [{i}/{len(rows)}] {os.path.basename(r[0])} {r[1]} {r[2]:.0f}s", flush=True)
    nok = sum(1 for _, s, _ in res if s == "ok"); nsk = sum(1 for _, s, _ in res if s == "skip")
    bad = [(o, s) for o, s, _ in res if s not in ("ok", "skip")]
    print(f"\nDONE ok={nok} skip={nsk} fail={len(bad)}")
    for o, s in bad:
        print("  FAIL", s, os.path.basename(o))


if __name__ == "__main__":
    main()
