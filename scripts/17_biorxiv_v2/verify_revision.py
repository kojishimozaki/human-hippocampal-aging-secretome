#!/usr/bin/env python3
"""Mechanical checks over the bioRxiv v2 manuscript surfaces.

Every check re-derives its expectation from data or from a second surface, never from the
manuscript restating itself. Run before calling the revision done.

Usage:
    PROJ=$(pwd) python scripts/17_biorxiv_v2/verify_revision.py
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

PROJ = Path(os.environ.get("PROJ", Path(__file__).resolve().parents[2]))
EN = PROJ / "manuscript/biorxiv_v2/manuscript_biorxiv_v2.md"
JA = PROJ / "manuscript/biorxiv_v2/manuscript_biorxiv_v2_ja.md"
DOCX = PROJ / "submission/v2/shimozaki-aging-secretome-v2.docx"
FIGS = PROJ / "submission/v2/figures"
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"

failures: list[str] = []
notes: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" — {detail}" if detail else ""))
    if not ok:
        failures.append(f"{name}: {detail}")


def docx_text() -> str:
    with zipfile.ZipFile(DOCX) as z:
        root = ET.fromstring(z.read("word/document.xml"))
    return "\n".join(
        "".join(n.text or "" for n in p.iter(W + "t"))
        for p in root.find(W + "body").findall(W + "p")
    )


def main() -> int:
    en = EN.read_text(encoding="utf-8")
    # The submitted Word file and its figures live in the working repository; the public
    # release carries the text of record instead, so the checks that read them report as
    # skipped rather than failing on a checkout that does not have them.
    dx = docx_text() if DOCX.exists() else None
    # The Japanese mirror is maintained in the working repository and is not part of the
    # public release, so the two checks that read it degrade to a note when it is absent.
    ja = JA.read_text(encoding="utf-8") if JA.exists() else None

    print("\nA. Approximation language (iron rule 13)")
    pat = r"≈|~\d|\bapproximately\b|\balmost exactly\b|\broughly\b|\babout \d"
    surfaces = [("English markdown", en)] + ([("DOCX", dx)] if dx is not None else [])
    for label, text in surfaces:
        hits = re.findall(pat, text)
        check(f"{label} free of approximations", not hits, ", ".join(sorted(set(hits))))
    if ja is None:
        notes.append("Japanese-mirror checks skipped: the mirror is not part of this "
                     "checkout (it is maintained in the working repository)")
        print("  [skip] Japanese mirror not present in this checkout")
    else:
        ja_body = re.sub(r"\A(?:<!--.*?-->\n|>.*\n|\n)+", "", ja)
        # "CLDN5 が内皮 subcluster の約 1/12" is a qualitative marker level, not a test
        # statistic; audit_log/2026-08-28_rho_notation/RESOLUTION.md ruled on it explicitly.
        ja_hits = [h for h in re.findall(r"≈|約\d+分の\d+|約\d|およそ", ja_body)
                   if "分の" not in h]
        check("Japanese mirror free of approximations", not ja_hits,
              ", ".join(sorted(set(ja_hits))))

    print("\nA2. Statistical symbols are italic and P-value is hyphenated")
    if dx is None:
        notes.append("DOCX checks skipped: the submitted Word file is not part of this "
                     "checkout (it is built in the working repository)")
        print("  [skip] submitted Word file not present in this checkout")
    else:
        with zipfile.ZipFile(DOCX) as z:
            root = ET.fromstring(z.read("word/document.xml"))
        plain, italic = [], []
        for para in root.find(W + "body").findall(W + "p"):
            st = para.find(f"{W}pPr/{W}pStyle")
            if st is not None and st.get(W + "val") == "ListParagraph":
                continue  # the reference list keeps its published typography
            for run in para.findall(W + "r"):
                txt = "".join(t.text or "" for t in run.findall(W + "t"))
                rpr = run.find(W + "rPr")
                is_it = rpr is not None and rpr.find(W + "i") is not None
                (italic if is_it else plain).append(txt)
        bare = re.findall(r"(?<![A-Za-z0-9])(?:[PQ](?![A-Za-z0-9])|n(?=\s*=))", "".join(plain))
        check("no statistical symbol left un-italicised", not bare,
              f"{len(bare)} found: {sorted(set(bare))}")
        check("italic runs carry only P, Q and n", set(italic) <= {"P", "Q", "n"},
              str(sorted(set(italic) - {"P", "Q", "n"})[:6]))
    for label, text in surfaces:
        check(f"{label} writes P-value, not P value",
              not re.search(r"\b[PQ]\s+values?\b", text))
    if ja is not None:
        ja_body = re.sub(r"\A(?:<!--.*?-->\n|>.*\n|\n)+", "", ja).split("## 参考文献")[0]
        check("Japanese mirror uses an uppercase italic P",
              not re.search(r"(?<![A-Za-z0-9])p\s*[=<]", ja_body)
              and not re.search(r"(?<![A-Za-z0-9*])P(?![A-Za-z0-9*\u5024])", ja_body))

    print("\nB. Calibration numbers against Table S1")
    t = pd.read_csv(PROJ / "manuscript/supplementary_tables/Table_S1.csv")
    p = t[t.cohort == "GSE268609"].set_index("celltype")
    q = t[(t.cohort == "GSE278576") & (t.status == "OK")]
    expect = {
        "Micro secretome P": (round(p.loc["Micro", "emp_p_secretome_fdr01"], 3), 0.020),
        "Astro secretome P": (round(p.loc["Astro", "emp_p_secretome_fdr01"], 3), 0.040),
        "Endo secretome P": (round(p.loc["Endo", "emp_p_secretome_fdr01"], 3), 0.040),
        "OPC secretome P": (round(p.loc["OPC", "emp_p_secretome_fdr01"], 3), 0.050),
        "Oligo secretome P": (round(p.loc["Oligo", "emp_p_secretome_fdr01"], 3), 0.109),
        "Astro genome-wide P": (round(p.loc["Astro", "emp_p_genomewide_fdr005"], 3), 0.030),
        "Micro genome-wide P": (round(p.loc["Micro", "emp_p_genomewide_fdr005"], 3), 0.030),
        "OPC genome-wide P": (round(p.loc["OPC", "emp_p_genomewide_fdr005"], 3), 0.050),
        "Oligo genome-wide P": (round(p.loc["Oligo", "emp_p_genomewide_fdr005"], 3), 0.139),
        "Endo genome-wide P": (round(p.loc["Endo", "emp_p_genomewide_fdr005"], 3), 0.069),
    }
    en_plain = en.replace("*", "")   # the markdown italicises the symbol
    for name, (got, want) in expect.items():
        txt = f"P = {want:.3f}"
        check(name, got == want and txt in en_plain,
              f"table {got}, manuscript text should carry {want:.3f}"
              + ("" if txt in en_plain else " — NOT FOUND IN TEXT"))
    check("null median zero in all five compartments",
          bool((p["null_median_secretome_fdr01"] == 0).all()))
    check("primary >=1 BH-significant rate is 29-42%",
          (round(p["null_frac_ge1_genomewide_fdr005"].min() * 100) == 29
           and round(p["null_frac_ge1_genomewide_fdr005"].max() * 100) == 42))
    check("primary >=10 rate is 8-16%",
          (round(p["null_frac_ge10_genomewide_fdr005"].min() * 100) == 8
           and round(p["null_frac_ge10_genomewide_fdr005"].max() * 100) == 16))
    check("external >=1 rate is 29-51%",
          (round(q["null_frac_ge1_genomewide_fdr005"].min() * 100) == 29
           and round(q["null_frac_ge1_genomewide_fdr005"].max() * 100) == 51))
    check("external >=10 rate is 3-17%",
          (round(q["null_frac_ge10_genomewide_fdr005"].min() * 100) == 3
           and round(q["null_frac_ge10_genomewide_fdr005"].max() * 100) == 17))
    check("24,310 / 11,440 label assignments",
          (p.loc["Astro", "n_distinct_label_splits"] == 24310
           and p.loc["Micro", "n_distinct_label_splits"] == 11440))
    check("microglial and vascular arms are 7 versus 9",
          bool((p.loc[["Micro", "Endo"], "n_donors_YA"] == 7).all()
               and (p.loc[["Micro", "Endo"], "n_donors_HA"] == 9).all()))

    print("\nC. Frozen-60 composition against the DE table")
    de = pd.read_csv(PROJ / "results/de/de_GSE268609_per_celltype.csv")
    sec = set(pd.read_csv(PROJ / "refs/secretome_union.csv")["gene"])
    hits = de[(de.padj < 0.1) & de.gene.isin(sec)]
    counts = hits.groupby("celltype").size().to_dict()
    check("60 pairs, 26/12/11/6/5 by compartment",
          len(hits) == 60 and counts == {"Astro": 26, "Micro": 12, "Endo": 11, "OPC": 6, "Oligo": 5},
          str(counts))
    check("57 distinct genes", hits.gene.nunique() == 57)
    check("36 increased and 24 decreased",
          int((hits.log2FoldChange > 0).sum()) == 36 and int((hits.log2FoldChange < 0).sum()) == 24)
    check("99,544 tests in the primary DE table", len(de) == 99544)
    check("Table S1 test counts sum to 99,543 (LINC01238, disclosed in Methods)",
          int(t[t.cohort == "GSE268609"].n_genes_tested.sum()) == 99543)

    print("\nD. Multiple-testing-scope claim")
    scope = pd.read_csv(PROJ / "results/de/global_bh_scope_GSE268609.csv")
    check("49 of 60 retained under a single genome-wide family",
          int((scope.scope_status == "retained_from_frozen60").sum()) == 49)
    new = scope[scope.scope_status == "new_under_global_bh"]
    check("4 new pairs, all oligodendrocyte, named in the text",
          len(new) == 4 and set(new.celltype) == {"Oligo"}
          and all(g in en_plain for g in new.gene), ", ".join(sorted(new.gene)))

    print("\nE. Nucleus counts against the deposited metadata")
    meta_path = PROJ / "processed/per_dataset/GSE268609_metadata.csv"
    if meta_path.exists():
        meta = pd.read_csv(meta_path, usecols=["Cluster"])
        vc = meta.Cluster.value_counts()
        for label, cluster in [("25,218 astrocyte", "Astrocytes"),
                               ("12,165 microglial", "Microglia"),
                               ("59,100 oligodendrocyte", "mOli"),
                               ("12,597 oligodendrocyte precursor", "OPCs")]:
            n = int(label.split()[0].replace(",", ""))
            check(f"{label} nuclei", vc[cluster] == n and label in en, f"metadata {vc[cluster]}")
        check("153,530 nuclei total", len(meta) == 153530 and "153,530" in en_plain)
    else:
        # processed/ is regenerable from GSE268609 but is not redistributed, so this
        # check runs in the working repository and is reported as skipped elsewhere.
        notes.append("nucleus-count check skipped: processed/ is not part of this "
                     "checkout (regenerate it from GSE268609 to run it)")
        print("  [skip] deposited metadata not present in this checkout")

    print("\nF. Supplementary figure numbering")
    legend = dict(re.findall(r"Fig\. (S\d) \| ([^.]{10,90})", en))
    if not FIGS.is_dir():
        notes.append("figure checks skipped: the submitted figures are not part of this "
                     "checkout (they are held in the working repository)")
        print("  [skip] submitted figures not present in this checkout")
        for n in range(1, 10):
            check(f"Fig. S{n} has a legend", f"S{n}" in legend, legend.get(f"S{n}", "")[:50])
    else:
        for n in range(1, 10):
            f = FIGS / f"Figure_S{n}.pdf"
            if not f.exists():
                check(f"Figure_S{n}.pdf present", False)
                continue
            out = subprocess.run(["pdftotext", "-layout", str(f), "-"],
                                 capture_output=True, text=True).stdout
            m = re.search(r"Shimozaki, Fig(?:ure)? (S\d+)", out)
            check(f"Figure_S{n}.pdf footer says S{n}", bool(m) and m.group(1) == f"S{n}",
                  m.group(1) if m else "no footer found")
            check(f"Fig. S{n} has a legend", f"S{n}" in legend,
                  legend.get(f"S{n}", "")[:50])

    print("\nF2. Figure bodies do not claim more than the text")
    # Editorial review (2026-09-07) found three claims printed inside figures that the
    # text contradicts. The figures under submission/ are composed by hand, so fixing the
    # scripts does not fix them; this check keeps the gap visible until it is closed.
    banned = {
        "survives genome-wide ATAC FDR":
            "text reports 0 of 369,393 genome-wide peaks significant; COL21A1 reaches "
            "significance only within the RNA-significant family",
        "stays control-like":
            "2 versus 37 differentially expressed genes does not establish equivalence",
        "Age-matched":
            "resilient donors are older (86 y) than severe-AD donors (81 y), unadjusted",
    }
    outlined = []
    for pdf in sorted(FIGS.glob("Figure_*.pdf")) if FIGS.is_dir() else []:
        out = subprocess.run(["pdftotext", "-layout", str(pdf), "-"],
                             capture_output=True, text=True).stdout
        if not out.strip():
            # Type converted to outlines: there is no text to read, so passing this check
            # would mean nothing. Say so rather than let it pass on an empty string.
            outlined.append(pdf.name)
            print(f"  [skip] {pdf.name} has no text layer; claims cannot be read")
            continue
        hits = [b for b in banned if b in out]
        check(f"{pdf.name} free of claims the text contradicts", not hits,
              "; ".join(f"{h} — {banned[h]}" for h in hits))
    if outlined:
        notes.append(
            "figures with type converted to outlines are exempt from the claim check and "
            "must be read by eye: " + ", ".join(outlined)
        )

    print("\nG. Supplementary figures are cited in numeric order")
    body = en.split("## Figure legends")[0]
    order, seen = [], set()
    for m in re.finditer(r"Fig\. (S\d+)", body):
        if m.group(1) not in seen:
            seen.add(m.group(1))
            order.append(m.group(1))
    check("first-citation order is S1..S9", order == [f"S{i}" for i in range(1, 10)], " ".join(order))

    print("\nH. Surfaces agree")
    def flat(s):
        s = re.sub(r"\A(?:<!--.*?-->\n|>.*\n|\n)+", "", s)
        # the markdown numbers the reference list explicitly; Word numbers it automatically
        s = re.sub(r"^\d{1,2}\. ", "", s, flags=re.M)
        return re.sub(r"[\s*`#>]", "", s)
    if dx is not None:
        en_flat, dx_flat = flat(en), flat(dx)
        check("DOCX body text equals the English markdown", en_flat == dx_flat,
              f"markdown {len(en_flat)} chars vs docx {len(dx_flat)}")
        check("Code availability present in both",
              "Code availability" in en and "Code availability" in dx)
    else:
        check("Code availability present", "Code availability" in en)
    check("Methods Reproducibility subsection present", "### Reproducibility" in en)
    check("38 references", len(re.findall(r"^\d+\. ", en.split("## References")[1], re.M)) == 38)

    print("\nI. English/Japanese measured-value parity")
    r = subprocess.run([sys.executable, str(PROJ / "scripts/17_biorxiv_v2/check_numbers.py")],
                       capture_output=True, text=True, env={**os.environ, "PROJ": str(PROJ)})
    last = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else ""
    if ja is None:
        print(f"  [skip] {last}")
    else:
        check("check_numbers.py passes", r.returncode == 0, last)

    print("\nJ. Placeholders")
    ph = re.findall(r"\[(?:repository URL|Zenodo DOI)[^\]]*\]", en)
    if ph:
        notes.append(f"Code availability still carries {len(ph)} placeholder(s) — fill before submission.")
    print(f"  [{'note' if ph else 'PASS'}] availability placeholders: {len(ph)}")

    print()
    for n in notes:
        print(f"NOTE: {n}")
    if failures:
        print(f"\n{len(failures)} CHECK(S) FAILED:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("\nAll checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
