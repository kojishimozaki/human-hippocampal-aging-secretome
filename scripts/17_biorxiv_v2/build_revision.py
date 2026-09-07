#!/usr/bin/env python3
"""Apply the bioRxiv v2 revision plan to the submitted DOCX.

Emits, from one plan:
  submission/v2/shimozaki-aging-secretome-v2.docx        (Word, formatting preserved)
  manuscript/biorxiv_v2/manuscript_biorxiv_v2.md         (English markdown of record)

Every operation is anchored on a prefix of the base paragraph it targets, so the build
fails loudly rather than silently mis-editing if the base document changes.

Usage:
    PROJ=$(pwd) python scripts/17_biorxiv_v2/build_revision.py [--allow-placeholders]
"""

from __future__ import annotations

import copy
import os
import re
import shutil
import sys
import zipfile
from pathlib import Path
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parent))
import revision_edits as R  # noqa: E402

PROJ = Path(os.environ.get("PROJ", Path(__file__).resolve().parents[2]))
BASE_DOCX = PROJ / "submission/shimozaki-aging-secretome-final_draft.docx"
OUT_DOCX = PROJ / "submission/v2/shimozaki-aging-secretome-v2.docx"
OUT_MD = PROJ / "manuscript/biorxiv_v2/manuscript_biorxiv_v2.md"

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
Q = lambda tag: f"{{{W}}}{tag}"  # noqa: E731
ET.register_namespace("w", W)

ALLOW_PLACEHOLDERS = "--allow-placeholders" in sys.argv


def para_text(p: ET.Element) -> str:
    return "".join(n.text or "" for n in p.iter(Q("t")))


def para_style(p: ET.Element) -> str:
    st = p.find(f"{Q('pPr')}/{Q('pStyle')}")
    return st.get(Q("val")) if st is not None else "Normal"


def set_text(p: ET.Element, text: str) -> None:
    """Replace a paragraph's text, keeping its style and its first run's formatting."""
    runs = p.findall(Q("r"))
    keep = copy.deepcopy(runs[0]) if runs else None
    for r in runs:
        p.remove(r)
    # Drop bookmarks/hyperlinks that would otherwise carry stale text.
    for tag in ("hyperlink", "bookmarkStart", "bookmarkEnd"):
        for el in p.findall(Q(tag)):
            p.remove(el)
    if keep is None:
        keep = ET.SubElement(p, Q("r"))
    for child in list(keep):
        if child.tag in (Q("t"), Q("br"), Q("tab")):
            keep.remove(child)
    t = ET.SubElement(keep, Q("t"))
    t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    t.text = text
    p.append(keep)


def clone_with_text(template: ET.Element, text: str) -> ET.Element:
    p = copy.deepcopy(template)
    set_text(p, text)
    return p


# --------------------------------------------------------------------------------------
# Statistical symbols are set in italic, and "P value" is written "P-value".
#
# The base document already italicised most of its P symbols; rewriting a paragraph
# through set_text() collapses it to one run and loses that, so this pass runs over every
# body paragraph at the end and restores it uniformly. It works run by run, cloning each
# run's own properties into the pieces it splits, so bold, superscript, size and colour
# survive. The hyphenation replaces a space with a hyphen and so does not move any
# character, which keeps the per-character property map valid.
# --------------------------------------------------------------------------------------

HYPHENATE = re.compile(r"(?<=[PQ]) (?=values?\b)")
# A standalone P or Q is always the statistic here (139 and 9 occurrences, all verified);
# n is only italicised where it introduces a count.
ITALIC_TOKEN = re.compile(r"(?<![A-Za-z0-9])(?:[PQ](?![A-Za-z0-9])|n(?=\s*=))")
# Exponents are superscript. Requiring the multiplication sign and the true minus keeps
# "10-nucleus" and "10-versus-20" out of it.
SUPERSCRIPT = re.compile("(?<=\u00d710)\u2212\\d+")


def rpr_signature(rpr) -> str:
    if rpr is None:
        return ""
    clone = copy.deepcopy(rpr)
    for it in clone.findall(Q("i")):
        clone.remove(it)
    return ET.tostring(clone).decode()


def apply_statistical_italics(p: ET.Element) -> None:
    runs = p.findall(Q("r"))
    if not runs:
        return
    chars: list[str] = []
    props: list[ET.Element | None] = []
    for r in runs:
        rpr = r.find(Q("rPr"))
        for t in r.findall(Q("t")):
            for ch in (t.text or ""):
                chars.append(ch)
                props.append(rpr)
    if not chars:
        return
    text = HYPHENATE.sub("-", "".join(chars))
    assert len(text) == len(chars), "hyphenation must not change the character count"

    italic = [False] * len(text)
    for m in ITALIC_TOKEN.finditer(text):
        italic[m.start()] = True
    superscript = [False] * len(text)
    for m in SUPERSCRIPT.finditer(text):
        for k in range(m.start(), m.end()):
            superscript[k] = True

    for r in runs:
        p.remove(r)
    i = 0
    while i < len(text):
        j = i + 1
        key = (rpr_signature(props[i]), italic[i], superscript[i])
        while j < len(text) and (rpr_signature(props[j]), italic[j], superscript[j]) == key:
            j += 1
        run = ET.SubElement(p, Q("r"))
        if props[i] is not None:
            rpr = copy.deepcopy(props[i])
            for it in rpr.findall(Q("i")):
                rpr.remove(it)
            if italic[i]:
                rpr.insert(0, ET.Element(Q("i")))
            if superscript[i] and rpr.find(Q("vertAlign")) is None:
                el = ET.SubElement(rpr, Q("vertAlign"))
                el.set(Q("val"), "superscript")
            run.append(rpr)
        elif italic[i] or superscript[i]:
            rpr = ET.SubElement(run, Q("rPr"))
            if italic[i]:
                ET.SubElement(rpr, Q("i"))
            if superscript[i]:
                ET.SubElement(rpr, Q("vertAlign")).set(Q("val"), "superscript")
        t = ET.SubElement(run, Q("t"))
        t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
        t.text = text[i:j]
        i = j


def para_markdown(p: ET.Element) -> str:
    """Paragraph text with italic runs wrapped for markdown."""
    out = []
    for r in p.findall(Q("r")):
        txt = "".join(t.text or "" for t in r.findall(Q("t")))
        if not txt:
            continue
        rpr = r.find(Q("rPr"))
        out.append(f"*{txt}*" if rpr is not None and rpr.find(Q("i")) is not None else txt)
    return "".join(out)


def build() -> None:
    if not BASE_DOCX.exists():
        raise SystemExit(
            f"base document not found: {BASE_DOCX}\n"
            "This target edits the superseded submission DOCX, which is deliberately not\n"
            "published (see .gitignore), so it is absent from a fresh clone. The manuscript\n"
            "text of record is manuscript/biorxiv_v2/manuscript_biorxiv_v2.md, and every\n"
            "data-driven target (figures, validation, supplementary-tables, bh-scope) runs\n"
            "without it."
        )
    with zipfile.ZipFile(BASE_DOCX) as z:
        names = z.namelist()
        blobs = {n: z.read(n) for n in names}
        # Carry each entry's original metadata so a rebuild is byte-identical: writing a
        # fresh archive would stamp the current time and leave the DOCX permanently dirty
        # in git even when nothing about the document changed.
        infos = {n: z.getinfo(n) for n in names}

    root = ET.fromstring(blobs["word/document.xml"])
    body = root.find(Q("body"))
    paras = body.findall(Q("p"))

    def check(idx: int, prefix: str) -> ET.Element:
        got = para_text(paras[idx])
        if not got.startswith(prefix):
            raise SystemExit(
                f"anchor failed at paragraph {idx}\n  expected prefix: {prefix!r}\n"
                f"  found          : {got[:len(prefix) + 20]!r}"
            )
        return paras[idx]

    # ---- anchors -----------------------------------------------------------------
    p_title = check(0, "A support-cell-enriched secretome-related program marks aging")
    p_abstract = check(4, "Research on brain aging has centred on neurons")
    p_keywords = check(5, "Keywords: hippocampal aging")
    p_statproc = check(57, "Analyses used pyDESeq2 0.5.4, scanpy 1.11.5")
    p_results_h = check(58, "Results")
    p_calib_old = check(65, "The pattern was stronger than expected after donor labels")
    p_disc_h1 = check(115, "Aging is accompanied by a change in the hippocampal support")
    check(122, "Convergent evidence supports a broader program")
    p_conv_body = check(123, "The second human hippocampal cohort provided the strongest")
    check(124, "The senescence and proteomic comparisons supply")
    check(125, "This distinction matters for intervention design")
    p_concl_h = check(146, "Conclusions and future directions")
    p_concl_1 = check(147, "Human hippocampal aging is associated with a support-cell")
    p_concl_2 = check(148, "The program remains a hypothesis-generating map")
    p_data_h = check(154, "Data availability")
    p_data_body = check(155, "All primary and secondary datasets are publicly available")
    p_atac_last = check(91, "Two further checks narrowed the interpretation")
    p_tables1 = check(174, "Table S1 | Donor-label permutation calibration")

    tmpl_body = paras[147]        # Normal body paragraph
    tmpl_h2 = paras[115]          # Heading2
    tmpl_datah = paras[154]       # the bold "Data availability" pseudo-heading

    # ---- edits -------------------------------------------------------------------
    set_text(p_title, R.TITLE)
    set_text(p_abstract, R.ABSTRACT)
    set_text(p_keywords, R.KEYWORDS)
    set_text(p_disc_h1, R.DISCUSSION_HEADING_SOFTENED)
    set_text(p_conv_body, R.CONVERGENCE_BRIDGE + para_text(p_conv_body))
    set_text(p_concl_h, R.CONCLUSIONS_HEADING)
    set_text(p_concl_1, R.CONCLUSIONS_1)
    set_text(p_concl_2, R.CONCLUSIONS_2)
    set_text(p_data_body, para_text(p_data_body) + R.DATA_AVAIL_ADDENDUM)
    set_text(p_tables1, para_text(p_tables1) + R.TABLE_S1_ADDENDUM)

    code_body = R.CODE_AVAIL_BODY.format(**R.PLACEHOLDERS)
    if not ALLOW_PLACEHOLDERS and any(v.startswith("[") for v in R.PLACEHOLDERS.values()):
        print(
            "note: Code availability still carries placeholders; "
            "fill PLACEHOLDERS in revision_edits.py before final submission.",
            file=sys.stderr,
        )

    insertions = {
        # after paragraph index -> list of (template, text)
        57: [(tmpl_h2, R.METHODS_REPRO_HEADING),
             (tmpl_body, R.METHODS_REPRO_1),
             (tmpl_body, R.METHODS_REPRO_2)],
        91: [(tmpl_h2, R.RESULTS_CALIB_HEADING),
             (tmpl_body, R.RESULTS_CALIB_1),
             (tmpl_body, R.RESULTS_CALIB_2),
             (tmpl_body, R.RESULTS_CALIB_3),
             (tmpl_body, R.RESULTS_CALIB_4),
             (tmpl_body, R.RESULTS_CALIB_4B),
             (tmpl_body, R.RESULTS_CALIB_5)],
        148: [(tmpl_body, R.CONCLUSIONS_3),
              (tmpl_body, R.CONCLUSIONS_4)],
        155: [(tmpl_datah, R.CODE_AVAIL_HEADING),
              (tmpl_body, code_body)],
    }
    deletions = {65}
    move_block = list(range(122, 126))   # convergence subsection
    move_before = 115                    # to the head of the Discussion

    # ---- rule-13 exact-value corrections -------------------------------------------
    for old_text, new_text in R.EXACT_VALUE_FIXES + R.REVIEW_FIXES:
        hits = [p for p in paras if old_text in para_text(p)]
        if len(hits) != 1:
            raise SystemExit(
                f"exact-value fix matched {len(hits)} paragraphs, expected 1: {old_text[:60]!r}"
            )
        set_text(hits[0], para_text(hits[0]).replace(old_text, new_text))

    # ---- assemble ----------------------------------------------------------------
    moved = [paras[i] for i in move_block]
    new_order: list[ET.Element] = []
    for i, p in enumerate(paras):
        if i in deletions or i in move_block:
            continue
        if i == move_before:
            new_order.extend(moved)
        new_order.append(p)
        for tmpl, text in insertions.get(i, []):
            new_order.append(clone_with_text(tmpl, text))

    assert len(new_order) == len(paras) - len(deletions) + sum(
        len(v) for v in insertions.values()
    ), "paragraph accounting failed"

    for p in new_order:
        if para_style(p) != "ListParagraph":
            apply_statistical_italics(p)

    for p in paras:
        body.remove(p)
    sect = body.find(Q("sectPr"))
    for p in new_order:
        if sect is not None:
            body.insert(list(body).index(sect), p)
        else:
            body.append(p)

    # ---- write DOCX --------------------------------------------------------------
    OUT_DOCX.parent.mkdir(parents=True, exist_ok=True)
    blobs["word/document.xml"] = ET.tostring(root, encoding="UTF-8", xml_declaration=True)
    with zipfile.ZipFile(OUT_DOCX, "w", zipfile.ZIP_DEFLATED) as z:
        for n in names:
            src = infos[n]
            info = zipfile.ZipInfo(n, date_time=src.date_time)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = src.external_attr
            info.internal_attr = src.internal_attr
            info.create_system = src.create_system
            z.writestr(info, blobs[n])

    # ---- write markdown of record ------------------------------------------------
    lines: list[str] = []
    in_refs = False
    ref_n = 0
    for p in new_order:
        style, text = para_style(p), para_markdown(p).strip()
        if not text:
            continue
        if style == "Title":
            lines += [f"# {text}", ""]
        elif style == "Heading1":
            in_refs = text.strip().lower().startswith("references")
            lines += [f"## {text}", ""]
        elif style == "Heading2":
            lines += [f"### {text}", ""]
        elif style == "ListParagraph" and in_refs:
            ref_n += 1
            lines += [f"{ref_n}. {text}", ""]
        else:
            if text in ("Data availability", R.CODE_AVAIL_HEADING):
                lines += [f"### {text}", ""]
            else:
                lines += [text, ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")

    print(f"base paragraphs : {len(paras)}")
    print(f"revised         : {len(new_order)}")
    print(f"references found: {ref_n}")
    print(f"wrote {OUT_DOCX.relative_to(PROJ)}")
    print(f"wrote {OUT_MD.relative_to(PROJ)}")


if __name__ == "__main__":
    build()
