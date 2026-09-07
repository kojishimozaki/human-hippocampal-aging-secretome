#!/usr/bin/env python3
"""Create a citation-corrected redline DOCX while preserving existing revisions."""

from __future__ import annotations

import copy
import io
import os
import sys
import tempfile
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path


W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
W14_NS = "http://schemas.microsoft.com/office/word/2010/wordml"
XML_NS = "http://www.w3.org/XML/1998/namespace"


def qn(namespace: str, local: str) -> str:
    return f"{{{namespace}}}{local}"


W_P = qn(W_NS, "p")
W_PPR = qn(W_NS, "pPr")
W_R = qn(W_NS, "r")
W_RPR = qn(W_NS, "rPr")
W_T = qn(W_NS, "t")
W_DEL_TEXT = qn(W_NS, "delText")
W_COLOR = qn(W_NS, "color")
W_STRIKE = qn(W_NS, "strike")
W_VAL = qn(W_NS, "val")
XML_SPACE = qn(XML_NS, "space")
W14_PARA_ID = qn(W14_NS, "paraId")
W14_TEXT_ID = qn(W14_NS, "textId")


def register_namespaces(xml_bytes: bytes) -> None:
    seen: set[tuple[str, str]] = set()
    for _event, item in ET.iterparse(io.BytesIO(xml_bytes), events=("start-ns",)):
        prefix, uri = item
        key = (prefix or "", uri)
        if key in seen or prefix == "xml":
            continue
        seen.add(key)
        try:
            ET.register_namespace(prefix or "", uri)
        except ValueError:
            # ElementTree reserves ns\d+ prefixes. Word accepts generated prefixes
            # for those rare namespaces, while the important w/w14/mc prefixes are
            # registered normally.
            pass


def paragraph_text(paragraph: ET.Element) -> str:
    return "".join(
        node.text or ""
        for node in paragraph.iter()
        if node.tag in {W_T, W_DEL_TEXT}
    )


def find_paragraph(root: ET.Element, anchor: str) -> ET.Element:
    matches = [p for p in root.iter(W_P) if anchor in paragraph_text(p)]
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one paragraph for anchor {anchor!r}; found {len(matches)}"
        )
    return matches[0]


def make_run(template: ET.Element, text: str, *, red: bool = False, strike: bool = False) -> ET.Element:
    run = ET.Element(W_R, dict(template.attrib))
    template_rpr = template.find(W_RPR)
    if template_rpr is not None:
        rpr = copy.deepcopy(template_rpr)
        run.append(rpr)
    elif red or strike:
        rpr = ET.SubElement(run, W_RPR)
    else:
        rpr = None

    if red:
        assert rpr is not None
        for existing in list(rpr.findall(W_COLOR)):
            rpr.remove(existing)
        color = ET.SubElement(rpr, W_COLOR)
        color.set(W_VAL, "FF0000")

    if strike:
        assert rpr is not None
        if rpr.find(W_STRIKE) is None:
            ET.SubElement(rpr, W_STRIKE)

    text_node = ET.SubElement(run, W_T)
    text_node.text = text
    if text.startswith(" ") or text.endswith(" "):
        text_node.set(XML_SPACE, "preserve")
    return run


def split_text_run(
    paragraph: ET.Element,
    target: str,
    *,
    inserted_after: str | None = None,
    style_target_red_strike: bool = False,
) -> None:
    parent_map = {child: parent for parent in paragraph.iter() for child in parent}
    candidates = [node for node in paragraph.iter(W_T) if target in (node.text or "")]
    if len(candidates) != 1:
        raise RuntimeError(
            f"Expected one text node containing {target!r} in paragraph; found {len(candidates)}"
        )

    text_node = candidates[0]
    run = parent_map[text_node]
    if run.tag != W_R:
        raise RuntimeError(f"Target {target!r} is not inside a standard Word run")
    container = parent_map[run]

    allowed = {W_RPR, W_T}
    unexpected = [child.tag for child in run if child.tag not in allowed]
    if unexpected:
        raise RuntimeError(f"Run for {target!r} has unexpected children: {unexpected}")

    original = text_node.text or ""
    before, after = original.split(target, 1)
    replacement: list[ET.Element] = []
    if before:
        replacement.append(make_run(run, before))
    replacement.append(
        make_run(run, target, red=style_target_red_strike, strike=style_target_red_strike)
    )
    if inserted_after is not None:
        replacement.append(make_run(run, inserted_after, red=True))
    if after:
        replacement.append(make_run(run, after))

    index = list(container).index(run)
    container.remove(run)
    for offset, new_run in enumerate(replacement):
        container.insert(index + offset, new_run)


def add_red_citation(root: ET.Element, anchor: str, accession: str, citation: str) -> None:
    paragraph = find_paragraph(root, anchor)
    split_text_run(paragraph, accession, inserted_after=f" [{citation}]")


def add_reference_after(
    root: ET.Element,
    after_anchor: str,
    text: str,
    para_id: str,
) -> ET.Element:
    paragraph = find_paragraph(root, after_anchor)
    parent_map = {child: parent for parent in root.iter() for child in parent}
    container = parent_map[paragraph]
    index = list(container).index(paragraph)

    new_p = ET.Element(W_P)
    new_p.set(W14_PARA_ID, para_id)
    new_p.set(W14_TEXT_ID, "77777777")
    ppr = paragraph.find(W_PPR)
    if ppr is not None:
        new_p.append(copy.deepcopy(ppr))

    template_run = paragraph.find(W_R)
    if template_run is None:
        raise RuntimeError("Reference template paragraph has no run")
    new_p.append(make_run(template_run, text, red=True))
    container.insert(index + 1, new_p)
    return new_p


def transform_document(xml_bytes: bytes) -> bytes:
    register_namespaces(xml_bytes)
    root = ET.fromstring(xml_bytes)

    # Minor bibliographic correction: only the newly added suffix is red.
    ref3 = find_paragraph(root, "Boldrini M, Fulmore CA, Tartt AN")
    split_text_run(ref3, "589–599", inserted_after=".e5")

    # Red strikethrough marks deletion candidates while keeping the redline legible.
    ref18 = find_paragraph(root, "Zemke NR, Lee S, Mamde S")
    split_text_run(
        ref18,
        " (preprint: bioRxiv 2024.10.14.618338, doi:10.1101/2024.10.14.618338)",
        style_target_red_strike=True,
    )
    ref22 = find_paragraph(root, "Uhlén M, Karlsson MJ, Hober A")
    split_text_run(
        ref22,
        " Tissue-based map: Uhlén M, Fagerberg L, Hallström BM, et al. Science. 2015;347(6220):1260419. doi:10.1126/science.1260419",
        style_target_red_strike=True,
    )

    # Add the two missing primary-paper citations wherever the datasets are
    # introduced or used as a named source.
    citation_insertions = [
        (
            "Spatial analyses used GSE264692 (Visium), and the cognitive-resilience analysis used GSE325391.",
            "GSE264692",
            "37",
        ),
        (
            "the cognitive-resilience analysis used GSE325391.",
            "GSE325391",
            "38",
        ),
        ("For GSE268609 and GSE325391, we retained", "GSE325391", "38"),
        ("Spatial localisation was tested in GSE264692", "GSE264692", "37"),
        ("In GSE325391, which profiles the subgranular-zone", "GSE325391", "38"),
        ("GSE325391 offered a separate, descriptive view", "GSE325391", "38"),
        ("receiver-availability projection in GSE325391", "GSE325391", "38"),
        (
            "All primary and secondary datasets are publicly available under GSE268609",
            "GSE264692",
            "37",
        ),
        (
            "All primary and secondary datasets are publicly available under GSE268609",
            "GSE325391",
            "38",
        ),
        ("Visium spatial transcriptomics from 34 quality-controlled GSE264692", "GSE264692", "37"),
        ("Fig. 6 | Granule-lineage transcription in a cognitive-resilience cohort", "GSE325391", "38"),
        ("separate GSE325391 granule-lineage cohort", "GSE325391", "38"),
    ]
    for anchor, accession, citation in citation_insertions:
        add_red_citation(root, anchor, accession, citation)

    ref37 = (
        "Thompson JR, Nelson ED, Tippani M, et al. An integrated single-nucleus and "
        "spatial transcriptomics atlas reveals the molecular landscape of the human "
        "hippocampus. Nat Neurosci. 2025;28(9):1990–2004. "
        "doi:10.1038/s41593-025-02022-0"
    )
    ref38 = (
        "Tosoni G, Ayyildiz D, Snoeck S, et al. Transcriptional profiles of immature "
        "neurons in aged human hippocampus track Alzheimer’s pathology and cognitive "
        "resilience. Cell Stem Cell. 2026;33(5):763–783.e9. "
        "doi:10.1016/j.stem.2026.04.002"
    )
    inserted_ref37 = add_reference_after(
        root,
        "Eng CL, Lawson M, Zhu Q, et al.",
        ref37,
        "C17A0037",
    )

    # Insert reference 38 immediately after the newly inserted reference 37.
    parent_map = {child: parent for parent in root.iter() for child in parent}
    container = parent_map[inserted_ref37]
    index = list(container).index(inserted_ref37)
    new_p = ET.Element(W_P)
    new_p.set(W14_PARA_ID, "C17A0038")
    new_p.set(W14_TEXT_ID, "77777777")
    ppr = inserted_ref37.find(W_PPR)
    if ppr is not None:
        new_p.append(copy.deepcopy(ppr))
    template_run = inserted_ref37.find(W_R)
    if template_run is None:
        raise RuntimeError("Inserted reference 37 has no run")
    new_p.append(make_run(template_run, ref38, red=True))
    container.insert(index + 1, new_p)

    return ET.tostring(root, encoding="utf-8", xml_declaration=True)


def build_docx(source: Path, output: Path) -> None:
    if source.resolve() == output.resolve():
        raise ValueError("Source and output DOCX paths must differ")
    if not source.is_file():
        raise FileNotFoundError(source)

    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(source, "r") as src_zip:
        document_xml = src_zip.read("word/document.xml")
        transformed = transform_document(document_xml)

        fd, tmp_name = tempfile.mkstemp(
            prefix=f".{output.stem}.", suffix=".docx", dir=output.parent
        )
        os.close(fd)
        tmp_path = Path(tmp_name)
        try:
            with zipfile.ZipFile(tmp_path, "w") as out_zip:
                for info in src_zip.infolist():
                    payload = transformed if info.filename == "word/document.xml" else src_zip.read(info.filename)
                    out_zip.writestr(info, payload)
            os.replace(tmp_path, output)
        finally:
            if tmp_path.exists():
                tmp_path.unlink()


def main() -> int:
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} SOURCE.docx OUTPUT.docx", file=sys.stderr)
        return 2
    source = Path(sys.argv[1])
    output = Path(sys.argv[2])
    build_docx(source, output)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
