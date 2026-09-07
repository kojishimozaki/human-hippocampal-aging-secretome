#!/usr/bin/env python3
"""Renumber four supplementary-figure footers to match the manuscript legends.

The figure PDFs were exported while the supplementary figures still carried mnemonic
labels (S(N), S(E1), S(R1)); the manuscript later assigned final numbers in Word and the
footers were never re-exported. Four files therefore carry a number that sends the reader
to a different figure:

    file (content)                         footer said   legend requires
    Figure _S3.pdf  (CSF proteome)         Fig S3        S4
    Figure _S4.pdf  (external cohorts)     Figure S4     S7
    Figure_ S6.pdf  (neuron specificity)   Figure S6     S3
    Figure_ S7.pdf  (TF RNA vs motif)      Figure S7     S6

Simply swapping the digit is not safe: the embedded LiberationSans-Bold subsets are missing
the digits two of the four files need ('4' and '6'), so the character maps through
ToUnicode -- pdftotext reports the right string -- but no glyph is drawn. The footers are
therefore re-typeset in Helvetica-Bold, a standard-14 font requiring no embedding and
metrically compatible with the Liberation Sans Bold already in use. The one footer whose
text length changes ("Fig S3" -> "Figure S4") also has its text-matrix origin shifted left
so that its right edge stays where it was.

Each change is a standards-compliant incremental update: new revisions of the content
stream and the page dictionary plus one new font object, and an appended xref section with
/Prev pointing at the previous table. The source files are not modified.

Usage:
    PROJ=$(pwd) python scripts/17_biorxiv_v2/fix_figure_footers.py
"""

from __future__ import annotations

import os
import re
import shutil
import zlib
from pathlib import Path

PROJ = Path(os.environ.get("PROJ", Path(__file__).resolve().parents[2]))
SRC = PROJ / "submission/figures"
DST = PROJ / "submission/v2/figures"

FONT_OBJ = b"<</Type/Font/Subtype/Type1/BaseFont/Helvetica-Bold/Encoding/WinAnsiEncoding>>"
FONT_NAME = b"FIXF"

# Helvetica-Bold advance widths (units/1000) for the characters used in a footer.
HB_WIDTHS = {
    "K": 722, ".": 278, " ": 278, "S": 667, "h": 611, "i": 278, "m": 889, "o": 611,
    "z": 556, "a": 556, "k": 611, ",": 278, "F": 611, "g": 611, "u": 611, "r": 389,
    "e": 556, "1": 556, "2": 556, "3": 556, "4": 556, "5": 556, "6": 556, "7": 556,
    "8": 556, "9": 556, "0": 556,
}

# source filename -> (old footer, new footer, output filename)
RENUMBER = {
    "Figure _S3.pdf": (b"K. Shimozaki, Fig S3", b"K. Shimozaki, Figure S4", "Figure_S4.pdf"),
    "Figure _S4.pdf": (b"K. Shimozaki, Figure S4", b"K. Shimozaki, Figure S7", "Figure_S7.pdf"),
    "Figure_ S6.pdf": (b"K. Shimozaki, Figure S6", b"K. Shimozaki, Figure S3", "Figure_S3.pdf"),
    "Figure_ S7.pdf": (b"K. Shimozaki, Figure S7", b"K. Shimozaki, Figure S6", "Figure_S6.pdf"),
}
RENAME_ONLY = {
    "Figure _1.pdf": "Figure_1.pdf", "Figure _2.pdf": "Figure_2.pdf",
    "Figure _3.pdf": "Figure_3.pdf", "Figure _4.pdf": "Figure_4.pdf",
    "Figure _5.pdf": "Figure_5.pdf", "Figure _6.pdf": "Figure_6.pdf",
    "Figure _S1.pdf": "Figure_S1.pdf", "Figure _S2.pdf": "Figure_S2.pdf",
    "Figure_ S5.pdf": "Figure_S5.pdf", "Figure_ S8.pdf": "Figure_S8.pdf",
    "Figure_ S9.pdf": "Figure_S9.pdf",
}


def text_width(s: str, size: float) -> float:
    return sum(HB_WIDTHS[c] for c in s) / 1000.0 * size


def index_objects(raw: bytes) -> dict[int, tuple[int, int]]:
    idx: dict[int, tuple[int, int]] = {}
    for m in re.finditer(rb"(?<![0-9])(\d+)\s+0\s+obj", raw):
        end = raw.find(b"endobj", m.end())
        idx[int(m.group(1))] = (m.end(), end if end != -1 else len(raw))
    return idx


def find_footer_stream(raw: bytes, idx: dict, needle: bytes):
    for num, (a, b) in sorted(idx.items()):
        s = raw.find(b"stream", a)
        if s == -1 or s > b:
            continue
        s2 = raw.find(b"\n", s) + 1
        es = raw.find(b"endstream", s2)
        if es == -1 or es > b:
            continue
        try:
            data = zlib.decompress(raw[s2:es])
        except zlib.error:
            continue
        if needle in data:
            return num, raw[a:s], data
    raise SystemExit(f"footer stream not found for {needle!r}")


def find_page_for(raw: bytes, idx: dict, contents_obj: int) -> int:
    want = f"/Contents {contents_obj} 0 R".encode()
    for num, (a, b) in idx.items():
        chunk = raw[a:b]
        if want in chunk and re.search(rb"/Type\s*/Page\b", chunk):
            return num
    raise SystemExit(f"page object referencing /Contents {contents_obj} not found")


def rebuild_footer(stream: bytes, old: bytes, new: bytes) -> bytes:
    """Point the footer at our own font, and re-origin it if the string got longer."""
    show_old = b"(" + old + b")Tj"
    if stream.count(show_old) != 1:
        raise SystemExit(f"expected exactly one {show_old!r}, found {stream.count(show_old)}")
    i = stream.index(show_old)
    tm = re.search(
        rb"([\d.]+) 0 0 ([\d.]+) ([\d.]+) ([\d.]+) Tm\s*$", stream[max(0, i - 120):i]
    )
    if not tm:
        raise SystemExit("could not read the footer text matrix")
    size = float(tm.group(1))
    x, y = float(tm.group(3)), float(tm.group(4))
    dx = text_width(new.decode(), size) - text_width(old.decode(), size)
    replacement = b"/" + FONT_NAME + b" 1 Tf\n(" + new + b")Tj"
    out = stream[:i] + replacement + stream[i + len(show_old):]
    if abs(dx) > 0.5:  # keep the right edge where it was
        old_tm = tm.group(0).rstrip()
        new_tm = f"{tm.group(1).decode()} 0 0 {tm.group(2).decode()} {x - dx:.4f} {y} Tm".encode()
        head, sep, tail = out.rpartition(old_tm)
        if not sep:
            raise SystemExit("could not rewrite the footer text matrix")
        out = head + new_tm + tail
    return out


def add_font_resource(page_body: bytes, font_obj: int) -> bytes:
    m = re.search(rb"/Font\s*<<", page_body)
    if not m:
        raise SystemExit("page /Resources has no inline /Font dictionary")
    ins = b"/" + FONT_NAME + f" {font_obj} 0 R".encode()
    return page_body[:m.end()] + ins + page_body[m.end():]


def incremental_update(raw: bytes, new_objects: dict[int, bytes], new_size: int) -> bytes:
    tail = raw[-3072:]
    m_start = re.search(rb"startxref\s+(\d+)\s*%%EOF\s*$", tail)
    m_tr = re.search(rb"trailer\s*<<(.*?)>>\s*startxref", tail, re.S)
    if not (m_start and m_tr):
        raise SystemExit("could not read the existing xref trailer")
    prev_xref = int(m_start.group(1))
    trailer = re.sub(rb"/Prev\s+\d+", b"", m_tr.group(1))
    trailer = re.sub(rb"/Size\s+\d+", f"/Size {new_size}".encode(), trailer, count=1)

    out = bytearray(raw)
    if not out.endswith(b"\n"):
        out += b"\n"
    offsets: dict[int, int] = {}
    for num in sorted(new_objects):
        offsets[num] = len(out)
        out += f"{num} 0 obj\n".encode() + new_objects[num] + b"\nendobj\n"

    xref_off = len(out)
    out += b"xref\n0 1\n0000000000 65535 f \n"
    for num in sorted(offsets):  # one subsection per object keeps this simple and valid
        out += f"{num} 1\n".encode() + f"{offsets[num]:010d} 00000 n \n".encode()
    out += b"trailer\n<<" + trailer.strip() + f"/Prev {prev_xref}".encode() + b">>\n"
    out += f"startxref\n{xref_off}\n%%EOF\n".encode()
    return bytes(out)


def main() -> None:
    if not SRC.exists():
        raise SystemExit(
            f"source figures not found: {SRC}\n"
            "This was a one-off renumbering of the superseded submission figures, which\n"
            "are not redistributed. The corrected figures it produced are in\n"
            "submission/v2/figures/ and are part of this repository."
        )
    DST.mkdir(parents=True, exist_ok=True)
    for src_name, dst_name in RENAME_ONLY.items():
        shutil.copy2(SRC / src_name, DST / dst_name)
        print(f"copied   {src_name!r:20s} -> {dst_name}")

    for src_name, (old, new, dst_name) in RENUMBER.items():
        raw = (SRC / src_name).read_bytes()
        idx = index_objects(raw)
        cnum, chead, cdata = find_footer_stream(raw, idx, old)
        pnum = find_page_for(raw, idx, cnum)
        size = int(re.search(rb"/Size\s+(\d+)", raw[-3072:]).group(1))
        fnum = size  # first free object number

        patched_stream = rebuild_footer(cdata, old, new)
        comp = zlib.compress(patched_stream, 9)
        chead = re.sub(rb"/Length\s+\d+", b"/Length " + str(len(comp)).encode(), chead, count=1)
        if b"/Length" not in chead:
            raise SystemExit(f"{src_name}: content stream dictionary has no /Length")

        a, b = idx[pnum]
        page_body = add_font_resource(raw[a:b], fnum)

        new_objects = {
            cnum: chead.strip() + b"\nstream\n" + comp + b"\nendstream",
            pnum: page_body.strip(),
            fnum: FONT_OBJ,
        }
        (DST / dst_name).write_bytes(incremental_update(raw, new_objects, size + 1))
        print(
            f"patched  {src_name!r:20s} -> {dst_name}  "
            f"(content {cnum}, page {pnum}, font {fnum}): "
            f"{old.decode()!r} -> {new.decode()!r}"
        )


if __name__ == "__main__":
    main()
