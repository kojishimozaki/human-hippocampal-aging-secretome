#!/usr/bin/env python
"""Correct the Zenodo DOI inside the merged manuscript PDF.

The merged PDF (manuscript text pages plus the fifteen composed figures) is assembled by
hand, so it cannot be regenerated from `build_revision.py`.  When the manuscript moved
from citing the first Zenodo version (10.5281/zenodo.22640703, code and tables only) to
the concept DOI (10.5281/zenodo.22640702, which always resolves to the newest version),
the PDF had to be corrected in place instead.

Why this is safe rather than a hack:

  * The text is drawn with a subset of Times New Roman under a custom single-byte
    encoding, so the digits are not ASCII in the content stream.  This script derives the
    byte codes from the font's own ToUnicode CMap rather than guessing them.
  * It then asserts that the glyphs for "2" and "3" have identical advance widths in that
    font's /Widths array.  They do (500 each), so exchanging one for the other moves no
    following glyph: the line is typeset exactly as before.
  * The patched content stream recompresses smaller than the original, so the replacement
    object is padded with dictionary whitespace to occupy the *same byte span* as the one
    it replaces.  Every other file offset, and the cross-reference stream, stay valid.

The script refuses to run if the string it expects is absent, so re-running it on an
already-corrected file fails loudly instead of quietly doing nothing.

Usage:
    python scripts/17_biorxiv_v2/patch_merged_pdf_doi.py --in SRC.pdf --out DST.pdf
"""
from __future__ import annotations

import argparse
import re
import sys
import zlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from revision_edits import PLACEHOLDERS

OLD_ID = "22640703"
NEW_ID = PLACEHOLDERS["ZENODO_DOI"].rsplit(".", 1)[-1]


def inflate(body: bytes) -> bytes | None:
    m = re.search(rb"stream\r?\n", body)
    if not m:
        return None
    data = body[m.end():]
    i = data.rfind(b"endstream")
    if i >= 0:
        data = data[:i]
    if b"FlateDecode" not in body[:m.start()]:
        return data
    try:
        return zlib.decompress(data)
    except zlib.error:
        return None


def objects(raw: bytes) -> dict[int, tuple[int, int, int, bytes]]:
    out = {}
    for m in re.finditer(rb"(?m)^(\d+)\s+(\d+)\s+obj\b", raw):
        end = raw.find(b"endobj", m.end())
        out[int(m.group(1))] = (m.start(), m.end(), end, raw[m.end():end])
    return out


def expand(raw: bytes, top) -> dict[int, bytes]:
    allobj = {n: b[3] for n, b in top.items()}
    for _, (_, _, _, body) in list(top.items()):
        if b"/ObjStm" not in body[:500]:
            continue
        d = inflate(body)
        if d is None:
            continue
        n = int(re.search(rb"/N\s+(\d+)", body).group(1))
        first = int(re.search(rb"/First\s+(\d+)", body).group(1))
        nums = d[:first].split()
        for i in range(n):
            onum, off = int(nums[2 * i]), int(nums[2 * i + 1])
            nxt = int(nums[2 * i + 3]) + first if i + 1 < n else len(d)
            allobj.setdefault(onum, d[first + off:nxt])
    return allobj


def tounicode(cmap: bytes) -> dict[int, str]:
    m: dict[int, str] = {}
    for blk in re.findall(rb"beginbfchar(.*?)endbfchar", cmap, re.S):
        for a, b in re.findall(rb"<([0-9A-Fa-f]+)>\s*<([0-9A-Fa-f]+)>", blk):
            m[int(a, 16)] = bytes.fromhex(b.decode()).decode("utf-16-be", "replace")
    for blk in re.findall(rb"beginbfrange(.*?)endbfrange", cmap, re.S):
        for a, b, c in re.findall(rb"<([0-9A-Fa-f]+)>\s*<([0-9A-Fa-f]+)>\s*<([0-9A-Fa-f]+)>", blk):
            lo, hi, st = int(a, 16), int(b, 16), int(c, 16)
            for i in range(lo, hi + 1):
                m[i] = chr(st + (i - lo))
    return m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", required=True)
    ap.add_argument("--out", dest="dst", required=True)
    args = ap.parse_args()

    raw = Path(args.src).read_bytes()
    top = objects(raw)
    allobj = expand(raw, top)

    # Every page whose text contains the old identifier, found through the font encoding.
    cat = next(n for n, b in allobj.items() if re.search(rb"/Type\s*/Catalog", b))
    root = int(re.search(rb"/Pages\s+(\d+)\s+0\s+R", allobj[cat]).group(1))
    order: list[int] = []

    def walk(n: int) -> None:
        b = allobj[n]
        if re.search(rb"/Type\s*/Page\b", b) and b"/Kids" not in b:
            order.append(n)
            return
        kids = re.search(rb"/Kids\s*\[(.*?)\]", b, re.S)
        for k in re.finditer(rb"(\d+)\s+0\s+R", kids.group(1)):
            walk(int(k.group(1)))

    walk(root)

    patched_objects = 0
    new_raw = raw
    for page_no, pnum in enumerate(order, start=1):
        pb = allobj[pnum]
        cm = re.search(rb"/Contents\s+(\d+)\s+0\s+R", pb)
        rm = re.search(rb"/Resources\s+(\d+)\s+0\s+R", pb)
        if not cm or not rm:
            continue
        cnum = int(cm.group(1))
        data = inflate(top[cnum][3])
        if data is None:
            continue

        fonts = re.search(rb"/Font\s*<<(.*?)>>", allobj[int(rm.group(1))], re.S)
        if not fonts:
            continue

        replacements: list[tuple[bytes, bytes]] = []
        for fname, fnum in re.findall(rb"/(\w+)\s+(\d+)\s+0\s+R", fonts.group(1)):
            fb = allobj[int(fnum)]
            tu = re.search(rb"/ToUnicode\s+(\d+)\s+0\s+R", fb)
            if not tu:
                continue
            umap = tounicode(inflate(top[int(tu.group(1))][3]) or b"")
            inv: dict[str, int] = {}
            for code, ch in umap.items():
                inv.setdefault(ch, code)
            if not all(d in inv for d in set(OLD_ID + NEW_ID)):
                continue
            old_bytes = bytes(inv[d] for d in OLD_ID)
            new_bytes = bytes(inv[d] for d in NEW_ID)
            if old_bytes not in data:
                continue

            # Advance widths must match, or the rest of the line would shift.
            fc = int(re.search(rb"/FirstChar\s+(\d+)", fb).group(1))
            wm = re.search(rb"/Widths\s*(?:\[(.*?)\]|(\d+)\s+0\s+R)", fb, re.S)
            arr = wm.group(1) if wm.group(1) else allobj[int(wm.group(2))]
            widths = [int(float(x)) for x in re.findall(rb"[-\d.]+", arr)]
            for a, b in zip(old_bytes, new_bytes):
                wa, wb = widths[a - fc], widths[b - fc]
                if wa != wb:
                    raise SystemExit(
                        f"page {page_no} font {fname.decode()}: replacing code {a} with "
                        f"{b} changes the advance width ({wa} -> {wb}). Refusing: the "
                        "line would be re-spaced. Re-export the PDF instead."
                    )
            replacements.append((old_bytes, new_bytes))

        if not replacements:
            continue

        patched = data
        n_here = 0
        for old_bytes, new_bytes in replacements:
            n_here += patched.count(old_bytes)
            patched = patched.replace(old_bytes, new_bytes)
        print(f"  page {page_no}: {n_here} occurrence(s) of {OLD_ID} -> {NEW_ID}")

        start, _, end, body = top[cnum]
        span = raw[start:end + len(b"endobj")]
        comp = zlib.compress(patched, 9)
        head = b"%d 0 obj\n<<\n/Filter /FlateDecode\n/Length %d" % (cnum, len(comp))
        tail = b"\n>>\nstream\n" + comp + b"\nendstream\nendobj"
        pad = len(span) - (len(head) + len(tail))
        if pad < 0:
            raise SystemExit(
                f"object {cnum}: the replacement is {-pad} bytes longer than the original. "
                "Rewriting it in place would invalidate every later offset. Refusing."
            )
        # Dictionary whitespace is free, so the object keeps its exact byte span and the
        # cross-reference stream that points at it stays correct.
        replacement = head + b" " * pad + tail
        assert len(replacement) == len(span)
        new_raw = new_raw[:start] + replacement + new_raw[start + len(span):]
        patched_objects += 1

    if not patched_objects:
        raise SystemExit(f"no page carries {OLD_ID}; nothing to correct in {args.src}")
    assert len(new_raw) == len(raw), "file length changed"
    Path(args.dst).write_bytes(new_raw)
    print(f"wrote {args.dst} ({patched_objects} content stream(s), {len(new_raw)} bytes, "
          f"file length unchanged)")


if __name__ == "__main__":
    main()
