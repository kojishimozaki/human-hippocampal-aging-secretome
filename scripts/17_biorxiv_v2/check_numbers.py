#!/usr/bin/env python3
"""Compare the numeric content of the English manuscript and its Japanese mirror.

Iron rule 2: numbers must not drift between surfaces. This extracts every numeric token
from both documents and reports any that appear a different number of times, so a
translation slip or a stale figure cannot hide in the Japanese file.

Usage:
    PROJ=$(pwd) python scripts/17_biorxiv_v2/check_numbers.py
"""

from __future__ import annotations

import os
import re
import sys
from collections import Counter
from pathlib import Path

PROJ = Path(os.environ.get("PROJ", Path(__file__).resolve().parents[2]))
EN = PROJ / "manuscript/biorxiv_v2/manuscript_biorxiv_v2.md"
JA = PROJ / "manuscript/biorxiv_v2/manuscript_biorxiv_v2_ja.md"

TOKEN = re.compile(r"\d[\d,]*(?:\.\d+)?")


def tokens(path: Path) -> Counter:
    text = path.read_text(encoding="utf-8")
    # drop the Japanese mirror's provenance header (paragraph indices and match ratios)
    text = re.sub(r"\A(?:<!--.*?-->\n|>.*\n|\n)+", "", text)
    # drop the shared reference list: identical English in both files
    text = re.split(r"^##+ (?:References|参考文献)\s*$", text, flags=re.M)[0]
    return Counter(m.group(0).rstrip(".").replace(",", "") for m in TOKEN.finditer(text))


def measured(tok: str) -> bool:
    """Tokens that carry a measured value, as opposed to prose counts and list indices.

    Small bare integers legitimately differ in frequency between the two languages
    ("the five niche classes" versus "5 つのニッチ細胞クラス", "Fig. 2A,B" versus
    "図 2A、B"), and scientific notation tokenises into its parts. Decimals, accessions,
    and four-digit-plus counts are the ones a translation slip would corrupt silently.
    """
    if re.fullmatch(r"[12]\d{5}", tok):   # GEO accession digits: an identifier, not a value
        return False
    return "." in tok or len(tok) >= 4


def main() -> int:
    if not JA.exists():
        # The Japanese mirror is maintained in the working repository and is not part of
        # the public release, so there is nothing to compare against here.
        print(f"Japanese mirror not present ({JA.name}); parity check skipped.")
        return 0
    en, ja = tokens(EN), tokens(JA)
    hard, soft = [], []
    for tok in sorted(set(en) | set(ja), key=lambda t: (-max(en[t], ja[t]), t)):
        if en[tok] != ja[tok]:
            (hard if measured(tok) else soft).append((tok, en[tok], ja[tok]))
    print(f"EN numeric tokens: {sum(en.values())} ({len(en)} distinct)")
    print(f"JA numeric tokens: {sum(ja.values())} ({len(ja)} distinct)")
    acc = [(t, a, b) for t, a, b in list(soft) if re.fullmatch(r"[12]\d{5}", t)]
    soft = [x for x in soft if x not in acc]
    if acc:
        print("\naccession mentions differ in frequency (identifier, not a value): "
              + ", ".join(f"GSE{t}({a}/{b})" for t, a, b in acc))
    if soft:
        print(f"\n{len(soft)} bare small integer(s) differ in frequency "
              f"(prose counting, not a value mismatch):")
        print("  " + ", ".join(f"{t}({a}/{b})" for t, a, b in soft))
    if not hard:
        print("\nOK: every measured value appears the same number of times in both files.")
        return 0
    print(f"\n{len(hard)} MEASURED value(s) differ in count (EN vs JA):")
    for tok, a, b in hard:
        print(f"  {tok:>12s}  EN={a:<3d} JA={b}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
