# bioRxiv v2 — the current manuscript

Generated, not hand-edited. Everything here is built from one definition so that the
English text, the Japanese mirror and the Word file cannot drift apart.

| File | What it is |
|---|---|
| `manuscript_biorxiv_v2.md` | **English text of record.** |
| `manuscript_biorxiv_v2_ja.md` | Japanese mirror. Its header lists the paragraphs whose wording still tracks the pre-Word English (no number or claim differs). |
| `abstract_250w.md` | 248-word abstract for journals that cap the abstract. Not the bioRxiv version. |
| `biorxiv_resubmission_note.md` | Draft note to the bioRxiv screening team. Three bracketed items to fill. |
| `../../submission/v2/shimozaki-aging-secretome-v2.docx` | Word file for submission. |
| `../../submission/v2/figures/` | Figures with corrected supplementary numbering. |

## To change the manuscript

Edit `scripts/17_biorxiv_v2/revision_edits.py`, never these files, then:

```bash
conda activate bio
export PROJ=$(pwd)
python scripts/17_biorxiv_v2/build_revision.py         # -> .docx + English markdown
python scripts/17_biorxiv_v2/build_japanese_mirror.py  # -> Japanese mirror
python scripts/17_biorxiv_v2/check_numbers.py          # EN/JA measured-value parity
```

Each builder anchors its edits on a prefix of the paragraph it targets and fails loudly if
the base document has moved, so a silent mis-edit is not possible.

`build_revision.py` edits the superseded submission DOCX, which is deliberately not
published, so `make biorxiv-v2` runs only where that file is still on disk. It says so
if it is missing. The text of record here is committed either way, and every
data-driven target runs from a bare clone.

## Before submitting

1. Fill `PLACEHOLDERS` in `scripts/17_biorxiv_v2/revision_edits.py` with the repository
   URL, the Zenodo DOI and the commit hash, then rebuild. The build prints a warning
   while they are still placeholders.
2. Export the DOCX to PDF and append the figures in the order
   `Figure_1` … `Figure_6`, `Figure_S1` … `Figure_S9`.
3. Confirm the licence in `LICENSE` (MIT was chosen as a default, not by the author).

The revision itself is recorded in `audit_log/2026-09-04_biorxiv_v2_reframe/RESOLUTION.md`.
