# The manuscript

Generated, not hand-edited. The English text and the Word file are built from one
definition — `scripts/17_biorxiv_v2/revision_edits.py` — so they cannot drift apart.

| File | What it is |
|---|---|
| `shimozaki-aging-secretome-v2.pdf` | **The manuscript as it is meant to be read**: 54 pages, 39 of text followed by Figures 1–6 and S1–S9. |
| `manuscript_biorxiv_v2.md` | The same text as markdown. This is the text of record: every number in it is checked against `results/` by `make verify-manuscript`. |
| `shimozaki-aging-secretome-v2.docx` | The same text as a Word file, byte-checked against the markdown by the same target. |
| `figures/Figure_*.pdf` | The fifteen composed figures, one vector PDF each, in the order they appear in the merged PDF. |
| `../supplementary_tables/` | Tables S1, S2 and SR1, and the workbook that holds all three. |

The manuscript has not been through peer review. It reached this form through six
adversarial audit rounds; the pre-specification each round was held to is in
`../../preregistration/`.

## To change the text

Edit `scripts/17_biorxiv_v2/revision_edits.py`, never the files here, then:

```bash
conda activate bio
export PROJ=$(pwd)
python scripts/17_biorxiv_v2/build_revision.py         # -> .docx + English markdown
make verify-manuscript                                 # 76 checks against the tables
```

Each builder anchors every edit on a prefix of the paragraph it targets and fails loudly
if the base document has moved, so a silent mis-edit is not possible. `build_revision.py`
starts from a superseded Word file that is deliberately not published, so it runs only in
the author's working repository; the built text of record is committed here either way,
and every data-driven target runs from a bare clone.

The merged PDF is assembled by hand from the Word file and the fifteen figure PDFs, so it
is the one artefact here that no script regenerates. `make verify-manuscript` therefore
reads it back and confirms it has 54 pages and sends readers to the same three addresses
the text of record does.
