#!/usr/bin/env bash
# Build manuscript PDFs without gstack.
#
# The repo's documented route (INTEGRATION_MANIFEST.md) uses the gstack `make-pdf`
# binary, which is a macOS arm64 build and cannot run on this Linux host. This script
# uses tools that are present: pandoc (in the `sgz_r` env) for markdown -> HTML, and the
# Playwright Chromium already installed under ~/.cache/ms-playwright for HTML -> PDF.
#
# One preprocessing step is needed and is deliberately NOT applied to the manuscript
# itself: heading lines in manuscript.md carry 0-2 leading spaces inconsistently, and
# pandoc renders the indented ones as literal "# ..." paragraphs instead of headings.
# The leading spaces are stripped in a temporary copy only; the manuscript is untouched.
#
# It also assembles manuscript/manuscript_with_figures.pdf. That bundle used to be built by
# hand from a pdfunite command pasted in INTEGRATION_MANIFEST.md, and it silently went stale:
# the 2026-08-26 codex audit found Fig. S1 and S2 missing from it (13 of 15 figures, 42 pp)
# even though both legends were in the manuscript. The figure list now lives here, next to the
# build, and the script REFUSES to write the bundle unless the list matches the figure legends
# in manuscript.md one-for-one and in order.
#
# Usage:  bash scripts/15_manuscript_pdf/build_manuscript_pdf.sh
# Output: manuscript/manuscript{,_ja,_with_figures}.pdf
set -euo pipefail
: "${PROJ:=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
PANDOC="${PANDOC:-$HOME/miniforge3/envs/sgz_r/bin/pandoc}"
CHROME="${CHROME:-$HOME/.cache/ms-playwright/chromium-1208/chrome-linux64/chrome}"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT

[ -x "$PANDOC" ] || { echo "pandoc not found at $PANDOC" >&2; exit 1; }
[ -x "$CHROME" ] || { echo "chromium not found at $CHROME" >&2; exit 1; }

cat > "$WORK/print.css" <<'CSS'
@page { size: A4; margin: 20mm 18mm 20mm 18mm; }
body { font-family: "Source Serif 4","Noto Serif","Noto Serif CJK JP",Georgia,serif;
       font-size: 10.5pt; line-height: 1.55; color: #111; max-width: none; }
h1 { font-size: 15pt; line-height: 1.3; margin: 0 0 0.7em; }
h2 { font-size: 12.5pt; margin: 1.4em 0 0.5em; border-bottom: 1px solid #bbb; padding-bottom: 2px; }
h3 { font-size: 11pt; margin: 1.1em 0 0.4em; }
h2, h3 { page-break-after: avoid; }
p { margin: 0 0 0.55em; text-align: justify; }
strong { font-weight: 650; }
code { font-family: "DejaVu Sans Mono",monospace; font-size: 0.88em; background: #f3f3f3; padding: 0 2px; }
pre { background: #f6f6f6; padding: 6px 8px; font-size: 8.5pt; white-space: pre-wrap; overflow-wrap: break-word; }
table { border-collapse: collapse; font-size: 8.8pt; margin: 0.6em 0; width: 100%; }
th, td { border: 1px solid #bbb; padding: 2px 5px; vertical-align: top; }
th { background: #eee; }
ol, ul { margin: 0.3em 0 0.6em; padding-left: 1.5em; }
li { margin-bottom: 0.35em; text-align: justify; }
blockquote { border-left: 3px solid #999; margin: 0.6em 0; padding: 0.2em 0 0.2em 0.8em; color: #333; }
CSS

# cover_letter is built here too: its PDF had gone stale against cover_letter.md and still carried
# the withdrawn "independently replicated" claim in the title (codex audit 2026-08-26).
for f in ${DOCS:-manuscript manuscript_ja cover_letter cover_letter_ja manuscript_v2_niche_focus manuscript_v2_niche_focus_ja}; do
  sed -E 's/^[[:space:]]{1,3}(#{1,6}[[:space:]])/\1/' "$PROJ/manuscript/$f.md" > "$WORK/$f.md"
  "$PANDOC" "$WORK/$f.md" -f markdown+pipe_tables+east_asian_line_breaks -t html5 -s --metadata title="" \
            -c print.css -o "$WORK/$f.html"
  "$CHROME" --headless --disable-gpu --no-sandbox --hide-scrollbars \
            --print-to-pdf="$WORK/$f.pdf" --no-pdf-header-footer \
            --virtual-time-budget=20000 "file://$WORK/$f.html" >/dev/null 2>&1
  [ -s "$WORK/$f.pdf" ] || { echo "failed to render $f" >&2; exit 1; }
  cp "$WORK/$f.pdf" "$PROJ/manuscript/$f.pdf"
  echo "wrote manuscript/$f.pdf ($(pdfinfo "$PROJ/manuscript/$f.pdf" | awk '/^Pages/{print $2}') pages)"
done

# ---- figure-embedded bundle -------------------------------------------------------------
# One entry per figure, in manuscript-legend order. Keep in step with "## Figure legends".
FIG_LABELS=(
  "Fig. 1"      "Fig. S1"      "Fig. S2"     "Fig. 2"        "Fig. S3"
  "Fig. 3"      "Fig. S(E1)"   "Fig. 4"      "Fig. S(M1)"    "Fig. 5"
  "Fig. 6"      "Fig. S(N)"    "Fig. S(R1)"  "Fig. S(Rcv)"   "Fig. S(RT1)"
)
FIG_FILES=(
  figures/figure1/figure1.pdf
  figures/figure1/supp_S1_QC.pdf
  figures/figure1/supp_S2_composition_allgroups.pdf
  figures/figure2/figure2.pdf
  figures/figure2/supp_S3_csf_direction.pdf
  figures/figure3/replication_GSE278576.pdf
  figures/figureSE1/figureSE1.pdf
  figures/figure4/figure4.pdf
  figures/figureM1/figure_SM1_methylation.pdf
  figures/figure5/figure5.pdf
  figures/figure6/figure6.pdf
  figures/figureS_neuron/figureS_neuron.pdf
  figures/figureSR1/figureSR1.pdf
  figures/figureS_receiver/figureS_receiver.pdf
  figures/figureSRT1/figureSRT1.pdf
)

# Gate: the list above must equal the figure legends in manuscript.md, in order.
python3 - "$PROJ" "${FIG_LABELS[@]}" <<'PYGATE'
import io, re, sys
proj, expected = sys.argv[1], sys.argv[2:]
sec = (io.open(f"{proj}/manuscript/manuscript.md", encoding="utf-8").read()
       .split("## Figure legends", 1)[1].split("## Supplementary tables", 1)[0])
found = re.findall(r"\*\*(Fig\.\s*[^\s\u2014]+)\s*\u2014", sec)
found = [re.sub(r"\s+", " ", f).strip() for f in found]
if found != expected:
    print("FIGURE LIST OUT OF STEP WITH THE MANUSCRIPT LEGENDS", file=sys.stderr)
    print("  manuscript.md legends :", found, file=sys.stderr)
    print("  build-script list     :", expected, file=sys.stderr)
    sys.exit(1)
print(f"figure list matches the {len(found)} legends in manuscript.md")
PYGATE

for f in "${FIG_FILES[@]}"; do
  [ -s "$PROJ/$f" ] || { echo "missing figure: $f" >&2; exit 1; }
done

# The v2 (niche-focus) version carries the same figure legends verbatim, so it takes the same
# bundle. Both are built when the default DOCS list is used.
for base in manuscript manuscript_v2_niche_focus; do
  [ -s "$PROJ/manuscript/$base.pdf" ] || continue
  ( cd "$PROJ" && pdfunite "manuscript/$base.pdf" "${FIG_FILES[@]}" "manuscript/${base}_with_figures.pdf" )
  echo "wrote manuscript/${base}_with_figures.pdf ($(pdfinfo "$PROJ/manuscript/${base}_with_figures.pdf" | awk '/^Pages/{print $2}') pages = text + ${#FIG_FILES[@]} figures)"
done
