#!/usr/bin/env python
"""Conclusion model figure (graphical conclusion / "ponchi-e") for the paper.

A conceptual schematic of the study's conclusion — NOT a data plot:
  TOP    biological model: normal aging shifts the niche sender cells' secreted output
         (astrocyte matrix ↑, microglial inflammation ↑, endothelial vascular ↓)
         into a senescence-aligned (SASP-like) secretome.
  BOTTOM the bounded conclusion: what is ESTABLISHED (validated + reproduced + bounded)
         vs what remains OPEN (the upstream regulator) — the paper's honesty thesis.

NOTE (2026-08-26, codex audit): this figure must not claim multi-modal "coordination". The
DNA-methylation layer failed four pre-submission controls and its RNA-directed coordination
reading is WITHDRAWN (manuscript Results / Fig S(M1)); the surviving regulatory observation is a
gene-level, matched-null RNA<->ATAC direction concordance in ASTROCYTES ALONE. The reproduction
in the second cohort is likewise not demonstrably specific against an effect-size-matched
background, so the manuscript says "reproduces", not "replicated".

Quantitative claims are conceptual summaries that trace to the data figures / manuscript
(authoritative); like the overview infographics this figure restates rather than recomputes them.
English (matches Fig 1–6). Shared style: scripts/07_figures/_pubstyle.py. Renders PDF + PNG + SVG.
Run:  conda activate bio && python scripts/07_figures/11_conclusion_model.py
"""
import os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _pubstyle as PS

PS.apply()

def _tick(ax, x, y, size, color, zorder=6):
    """Draw a check mark as line segments. Liberation Sans has no U+2713 and
    mathtext \\checkmark pulls in STIXGeneral, which would put a second
    typeface in the figure (see skill: manuscript-production.md §12)."""
    h = size / 900.0
    ax.plot([x - 1.1 * h, x - 0.25 * h, x + 1.2 * h],
            [y, y - 1.15 * h, y + 1.25 * h],
            color=color, lw=1.7, solid_capstyle="round",
            solid_joinstyle="miter", zorder=zorder, clip_on=False)

P = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
OUT = f"{P}/figures/overview"; os.makedirs(OUT, exist_ok=True)

INK, SUB = PS.INK, PS.SUB
AST, MIC, END = PS.CT["Astro"], PS.CT["Micro"], PS.CT["Endo"]
GREEN, GREY, MUST = PS.GREEN, "#7c8696", PS.MUSTARD
NAVY = "#202d4a"
W, H = 8.8, 6.95
ASPECT = H / W

fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(-0.03, 1); ax.axis("off")


def card(x, y, w, h, fc, ec, lw=1.4, r=0.020, alpha=1.0, z=1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                fc=fc, ec=ec, lw=lw, alpha=alpha, zorder=z,
                                mutation_aspect=ASPECT))


def T(x, y, s, size=9, color=INK, weight="normal", ha="left", va="center", style="normal"):
    ax.text(x, y, s, fontsize=size, color=color, fontweight=weight, ha=ha, va=va,
            style=style, zorder=6)


# ───────────────────────── title
T(0.5, 0.965, "Aging reprograms the secreted-protein output of the human hippocampal niche",
  size=15, weight="bold", ha="center")
T(0.5, 0.927, "An externally validated signature that reproduces in a second cohort",
  size=10.5, color=SUB, ha="center")

# ───────────────────────── TOP: biological model
# young-adult niche
card(0.035, 0.565, 0.140, 0.255, "#f3f5f9", "#c8d0dc", 1.2)
T(0.105, 0.785, "Young-adult", size=10, weight="bold", ha="center")
T(0.105, 0.762, "niche", size=10, weight="bold", ha="center")
T(0.105, 0.715, "support cells —", size=8, color=SUB, ha="center")
T(0.105, 0.695, "homeostatic", size=8, color=SUB, ha="center", style="italic")
T(0.105, 0.677, "secreted output", size=8, color=SUB, ha="center", style="italic")
for i, c in enumerate([AST, MIC, END]):
    ax.plot([0.082 + i * 0.024], [0.612], "o", ms=8, color=c, alpha=0.85, zorder=5)

# aging arrow
ax.annotate("", xy=(0.252, 0.693), xytext=(0.180, 0.693),
            arrowprops=dict(arrowstyle="-|>", lw=2.6, color=INK))
T(0.216, 0.726, "normal aging", size=8.5, weight="bold", ha="center")
T(0.216, 0.662, "YA 8  vs  HA 9", size=7.5, color=SUB, ha="center")

# sender-cell rows
NX = 0.300
rows = [(0.792, AST, "Astrocyte", "extracellular matrix / matricellular  ↑",
         "TNC · COL21A1 · CCN2 · IGFBP5", None),
        (0.690, MIC, "Microglia", "inflammatory cytokines  ↑",
         "IL15 · SERPINE1 · APOE", "+ microgliosis"),
        (0.588, END, "Endothelium", "vascular / angiogenic factors  ↓",
         "KDR · PLAT · LAMA4", None)]
for (y, c, name, head, genes, extra) in rows:
    ax.plot([NX], [y], "o", ms=15, color=c, zorder=5)
    T(NX + 0.026, y + 0.016, name, size=9.5, weight="bold", color=c)
    T(NX + 0.026, y - 0.010, head, size=8.6, color=INK)
    T(NX + 0.026, y - 0.036, genes, size=7.5, color=SUB)
    if extra:
        T(NX + 0.026, y - 0.058, extra, size=7.2, color=MIC, style="italic")
    ax.annotate("", xy=(0.726, 0.693), xytext=(0.610, y),
                arrowprops=dict(arrowstyle="-|>", lw=0.9, color="#bcc4d0",
                                mutation_scale=9, shrinkA=2, shrinkB=3))

# convergence: SASP-like secretome
card(0.728, 0.572, 0.237, 0.243, "#fbf3ec", MUST, 1.5)
T(0.846, 0.778, "Senescence-aligned", size=9.6, weight="bold", color="#9a6a12", ha="center")
T(0.846, 0.751, "(SASP-like) secretome", size=9.6, weight="bold", color="#9a6a12", ha="center")
T(0.846, 0.706, "the niche's altered", size=8, color=SUB, ha="center")
T(0.846, 0.686, "signalling environment", size=8, color=SUB, ha="center")
T(0.846, 0.636, "→ the up-arm tracks", size=7.7, color=INK, ha="center")
T(0.846, 0.617, "cognitive decline (brain proteome)", size=7.7, color=INK, ha="center")

# ───────────────────────── divider
ax.plot([0.035, 0.965], [0.527, 0.527], color="#d6dbe4", lw=1.0)
T(0.5, 0.500, "What is established  —  and what remains open", size=10.5, weight="bold", ha="center")

# ───────────────────────── BOTTOM-LEFT: established
card(0.035, 0.070, 0.452, 0.400, "#eef7f1", GREEN, 1.6)
_tick(ax, 0.064, 0.438, 13, GREEN)
T(0.082, 0.438, "ESTABLISHED", size=11, weight="bold", color=GREEN)
T(0.300, 0.438, "validated · reproduced · bounded", size=7.8, color=SUB)
solid = [("Senescence-aligned & protein-validated",
          "SASP-Atlas secretome + brain proteome→cognition (Wingo);\nCSF concordance suggestive"),
         ("Reproduces in a 2nd human cohort",
          "40 donors, 30 in the aging contrast: ρ = 0.24, robust to sex &\n"
          "continuous age — reproduction, not demonstrable specificity"),
         ("Chromatin concordance — astrocytes only",
          "gene-level, matched-null RNA↔ATAC direction concordance in astrocytes\n"
          "— supporting. DNA-methylation failed 4 controls → withdrawn"),
         ("Niche-enriched; no DG concentration detected",
          "no detectable dentate-gyrus concentration; not reproduced in mature neurons; differs\n"
          "from AD's broad chromatin dysregulation")]
y = 0.392
for lead, sub in solid:
    _tick(ax, 0.067, y, 10.5, GREEN)
    T(0.086, y, lead, size=8.8, weight="bold")
    T(0.086, y - 0.030, sub, size=7.6, color=SUB, va="top")
    y -= 0.083

# ───────────────────────── BOTTOM-RIGHT: open
card(0.513, 0.070, 0.452, 0.400, "#f1f3f6", GREY, 1.6)
T(0.538, 0.438, "?", size=14, weight="bold", color="#5f6b7d")
T(0.560, 0.438, "OPEN", size=11, weight="bold", color="#5f6b7d")
T(0.640, 0.438, "the explicit question for future work", size=7.8, color=SUB)
opn = [("Upstream regulator unresolved",
        "cis-TF footprinting FDR-null (nominally reversed);\nautocrine ligand–receptor analysis null"),
       ("“No evidence for,” not a disproof",
        "under-powered at n = 8/9 donors; modest effect sizes;\nsingle primary multiome cohort"),
       ("Requires orthogonal data",
        "perturbation · CUT&RUN / ChIP for NF-κB / IRF ·\ndepth-matched bulk ATAC")]
y = 0.392
for lead, sub in opn:
    T(0.540, y, "?", size=10.5, weight="bold", color="#5f6b7d", ha="center")
    T(0.560, y, lead, size=8.8, weight="bold")
    T(0.560, y - 0.030, sub, size=7.6, color=SUB, va="top")
    y -= 0.108

# ───────────────────────── takeaway band
card(0.035, -0.012, 0.930, 0.058, NAVY, NAVY, 0, r=0.014)
T(0.5, 0.017, "A rigorously bounded, externally validated niche-secretome "
  "aging signature — its regulatory origin an explicit open question.",
  size=9, weight="bold", color="white", ha="center")

# footer
T(0.965, -0.027, "Conceptual conclusion figure · quantitative claims trace to the data figures / "
  "manuscript (authoritative) · scripts/07_figures/11_conclusion_model.py",
  size=6.2, color="#9aa3b0", ha="right")

for ext in ("pdf", "png", "svg"):
    fig.savefig(f"{OUT}/conclusion_model.{ext}", dpi=300 if ext == "png" else None,
                bbox_inches="tight", facecolor="white")
# strip trailing whitespace matplotlib emits on some SVG <path> lines (keeps `git diff --check` clean;
# whitespace/newlines are interchangeable separators in SVG path data, so per-line rstrip is safe)
_svg = f"{OUT}/conclusion_model.svg"
with open(_svg) as _fh: _lines = _fh.read().split("\n")
with open(_svg, "w") as _fh: _fh.write("\n".join(s.rstrip() for s in _lines))
plt.close(fig)
print(f"wrote {OUT}/conclusion_model.{{pdf,png,svg}}")
