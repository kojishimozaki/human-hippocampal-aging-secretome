#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Concise graphical abstract (English, at-a-glance).

English sibling of 09_graphical_abstract.py — same 3-act layout
(Question -> Finding (3 arms) -> Evidence status), translated text, Latin font
stack, *_en output filenames so the Japanese version is preserved. Numbers per
manuscript/manuscript.md. Renders SVG -> PNG + PDF.

Run:  conda activate bio && python scripts/07_figures/09b_graphical_abstract_en.py
"""
import os, html

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
OUTDIR = os.path.join(PROJ, "figures", "overview")
os.makedirs(OUTDIR, exist_ok=True)

W, H = 1640, 900
FONT = "Arial, Helvetica, Liberation Sans, DejaVu Sans, sans-serif"

INK, SUB, LINE = "#1f2733", "#5b6573", "#d6dbe4"
PAGEBG, NAVY, NAVYSUB = "#ffffff", "#202d4a", "#aeb8cc"
GREEN, TEAL, MUSTARD, GREY, SLATE = "#2c8a5b", "#1f8a8a", "#c0941f", "#8a94a6", "#2f3e5c"

C_AST, C_MIC, C_END = "#e0962f", "#c5453b", "#3b72c9"
T_AST, T_MIC, T_END = "#fbeed7", "#f7dcd9", "#dde9fa"   # light tints
UP, DN = "#c5453b", "#3b72c9"

svg = []
def add(s): svg.append(s)
def esc(s): return html.escape(str(s), quote=True)

def rect(x, y, w, h, fill="#ffffff", stroke=LINE, sw=1.2, rx=12):
    add(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" ry="{rx}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')

def line(x1, y1, x2, y2, stroke=LINE, sw=1.2):
    add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" stroke-width="{sw}"/>')

def text(x, y, s, size=15, fill=INK, weight="normal", anchor="start", spacing=None):
    ls = f' letter-spacing="{spacing}"' if spacing else ""
    add(f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
        f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}"{ls}>{esc(s)}</text>')

def lines(x, y, rows, size=15, fill=INK, lh=22, weight="normal", anchor="start"):
    for i, r in enumerate(rows):
        text(x, y + i*lh, r, size=size, fill=fill, weight=weight, anchor=anchor)

def eyebrow(x, y, s, color):
    add(f'<rect x="{x:.1f}" y="{y-12:.1f}" width="5" height="16" rx="2.5" fill="{color}"/>')
    text(x+13, y+2, s, size=14.5, fill=color, weight="bold", spacing="1")

def rarrow(cx, cy, color=SLATE, w=22, h=26):
    add(f'<path d="M{cx-w/2:.1f} {cy-h/2:.1f} L{cx+w/2:.1f} {cy:.1f} L{cx-w/2:.1f} {cy+h/2:.1f}" '
        f'fill="none" stroke="{color}" stroke-width="4" stroke-linecap="round" stroke-linejoin="round"/>')

def darrow(cx, y, color=SLATE, w=26, h=12):
    add(f'<path d="M{cx-w/2:.1f} {y:.1f} L{cx:.1f} {y+h:.1f} L{cx+w/2:.1f} {y:.1f}" '
        f'fill="none" stroke="{color}" stroke-width="3.4" stroke-linecap="round" stroke-linejoin="round"/>')

def cell_icon(cx, cy, color, tint):
    add(f'<ellipse cx="{cx:.1f}" cy="{cy:.1f}" rx="27" ry="21" fill="{tint}" stroke="{color}" stroke-width="2.2"/>')
    add(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="8.5" fill="{color}" opacity="0.85"/>')

def updown(x, y, color, up=True, size=22):
    text(x, y, "↑" if up else "↓", size=size, fill=color, weight="bold", anchor="middle")

def verdict(x, y, kind, r=12):
    col = {"ok": GREEN, "weak": MUSTARD, "null": GREY}[kind]
    add(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{col}"/>')
    if kind == "ok":
        add(f'<path d="M{x-5.5:.1f} {y:.1f} l4 4.5 l7 -9" stroke="#fff" stroke-width="2.6" '
            f'fill="none" stroke-linecap="round" stroke-linejoin="round"/>')
    elif kind == "weak":
        text(x, y+5, "?", size=15, fill="#fff", weight="bold", anchor="middle")
    else:
        add(f'<path d="M{x-5:.1f} {y-5:.1f} l10 10 M{x+5:.1f} {y-5:.1f} l-10 10" '
            f'stroke="#fff" stroke-width="2.6" stroke-linecap="round"/>')

# ================================================================= BUILD
add(f'<rect x="0" y="0" width="{W}" height="{H}" fill="{PAGEBG}"/>')

# ---- title ---------------------------------------------------------------
text(46, 50, "GRAPHICAL ABSTRACT  —  Aging of the human hippocampal niche secretome", size=14.5,
     fill=SUB, weight="bold", spacing="0.5")
text(44, 92, "Aging reprograms the secreted-protein output of the human hippocampal niche",
     size=29, fill=INK, weight="bold")
text(46, 122, "An externally validated SASP-like signature reproduced in a second cohort — mechanism unresolved",
     size=17, fill=SUB)
line(46, 138, W-46, 138, stroke=LINE, sw=1.4)

# zone geometry
ZA_x, ZA_w = 46, 432
ZB_x, ZB_w = 566, 540
ZC_x, ZC_w = 1166, 432
top, bot = 168, 792
ar1 = (ZA_x+ZA_w + ZB_x)/2
ar2 = (ZB_x+ZB_w + ZC_x)/2
midy = (top+bot)/2

# inter-zone arrows
rarrow(ar1, midy, SLATE); rarrow(ar2, midy, SLATE)

# ======================= ZONE A : question + data =========================
eyebrow(ZA_x, top+6, "QUESTION & DATA", SLATE)
ry = top + 34
# reframe: old (muted)
rect(ZA_x, ry, ZA_w, 64, fill="#f0f1f4", stroke="#cfd5df", rx=12)
verdict(ZA_x+28, ry+32, "null")
text(ZA_x+50, ry+27, "Conventional", size=14, fill=SUB, weight="bold")
text(ZA_x+50, ry+49, "Neural stem cells & neurogenesis (existence debated)", size=14, fill=SUB)
darrow(ZA_x+ZA_w/2, ry+70, SLATE)
ry += 88
# new (green)
rect(ZA_x, ry, ZA_w, 64, fill="#eef7f1", stroke=GREEN, rx=12)
verdict(ZA_x+28, ry+32, "ok")
text(ZA_x+50, ry+27, "This study  (neurogenesis-agnostic)", size=14, fill=GREEN, weight="bold")
text(ZA_x+50, ry+49, "Ask what the support cells secrete", size=14.5, fill=INK, weight="bold")
ry += 86
# data card
rect(ZA_x, ry, ZA_w, 168, fill="#f7f8fb", stroke=LINE, rx=12)
text(ZA_x+24, ry+34, "Human hippocampus single-nucleus multiome", size=16, fill=SLATE, weight="bold")
line(ZA_x+24, ry+44, ZA_x+ZA_w-24, ry+44, stroke=LINE)
lines(ZA_x+24, ry+72, [
    "GSE268609  ·  snRNA + ATAC (same nucleus)",
    "8 young (21–38 y)  vs  9 aged (60–93 y)",
    "neurotypical (AD separated)  ·  donor-level stats",
    "read through a 2,224-gene secretome lens",
], size=14.5, fill=INK, lh=30)

# ======================= ZONE B : the finding (hero) ======================
eyebrow(ZB_x, top+6, "WHAT WE FOUND", NAVY)
text(ZB_x, top+44, "The niche's secreted output is reprogrammed", size=21, fill=INK, weight="bold")
text(ZB_x+ZB_w, top+44, "60 significant hits", size=14, fill=SUB, anchor="end")

arms = [
    (C_AST, T_AST, "Astrocyte",              "Extracellular matrix (ECM)", True,  "TNC · IGFBP5 · CCN2"),
    (C_MIC, T_MIC, "Microglia",              "Inflammatory cytokines",     True,  "IL15 · SERPINE1 · APOE"),
    (C_END, T_END, "Endothelium / vascular", "Vascular / angiogenic factors", False, "KDR(VEGFR2) · PLAT"),
]
ay = top + 66
rh = 86
for col, tint, name, func, up, genes in arms:
    rect(ZB_x, ay, ZB_w, rh-12, fill="#ffffff", stroke=col, sw=1.5, rx=12)
    add(f'<rect x="{ZB_x:.1f}" y="{ay:.1f}" width="6" height="{rh-12:.1f}" rx="3" fill="{col}"/>')
    cell_icon(ZB_x+52, ay+(rh-12)/2, col, tint)
    text(ZB_x+96, ay+30, name, size=17, fill=col, weight="bold")
    text(ZB_x+96, ay+54, func, size=15, fill=INK)
    updown(ZB_x+ZB_w-220, ay+(rh-12)/2+8, UP if up else DN, up=up, size=30)
    text(ZB_x+ZB_w-22, ay+(rh-12)/2+6, genes, size=13.5, fill=SUB, anchor="end")
    ay += rh
# converge -> SASP banner
darrow(ZB_x+ZB_w/2, ay+2, NAVY); ay += 22
rect(ZB_x, ay, ZB_w, 60, fill=NAVY, stroke=NAVY, rx=12)
text(ZB_x+ZB_w/2, ay+27, "= a senescence-associated secretory phenotype (SASP)-like signature", size=15,
     fill="#ffffff", weight="bold", anchor="middle")
text(ZB_x+ZB_w/2, ay+49, "Broad inflammation + astrocyte-specific matrix (modest effect size)",
     size=13, fill=NAVYSUB, anchor="middle")

# ======================= ZONE C : evidence status =========================
eyebrow(ZC_x, top+6, "EVIDENCE STATUS", TEAL)
cy = top + 34
cards = [
    ("ok",   GREEN, "External validation (protein level)",
     ["senescence proteome (SASP Atlas)", "brain proteome → cognitive decline"]),
    ("ok",   TEAL,  "Reproduction in a 2nd cohort",
     ["GSE278576 · 40 donors · age 20–95", "effect ρ = 0.24 (positive in all 4 types)"]),
    ("null", GREY,  "Mechanism unresolved",
     ["neither cis-TF nor autocrine supported", "→ open question for perturbation studies"]),
]
ch = 138
for kind, col, ttl, rows in cards:
    rect(ZC_x, cy, ZC_w, ch-16, fill="#ffffff", stroke=col, sw=1.5, rx=12)
    verdict(ZC_x+28, cy+30, kind)
    text(ZC_x+50, cy+27, ttl, size=16, fill=col, weight="bold")
    line(ZC_x+24, cy+44, ZC_x+ZC_w-24, cy+44, stroke=LINE)
    lines(ZC_x+24, cy+72, rows, size=14, fill=INK, lh=26)
    cy += ch

# ---- bottom takeaway band ------------------------------------------------
by = bot + 12
rect(46, by, W-92, 60, fill=NAVY, stroke=NAVY, rx=14)
text(70, by+38, "Takeaway",  size=17, fill="#fff", weight="bold")
line(70, by+47, 152, by+47, stroke=TEAL, sw=3)
text(174, by+38,
     "A rigorously bounded, externally validated aging signature reproduced in a second cohort. Mechanism is an explicit open question; differs from AD broad chromatin dysregulation (no broad amplification in the AD group).",
     size=15, fill="#eaf0fb")

text(48, H-16,
     "Source: scripts/07_figures/09b_graphical_abstract_en.py  ·  numbers per manuscript/manuscript.md  ·  detailed: figures/overview/project_overview.*",
     size=11.5, fill=SUB)

# ================================================================= EMIT
doc = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
       f'viewBox="0 0 {W} {H}">\n' + "\n".join(svg) + "\n</svg>\n")
svg_path = os.path.join(OUTDIR, "graphical_abstract_en.svg")
with open(svg_path, "w", encoding="utf-8") as f:
    f.write(doc)
print("wrote", svg_path, f"({W}x{H})")
try:
    import cairosvg
    cairosvg.svg2png(bytestring=doc.encode("utf-8"),
                     write_to=os.path.join(OUTDIR, "graphical_abstract_en.png"),
                     output_width=2*W, output_height=2*H, background_color="white")
    cairosvg.svg2pdf(bytestring=doc.encode("utf-8"),
                     write_to=os.path.join(OUTDIR, "graphical_abstract_en.pdf"))
    print("wrote graphical_abstract_en.png / .pdf")
except Exception as e:
    print("cairosvg skipped:", e)
