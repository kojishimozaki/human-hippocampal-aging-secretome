#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Where this analysis sits in neuroscience (Japanese positioning illustration).

Not a methods/abstract figure — a *map*. It places the study inside the field of
hippocampal-aging neuroscience: (1) the reframe (neurogenesis/stem-cell-centric ->
niche support-cell secretome output), (2) the adjacent neuroscience domains it
bridges (brain senescence/SASP, neuroinflammation, gliovascular, cognitive aging,
fluid biomarkers, epigenetic regulation), (3) its boundaries (distinct from AD;
sidesteps the adult-neurogenesis debate; mechanism is the open frontier).

Shares the cell-type colour code + verdict semantics (ok/weak/open) with the
08/09 overview figures. Renders SVG -> PNG + PDF.

Run:  conda activate bio && python scripts/07_figures/10_neuroscience_positioning.py
"""
import os, math, html

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
OUTDIR = os.path.join(PROJ, "figures", "overview")
os.makedirs(OUTDIR, exist_ok=True)

W, H = 1680, 1130
FONT = "Noto Sans CJK JP, DejaVu Sans, sans-serif"

INK, SUB, LINE = "#1f2733", "#5b6573", "#d6dbe4"
PAGEBG, NAVY, NAVYSUB = "#ffffff", "#202d4a", "#aeb8cc"
GREEN, TEAL, MUSTARD, GREY, SLATE, PURPLE = "#2c8a5b", "#1f8a8a", "#c0941f", "#8a94a6", "#2f3e5c", "#5b6ba8"
C_AST, C_MIC, C_END = "#e0962f", "#c5453b", "#3b72c9"
T_AST, T_MIC, T_END = "#fbeed7", "#f7dcd9", "#dde9fa"
BRAIN_F, BRAIN_S = "#f7dbe6", "#d68aa9"
HIPPO = "#d9962a"

svg = []
def add(s): svg.append(s)
def esc(s): return html.escape(str(s), quote=True)

def rect(x, y, w, h, fill="#ffffff", stroke=LINE, sw=1.2, rx=12, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    add(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" ry="{rx}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>')

def line(x1, y1, x2, y2, stroke=LINE, sw=1.2):
    add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" stroke-width="{sw}"/>')

def text(x, y, s, size=15, fill=INK, weight="normal", anchor="start", spacing=None):
    ls = f' letter-spacing="{spacing}"' if spacing else ""
    add(f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
        f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}"{ls}>{esc(s)}</text>')

def lines(x, y, rows, size=15, fill=INK, lh=22, weight="normal", anchor="start"):
    for i, r in enumerate(rows):
        text(x, y+i*lh, r, size=size, fill=fill, weight=weight, anchor=anchor)

def circle(cx, cy, r, fill, stroke="none", sw=0, op=1.0):
    st = f' stroke="{stroke}" stroke-width="{sw}"' if stroke != "none" else ""
    add(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="{fill}"{st} opacity="{op}"/>')

def verdict(x, y, kind, r=11):
    col = {"ok": GREEN, "weak": MUSTARD, "open": GREY}[kind]
    add(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{col}"/>')
    if kind == "ok":
        add(f'<path d="M{x-5:.1f} {y:.1f} l3.6 4 l6.4 -8" stroke="#fff" stroke-width="2.4" '
            f'fill="none" stroke-linecap="round" stroke-linejoin="round"/>')
    elif kind == "weak":
        text(x, y+5, "△", size=13, fill="#fff", weight="bold", anchor="middle")
    else:
        text(x, y+5, "?", size=14, fill="#fff", weight="bold", anchor="middle")

def cell_icon(cx, cy, color, tint, rx=20, ry=15.5):
    add(f'<ellipse cx="{cx:.1f}" cy="{cy:.1f}" rx="{rx}" ry="{ry}" fill="{tint}" stroke="{color}" stroke-width="2"/>')
    add(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{rx*0.32:.1f}" fill="{color}" opacity="0.85"/>')

def updown(x, y, color, up=True, size=16):
    text(x, y, "↑" if up else "↓", size=size, fill=color, weight="bold", anchor="middle")

# ---- stylised brain with highlighted hippocampus ------------------------
def brain(cx, cy, s=1.0):
    add(f'<g transform="translate({cx},{cy}) scale({s})">')
    # cerebrum blob with a few bumpy gyri on top
    path = ("M -86 14 "
            "C -94 -8 -78 -24 -60 -24 "
            "C -60 -44 -30 -48 -20 -30 "
            "C -12 -50 20 -48 28 -28 "
            "C 40 -48 76 -36 70 -10 "
            "C 90 -6 88 20 68 26 "
            "C 72 46 40 52 26 40 "
            "C 14 54 -16 52 -28 38 "
            "C -46 50 -78 38 -70 14 Z")
    add(f'<path d="{path}" fill="{BRAIN_F}" stroke="{BRAIN_S}" stroke-width="2.4" stroke-linejoin="round"/>')
    # gyri squiggles
    for d in ["M -54 -8 C -44 -18 -30 -2 -18 -12",
              "M -8 -6 C 4 -18 18 -4 30 -14",
              "M -50 12 C -36 4 -24 18 -10 8",
              "M 2 12 C 16 4 30 18 46 8"]:
        add(f'<path d="{d}" fill="none" stroke="{BRAIN_S}" stroke-width="1.6" opacity="0.7"/>')
    # brainstem
    add(f'<path d="M -12 40 C -14 56 -8 64 0 70" fill="none" stroke="{BRAIN_S}" stroke-width="6" stroke-linecap="round"/>')
    # hippocampus (curled, highlighted) in the temporal region
    add(f'<path d="M -36 22 C -20 36 6 36 18 22 C 26 14 22 2 12 6 C 22 16 6 26 -6 22 C -18 18 -26 14 -36 22 Z" '
        f'fill="{HIPPO}" stroke="#a9711a" stroke-width="1.6" opacity="0.95"/>')
    add('</g>')

MARG = 40

# ================================================================= BUILD
add(f'<rect x="0" y="0" width="{W}" height="{H}" fill="{PAGEBG}"/>')

# ---- title ---------------------------------------------------------------
text(MARG+4, 56, "この解析の神経科学における立ち位置", size=32, fill=INK, weight="bold")
text(MARG+6, 88,
     "海馬の加齢研究を「神経新生」中心から「ニッチ支持細胞の分泌 output」へ転換し — 老化・炎症・血管・認知を架橋する",
     size=16.5, fill=SUB)
line(MARG+4, 104, W-MARG, 104, stroke=LINE, sw=1.4)

# ---- reframe ribbon (the core pivot) ------------------------------------
ry = 122; rh = 92
half = (W - 2*MARG - 96)/2
# old framing (faded, avoided)
rect(MARG, ry, half, rh, fill="#f0f1f4", stroke="#cfd5df", rx=14, dash="2 5")
add(f'<circle cx="{MARG+30}" cy="{ry+34}" r="13" fill="{GREY}"/>')
add(f'<path d="M{MARG+30-7} {ry+34} h14" stroke="#fff" stroke-width="2.6" stroke-linecap="round"/>')
text(MARG+52, ry+30, "従来の中心軸（本研究は回避）", size=15, fill=SUB, weight="bold")
text(MARG+52, ry+56, "神経幹細胞・神経新生に注目 — 成体ヒトでは存在自体が論争中", size=14.5, fill=SUB)
text(MARG+52, ry+78, "neurogenesis-centric ・ cell-intrinsic", size=12.5, fill=GREY)
# arrow
text(W/2, ry+rh/2+14, "→", size=52, fill=SLATE, weight="bold", anchor="middle")
# new lens
nx = MARG + half + 96
rect(nx, ry, half, rh, fill="#eef7f1", stroke=GREEN, rx=14)
verdict(nx+30, ry+34, "ok", r=13)
text(nx+52, ry+30, "本研究のレンズ（neurogenesis-agnostic）", size=15, fill=GREEN, weight="bold")
text(nx+52, ry+56, "支持細胞が「何を分泌するか」= ニッチ環境の output を読む", size=14.5, fill=INK, weight="bold")
text(nx+52, ry+78, "機能的 ・ プロテオーム/CSF と直接比較可能 ・ 論争を迂回", size=12.5, fill=SUB)

# ---- field backdrop ------------------------------------------------------
BX, BY, BW, BH = MARG, 244, W-2*MARG, 706
rect(BX, BY, BW, BH, fill="#f9fbfd", stroke="#e2e7ef", sw=1.4, rx=20, dash="3 6")
add(f'<rect x="{BX+24}" y="{BY-14}" width="360" height="28" rx="14" fill="#eef2f8" stroke="#e2e7ef"/>')
text(BX+44, BY+5, "海馬の加齢神経科学 という土俵での位置", size=14.5, fill=SLATE, weight="bold")

# geometry: hub + 6 radial domain bubbles
HX, HY, HR = 858, 604, 138
bub = [   # (cx, cy, title, role, verdict)
    (858,  366, "細胞老化・SASP", "老年学：SASP 様の分泌表現型", "ok"),
    (1216, 500, "神経炎症",       "ミクログリア活性化・炎症 ↑",   "ok"),
    (1216, 772, "認知老化",       "脳プロテオーム → 認知低下と連動", "ok"),
    (858,  842, "体液バイオマーカー", "CSF 加齢プロテオームと方向一致", "weak"),
    (500,  772, "グリア–血管ユニット", "アストロ ECM ↑ ・ 内皮 血管 ↓", "ok"),
    (500,  500, "エピゲノム制御",   "chromatin/メチル化 = 未解決の機構", "open"),
]
# spokes (behind)
for cx, cy, *_ in bub:
    add(f'<line x1="{HX}" y1="{HY}" x2="{cx}" y2="{cy}" stroke="{TEAL}" stroke-width="2.2" opacity="0.45"/>')

# anatomical inset (brain -> niche) at left, with a guiding arrow to hub
brain(196, 470, 1.18)
text(196, 388, "ヒト海馬", size=15, fill=INK, weight="bold", anchor="middle")
text(196, 556, "GSE268609 ほか", size=12.5, fill=SUB, anchor="middle")
text(196, 574, "公開単一核データ", size=12.5, fill=SUB, anchor="middle")
add(f'<path d="M 252 548 C 400 656 595 696 745 678" fill="none" '
    f'stroke="{SLATE}" stroke-width="2.4" stroke-dasharray="2 6" marker-end="url(#ah)"/>')
text(338, 640, "ニッチを拡大して読む", size=12.5, fill=SLATE)

# domain bubbles
bw, bh = 322, 86
for cx, cy, title, role, vk in bub:
    border = GREY if vk == "open" else TEAL
    dash = "3 5" if vk == "open" else None
    rect(cx-bw/2, cy-bh/2, bw, bh, fill="#ffffff", stroke=border, sw=1.6, rx=14, dash=dash)
    verdict(cx-bw/2+26, cy-bh/2+24, vk)
    text(cx-bw/2+48, cy-bh/2+29, title, size=16.5, fill=INK, weight="bold")
    text(cx-bw/2+22, cy+bh/2-18, role, size=13.5, fill=SUB)

# ---- hub (the study) -----------------------------------------------------
circle(HX, HY, HR+6, "#ffffff", stroke=NAVY, sw=0)
circle(HX, HY, HR, "#f3f6fc", stroke=NAVY, sw=3)
# niche emblem: three sender cells secreting (output)
cell_icon(HX-58, HY-64, C_AST, T_AST); updown(HX-26, HY-60, C_AST, True, 15)
cell_icon(HX+50, HY-66, C_MIC, T_MIC); updown(HX+82, HY-62, C_MIC, True, 15)
cell_icon(HX-2,  HY-22, C_END, T_END); updown(HX+30, HY-18, C_END, False, 15)
for dx, dy, c in [(-30,-30,C_AST),(20,-34,C_MIC),(-8,-2,C_END),(28,-8,C_MIC),
                  (-46,-40,C_AST),(8,-50,C_MIC),(-30,-6,C_END)]:
    circle(HX+dx, HY+dy, 3.0, c, op=0.8)
text(HX, HY+36, "海馬ニッチの分泌 output", size=15.5, fill=NAVY, weight="bold", anchor="middle")
text(HX, HY+58, "加齢シグネチャ ＝ 本研究", size=14.5, fill=INK, weight="bold", anchor="middle")
# boundary tag: distinct from AD
add(f'<rect x="{HX-74}" y="{HY+72}" width="148" height="26" rx="13" fill="#fdeceb" stroke="{C_MIC}" stroke-width="1.4"/>')
text(HX, HY+90, "正常加齢に限定（≠ AD）", size=12.5, fill="#a23a32", weight="bold", anchor="middle")

# ---- bottom: positioning coordinates ------------------------------------
cy0 = 978; chh = 92
gap = 22; cwid = (W - 2*MARG - 2*gap)/3
coords = [
    (SLATE, "扱う軸",   ["正常加齢の生物学に位置づけ", "神経変性（AD）とは別の軸として分離"]),
    (TEAL,  "つなぐ階層", ["単一細胞 RNA → タンパク質 → 認知へ接続", "外部裏付け・第2コホートで再現"]),
    (GREY,  "未踏の縁",  ["記述（output シグネチャ）は確立", "上流の制御機構は未解決フロンティア"]),
]
for i, (col, tag, rows) in enumerate(coords):
    x = MARG + i*(cwid+gap)
    rect(x, cy0, cwid, chh, fill="#f7f8fb", stroke=col, sw=1.5, rx=14)
    add(f'<path d="M{x+30} {cy0+30} l11 11 l-11 11 l-11 -11 Z" fill="{col}"/>')  # compass diamond
    text(x+56, cy0+36, tag, size=15.5, fill=col, weight="bold")
    lines(x+26, cy0+60, rows, size=14, fill=INK, lh=22)

text(MARG+2, H-16,
     "生成: scripts/07_figures/10_neuroscience_positioning.py  ・  位置づけは manuscript/manuscript.md（.ja）の Introduction/Discussion に準拠",
     size=11.5, fill=SUB)

# defs (arrowhead)
defs = (f'<defs><marker id="ah" markerWidth="10" markerHeight="10" refX="7" refY="3" '
        f'orient="auto" markerUnits="strokeWidth"><path d="M0 0 L7 3 L0 6 Z" fill="{SLATE}"/></marker></defs>')

doc = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">\n'
       + defs + "\n" + "\n".join(svg) + "\n</svg>\n")
svg_path = os.path.join(OUTDIR, "neuroscience_positioning.svg")
with open(svg_path, "w", encoding="utf-8") as f:
    f.write(doc)
print("wrote", svg_path, f"({W}x{H})")
try:
    import cairosvg
    cairosvg.svg2png(bytestring=doc.encode("utf-8"),
                     write_to=os.path.join(OUTDIR, "neuroscience_positioning.png"),
                     output_width=2*W, output_height=2*H, background_color="white")
    cairosvg.svg2pdf(bytestring=doc.encode("utf-8"),
                     write_to=os.path.join(OUTDIR, "neuroscience_positioning.pdf"))
    print("wrote neuroscience_positioning.png / .pdf")
except Exception as e:
    print("cairosvg skipped:", e)
