#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Project overview / graphical abstract (one-page, Japanese-primary).

Distils the whole study — reframe -> primary signature -> external validation ->
reproduction in a 2nd cohort -> the regulatory layer -> honest dual-null on
mechanism -> context -> bounded conclusion — into a single self-contained figure.

NOTE (2026-08-26, codex audit): card 5 must not claim methylation coordination. The
DNA-methylation layer failed four pre-submission controls and its RNA-directed
coordination reading is WITHDRAWN (manuscript Results / Fig S(M1)); what survives is a
gene-level, matched-null RNA<->ATAC direction concordance in ASTROCYTES ALONE.

All numbers are taken from manuscript/manuscript.md (.ja) (the authoritative
claims). Renders to SVG, then PNG + PDF via cairosvg.

Run:  conda activate bio && python scripts/07_figures/08_project_overview.py
"""
import os
import html

PROJ = os.environ.get("PROJ", "/home/neurofuture/Bioanalysis/01.SGZ_aging_project")
OUTDIR = os.path.join(PROJ, "figures", "overview")
os.makedirs(OUTDIR, exist_ok=True)

# ---------------------------------------------------------------- canvas / theme
W = 1720
FONT = "Noto Sans CJK JP, DejaVu Sans, sans-serif"

# palette
INK      = "#1f2733"   # primary text
SUB      = "#5b6573"   # secondary text
LINE     = "#d6dbe4"   # card border
CARDBG   = "#f7f8fb"   # card fill
PAGEBG   = "#ffffff"
NAVY     = "#202d4a"   # header / conclusion band
NAVYSUB  = "#aeb8cc"

# semantic accents (what is solid vs open)
GREEN    = "#2c8a5b"   # validated / positive
TEAL     = "#1f8a8a"   # replication
MUSTARD  = "#c0941f"   # supporting / partial / weak
GREY     = "#8a94a6"   # null / unresolved (honest negative, deliberately muted)
SLATE    = "#2f3e5c"   # data / structural
PURPLE   = "#5b6ba8"   # context

# cell-type colour code (Okabe-Ito-ish, muted)
C_AST  = "#e0962f"     # astrocyte   (ECM up)
C_MIC  = "#c5453b"     # microglia   (inflammation up)
C_END  = "#3b72c9"     # endothelium (vascular down)
C_OPC  = "#2b8c7a"     # OPC
C_OLI  = "#8e6fb3"     # oligodendrocyte
C_NSC  = "#9aa3b0"     # neurogenic (de-emphasised)

UP, DN = "#c5453b", "#3b72c9"   # up / down arrows

svg = []
def add(s): svg.append(s)
def esc(s): return html.escape(str(s), quote=True)

def rect(x, y, w, h, fill=CARDBG, stroke=LINE, sw=1.2, rx=12):
    add(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
        f'rx="{rx}" ry="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')

def line(x1, y1, x2, y2, stroke=LINE, sw=1.2, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
        f'stroke="{stroke}" stroke-width="{sw}"{d}/>')

def text(x, y, s, size=15, fill=INK, weight="normal", anchor="start", spacing=None):
    ls = f' letter-spacing="{spacing}"' if spacing else ""
    add(f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
        f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}"{ls}>{esc(s)}</text>')

def lines(x, y, rows, size=15, fill=INK, lh=21, weight="normal", anchor="start"):
    for i, r in enumerate(rows):
        text(x, y + i*lh, r, size=size, fill=fill, weight=weight, anchor=anchor)

def chevron(cx, y, color=SLATE, w=26, h=12):
    add(f'<path d="M{cx-w/2:.1f} {y:.1f} L{cx:.1f} {y+h:.1f} L{cx+w/2:.1f} {y:.1f}" '
        f'fill="none" stroke="{color}" stroke-width="3.2" stroke-linecap="round" '
        f'stroke-linejoin="round"/>')

def numbadge(x, y, n, color):
    add(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="16" fill="{color}"/>')
    text(x, y+6, n, size=18, fill="#fff", weight="bold", anchor="middle")

def chip(x, y, label, color, r=9):
    add(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{color}"/>')
    text(x+r+6, y+5, label, size=14, fill=INK)

# verdict glyph: ok / weak / null  (the "what's solid vs open" language)
def verdict(x, y, kind):
    if kind == "ok":
        add(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="11" fill="{GREEN}"/>')
        add(f'<path d="M{x-5:.1f} {y:.1f} l3.5 4 l6.5 -8" stroke="#fff" '
            f'stroke-width="2.4" fill="none" stroke-linecap="round" stroke-linejoin="round"/>')
    elif kind == "weak":
        add(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="11" fill="{MUSTARD}"/>')
        text(x, y+5, "△", size=14, fill="#fff", weight="bold", anchor="middle")
    else:  # null
        add(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="11" fill="{GREY}"/>')
        add(f'<path d="M{x-4.5:.1f} {y-4.5:.1f} l9 9 M{x+4.5:.1f} {y-4.5:.1f} l-9 9" '
            f'stroke="#fff" stroke-width="2.4" stroke-linecap="round"/>')

def arrow_glyph(x, y, color, up=True, size=15):
    text(x, y, "↑" if up else "↓", size=size, fill=color, weight="bold", anchor="middle")

MARG = 40
IW = W - 2*MARG          # inner content width
X0 = MARG

# running vertical cursor; total height computed at the end
parts_height_cursor = [0.0]
def card_header(x, y, w, n, title, color, subtitle=None):
    """Draw number badge + title + accent underline. Returns y of content start."""
    numbadge(x+26, y+26, n, color)
    text(x+54, y+33, title, size=20, fill=INK, weight="bold")
    line(x+54, y+44, x+54+min(360, w-80), y+44, stroke=color, sw=3)
    if subtitle:
        text(x+w-22, y+30, subtitle, size=13.5, fill=SUB, anchor="end")
    return y + 64

# ================================================================= BUILD
y = 34

# ---- HEADER (navy) -------------------------------------------------------
HH = 120
rect(X0, y, IW, HH, fill=NAVY, stroke=NAVY, rx=16)
text(X0+34, y+46, "加齢はヒト海馬ニッチの「分泌タンパク質 output」を再編する",
     size=31, fill="#ffffff", weight="bold")
text(X0+34, y+78, "外部裏付けを備え第2のヒトコホートで再現する、制御機構が未解決の老化シグネチャ",
     size=18, fill="#e7ecf6")
text(X0+34, y+103,
     "公開単一核データの再解析  —  primary: GSE268609 (Lazarov 2024) ヒト海馬 同一核マルチオーム snRNA + ATAC  ／  一貫してドナー単位の統計 (pseudoreplication-safe, seed 42)",
     size=13.5, fill=NAVYSUB)
y += HH + 18

# ---- REFRAME row (the core idea) ----------------------------------------
RH = 116
half = (IW - 90) / 2
# left: old framing (muted)
rect(X0, y, half, RH, fill="#f0f1f4", stroke="#cfd5df", rx=14)
verdict(X0+30, y+30, "null")
text(X0+48, y+35, "従来の枠組み", size=16, fill=SUB, weight="bold")
lines(X0+24, y+66, [
    "神経幹細胞そのものと「神経新生」に注目",
    "— 成体ヒトで存在するか自体が論争中の現象に依存",
], size=14.5, fill=SUB, lh=23)
# arrow
text(X0+half+45, y+RH/2+10, "→", size=46, fill=SLATE, weight="bold", anchor="middle")
# right: this study (green)
rx0 = X0 + half + 90
rect(rx0, y, half, RH, fill="#eef7f1", stroke=GREEN, rx=14)
verdict(rx0+30, y+30, "ok")
text(rx0+48, y+35, "本研究の問い  (neurogenesis-agnostic)", size=16, fill=GREEN, weight="bold")
lines(rx0+24, y+66, [
    "ニッチ支持細胞が「何を分泌するか」= secretome output を問う",
    "→ 神経新生論争に非依存・機能的・プロテオームと直接照合可能",
], size=14.5, fill=INK, lh=23)
y += RH + 8
chevron(W/2, y, SLATE); y += 24

# ---- CARD 1 : input + lens ----------------------------------------------
H1 = 168
rect(X0, y, IW, H1)
cy = card_header(X0, y, IW, "1", "入力データと解析レンズ", SLATE,
                 subtitle="primary cohort × secretome lens")
lines(X0+28, cy+4, [
    "対比:  若齢 YA 8名 (21–38歳)  vs  高齢 HA 9名 (60–93歳)  ・ neurotypical（AD は分離扱い）",
    "レンズ:  secretome 2,224 遺伝子（Human Protein Atlas predicted ∪ UniProt reviewed-secreted）",
    "解析:  ドナー疑似バルク → pyDESeq2（design ~group）で各 niche 細胞型の差分発現 × secretome",
], size=15, fill=INK, lh=27)
# cell-type colour legend (sender cells)
ly = cy + 92
text(X0+28, ly, "niche sender 細胞（色コード）:", size=14, fill=SUB, weight="bold")
cx = X0 + 250
for lab, col in [("アストロサイト", C_AST), ("ミクログリア", C_MIC),
                 ("オリゴ", C_OLI), ("OPC", C_OPC), ("内皮", C_END),
                 ("(神経新生系列)", C_NSC)]:
    chip(cx, ly-5, lab, col)
    cx += 60 + len(lab)*15
y += H1 + 6
chevron(W/2, y, SLATE); y += 24

# ---- CARD 2 : primary signature -----------------------------------------
H2 = 250
rect(X0, y, IW, H2)
cy = card_header(X0, y, IW, "2", "主軸の結果 — ニッチ secretome 加齢シグネチャ", NAVY,
                 subtitle="60 有意 hit ・ SASP 様 ・ 効果量は中程度")
# three arm sub-boxes
gut = 22
aw = (IW - 56 - 2*gut) / 3
ax = X0 + 28
arms = [
    (C_AST, "アストロサイト", True,  "細胞外マトリックス (ECM) ↑",
     ["TNC ・ COL21A1 ・ COL12A1", "CCN2(CTGF) ・ IGFBP5 ・ APOD", "（Wnt成分 WNT5B は ↓）"], "26 hit"),
    (C_MIC, "ミクログリア", True,  "炎症性サイトカイン ↑",
     ["IL15 (+3.1) ・ SERPINE1/PAI-1", "(+4.1) ・ APOE ・ LGALS9 ・ APOC1", "＋ ミクログリア拡大 (microgliosis)"], "12 hit"),
    (C_END, "内皮 / 血管", False, "血管・血管新生因子 ↓",
     ["KDR (VEGFR2) ・ PLAT", "LAMA4 ・ KITLG", ""], "11 hit"),
]
for i, (col, name, up, head, genes, nhit) in enumerate(arms):
    bx = ax + i*(aw+gut)
    rect(bx, cy, aw, 138, fill="#ffffff", stroke=col, sw=1.6, rx=10)
    add(f'<rect x="{bx:.1f}" y="{cy:.1f}" width="{aw:.1f}" height="6" rx="3" fill="{col}"/>')
    text(bx+18, cy+34, name, size=17, fill=col, weight="bold")
    text(bx+aw-16, cy+33, nhit, size=13, fill=SUB, anchor="end")
    arrow_glyph(bx+aw-16, cy+58, UP if up else DN, up=up, size=20)
    text(bx+18, cy+60, head, size=14.5, fill=INK, weight="bold")
    lines(bx+18, cy+85, genes, size=14, fill=INK, lh=21)
text(X0+28, cy+170,
     "→ niche 細胞型を横断して「広汎な炎症 ＋ アストロサイト特異的マトリックス」 = 細胞老化分泌表現型 (SASP) の形。"
     "  古典的 SASP マーカーはドナー異質で個別非有意（セット水準で実在）。",
     size=14, fill=SUB)
y += H2 + 6
chevron(W/2, y, GREEN); y += 24

# ---- CARD 3 : external validation ---------------------------------------
H3 = 178
rect(X0, y, IW, H3)
cy = card_header(X0, y, IW, "3", "外部検証 — 4 リファレンス（うち 3 つはタンパク質水準）", GREEN,
                 subtitle="核内RNA は分泌タンパクの代理 → 直交参照で照合")
gut = 18
cw = (IW - 56 - 3*gut) / 4
cx0 = X0 + 28
refs = [
    ("ok", "SenMayo", "転写体 SASP 遺伝子セット",
     ["astro 74%↑ (p=4.9e-3)", "micro 60%↑ (p=6.8e-3)"]),
    ("ok", "SASP Atlas", "質量分析 senescence proteome",
     ["astro 76%↑ (p=3.2e-4)", "micro 67%↑ (p=1.7e-3)"]),
    ("ok", "Wingo 脳プロテオーム", "→ 認知トラジェクトリ",
     ["UP アームが認知低下と連動", "(pooled p=0.023, 閾値に頑健)"]),
    ("weak", "CSF 加齢プロテオーム", "SomaScan, 994例",
     ["閾値依存・示唆的", "拡張セットで有意 / 核心は trend"]),
]
for i, (kind, title, sub, rows) in enumerate(refs):
    bx = cx0 + i*(cw+gut)
    edge = GREEN if kind == "ok" else MUSTARD
    rect(bx, cy, cw, 96, fill="#ffffff", stroke=edge, sw=1.4, rx=10)
    verdict(bx+22, cy+24, kind)
    text(bx+40, cy+22, title, size=14.5, fill=INK, weight="bold")
    text(bx+18, cy+46, sub, size=12.5, fill=SUB)
    lines(bx+18, cy+66, rows, size=12.5, fill=INK, lh=18)
text(X0+28, cy+128,
     "主検定は「方向性」（背景頑健）。3/4 がタンパク水準で一致し、特に加齢-UP (SASP) アームが情報的。",
     size=14, fill=SUB)
y += H3 + 6
chevron(W/2, y, TEAL); y += 24

# ---- CARD 4 : independent replication ------------------------------------
H4 = 158
rect(X0, y, IW, H4)
cy = card_header(X0, y, IW, "4", "第2のヒト海馬コホートでの再現", TEAL,
                 subtitle="凍結シグネチャを「再発見せず検定」")
# left text block
lines(X0+28, cy+6, [
    "GSE278576 (Zemke 2024)  40 donor ・ 年齢 20–95 ・ 性別均衡",
    "同一の donor-pseudobulk pipeline（HA≥60 n=20 vs YA<40 n=10）",
    "性別調整・連続年齢モデルに頑健 ／ 2,000× permutation 合格",
    "補助: GSE186538 が内皮アーム 11/11 を方向支持（leverage 限定）",
], size=14.5, fill=INK, lh=24)
# right metric callouts
mx = X0 + IW - 470
for j, (big, small) in enumerate([("ρ = 0.24", "遺伝子別 log2FC 相関\n(4 細胞型すべて正・中程度)"),
                                   ("81% vs 58%", "符号一致 48/60 hit\nvs 背景 (p = 5e-4)")]):
    bx = mx + j*235
    rect(bx, cy, 215, 90, fill="#eaf6f6", stroke=TEAL, sw=1.5, rx=10)
    text(bx+107, cy+40, big, size=26, fill=TEAL, weight="bold", anchor="middle")
    for k, sl in enumerate(small.split("\n")):
        text(bx+107, cy+62+k*17, sl, size=12, fill=SUB, anchor="middle")
y += H4 + 6
chevron(W/2, y, MUSTARD); y += 24

# ---- CARD 5 : regulatory coordination (supporting) ----------------------
H5 = 188
rect(X0, y, IW, H5)
cy = card_header(X0, y, IW, "5", "制御層 — 補助・方向性のみ（座位解像ではない）", MUSTARD,
                 subtitle="multiome を活かした 2 モダリティ（うち1つは撤回）")
gut = 22
cw = (IW - 56 - gut) / 2
bx1 = X0 + 28
bx2 = bx1 + cw + gut
# ATAC
rect(bx1, cy, cw, 96, fill="#ffffff", stroke=MUSTARD, sw=1.4, rx=10)
verdict(bx1+22, cy+24, "weak")
text(bx1+40, cy+22, "クロマチン accessibility (ATAC)", size=15, fill=INK, weight="bold")
lines(bx1+18, cy+48, [
    "RNA と同方向 68% (41/60)。遺伝子単位の matched-null で背景超えは",
    "アストロサイトのみ (p = 5.0×10⁻⁴)。0 / 369,393 peak が FDR を通過",
], size=13, fill=INK, lh=20)
# methylation
rect(bx2, cy, cw, 96, fill="#ffffff", stroke=GREY, sw=1.4, rx=10)
verdict(bx2+22, cy+24, "null")
text(bx2+40, cy+22, "DNA メチル化 (snm3C, GSE299139)", size=15, fill=INK, weight="bold")
lines(bx2+18, cy+48, [
    "40 donor ・ mCG ~ age + sex + coverage（gene-level 検定）",
    "唯一の名目陽性が4統制に不合格 → 協調の読みは撤回",
], size=13, fill=INK, lh=20)
text(X0+28, cy+128,
     "※ メチル化層は方向非特異的な mCG 偏りとして報告し、アクセシビリティ結果の支持には数えない"
     "（サブステート合算・他細胞型RNAでの向き付け・非セクレトーム遺伝子・撹乱帰無の4統制に不合格）。",
     size=13, fill=SUB)
y += H5 + 6
chevron(W/2, y, GREY); y += 24

# ---- CARD 6 : mechanism dual-null ---------------------------------------
H6 = 150
rect(X0, y, IW, H6)
cy = card_header(X0, y, IW, "6", "上流の制御機構 — 正直な「二重 null」", GREY,
                 subtitle="事前指定・locked された探索")
gut = 22
cw = (IW - 56 - gut) / 2
bx1 = X0 + 28; bx2 = bx1 + cw + gut
for bx, ttl, rows in [
    (bx1, "cis-TF プログラム", ["footprinting 0/21 が FDR 非有意",
                              "名目上の 5 hit はむしろ逆向き（加齢で footprint が浅く）"]),
    (bx2, "autocrine ligand–receptor", ["permutation null（ligand 82–91% 通過で非識別的）",
                                         "トップ候補 CCN2 / SPON1 も L–R 検証で脱落"]),
]:
    rect(bx, cy, cw, 78, fill="#f1f3f6", stroke=GREY, sw=1.4, rx=10)
    verdict(bx+22, cy+24, "null")
    text(bx+40, cy+22, ttl, size=15, fill=INK, weight="bold")
    lines(bx+18, cy+46, rows, size=12.8, fill=SUB, lh=18)
text(X0+28, cy+108,
     "→ n = 8/9 では「支持されない」≠ 積極的否定。上流制御因子は perturbation 研究への明示的な open question。",
     size=14, fill=INK, weight="bold")
y += H6 + 6
chevron(W/2, y, PURPLE); y += 24

# ---- CARD 7 : context ----------------------------------------------------
H7 = 168
rect(X0, y, IW, H7)
cy = card_header(X0, y, IW, "7", "文脈 — 空間 ・ クロスコホート ・ 認知 resilience", PURPLE)
gut = 18
cw = (IW - 56 - 2*gut) / 3
cx0 = X0 + 28
ctx = [
    ("空間 (Visium, GSE264692)", ["DG マーカー PROX1 と無相関", "(|ρ| 中央値 0.022/0.021/0.040)",
                                "→ DG への集中は検出されない"]),
    ("クロスコホート (GSE199243)", ["signed 再現性 ~0（WNT5B のみ）", "",
                                 "→ 方向再現であり meta 再現ではない"]),
    ("認知 resilience (GSE325391)", ["resilient granule は control 様 (2 DEG)", "vs 重症 AD (37 DEG)",
                                   "→ 病理下でも転写を保存"]),
]
for i, (ttl, rows) in enumerate(ctx):
    bx = cx0 + i*(cw+gut)
    rect(bx, cy, cw, 96, fill="#ffffff", stroke=PURPLE, sw=1.3, rx=10)
    text(bx+18, cy+26, ttl, size=14.5, fill=PURPLE, weight="bold")
    lines(bx+18, cy+50, rows, size=13, fill=INK, lh=19)
y += H7 + 14

# ---- CONCLUSION band -----------------------------------------------------
HC = 128
rect(X0, y, IW, HC, fill=NAVY, stroke=NAVY, rx=16)
text(X0+34, y+40, "結論", size=20, fill="#fff", weight="bold")
line(X0+34, y+50, X0+90, y+50, stroke=TEAL, sw=3)
lines(X0+110, y+30, [
    "厳密に境界づけられ・外部裏付けを備え・第2のヒトコホートで再現する niche-secretome 加齢シグネチャ（効果量は中程度）。",
    "上流の制御機構は未解決 = perturbation／直交アッセイ研究への明示的 open question。ADの広範なクロマチン異常とは異なる（同一コホートのAD群でも広範増幅は見えない）。",
    "「正直に境界を示すこと」自体が、この positive な主張を信頼できるものにしている。",
], size=15.5, fill="#eaf0fb", lh=30)
y += HC + 26

# ---- footer caption ------------------------------------------------------
text(X0+2, y, "生成: scripts/07_figures/08_project_overview.py  ・  数値は manuscript/manuscript.md（.ja, authoritative）に準拠  ・  raw データは未コミット (.gitignore)",
     size=12, fill=SUB)
y += 22

H = y + 14
# ================================================================= EMIT
svg_doc = (
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H:.0f}" '
    f'viewBox="0 0 {W} {H:.0f}">\n'
    f'<rect x="0" y="0" width="{W}" height="{H:.0f}" fill="{PAGEBG}"/>\n'
    + "\n".join(svg) + "\n</svg>\n"
)

svg_path = os.path.join(OUTDIR, "project_overview.svg")
with open(svg_path, "w", encoding="utf-8") as f:
    f.write(svg_doc)
print("wrote", svg_path, f"({W}x{int(H)})")

try:
    import cairosvg
    png_path = os.path.join(OUTDIR, "project_overview.png")
    pdf_path = os.path.join(OUTDIR, "project_overview.pdf")
    cairosvg.svg2png(bytestring=svg_doc.encode("utf-8"), write_to=png_path,
                     output_width=2*W, output_height=2*int(H), background_color="white")
    cairosvg.svg2pdf(bytestring=svg_doc.encode("utf-8"), write_to=pdf_path)
    print("wrote", png_path)
    print("wrote", pdf_path)
except Exception as e:
    print("cairosvg conversion skipped:", e)
