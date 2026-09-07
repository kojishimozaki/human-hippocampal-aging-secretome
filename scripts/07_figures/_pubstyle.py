#!/usr/bin/env python
"""Shared publication style for the SGZ niche-secretome figure set.

ONE source of truth for fonts / palette / panel letters / saving, so every
Figure (1-5, S*, A1, methyl) is chromatically and typographically identical.
Colours are inherited from the overview infographics
(`figures/overview/`, scripts 08-10) so the whole package reads as one piece.

Usage
-----
    import _pubstyle as PS
    PS.apply()                      # set rcParams once, at top of each fig script
    ... draw into Axes / SubFigure ...
    PS.letter(subfig, "A")          # stamp a panel letter (composites only)
    PS.save(fig, f"{FIG}/panel_X.pdf")

Design notes
------------
* Liberation Sans (Arial-metric-compatible) is the journal-standard look that is
  actually installed here; we list a stack so a bare clone with only DejaVu still
  renders. pdf.fonttype=42 keeps text as *editable* TrueType in the vector PDF
  (so panels can be re-typeset in Illustrator without outlining).
* No randomness lives here; callers own `np.random.seed(42)` (rule 10).
* Panel functions take a *container* (a Figure or a SubFigure) and build their own
  axes inside it, so the identical drawing code serves both the standalone panel
  PDF (legend-referenced filename, unchanged) and the assembled composite.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt

# ----------------------------------------------------------------- palette
INK   = "#1f2733"   # primary text / axes
SUB   = "#5b6573"   # secondary text
GRIDC = "#e6e9ef"   # faint gridlines

# niche cell types (Okabe-Ito-ish; identical to figures/overview)
CT = {
    "Astro": "#e0962f",   # astrocyte    (ECM up)
    "Micro": "#c5453b",   # microglia    (inflammation up)
    "Endo":  "#3b72c9",   # endothelium  (vascular down)
    "OPC":   "#2b8c7a",   # OPC
    "Oligo": "#8e6fb3",   # oligodendrocyte
    "NSC":   "#9aa3b0",   # neurogenic   (de-emphasised)
}
NICHE = ["Astro", "Micro", "Endo", "OPC", "Oligo"]
# spelled-out labels (legends spell them; axes use the short canonical tag)
CT_FULL = {"Astro": "Astrocyte", "Micro": "Microglia", "Endo": "Endothelium",
           "OPC": "OPC", "Oligo": "Oligodendrocyte", "NSC": "Neural stem cell"}

# young- vs aged-donor contrast — fixed across every figure (young = cool, aged = warm)
YA = "#5891c0"
HA = "#cf5b46"

# semantic accents
UP      = "#c5453b"   # up with age / SASP
DN      = "#3b72c9"   # down with age
SECR    = "#3b72c9"   # secretome highlight
SASP    = "#c5453b"   # senescence (SASP) highlight
NEUTRAL = "#c9ced6"   # non-significant / background scatter
GREEN   = "#2c8a5b"   # validated / positive  (verdict)
MUSTARD = "#c0941f"   # supporting / partial  (verdict)
GREY    = "#8a94a6"   # null / unresolved     (verdict)


def _resolve_sans():
    """Return the first sans family actually installed. Naming an absent font
    (e.g. Arial on Linux) makes matplotlib fall back silently -- for text to
    DejaVu, and for mathtext to DejaVu as well, mixing typefaces in one figure."""
    from matplotlib.font_manager import fontManager
    have = {f.name for f in fontManager.ttflist}
    for cand in ("Arial", "Liberation Sans", "Nimbus Sans", "Helvetica", "DejaVu Sans"):
        if cand in have:
            return cand
    return "DejaVu Sans"


_SANS = _resolve_sans()


def apply():
    """Install the shared rcParams. Idempotent; call at the top of each script."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "Nimbus Sans",
                            "Helvetica", "DejaVu Sans"],
        # mathtext must resolve to the SAME family, or ρ/$...$ embed a second
# typeface and the figure carries two fonts (verified with pdffonts).
        "mathtext.fontset": "custom",
        # Resolved at runtime: naming a font that is NOT installed makes mathtext
        # fall back to DejaVu and the figure ends up carrying two typefaces.
        "mathtext.rm": _SANS, "mathtext.it": f"{_SANS}:italic",
        "mathtext.bf": f"{_SANS}:bold",
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "figure.dpi": 150, "savefig.dpi": 300,
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
        "axes.linewidth": 0.8, "axes.edgecolor": INK,
        "axes.labelcolor": INK, "axes.titlecolor": INK,
        "axes.titlesize": 10.5, "axes.titleweight": "bold", "axes.labelsize": 9.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "grid.color": GRIDC, "grid.linewidth": 0.6,
        "text.color": INK, "xtick.color": INK, "ytick.color": INK,
        "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "legend.fontsize": 8, "legend.frameon": False, "legend.handletextpad": 0.5,
        "lines.linewidth": 1.2, "lines.solid_capstyle": "round",
    })


def letter(container, s, x=0.005, y=0.995, size=15):
    """Stamp a bold panel letter at the top-left of a Figure or SubFigure."""
    container.text(x, y, s, fontsize=size, fontweight="bold", color=INK,
                   ha="left", va="top")


def save(fig, path, png=False, dpi=300):
    """Save a vector PDF (and optionally a same-name PNG preview), then close."""
    fig.savefig(path)
    if png:
        fig.savefig(path[:-4] + ".png" if path.endswith(".pdf") else path + ".png",
                    dpi=dpi)
    plt.close(fig)


def fmt_p(p):
    """Journal-style p-value string: never '0.000'; values near 1e-3 read as '0.001'/'0.002'
    (matching the manuscript), scientific only below 1e-4."""
    if p is None or (isinstance(p, float) and (p != p)):
        return "p = n/a"
    if p >= 0.1:
        return f"p = {p:.2f}"
    if p >= 1e-4:
        return "p = " + f"{p:.4f}".rstrip("0").rstrip(".")
    # scientific notation via mathtext (font-independent superscripts; Liberation Sans
    # lacks the unicode superscript-minus/zero glyphs, so we render math instead)
    import math
    e = int(math.floor(math.log10(p)))
    m = p / 10**e
    return rf"$p = {m:.1f}\times10^{{{e}}}$"


def style_ax(ax, grid=False):
    """Light per-axes cleanup beyond rcParams (call where a faint y-grid helps bars)."""
    ax.tick_params(length=3)
    if grid:
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color=GRIDC, linewidth=0.6)


def load_script(path, name="loaded_fig"):
    """Import a sibling figure script by file path (its name has a digit prefix, so the
    normal import machinery can't see it). The loaded module must guard its standalone
    drawing under `if __name__ == '__main__'` so importing only defines its panel fns.
    Used to draw a panel that physically lives in another script into a composite here."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m
