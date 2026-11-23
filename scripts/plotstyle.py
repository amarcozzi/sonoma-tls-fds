"""Shared matplotlib style for the sonoma-tls-fds figures.

Journal conventions: no in-figure titles (the caption carries that), lowercase
bold panel letters, 8 pt Arial, black axes with outward ticks and no grid,
figure widths fixed to the standard single (89 mm) and double (183 mm) column,
vector PDF alongside a 300 dpi PNG.

Colour carries exactly one thing: the canopy group (c4 / c8), using the
Okabe-Ito blue/vermillion pair, which stays distinguishable under all common
colour-vision deficiencies and in greyscale. Individual simulations are
identified by direct labels, never by hue. Annotation (ROI, reference lines)
is drawn in greys so it never competes with the data.
"""

import matplotlib as mpl

# figure widths (inches)
SINGLE = 3.5     # 89 mm
DOUBLE = 7.2     # 183 mm

# Okabe-Ito
GROUP_COLORS = {"c4": "#0072B2", "c8": "#D55E00"}
NEUTRAL = "#000000"

INK = "#000000"
INK_2 = "#4d4d4d"     # secondary text, reference lines
MUTED = "#8c8c8c"
REGION = "#ebebeb"    # ROI shading
UNBURNED = "#d9d9d9"  # never-burned cells on the arrival map


def group_color(group, default=NEUTRAL):
    """Colour for a canopy group label ('c4', 'c8', ...)."""
    return GROUP_COLORS.get(group, default)


def apply():
    """Install the house style. Call once, after matplotlib.use('Agg')."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.edgecolor": INK,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "axes.axisbelow": True,
        "axes.labelpad": 3,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.pad": 2,
        "ytick.major.pad": 2,
        "lines.linewidth": 1.0,
        "lines.solid_capstyle": "butt",
        "patch.linewidth": 0.6,
        "legend.frameon": False,
        "legend.handlelength": 1.8,
        "legend.borderaxespad": 0.3,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,     # editable text in Illustrator/Inkscape
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def panel_label(ax, letter, x=-0.02, y=1.0):
    """Bold lowercase panel letter just outside the top-left of the axes."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", ha="right", va="bottom")


def save(fig, stem):
    """Write <stem>.pdf (vector) and <stem>.png (300 dpi)."""
    fig.savefig(f"{stem}.pdf")
    fig.savefig(f"{stem}.png")
