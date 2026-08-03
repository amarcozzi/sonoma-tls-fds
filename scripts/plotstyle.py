"""Shared matplotlib styling for the sonoma-tls-fds figures.

One place for the palette so per-simulation figures (postprocess.py) and the
cross-simulation figures (aggregate_metrics.py) read as one set.

Palette rules followed here:
  * Categorical hues are assigned in a fixed order and never cycled. There are
    21 simulations, which is far past the point where hue can carry identity,
    so simulations are NOT colour-coded individually: colour encodes the canopy
    group (c4 / c8, two hues) and identity comes from small multiples, sorted
    axes, and direct labels.
  * Grid/axis furniture is recessive; text never wears a series colour.
  * The two categorical hues pass the CVD checks against a white surface
    (worst adjacent pair dE 24.7 protan / 32.7 tritan, 33.6 normal vision).
"""

import matplotlib as mpl

# categorical slots, fixed order
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]

PRIMARY = SERIES[0]
ACCENT = SERIES[1]

# status / annotation
CRITICAL = "#d03b3b"   # ROI markers, extinction flags
GOOD = "#0ca30c"

# ink + furniture
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#ffffff"

GROUP_COLORS = {"c4": SERIES[0], "c8": SERIES[1]}


def group_color(group, default=MUTED):
    """Colour for a canopy group label ('c4', 'c8', ...)."""
    return GROUP_COLORS.get(group, default)


def apply():
    """Install the house style. Call once, after matplotlib.use('Agg')."""
    mpl.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "savefig.bbox": "tight",
        "figure.dpi": 160,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "axes.labelcolor": INK_2,
        "axes.titlecolor": INK,
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.prop_cycle": mpl.cycler(color=SERIES),
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK_2,
        "ytick.labelcolor": INK_2,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.frameon": False,
        "legend.fontsize": 8,
        "lines.linewidth": 1.6,
        "lines.solid_capstyle": "round",
    })


def roi_marker(ax, positions, axis="y", label_fmt=None):
    """Draw the recessive dashed ROI bounds without stealing the legend."""
    for pos in positions:
        line = ax.axhline if axis == "y" else ax.axvline
        line(pos, color=CRITICAL, ls="--", lw=1.0, alpha=0.8, zorder=0)
        if label_fmt:
            if axis == "y":
                ax.annotate(label_fmt.format(pos), (0.005, pos),
                            xycoords=("axes fraction", "data"), va="bottom",
                            fontsize=7, color=CRITICAL)
            else:
                ax.annotate(label_fmt.format(pos), (pos, 0.99),
                            xycoords=("data", "axes fraction"), va="top",
                            ha="right", rotation=90, fontsize=7,
                            color=CRITICAL)
