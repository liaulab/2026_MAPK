"""Option dataclasses for the figure types this manuscript adds.

be_scan.figure_plot drives every figure from frozen option dataclasses —
``AxisLabelOpts`` for labels and limits, a per-figure ``*StyleOpts`` for
appearance, ``DomainOpts`` for domain shading, ``LegendOpts`` for the separate
legend file, ``OutputOpts`` for size and destination — with a module-level
default instance in capitals for each. The figures here follow the same
convention, so a call reads the same whichever module it lands in, and a
static figure and its interactive twin take the identical objects.

Anything that is data (the values, the labels, the palette) stays a normal
argument; only appearance lives in these classes, matching be_scan.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

# be_scan's own option classes are re-exported so callers need one import.
from be_scan.figure_plot import (  # noqa: F401
    AXIS, BOXPLOTSTYLE, DOMAIN, LEGEND, LOLLIPOPSTYLE, OUTPUT, SCATTERSTYLE,
    AxisLabelOpts, BoxStyleOpts, DomainOpts, LegendOpts, LollipopStyleOpts,
    OutputOpts, ScatterStyleOpts,
)

import matplotlib as mpl

# Importing be_scan.figure_plot applies the shared figure rcParams. This adds
# the one it does not set. Matplotlib derives the clip-path ids in an SVG from
# a salt drawn at random once per process, so saving the same figure twice
# produces different bytes and no two runs of the notebook can be compared —
# every figure with axes differs, every figure without them (pie charts,
# legends) does not. Fixing the salt makes the SVGs reproducible, which is the
# only way to tell a real change from noise when re-running the analysis.
mpl.rcParams['svg.hashsalt'] = 'MAPK'


@dataclass(frozen=True)
class BarStyleOpts:
    """Stacked and mirrored bar charts."""

    linewidth: float = 0.5
    edgecolor: str = 'black'
    ylog: bool = False
    # Gap inserted after every pair of bars, so each structure's ABE and CBE
    # bars read as one group.
    group_gap: float = 0.5
    # Bar labels carry the structure as well as the gene; only the first word
    # is shown on the axis.
    short_labels: bool = True
    xlabel_rotation: float = 90
    tick_kws: Dict[str, Any] = field(default_factory=lambda:
        {'which': 'both', 'length': 2})
    # Axis break. ybreak is the value the lower panel stops at; ybreak_gap is
    # how much of the range is hidden; ybreak_ratio is (lower, upper) height.
    ybreak: Optional[float] = None
    ybreak_gap: float = 5
    ybreak_ratio: Tuple[float, float] = (2, 1)


@dataclass(frozen=True)
class HistogramStyleOpts:
    """Stacked histogram of hit-to-ligand distances."""

    linewidth: float = 0.5
    binwidth: Optional[float] = 1
    bins: int = 10
    multiple: str = 'stack'
    hist_kws: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class PieStyleOpts:
    """Composition pie charts."""

    linewidth: float = 0.5
    edgecolor: str = 'black'
    startangle: float = 0
    counterclock: bool = True
    wedge_kws: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class HeatmapStyleOpts:
    """Gene-by-condition heatmaps."""

    cmap: str = 'bwr'
    # Symmetric about zero, so depletion and enrichment are comparable.
    center: float = 0.0
    aspect: str = 'auto'
    linewidth: float = 0.5
    colorbar_figsize: Tuple[float, float] = (0.2, 1.0)
    xlabel_rotation: float = 90


@dataclass(frozen=True)
class TrajectoryStyleOpts:
    """GMM trajectory scatter panels and clustered heatmaps."""

    cluster_cmap: str = 'tab20'
    trajectory_cmap: str = 'Reds_r'
    gene_colors: List[str] = field(default_factory=lambda: ['red', 'purple', 'blue'])

    point_size: float = 3
    point_alpha: float = 0.75
    node_size: float = 30
    root_size: float = 50
    edge_kws: Dict[str, Any] = field(default_factory=lambda:
        {'color': 'black', 'lw': 0.5, 'alpha': 0.5})
    ellipse_alpha: float = 0.2
    draw_ellipses: bool = True
    # Points below this posterior confidence are drawn grey and excluded from
    # the coloured layers.
    prob_threshold: float = 0.0

    # Heatmap geometry.
    gray_gene_bars: bool = False
    gene_bar_gray: Tuple[float, float, float, float] = (0.30, 0.30, 0.30, 1.0)
    off_color: Tuple[float, float, float, float] = (0.95, 0.95, 0.95, 1.0)
    col_width: float = 0.12
    heatmap_width: float = 2.0
    row_height: float = 0.015
    hspace: float = 0.01
    spine_linewidth: float = 0.5
    border_linewidth: float = 0.4
    heatmap_cmap: str = 'bwr'

    scatter_figsize: Tuple[float, float] = (5 / 2.54, 5 / 2.54)
    heatmap_figsize: Tuple[float, float] = (10 / 2.54, 10 / 2.54)
    small_heatmap_figsize: Tuple[float, float] = (8 / 2.54, 4 / 2.54)


@dataclass(frozen=True)
class BracketStyleOpts:
    """Significance bracket drawn over a pair of categorical boxes.

    Heights are fractions of the y-axis range rather than data units, so one
    setting works across panels whose Z-score ranges differ.
    """

    # Height of the first bracket, and the gap to each bracket above it.
    # Clears the whisker caps on the control panels without running into the
    # title, at both the (-6, 4) and (-8, 4) ranges those panels use.
    start: float = 0.88
    step: float = 0.08
    # Length of the downward tick at each end of the bar.
    tick: float = 0.03
    # Gap between the bar and its label.
    text_pad: float = 0.01
    linewidth: float = 0.5
    color: str = 'black'


@dataclass(frozen=True)
class SplitLollipopStyleOpts(LollipopStyleOpts):
    """Lollipop with a pocket-distance panel above the stems.

    Extends be_scan's LollipopStyleOpts rather than replacing it, so the
    marker, stem and base settings carry over unchanged.
    """

    distance_line_kws: Dict[str, Any] = field(default_factory=lambda:
        {'color': 'black', 'linewidth': 0.5})
    distance_label: str = 'Distance'
    distance_yticks: Tuple[float, ...] = (0, 50)
    # Distance panel, hit panel, spacer, count panel.
    height_ratios: Tuple[float, ...] = (0.49, 1, 0.01, 1)
    editor_col: str = 'Editor'
    editor_markers: Dict[str, str] = field(default_factory=lambda:
        {'ABE': 'o', 'CBE': 'D'})
    # Lower-panel axis break, for genes with one dominant hotspot.
    ybreak: Optional[Tuple[float, float]] = None
    ybreak_ratio: Tuple[float, float] = (0.4, 0.6)


BARSTYLE = BarStyleOpts()
HISTOGRAMSTYLE = HistogramStyleOpts()
PIESTYLE = PieStyleOpts()
HEATMAPSTYLE = HeatmapStyleOpts()
TRAJECTORYSTYLE = TrajectoryStyleOpts()
SPLITLOLLIPOPSTYLE = SplitLollipopStyleOpts()
BRACKET = BracketStyleOpts()

__all__ = [
    'AXIS', 'BOXPLOTSTYLE', 'DOMAIN', 'LEGEND', 'LOLLIPOPSTYLE', 'OUTPUT',
    'SCATTERSTYLE', 'AxisLabelOpts', 'BoxStyleOpts', 'DomainOpts', 'LegendOpts',
    'LollipopStyleOpts', 'OutputOpts', 'ScatterStyleOpts',
    'BarStyleOpts', 'HistogramStyleOpts', 'PieStyleOpts', 'HeatmapStyleOpts',
    'TrajectoryStyleOpts', 'SplitLollipopStyleOpts', 'BracketStyleOpts',
    'BARSTYLE', 'HISTOGRAMSTYLE', 'PIESTYLE', 'HEATMAPSTYLE',
    'TRAJECTORYSTYLE', 'SPLITLOLLIPOPSTYLE', 'BRACKET',
]
