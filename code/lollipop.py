"""Split lollipop with a pocket-distance panel.

be_scan ships ``lollipop_split_figure``, but it has neither the distance panel
above the stems nor a breakable lower axis. Both are needed to put screen
hits, patient mutation frequency and distance-to-pocket on one position axis,
so that variant lives here, following the same call convention as its be_scan
counterpart: data positionally, appearance in ``SplitLollipopStyleOpts``,
domains in ``DomainOpts``, destination in ``OutputOpts``.

Layout, top to bottom:
    1. distance from the inhibitor to each residue
    2. screen hits, stems upward, coloured by condition
    3. patient mutation counts, stems downward, optionally axis-broken
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import gridspec

from be_scan.figure_plot import arial_font6, finish_figure, load_frame, resolve_domains

from . import config as cfg
from .figure_options import (
    AXIS, DOMAIN, OUTPUT, SPLITLOLLIPOPSTYLE,
    AxisLabelOpts, DomainOpts, OutputOpts, SplitLollipopStyleOpts,
)


def lollipop_split_figure(
    # DISTANCE LINE (top panel) #
    df_line, line_xcol, line_ycol,

    # POSITIVE, screen hits (middle panel) #
    df_pos, pos_xcol, pos_ycol, pos_threshold,
    pos_category_col, pos_categories, pos_colors,
    ylim_top, yticks_top,

    # NEGATIVE, patient counts (bottom panel) #
    df_neg, neg_xcol, neg_ycol, neg_threshold,
    neg_category_col, neg_colors,
    ylim_bot, yticks_bot,

    axis: AxisLabelOpts = AXIS,
    domain: DomainOpts = DOMAIN,
    style: SplitLollipopStyleOpts = SPLITLOLLIPOPSTYLE,
    output: OutputOpts = OUTPUT,
    ):
    """Draw the three-panel lollipop.

    ``pos_categories`` fixes the condition order, so a condition keeps its
    colour across genes even where it contributes no hits. Reading the
    categories off whichever subset happened to be present, as the original
    notebooks did, let colours shift between panels.

    Returns ``(fig, axes, counts)``.
    """
    df_pos = load_frame(df_pos)
    df_neg = load_frame(df_neg)
    df_line = load_frame(df_line)

    fig = plt.figure(figsize=output.figsize, **output.subplots_kws)
    grid = gridspec.GridSpec(4, 1, height_ratios=list(style.height_ratios), figure=fig)
    ax_line = fig.add_subplot(grid[0])
    ax_top = fig.add_subplot(grid[1])

    if style.ybreak is not None:
        grid_bot = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=grid[3],
            height_ratios=[style.ybreak_ratio[1], style.ybreak_ratio[0]], hspace=0.08)
        ax_bot_hi = fig.add_subplot(grid_bot[0])
        ax_bot_lo = fig.add_subplot(grid_bot[1])
        bot_axes = [ax_bot_hi, ax_bot_lo]
    else:
        ax_bot_hi = ax_bot_lo = fig.add_subplot(grid[3])
        bot_axes = [ax_bot_lo]

    pos_hits = df_pos.loc[df_pos[pos_ycol] > pos_threshold]
    neg_hits = df_neg.loc[df_neg[neg_ycol] < neg_threshold]

    # DOMAINS #
    for dom in resolve_domains(domain, df_pos):
        for ax in [ax_line, ax_top] + bot_axes:
            ax.axvspan(dom['start'], dom['end'],
                       facecolor=dom.get('color', domain.default_color),
                       alpha=dom.get('alpha', domain.default_alpha))

    # DISTANCE LINE #
    line_sorted = df_line.sort_values(by=line_xcol)
    ax_line.plot(line_sorted[line_xcol], line_sorted[line_ycol],
                 **style.distance_line_kws)

    # SCREEN HITS #
    for category, color in zip(pos_categories, pos_colors):
        subset = pos_hits[pos_hits[pos_category_col] == category]
        for editor, marker in style.editor_markers.items():
            points = subset[subset[style.editor_col] == editor]
            if points.empty:
                continue
            mark, stem, base = ax_top.stem(
                points[pos_xcol], points[pos_ycol],
                label=f"Pos {category}", basefmt=" ", markerfmt=marker)
            plt.setp(mark, markerfacecolor=color, **style.marker_kws)
            plt.setp(stem, **style.stem_kws)
            plt.setp(base, **style.base_kws)
            mark.set_rasterized(style.marker_rasterized)
            stem.set_rasterized(style.stem_rasterized)

    # PATIENT COUNTS #
    for category, color in zip(neg_hits[neg_category_col].unique(), neg_colors):
        subset = neg_hits[neg_hits[neg_category_col] == category]
        for ax in bot_axes:
            mark, stem, base = ax.stem(
                subset[neg_xcol], subset[neg_ycol],
                label=f"Neg {category}", basefmt=" ")
            plt.setp(mark, markerfacecolor=color, **style.marker_kws)
            plt.setp(stem, **style.stem_kws)
            plt.setp(base, **style.base_kws)
            mark.set_rasterized(style.marker_rasterized)
            stem.set_rasterized(style.stem_rasterized)

    # SHARED AXIS CHROME #
    all_axes = [ax_line, ax_top] + bot_axes
    for ax in all_axes:
        for spine in ax.spines.values():
            spine.set_linewidth(axis.linewidth)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_position(('outward', axis.space))
        if axis.xticks is not None:
            ax.set_xticks(axis.xticks)
        if axis.xlim is not None:
            ax.set_xlim(axis.xlim)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontproperties(arial_font6)

    ax_line.set_xticklabels([])
    if axis.title:
        ax_line.set_title(axis.title, fontproperties=arial_font6)
    ax_line.set_ylim(style.distance_yticks[0], style.distance_yticks[-1])
    ax_line.set_yticks(list(style.distance_yticks))
    ax_line.set_ylabel(style.distance_label, fontproperties=arial_font6)

    ax_top.spines['bottom'].set_visible(False)
    ax_top.set_xticklabels([])
    ax_top.axhline(y=0, **style.line_kws)
    ax_top.set_ylim(*ylim_top)
    ax_top.set_yticks(list(yticks_top))
    if axis.ylabel is not None:
        ax_top.set_ylabel(axis.ylabel, fontproperties=arial_font6)

    if style.ybreak is not None:
        ax_bot_hi.set_ylim(style.ybreak[1], ylim_bot[1])
        ax_bot_hi.set_yticks([t for t in yticks_bot if t >= style.ybreak[1]])
        ax_bot_hi.spines['bottom'].set_visible(False)
        ax_bot_hi.set_xticklabels([])
        ax_bot_hi.axhline(y=0, **style.line_kws)

        ax_bot_lo.set_ylim(ylim_bot[0], style.ybreak[0])
        ax_bot_lo.set_yticks([t for t in yticks_bot if t <= style.ybreak[0]])
        ax_bot_lo.spines['top'].set_visible(False)
        ax_bot_lo.axhline(y=0, **style.line_kws)

        break_kws = dict(transform=fig.transFigure, color='k',
                         linewidth=axis.linewidth, clip_on=False)
        for fragment in (ax_bot_hi, ax_bot_lo):
            box = fragment.get_position()
            y = box.y0 if fragment is ax_bot_hi else box.y1
            dx, dy = 0.008, 0.004
            for x in (box.x0, box.x1):
                fig.add_artist(plt.Line2D([x - dx, x + dx], [y - dy, y + dy], **break_kws))
    else:
        ax_bot_lo.set_ylim(*ylim_bot)
        ax_bot_lo.set_yticks(list(yticks_bot))
        ax_bot_lo.axhline(y=0, **style.line_kws)

    if axis.xlabel is not None:
        ax_bot_lo.set_xlabel(axis.xlabel, fontproperties=arial_font6)

    counts = {'positive_hits': len(pos_hits), 'negative_sites': len(neg_hits)}
    finish_figure(fig, output)
    return fig, all_axes, counts


def load_cbioportal(path=None) -> dict:
    """Every gene's cBioPortal export, keyed by gene symbol.

    One sheet per gene in ``Inputs-cBioPortal.xlsx``, which replaced the
    twenty-two per-gene TSV downloads; the workbook's first sheet is the
    transcript summary and is not a gene.
    """
    sheets = pd.read_excel(path or cfg.CBIOPORTAL_XLSX, sheet_name=None)
    return {gene: sheets[gene] for gene in cfg.CBIOPORTAL_FILES}


def cbioportal_counts(df, mutation_types):
    """Per-residue count of coding mutations reported in cBioPortal.

    Counts are negated so they draw downward beneath the screen hits.
    """
    df = df[df['Mutation Type'].isin(mutation_types)]
    df = df[df['Protein Change'].notna()]
    positions = pd.to_numeric(df['Protein Change'].str.extract(r'(\d+)')[0],
                              errors='coerce')
    df = df.assign(Position=positions)
    df = df[df['Position'].notna()]
    df['Position'] = df['Position'].astype(int)

    counts = df.groupby('Position').size().rename('Count').reset_index()
    counts['pos'] = counts['Position']
    counts['Score'] = -counts['Count']
    counts['Database'] = 'cBioPortal'
    return counts
