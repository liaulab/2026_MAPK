"""Figure types this manuscript needs that be_scan does not provide.

Stacked bars, mirrored bars, broken-axis bars, distance histograms, pie charts
and the gene-by-condition heatmap. Grouped boxplots come from be_scan's
``boxplot_figure``; the split lollipop lives in ``lollipop.py``.

These follow be_scan.figure_plot's conventions exactly: data as positional
arguments, appearance in frozen option dataclasses, ``load_frame`` on any
frame argument, ``apply_axis_limits`` for limits and ticks, ``arial_font6``
for every label, and ``finish_figure`` to honour ``output.save`` and
``output.show``. Each returns ``(fig, ax)`` like its be_scan counterparts, so
the two are interchangeable at a call site.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch

from be_scan.figure_plot import (
    apply_axis_limits, arial_font6, finish_figure, load_frame,
)

from .figure_options import (
    AXIS, HEATMAPSTYLE, HISTOGRAMSTYLE, LEGEND, OUTPUT, PIESTYLE, BARSTYLE,
    AxisLabelOpts, BarStyleOpts, HeatmapStyleOpts, HistogramStyleOpts,
    LegendOpts, OutputOpts, PieStyleOpts,
)

# ---------------------------------------------------------------------------
# SHARED HELPERS
# ---------------------------------------------------------------------------


def _style_axes(ax, linewidth, hide=('top', 'right')):
    """Thin spines with the top and right removed, as in be_scan's figures."""
    for name in hide:
        ax.spines[name].set_visible(False)
    for spine in ax.spines.values():
        spine.set_linewidth(linewidth)
    ax.tick_params(axis='both', which='both', length=2, width=linewidth)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(arial_font6)


def _label_axes(ax, axis: AxisLabelOpts):
    if axis.xlabel is not None:
        ax.set_xlabel(axis.xlabel, fontproperties=arial_font6)
    if axis.ylabel is not None:
        ax.set_ylabel(axis.ylabel, fontproperties=arial_font6)
    if axis.title:
        ax.set_title(axis.title, fontproperties=arial_font6)


def _legend_path(legend: LegendOpts, output: OutputOpts):
    """Destination for a separately-written legend, from ``legend.path``.

    be_scan treats LegendOpts.path as a stem and appends the output format,
    the same way OutputOpts.path works; this follows that.
    """
    path = Path(f"{legend.path}.{output.out_type}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def save_legend(handles, labels, legend: LegendOpts, output: OutputOpts):
    """Write the legend to its own file, cropped to the legend box.

    be_scan writes legends separately so they can be placed by hand during
    figure assembly; this does the same, keyed off ``legend.path``.
    """
    if not (output.save and legend.path):
        return None
    fig, ax = plt.subplots(figsize=(1, 1))
    ax.axis('off')
    drawn = ax.legend(handles, labels, loc='center', prop=arial_font6,
                      frameon=legend.frameon, ncol=legend.ncol)
    if legend.title:
        drawn.set_title(legend.title, prop=arial_font6)
    fig.canvas.draw()
    bbox = drawn.get_window_extent().transformed(fig.dpi_scale_trans.inverted())

    path = _legend_path(legend, output)
    fig.savefig(path, dpi=output.dpi, format=output.out_type,
                transparent=output.transparent, bbox_inches=bbox)
    plt.close(fig)
    return path


def patch_legend(palette, legend: LegendOpts, output: OutputOpts, linewidth=0.5):
    """Legend of filled swatches, one per palette entry."""
    handles = [Patch(facecolor=color, edgecolor='black', linewidth=linewidth, label=label)
               for label, color in palette.items()]
    return save_legend(handles, list(palette), legend, output)


def rearrange(values):
    """Interleave the two halves of a list: [a1,a2,b1,b2] -> [a1,b1,a2,b2].

    Places each structure's ABE and CBE bars side by side after the two
    editors have been accumulated as separate series.
    """
    mid = len(values) // 2
    return [x for pair in zip(values[:mid], values[mid:]) for x in pair]


def _bar_positions(n, group_gap):
    """Bar x positions with a gap inserted after every pair."""
    xpos = np.arange(n, dtype=float)
    return xpos + (xpos // 2) * group_gap


# ---------------------------------------------------------------------------
# BAR CHARTS
# ---------------------------------------------------------------------------


def stacked_bar_chart(
    xlabels, data_list, colors, labels,

    axis: AxisLabelOpts = AXIS,
    style: BarStyleOpts = BARSTYLE,
    legend: LegendOpts = LEGEND,
    output: OutputOpts = OUTPUT,
    ):
    """Stacked bars, one series per entry in ``data_list``."""
    assert len(colors) == len(data_list), "one colour per series"

    xpos = _bar_positions(len(xlabels), style.group_gap)
    fig, ax = plt.subplots(1, 1, figsize=output.figsize, **output.subplots_kws)

    running = np.zeros(len(data_list[0]))
    for data, color, label in zip(data_list, colors, labels):
        ax.bar(xpos, data, bottom=running, color=color, label=label,
               edgecolor=style.edgecolor, linewidth=style.linewidth)
        running = running + np.asarray(data, dtype=float)

    _style_axes(ax, style.linewidth)
    _label_axes(ax, axis)

    ax.set_xticks(xpos)
    tick_labels = ([l.split(' ')[0] for l in xlabels] if style.short_labels
                   else list(xlabels))
    ax.set_xticklabels(tick_labels, rotation=style.xlabel_rotation, ha='center')
    apply_axis_limits(ax, axis)
    if style.ylog:
        ax.set_yscale('symlog')
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(arial_font6)

    handles, handle_labels = ax.get_legend_handles_labels()
    save_legend(handles, handle_labels, legend, output)
    finish_figure(fig, output)
    return fig, ax


def broken_stacked_bar_chart(
    xlabels, data_list, colors, labels,

    axis: AxisLabelOpts = AXIS,
    style: BarStyleOpts = BARSTYLE,
    legend: LegendOpts = LEGEND,
    output: OutputOpts = OUTPUT,
    ):
    """Stacked bars on a y-axis broken at ``style.ybreak``.

    One tall category would otherwise flatten every other bar. The same bars
    are drawn on both panels and each clips to its own range. Falls back to
    ``stacked_bar_chart`` when no break is configured.
    """
    if style.ybreak is None:
        return stacked_bar_chart(xlabels, data_list, colors, labels,
                                 axis=axis, style=style, legend=legend, output=output)

    xpos = _bar_positions(len(xlabels), style.group_gap)
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=output.figsize, sharex=True,
        gridspec_kw={'height_ratios': [style.ybreak_ratio[1], style.ybreak_ratio[0]],
                     'hspace': 0.08},
    )
    axes = [ax_bot, ax_top]

    running = np.zeros(len(data_list[0]))
    for data, color, label in zip(data_list, colors, labels):
        for ax in axes:
            ax.bar(xpos, data, bottom=running, color=color, label=label,
                   edgecolor=style.edgecolor, linewidth=style.linewidth)
        running = running + np.asarray(data, dtype=float)

    break_lo = style.ybreak
    break_hi = style.ybreak + style.ybreak_gap
    yticks = list(axis.yticks or [])
    ax_bot.set_ylim(axis.ylim[0] if axis.ylim else 0, break_lo)
    ax_top.set_ylim(break_hi, axis.ylim[1] if axis.ylim else running.max() * 1.05)
    ax_bot.set_yticks([t for t in yticks if t <= break_lo])
    ax_top.set_yticks([t for t in yticks if t >= break_hi])

    ax_top.spines['bottom'].set_visible(False)
    ax_bot.spines['top'].set_visible(False)
    ax_top.tick_params(bottom=False)

    # Diagonal marks across the break.
    d = 0.015
    kws = dict(transform=ax_top.transAxes, color='black',
               linewidth=style.linewidth, clip_on=False)
    ax_top.plot((-d, +d), (-d * 2, +d * 2), **kws)
    ax_top.plot((1 - d, 1 + d), (-d * 2, +d * 2), **kws)
    kws.update(transform=ax_bot.transAxes)
    ax_bot.plot((-d, +d), (1 - d * 2, 1 + d * 2), **kws)
    ax_bot.plot((1 - d, 1 + d), (1 - d * 2, 1 + d * 2), **kws)

    for ax in axes:
        _style_axes(ax, style.linewidth)

    if axis.xlabel is not None:
        ax_bot.set_xlabel(axis.xlabel, fontproperties=arial_font6)
    if axis.title:
        ax_top.set_title(axis.title, fontproperties=arial_font6)
    if axis.ylabel is not None:
        fig.text(0.02, 0.5, axis.ylabel, va='center', rotation='vertical',
                 fontproperties=arial_font6)

    ax_bot.set_xticks(xpos)
    ax_bot.set_xticklabels([l.split(' ')[0] for l in xlabels],
                           rotation=style.xlabel_rotation, ha='center')
    for label in ax_bot.get_xticklabels():
        label.set_fontproperties(arial_font6)

    # Each series was drawn on both panels; keep one handle per label.
    handles, handle_labels = ax_bot.get_legend_handles_labels()
    seen = dict(zip(handle_labels, handles))
    save_legend(list(seen.values()), list(seen), legend, output)

    finish_figure(fig, output)
    return fig, axes


def posneg_stacked_bar_chart(
    xlabels, pos_data_list, neg_data_list, colors, labels,
    ylim_top=None, yticks_top=None, ylim_bot=None, yticks_bot=None,

    axis: AxisLabelOpts = AXIS,
    style: BarStyleOpts = BARSTYLE,
    legend: LegendOpts = LEGEND,
    output: OutputOpts = OUTPUT,
    ):
    """Mirrored stacked bars: enriched guides above, depleted below.

    ``colors`` and ``labels`` cover the positive series first, then the
    negative ones. Negative values are supplied positive and drawn downward.
    """
    assert len(colors) == len(pos_data_list) + len(neg_data_list)

    n_pos = len(pos_data_list)
    fig, ax = plt.subplots(2, 1, figsize=output.figsize, **output.subplots_kws)

    running = np.zeros(len(pos_data_list[0]))
    for data, color, label in zip(pos_data_list, colors[:n_pos], labels[:n_pos]):
        ax[0].bar(xlabels, data, bottom=running, color=color, label=label,
                  edgecolor=style.edgecolor, linewidth=style.linewidth)
        running = running + np.asarray(data, dtype=float)

    running = np.zeros(len(neg_data_list[0]))
    for data, color, label in zip(neg_data_list, colors[n_pos:], labels[n_pos:]):
        flipped = -np.asarray(data, dtype=float)
        ax[1].bar(xlabels, flipped, bottom=running, color=color, label=label,
                  edgecolor=style.edgecolor, linewidth=style.linewidth)
        running = running + flipped

    _style_axes(ax[0], style.linewidth, hide=('top', 'right'))
    _style_axes(ax[1], style.linewidth, hide=('bottom', 'right'))

    ax[0].set_xlabel(None)
    if axis.xlabel is not None:
        ax[1].set_xlabel(axis.xlabel, fontproperties=arial_font6)
    if axis.ylabel is not None:
        for panel in ax:
            panel.set_ylabel(axis.ylabel, fontproperties=arial_font6)
    if axis.title:
        ax[0].set_title(axis.title, fontproperties=arial_font6)

    ax[0].set_xticks([])
    ax[1].set_xticks(range(len(xlabels)))
    ax[1].set_xticklabels(xlabels, rotation=style.xlabel_rotation, ha='center')
    if yticks_top is not None:
        ax[0].set_yticks(list(yticks_top))
    if ylim_top is not None:
        ax[0].set_ylim(ylim_top)
    if yticks_bot is not None:
        ax[1].set_yticks(list(yticks_bot))
    if ylim_bot is not None:
        ax[1].set_ylim(ylim_bot)
    for panel in ax:
        for label in panel.get_xticklabels() + panel.get_yticklabels():
            label.set_fontproperties(arial_font6)

    handles = ax[0].get_legend_handles_labels()[0] + ax[1].get_legend_handles_labels()[0]
    handle_labels = (ax[0].get_legend_handles_labels()[1]
                     + ax[1].get_legend_handles_labels()[1])
    save_legend(handles, handle_labels, legend, output)

    finish_figure(fig, output)
    return fig, ax


# ---------------------------------------------------------------------------
# DISTRIBUTIONS
# ---------------------------------------------------------------------------


def distance_histogram_figure(
    df, x_col, hue_col, palette,

    axis: AxisLabelOpts = AXIS,
    style: HistogramStyleOpts = HISTOGRAMSTYLE,
    legend: LegendOpts = LEGEND,
    output: OutputOpts = OUTPUT,
    ):
    """Stacked histogram of hit-to-ligand distances, split by editor."""
    df_data = load_frame(df)
    fig, ax = plt.subplots(1, 1, figsize=output.figsize, **output.subplots_kws)

    hist_kws = dict(data=df_data, x=x_col, hue=hue_col, palette=palette,
                    multiple=style.multiple, ax=ax, **style.hist_kws)
    try:
        sns.histplot(binwidth=style.binwidth, **hist_kws)
    except ValueError:
        # An empty or single-valued set cannot infer a bin width.
        sns.histplot(bins=style.bins, **hist_kws)

    _style_axes(ax, style.linewidth)
    _label_axes(ax, axis)
    apply_axis_limits(ax, axis)
    if ax.legend_ is not None:
        ax.legend_.remove()

    patch_legend(palette, legend, output, linewidth=style.linewidth)
    finish_figure(fig, output)
    return fig, ax


def pie_chart_figure(
    values, labels, colors,

    axis: AxisLabelOpts = AXIS,
    style: PieStyleOpts = PIESTYLE,
    legend: LegendOpts = LEGEND,
    output: OutputOpts = OUTPUT,
    ):
    """Composition pie. Returns ``(None, None)`` when every wedge is zero."""
    if sum(values) == 0:
        return None, None

    fig, ax = plt.subplots(1, 1, figsize=output.figsize, **output.subplots_kws)
    ax.pie(values, labels=labels, colors=colors,
           startangle=style.startangle, counterclock=style.counterclock,
           wedgeprops={'edgecolor': style.edgecolor,
                       'linewidth': style.linewidth, **style.wedge_kws})
    if axis.title:
        ax.set_title(axis.title, fontproperties=arial_font6)
    for text in ax.texts:
        text.set_fontproperties(arial_font6)

    patch_legend(dict(zip(labels, colors)), legend, output, linewidth=style.linewidth)
    finish_figure(fig, output)
    return fig, ax


# ---------------------------------------------------------------------------
# HEATMAP
# ---------------------------------------------------------------------------


def gene_condition_heatmap(
    df, conditions, genes, gene_col='Gene',

    axis: AxisLabelOpts = AXIS,
    style: HeatmapStyleOpts = HEATMAPSTYLE,
    legend: LegendOpts = LEGEND,
    output: OutputOpts = OUTPUT,
    ):
    """Mean score per gene and condition, on a scale symmetric about zero.

    Returns ``(fig, ax, pivot)``; the pivot is the table that was drawn, so
    the numbers behind the figure can be written out alongside it.
    """
    df_data = load_frame(df)
    pivot = (df_data[df_data[gene_col].isin(genes)]
             .groupby(gene_col)[conditions].mean()
             .reindex(index=genes).T
             .reindex(index=conditions, columns=genes))
    values = pivot.values

    abs_max = np.nanmax(np.abs(values))
    norm = mcolors.TwoSlopeNorm(vmin=-abs_max, vcenter=style.center, vmax=abs_max)
    cmap = plt.get_cmap(style.cmap)

    fig, ax = plt.subplots(figsize=output.figsize, **output.subplots_kws)
    ax.imshow(values, cmap=cmap, norm=norm, aspect=style.aspect)

    _label_axes(ax, axis)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=style.xlabel_rotation, ha='right')
    ax.set_yticks(range(len(conditions)))
    ax.set_yticklabels(conditions)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(arial_font6)

    # Colorbar as its own file, matching how legends are handled elsewhere.
    if output.save and legend.path:
        fig_cb, ax_cb = plt.subplots(figsize=style.colorbar_figsize)
        bar = fig_cb.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                              cax=ax_cb, orientation='vertical')
        for label in bar.ax.get_yticklabels():
            label.set_fontproperties(arial_font6)
        path = _legend_path(legend, output)
        fig_cb.savefig(path, dpi=output.dpi, format=output.out_type,
                       transparent=output.transparent, bbox_inches='tight')
        plt.close(fig_cb)

    finish_figure(fig, output)
    return fig, ax, pivot
