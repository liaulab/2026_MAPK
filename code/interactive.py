"""Interactive twins of the figures in ``charts.py``, ``lollipop.py`` and
``trajectory.py``.

be_scan already ships interactive counterparts of its own figures —
``interactive_scatterplot_figure``, ``interactive_correlation_scatterplot_figure``,
``interactive_boxplot_figure``, ``interactive_jitterbox_kdeplot_figure``,
``interactive_lollipop_figure`` — driven by the same option dataclasses as the
static versions and styled from the same matplotlib rcParams. Those are used
directly wherever they fit. This module covers only the figure types this
manuscript adds on top.

Each function here mirrors its static twin's signature exactly, with one extra
argument, ``scale``, exactly as be_scan's own interactive functions do. Swap
``charts.stacked_bar_chart`` for ``interactive.stacked_bar_chart`` at a call
site and nothing else changes.

Styling is inherited, not reimplemented: ``base_layout``, ``axis_style``,
``apply_axis``, ``add_domains``, ``marker_from_kws``, ``to_css`` and
``hover_text`` all come from ``be_scan.figure_plot.figure_interactive_style``,
so restyling the print figures restyles these too.

What deliberately does not carry across is absolute size. Print figures are
about 4 x 2.2 inches with 6 pt text, illegible in a browser, so ``scale``
multiplies pixel dimensions and font sizes together: proportions are preserved
exactly and the result is readable. ``scale=1`` reproduces print geometry.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from be_scan.figure_plot import load_frame, resolve_domains
from be_scan.figure_plot import figure_interactive_style as bstyle

from . import config as cfg
from .figure_options import (
    AXIS, BARSTYLE, BRACKET, DOMAIN, HEATMAPSTYLE, HISTOGRAMSTYLE, OUTPUT,
    PIESTYLE, SPLITLOLLIPOPSTYLE, TRAJECTORYSTYLE,
    AxisLabelOpts, BarStyleOpts, BracketStyleOpts, DomainOpts,
    HeatmapStyleOpts, HistogramStyleOpts, OutputOpts, PieStyleOpts,
    SplitLollipopStyleOpts, TrajectoryStyleOpts,
)

DEFAULT_SCALE = 2.0

# matplotlib colormaps that plotly does not ship. Defining them explicitly
# keeps the interactive heatmaps on exactly the static colour ramp rather than
# a lookalike.
MPL_COLORSCALES = {
    'bwr': [[0.0, 'rgb(0,0,255)'], [0.5, 'rgb(255,255,255)'], [1.0, 'rgb(255,0,0)']],
    'seismic': [[0.0, 'rgb(0,0,76)'], [0.25, 'rgb(0,0,255)'], [0.5, 'rgb(255,255,255)'],
                [0.75, 'rgb(255,0,0)'], [1.0, 'rgb(127,0,0)']],
}


def colorscale(name):
    """Plotly colorscale for a matplotlib colormap name."""
    return MPL_COLORSCALES.get(name, name)

# How plotly.js is delivered. 'cdn' keeps each figure a few tens of kilobytes
# rather than ~3.7 MB, which matters across hundreds of figures; the browser
# caches the library once and reuses it for all of them.
PLOTLYJS = 'cdn'


def configure(include_plotlyjs='cdn'):
    """Set plotly.js delivery for be_scan's interactive figures and these."""
    global PLOTLYJS
    PLOTLYJS = include_plotlyjs
    bstyle.INCLUDE_PLOTLYJS = include_plotlyjs


# ---------------------------------------------------------------------------
# SHARED HELPERS
# ---------------------------------------------------------------------------


def _layout(output: OutputOpts, scale=DEFAULT_SCALE, title=None,
            showlegend=False, n_rows=1):
    """be_scan's base layout, with the title and legend flag applied."""
    layout = bstyle.base_layout(output, scale, n_rows=n_rows)
    layout['showlegend'] = showlegend
    if title:
        layout['title'] = dict(text=title, font=dict(
            family=bstyle.FONT_FAMILY, size=bstyle.BASE_FONT * scale, color='black'))
    return layout


def _axis_kwargs(scale=DEFAULT_SCALE):
    return bstyle.axis_style(scale)


def _apply(fig, axis: AxisLabelOpts, scale, row=None, col=None, transparent=True):
    bstyle.apply_axis(fig, axis, scale, row=row, col=col, transparent=transparent)
    return repair_axis_labels(fig, axis, row=row, col=col)


def repair_axis_labels(fig, axis: AxisLabelOpts, row=None, col=None):
    """Re-apply the x and y titles after be_scan's ``apply_axis``.

    ``apply_axis`` builds its two axis dicts with ``dict(style)``, a shallow
    copy, so ``xaxis['title']`` and ``yaxis['title']`` are the same nested
    object and whichever label is assigned second wins. Any interactive figure
    that sets both ends up with the y label on both axes. This reasserts them
    independently; harmless once be_scan is fixed upstream.
    """
    if axis.xlabel is not None:
        fig.update_xaxes(title_text=axis.xlabel, row=row, col=col)
    if axis.ylabel is not None:
        fig.update_yaxes(title_text=axis.ylabel, row=row, col=col)
    return fig


def add_pvalue_brackets(fig, comparisons, lab_order,
                        bracket: BracketStyleOpts = BRACKET,
                        scale: float = DEFAULT_SCALE):
    """Draw significance brackets over a categorical plotly figure.

    The interactive twin of what ``boxplot.control_boxplot_figure`` draws, and
    takes the identical ``comparisons`` list, so the web panel and the print
    panel carry the same labels. Categories sit at integer x positions, and
    heights are fractions of the y range, read off the axis the caller has
    already set.
    """
    order = list(lab_order)
    bottom, top = fig.layout.yaxis.range or (0, 1)
    span = top - bottom

    for index, (left, right, text) in enumerate(comparisons):
        if left not in order or right not in order:
            continue
        x_left, x_right = order.index(left), order.index(right)
        y = bottom + span * (bracket.start + index * bracket.step)
        tick = span * bracket.tick
        fig.add_trace(go.Scatter(
            x=[x_left, x_left, x_right, x_right],
            y=[y - tick, y, y, y - tick],
            mode='lines', showlegend=False, hoverinfo='skip',
            line=dict(color=bracket.color, width=bracket.linewidth * scale),
        ))
        fig.add_annotation(
            x=(x_left + x_right) / 2, y=y + span * bracket.text_pad,
            text=text, showarrow=False, yanchor='bottom',
            font=dict(size=bstyle.BASE_FONT * scale, color=bracket.color),
        )
    return fig


def finish(fig, output: OutputOpts):
    """Write the HTML when ``output.save`` is set, and return the figure.

    Mirrors be_scan's ``finish_interactive``, but writes with the viewer
    config this site wants (SVG export, no plotly logo, responsive).
    """
    if output.save:
        path = Path(f"{output.path}.html")
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(
            str(path), include_plotlyjs=PLOTLYJS, full_html=True,
            config={'displaylogo': False, 'responsive': True,
                    'toImageButtonOptions': {'format': 'svg'}},
        )
    return fig


def save(fig, out_path, include_plotlyjs=None):
    """Write a figure as a standalone HTML page at an explicit path."""
    out_path = Path(out_path)
    if out_path.suffix != '.html':
        out_path = out_path.with_suffix('.html')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        str(out_path),
        include_plotlyjs=PLOTLYJS if include_plotlyjs is None else include_plotlyjs,
        full_html=True,
        config={'displaylogo': False, 'responsive': True,
                'toImageButtonOptions': {'format': 'svg'}},
    )
    return out_path


class Manifest:
    """Record of every interactive figure written, for the site builder.

    Each entry carries the section it belongs to, the group it is filed under
    on that section's page (a gene or a structure), a title, and the figure's
    path relative to the site root.
    """

    def __init__(self, site_root):
        self.site_root = Path(site_root)
        self.entries = []

    def add(self, path, section, group, title, caption='', height=None):
        self.entries.append({
            'section': section, 'group': group, 'title': title, 'caption': caption,
            'height': int(height) if height else None,
            'path': Path(path).relative_to(self.site_root).as_posix(),
        })
        return path

    def save_figure(self, fig, out_path, section, group, title, caption=''):
        # The figure knows its own pixel height; the site uses it so each
        # frame matches its figure rather than a per-section guess.
        height = getattr(fig.layout, 'height', None)
        return self.add(save(fig, out_path), section, group, title, caption,
                        height=height)

    def write(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.entries, indent=1))
        return path

    def __len__(self):
        return len(self.entries)


def _bar_positions(n, group_gap):
    xpos = np.arange(n, dtype=float)
    return xpos + (xpos // 2) * group_gap


# ---------------------------------------------------------------------------
# BAR CHARTS
# ---------------------------------------------------------------------------


def stacked_bar_chart(
    xlabels, data_list, colors, labels,

    axis: AxisLabelOpts = AXIS,
    style: BarStyleOpts = BARSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = DEFAULT_SCALE,
    ):
    """Interactive twin of ``charts.stacked_bar_chart``."""
    xpos = _bar_positions(len(xlabels), style.group_gap)
    tick_labels = ([l.split(' ')[0] for l in xlabels] if style.short_labels
                   else list(xlabels))

    fig = go.Figure()
    for data, color, label in zip(data_list, colors, labels):
        fig.add_trace(go.Bar(
            x=xpos, y=list(data), name=label, marker_color=color,
            marker_line=dict(color=style.edgecolor, width=style.linewidth * scale),
            customdata=list(xlabels),
            hovertemplate='%{customdata}<br>' + f'{label}: %{{y}}<extra></extra>',
        ))

    fig.update_layout(barmode='stack',
                      **_layout(output, scale, axis.title, showlegend=True))
    kwargs = _axis_kwargs(scale)
    fig.update_xaxes(**{**kwargs, 'tickmode': 'array', 'tickvals': list(xpos),
                        'ticktext': tick_labels, 'tickangle': -style.xlabel_rotation},
                     title_text=axis.xlabel)
    ykw = dict(kwargs)
    if axis.yticks is not None:
        ykw.update(tickmode='array', tickvals=list(axis.yticks))
    if axis.ylim is not None:
        ykw['range'] = list(axis.ylim)
    if style.ylog:
        ykw['type'] = 'log'
        ykw.pop('range', None)
    fig.update_yaxes(**ykw, title_text=axis.ylabel)
    return finish(fig, output)


def broken_stacked_bar_chart(
    xlabels, data_list, colors, labels,

    axis: AxisLabelOpts = AXIS,
    style: BarStyleOpts = BARSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = DEFAULT_SCALE,
    ):
    """Interactive twin of ``charts.broken_stacked_bar_chart``."""
    if style.ybreak is None:
        return stacked_bar_chart(xlabels, data_list, colors, labels,
                                 axis=axis, style=style, output=output, scale=scale)

    xpos = _bar_positions(len(xlabels), style.group_gap)
    break_hi = style.ybreak + style.ybreak_gap
    yticks = list(axis.yticks or [])

    fig = make_subplots(rows=2, cols=1, shared_xaxes=True,
                        row_heights=[style.ybreak_ratio[1], style.ybreak_ratio[0]],
                        vertical_spacing=0.06)
    for data, color, label in zip(data_list, colors, labels):
        for row in (1, 2):
            fig.add_trace(go.Bar(
                x=xpos, y=list(data), name=label, marker_color=color,
                marker_line=dict(color=style.edgecolor, width=style.linewidth * scale),
                legendgroup=label, showlegend=(row == 1),
                customdata=list(xlabels),
                hovertemplate='%{customdata}<br>' + f'{label}: %{{y}}<extra></extra>',
            ), row=row, col=1)

    fig.update_layout(barmode='stack',
                      **_layout(output, scale, axis.title, showlegend=True))
    kwargs = _axis_kwargs(scale)
    fig.update_yaxes(**{**kwargs, 'range': [break_hi, axis.ylim[1] if axis.ylim else None],
                        'tickmode': 'array',
                        'tickvals': [t for t in yticks if t >= break_hi]}, row=1, col=1)
    fig.update_yaxes(**{**kwargs, 'range': [axis.ylim[0] if axis.ylim else 0, style.ybreak],
                        'tickmode': 'array',
                        'tickvals': [t for t in yticks if t <= style.ybreak]},
                     title_text=axis.ylabel, row=2, col=1)
    fig.update_xaxes(**kwargs, row=1, col=1)
    fig.update_xaxes(**{**kwargs, 'tickmode': 'array', 'tickvals': list(xpos),
                        'ticktext': [l.split(' ')[0] for l in xlabels],
                        'tickangle': -style.xlabel_rotation},
                     title_text=axis.xlabel, row=2, col=1)
    return finish(fig, output)


def posneg_stacked_bar_chart(
    xlabels, pos_data_list, neg_data_list, colors, labels,
    ylim_top=None, yticks_top=None, ylim_bot=None, yticks_bot=None,

    axis: AxisLabelOpts = AXIS,
    style: BarStyleOpts = BARSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = DEFAULT_SCALE,
    ):
    """Interactive twin of ``charts.posneg_stacked_bar_chart``."""
    n_pos = len(pos_data_list)
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.05)

    for data, color, label in zip(pos_data_list, colors[:n_pos], labels[:n_pos]):
        fig.add_trace(go.Bar(
            x=list(xlabels), y=list(data), name=label, marker_color=color,
            marker_line=dict(color=style.edgecolor, width=style.linewidth * scale),
            hovertemplate='%{x}<br>' + f'{label}: %{{y}}<extra></extra>'), row=1, col=1)
    for data, color, label in zip(neg_data_list, colors[n_pos:], labels[n_pos:]):
        fig.add_trace(go.Bar(
            x=list(xlabels), y=[-v for v in data], name=label, marker_color=color,
            marker_line=dict(color=style.edgecolor, width=style.linewidth * scale),
            hovertemplate='%{x}<br>' + f'{label}: %{{y}}<extra></extra>'), row=2, col=1)

    fig.update_layout(barmode='relative',
                      **_layout(output, scale, axis.title, showlegend=True))
    kwargs = _axis_kwargs(scale)
    fig.update_yaxes(**{**kwargs, 'tickmode': 'array', 'tickvals': list(yticks_top or []),
                        'range': list(ylim_top) if ylim_top else None},
                     title_text=axis.ylabel, row=1, col=1)
    fig.update_yaxes(**{**kwargs, 'tickmode': 'array', 'tickvals': list(yticks_bot or []),
                        'range': list(ylim_bot) if ylim_bot else None},
                     title_text=axis.ylabel, row=2, col=1)
    fig.update_xaxes(**kwargs, row=1, col=1)
    fig.update_xaxes(**{**kwargs, 'tickangle': -style.xlabel_rotation},
                     title_text=axis.xlabel, row=2, col=1)
    return finish(fig, output)


# ---------------------------------------------------------------------------
# DISTRIBUTIONS
# ---------------------------------------------------------------------------


def distance_histogram_figure(
    df, x_col, hue_col, palette,

    axis: AxisLabelOpts = AXIS,
    style: HistogramStyleOpts = HISTOGRAMSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = DEFAULT_SCALE,
    ):
    """Interactive twin of ``charts.distance_histogram_figure``."""
    df_data = load_frame(df)
    fig = go.Figure()
    for value, color in palette.items():
        subset = df_data[df_data[hue_col] == value]
        fig.add_trace(go.Histogram(
            x=subset[x_col], name=str(value), marker_color=color,
            marker_line=dict(color='black', width=style.linewidth * scale),
            xbins=dict(size=style.binwidth) if style.binwidth else None,
            hovertemplate=f'{value}<br>{x_col}: %{{x}}<br>count: %{{y}}<extra></extra>',
        ))

    fig.update_layout(barmode='stack',
                      **_layout(output, scale, axis.title, showlegend=True))
    _apply(fig, axis, scale, transparent=output.transparent)
    return finish(fig, output)


def pie_chart_figure(
    values, labels, colors,

    axis: AxisLabelOpts = AXIS,
    style: PieStyleOpts = PIESTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = DEFAULT_SCALE,
    ):
    """Interactive twin of ``charts.pie_chart_figure``. None if all zero."""
    if sum(values) == 0:
        return None

    fig = go.Figure(go.Pie(
        values=list(values), labels=list(labels),
        marker=dict(colors=list(colors),
                    line=dict(color=style.edgecolor, width=style.linewidth * scale)),
        textinfo='label+value', sort=False,
        direction='counterclockwise' if style.counterclock else 'clockwise',
        hovertemplate='%{label}<br>%{value} guides (%{percent})<extra></extra>',
    ))
    fig.update_layout(**_layout(output, scale, axis.title, showlegend=True))
    return finish(fig, output)


# ---------------------------------------------------------------------------
# HEATMAP
# ---------------------------------------------------------------------------


def gene_condition_heatmap(
    df, conditions, genes, gene_col='Gene',

    axis: AxisLabelOpts = AXIS,
    style: HeatmapStyleOpts = HEATMAPSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = DEFAULT_SCALE,
    ):
    """Interactive twin of ``charts.gene_condition_heatmap``.

    Returns ``(fig, pivot)``, matching the static version's extra return.
    """
    df_data = load_frame(df)
    pivot = (df_data[df_data[gene_col].isin(genes)]
             .groupby(gene_col)[conditions].mean()
             .reindex(index=genes).T
             .reindex(index=conditions, columns=genes))
    values = pivot.values
    abs_max = float(np.nanmax(np.abs(values)))

    fig = go.Figure(go.Heatmap(
        z=values, x=list(pivot.columns), y=list(pivot.index),
        colorscale=colorscale(style.cmap), zmid=style.center, zmin=-abs_max, zmax=abs_max,
        colorbar=dict(title=dict(text='Mean Z',
                                 font=dict(size=bstyle.BASE_FONT * scale)),
                      tickfont=dict(size=bstyle.BASE_FONT * scale),
                      thickness=8 * scale, len=0.9),
        hovertemplate='%{x}<br>%{y}<br>mean Z: %{z:.2f}<extra></extra>',
    ))
    fig.update_layout(**_layout(output, scale, axis.title))
    kwargs = _axis_kwargs(scale)
    fig.update_xaxes(**{**kwargs, 'tickangle': -style.xlabel_rotation},
                     title_text=axis.xlabel or 'Gene')
    fig.update_yaxes(**kwargs, title_text=axis.ylabel or 'Condition')
    finish(fig, output)
    return fig, pivot


# ---------------------------------------------------------------------------
# LOLLIPOP
# ---------------------------------------------------------------------------


def lollipop_split_figure(
    df_line, line_xcol, line_ycol,

    df_pos, pos_xcol, pos_ycol, pos_threshold,
    pos_category_col, pos_categories, pos_colors,
    ylim_top, yticks_top,

    df_neg, neg_xcol, neg_ycol, neg_threshold,
    neg_category_col, neg_colors,
    ylim_bot, yticks_bot,

    axis: AxisLabelOpts = AXIS,
    domain: DomainOpts = DOMAIN,
    style: SplitLollipopStyleOpts = SPLITLOLLIPOPSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = DEFAULT_SCALE,
    hover_cols=('mutations',),
    ):
    """Interactive twin of ``lollipop.lollipop_split_figure``.

    Hovering a screen hit gives its position, condition, editor and the
    amino-acid edits the guide makes; hovering a lower stem gives the patient
    mutation count at that residue.
    """
    df_pos = load_frame(df_pos)
    df_neg = load_frame(df_neg)
    df_line = load_frame(df_line)

    pos_hits = df_pos.loc[df_pos[pos_ycol] > pos_threshold]
    neg_hits = df_neg.loc[df_neg[neg_ycol] < neg_threshold]

    # A single dominant hotspot (KRAS G12, BRAF V600) would otherwise flatten
    # every other stem, which is why the static figure breaks this axis.
    broken = style.ybreak is not None
    rows = 4 if broken else 3
    heights = [0.20, 0.36, 0.22, 0.22] if broken else [0.22, 0.42, 0.36]
    fig = make_subplots(rows=rows, cols=1, shared_xaxes=True,
                        row_heights=heights, vertical_spacing=0.04)
    neg_rows = (3, 4) if broken else (3,)

    for dom in resolve_domains(domain, df_pos):
        for row in range(1, rows + 1):
            fig.add_shape(row=row, col=1, type='rect', xref='x', yref='paper',
                          x0=dom['start'], x1=dom['end'], y0=0, y1=1,
                          fillcolor=bstyle.to_css(dom.get('color', domain.default_color)),
                          opacity=dom.get('alpha', domain.default_alpha),
                          line_width=0, layer='below')

    line_sorted = df_line.sort_values(by=line_xcol)
    fig.add_trace(go.Scatter(
        x=line_sorted[line_xcol], y=line_sorted[line_ycol], mode='lines',
        line=bstyle.line_from_kws(style.distance_line_kws, scale),
        name=style.distance_label,
        hovertemplate='residue %{x}<br>%{y:.1f} A from ligand<extra></extra>',
    ), row=1, col=1)

    marker_symbols = {'o': 'circle', 'D': 'diamond', 's': 'square', '^': 'triangle-up'}
    for category, color in zip(pos_categories, pos_colors):
        subset = pos_hits[pos_hits[pos_category_col] == category]
        for editor, marker in style.editor_markers.items():
            points = subset[subset[style.editor_col] == editor]
            if points.empty:
                continue
            stem_x, stem_y = [], []
            for x, y in zip(points[pos_xcol], points[pos_ycol]):
                stem_x += [x, x, None]
                stem_y += [0, y, None]
            fig.add_trace(go.Scatter(
                x=stem_x, y=stem_y, mode='lines',
                line=bstyle.line_from_kws(style.stem_kws, scale),
                showlegend=False, hoverinfo='skip', legendgroup=category,
            ), row=2, col=1)
            text = bstyle.hover_text(points, list(hover_cols))
            fig.add_trace(go.Scatter(
                x=points[pos_xcol], y=points[pos_ycol], mode='markers',
                name=f'{category} {editor}', legendgroup=category,
                marker=dict(symbol=marker_symbols.get(marker, marker),
                            size=style.marker_kws.get('markersize', 3) * scale,
                            color=color,
                            line=dict(color='black',
                                      width=style.marker_kws.get('markeredgewidth', 0.4) * scale)),
                hovertext=text,
                hovertemplate=('residue %{x}<br>Z: %{y:.2f}<br>' + f'{category}, {editor}'
                               + ('<br>%{hovertext}' if text is not None else '')
                               + '<extra></extra>'),
            ), row=2, col=1)

    if not neg_hits.empty:
        stem_x, stem_y = [], []
        for x, y in zip(neg_hits[neg_xcol], neg_hits[neg_ycol]):
            stem_x += [x, x, None]
            stem_y += [0, y, None]
        for index, row in enumerate(neg_rows):
            fig.add_trace(go.Scatter(
                x=stem_x, y=stem_y, mode='lines',
                line=bstyle.line_from_kws(style.stem_kws, scale),
                showlegend=False, hoverinfo='skip'), row=row, col=1)
            fig.add_trace(go.Scatter(
                x=neg_hits[neg_xcol], y=neg_hits[neg_ycol], mode='markers',
                name='cBioPortal', legendgroup='cBioPortal', showlegend=(index == 0),
                marker=dict(symbol='circle',
                            size=style.marker_kws.get('markersize', 3) * scale,
                            color=neg_colors[0],
                            line=dict(color='black', width=0.4 * scale)),
                customdata=(-neg_hits[neg_ycol]).tolist(),
                hovertemplate='residue %{x}<br>%{customdata} patient mutations<extra></extra>',
            ), row=row, col=1)

    fig.update_layout(**_layout(output, scale, axis.title, showlegend=True))
    kwargs = _axis_kwargs(scale)
    xkw = dict(kwargs)
    if axis.xlim is not None:
        xkw['range'] = list(axis.xlim)
    if axis.xticks is not None:
        xkw.update(tickmode='array', tickvals=list(axis.xticks))

    fig.update_yaxes(**{**kwargs, 'range': [style.distance_yticks[0], style.distance_yticks[-1]],
                        'tickmode': 'array', 'tickvals': list(style.distance_yticks)},
                     title_text=style.distance_label, row=1, col=1)
    fig.update_yaxes(**{**kwargs, 'range': list(ylim_top),
                        'tickmode': 'array', 'tickvals': list(yticks_top)},
                     title_text=axis.ylabel, row=2, col=1)
    if broken:
        fig.update_yaxes(**{**kwargs, 'range': [style.ybreak[1], ylim_bot[1]],
                            'tickmode': 'array',
                            'tickvals': [t for t in yticks_bot if t >= style.ybreak[1]]},
                         title_text='Patients', row=3, col=1)
        fig.update_yaxes(**{**kwargs, 'range': [ylim_bot[0], style.ybreak[0]],
                            'tickmode': 'array',
                            'tickvals': [t for t in yticks_bot if t <= style.ybreak[0]]},
                         row=4, col=1)
    else:
        fig.update_yaxes(**{**kwargs, 'range': list(ylim_bot),
                            'tickmode': 'array', 'tickvals': list(yticks_bot)},
                         title_text='Patients', row=3, col=1)
    for row in range(1, rows + 1):
        fig.update_xaxes(**xkw, row=row, col=1)
    fig.update_xaxes(title_text=axis.xlabel, row=rows, col=1)
    return finish(fig, output)


# ---------------------------------------------------------------------------
# TRAJECTORY
# ---------------------------------------------------------------------------


def trajectory_scatter(
    fit, genes_list, annotations=None, color_by='cluster',

    axis: AxisLabelOpts = AXIS,
    style: TrajectoryStyleOpts = TRAJECTORYSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = 3.0,
    ):
    """PCA projection coloured by cluster or pseudotime, with the MST drawn on.

    The interactive counterpart of the scatter panels ``trajectory.plot_trajectory``
    writes. Hovering a point gives its gene, cluster, pseudotime and, when
    supplied, the amino-acid edits the guide makes.
    """
    import matplotlib.pyplot as plt

    coords, labels = fit['coords'], fit['hard_labels']
    comp_to_label, comp_to_rank = fit['comp_to_label'], fit['comp_to_rank']
    means, mst = fit['means_pca'], fit['mst']
    pseudotime, point_pt = fit['pseudotime_component'], fit['pseudotime_points']
    n, ev = fit['n_components'], fit['explained_variance']

    cluster_cmap = plt.get_cmap(style.cluster_cmap)
    pt_cmap = plt.get_cmap(style.trajectory_cmap)

    frame = pd.DataFrame({
        'PC1': coords[:, 0], 'PC2': coords[:, 1],
        'gene': np.asarray(genes_list),
        'cluster': comp_to_label[labels], 'pseudotime': point_pt.round(3),
    })
    if annotations is not None:
        for column in annotations.columns:
            frame[column] = annotations[column].to_numpy()
    hover = bstyle.hover_text(frame, ['gene', 'cluster', 'pseudotime']
                             + (list(annotations.columns) if annotations is not None else []))

    fig = go.Figure()
    if color_by == 'cluster':
        for k in range(n):
            mask = labels == k
            if not mask.any():
                continue
            fig.add_trace(go.Scattergl(
                x=frame.loc[mask, 'PC1'], y=frame.loc[mask, 'PC2'], mode='markers',
                name=f"Cluster {comp_to_label[k]}",
                marker=dict(size=style.point_size * scale / 2,
                            color=bstyle.to_css(cluster_cmap(comp_to_rank[k] / max(n - 1, 1))),
                            opacity=style.point_alpha),
                hovertext=hover[mask], hovertemplate='%{hovertext}<extra></extra>',
            ))
    else:
        fig.add_trace(go.Scattergl(
            x=frame['PC1'], y=frame['PC2'], mode='markers',
            marker=dict(size=style.point_size * scale / 2, color=point_pt,
                        colorscale='Reds_r', cmin=0, cmax=1, showscale=True,
                        colorbar=dict(title=dict(text='Pseudotime'),
                                      thickness=8 * scale, len=0.8)),
            hovertext=hover, hovertemplate='%{hovertext}<extra></extra>', showlegend=False,
        ))

    edge_x, edge_y = [], []
    for i in range(n):
        for j in range(i + 1, n):
            if mst[i, j] > 0:
                edge_x += [means[i, 0], means[j, 0], None]
                edge_y += [means[i, 1], means[j, 1], None]
    fig.add_trace(go.Scatter(
        x=edge_x, y=edge_y, mode='lines',
        line=dict(color=style.edge_kws.get('color', 'black'),
                  width=style.edge_kws.get('lw', 0.5) * scale),
        opacity=style.edge_kws.get('alpha', 0.5), showlegend=False, hoverinfo='skip'))
    fig.add_trace(go.Scatter(
        x=means[:, 0], y=means[:, 1], mode='markers+text',
        text=[str(comp_to_label[k]) for k in range(n)], textposition='middle center',
        textfont=dict(size=bstyle.BASE_FONT * scale, color='black'),
        marker=dict(size=style.node_size * scale / 3,
                    color=[bstyle.to_css(pt_cmap(pseudotime[k])) for k in range(n)],
                    line=dict(color='black', width=0.4 * scale)),
        customdata=[[int(comp_to_label[k]), float(pseudotime[k])] for k in range(n)],
        hovertemplate='Cluster %{customdata[0]}<br>pseudotime %{customdata[1]:.3f}<extra></extra>',
        showlegend=False, name='clusters'))

    fig.update_layout(**_layout(output, scale, axis.title,
                                showlegend=(color_by == 'cluster')))
    kwargs = _axis_kwargs(scale)
    fig.update_xaxes(**{**kwargs, 'range': list(axis.xlim or (-12, 18)),
                        'tickmode': 'array',
                        'tickvals': list(axis.xticks or (-12, -6, 0, 6, 12, 18))},
                     title_text=f"PC-1 ({ev[0] * 100:.1f}% var)")
    fig.update_yaxes(**{**kwargs, 'range': list(axis.ylim or (-9, 9)),
                        'tickmode': 'array',
                        'tickvals': list(axis.yticks or (-9, -6, -3, 0, 3, 6, 9))},
                     title_text=f"PC-2 ({ev[1] * 100:.1f}% var)")
    return finish(fig, output)


def _discrete_colorscale(colors):
    """Plotly colorscale of ``len(colors)`` flat bands.

    Paired with ``zmin=-0.5``, ``zmax=len(colors) - 0.5``, band *i* covers
    exactly the integer *i*, so a categorical matrix maps to fixed colours
    instead of being interpolated.
    """
    n = len(colors)
    scale = []
    for index, color in enumerate(colors):
        scale.append([index / n, color])
        scale.append([(index + 1) / n, color])
    return scale


def trajectory_heatmap(
    fit, genes_list, components=None,

    axis: AxisLabelOpts = AXIS,
    style: TrajectoryStyleOpts = TRAJECTORYSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = 1.6,
    ):
    """Interactive twin of the clustered heatmap ``trajectory.plot_trajectory`` writes.

    Guides are rows and conditions are columns, blocked by cluster along the
    trajectory and ordered by pseudotime within each block. A strip of gene
    bars sits to the right, one column per gene, marking which gene each guide
    targets — the same panel the static figure carries, and the reason the
    figure is read as "which nodes populate which part of the landscape".

    ``components`` selects a subset of GMM components, in the order given, for
    the per-branch heatmaps; the default is every component in trajectory order.
    """
    values, labels = fit['values'], fit['hard_labels']
    point_pt = fit['pseudotime_points']
    genes_arr = np.asarray(genes_list)
    genes_ordered = list(cfg.COLOR_MAP)

    chosen = list(fit['trajectory_order'] if components is None else components)
    order, boundaries, tick_pos, tick_text = [], [], [], []
    for component in chosen:
        indices = np.where(labels == component)[0]
        if len(indices) == 0:
            continue
        indices = indices[np.argsort(point_pt[indices])]
        tick_pos.append(len(order) + len(indices) / 2)
        tick_text.append(f"Cluster {fit['comp_to_label'][component]} (n={len(indices)})")
        order.extend(indices.tolist())
        boundaries.append(len(order))
    if not order:
        return None

    ordered = np.round(values[order], 3)
    row_genes = genes_arr[order]
    # vmax spans the whole fit, not the selected subset, so every branch
    # heatmap shares one colour scale with the full one.
    vmax = float(np.abs(values).max())

    # Two panels sharing the row axis: scores on the left, gene identity on
    # the right. Widths follow the static figure's gridspec ratios.
    gene_width = style.col_width * len(genes_ordered)
    fig = make_subplots(
        rows=1, cols=2, shared_yaxes=True, horizontal_spacing=0.015,
        column_widths=[style.heatmap_width, gene_width],
    )

    hover = [[f"{row_genes[r]}<br>cluster {fit['comp_to_label'][labels[i]]}"
              f"<br>pseudotime {point_pt[i]:.3f}" for _ in fit['columns']]
             for r, i in enumerate(order)]
    fig.add_trace(go.Heatmap(
        z=ordered, x=fit['columns'], y=list(range(len(order))),
        colorscale=colorscale(style.heatmap_cmap), zmid=0, zmin=-vmax, zmax=vmax,
        customdata=hover,
        hovertemplate='%{x}<br>Z: %{z:.2f}<br>%{customdata}<extra></extra>',
        colorbar=dict(title=dict(text='Z-score',
                                 font=dict(size=bstyle.BASE_FONT * scale)),
                      tickfont=dict(size=bstyle.BASE_FONT * scale),
                      thickness=8 * scale, len=0.5, y=0.5, x=1.02),
    ), row=1, col=1)

    # Gene bars. One column per gene; a cell is "on" when that guide targets
    # that gene, so each row has exactly one mark. Encoded as a categorical
    # matrix with a banded colorscale rather than an image, so it stays
    # hoverable and rescales with the panel.
    gene_index = {gene: i for i, gene in enumerate(genes_ordered)}
    codes = np.zeros((len(order), len(genes_ordered)))
    for row, gene in enumerate(row_genes):
        column = gene_index.get(gene)
        if column is not None:
            codes[row, column] = column + 1

    off = bstyle.to_css(style.off_color)
    on_colors = ([bstyle.to_css(style.gene_bar_gray)] * len(genes_ordered)
                 if style.gray_gene_bars
                 else [cfg.COLOR_MAP[g] for g in genes_ordered])
    fig.add_trace(go.Heatmap(
        z=codes, x=genes_ordered, y=list(range(len(order))),
        colorscale=_discrete_colorscale([off] + on_colors),
        zmin=-0.5, zmax=len(genes_ordered) + 0.5, showscale=False,
        customdata=[[g] * len(genes_ordered) for g in row_genes],
        hovertemplate='%{customdata}<extra></extra>',
    ), row=1, col=2)

    for boundary in boundaries[:-1]:
        for column in (1, 2):
            fig.add_hline(y=boundary - 0.5, row=1, col=column,
                          line=dict(color='black', width=style.spine_linewidth * scale))

    fig.update_layout(**_layout(output, scale, axis.title))
    kwargs = _axis_kwargs(scale)
    fig.update_xaxes(**{**kwargs, 'tickangle': -90}, row=1, col=1)
    fig.update_xaxes(**{**kwargs, 'tickangle': -90}, row=1, col=2)
    fig.update_yaxes(**{**kwargs, 'tickmode': 'array', 'tickvals': tick_pos,
                        'ticktext': tick_text, 'autorange': 'reversed'},
                     row=1, col=1)
    fig.update_yaxes(**{**kwargs, 'showticklabels': False, 'ticks': '',
                        'autorange': 'reversed'}, row=1, col=2)
    return finish(fig, output)
