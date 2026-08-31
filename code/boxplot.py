"""Control boxplot carrying its significance bracket.

be_scan's ``boxplot_figure`` draws the boxes, and every per-gene panel uses it
unchanged. The control panel is the one boxplot that exists to make a claim
rather than to show a distribution: that the non-targeting guides and the
essential-splice-site guides are drawn from different populations, separately
for each editor. That claim needs a bracket and a significance label over each
editor's pair, which be_scan has no notion of, so this wraps it — the same
arrangement as `lollipop.py`, and the same call convention: data positionally,
appearance in the frozen option classes, destination in ``OutputOpts``.

The labels are passed in rather than computed here. Which test to run and how
to render its p-value are analysis decisions and live in `analysis.py`; this
module only knows where to put the bar.
"""

from __future__ import annotations

from dataclasses import replace

from be_scan.figure_plot import arial_font6, boxplot_figure, finish_figure

from .figure_options import (
    AXIS, BOXPLOTSTYLE, BRACKET, OUTPUT,
    AxisLabelOpts, BoxStyleOpts, BracketStyleOpts, OutputOpts,
)


def control_boxplot_figure(
    df, x_col, y_col, palette, comparisons, lab_order=None,

    axis: AxisLabelOpts = AXIS,
    style: BoxStyleOpts = BOXPLOTSTYLE,
    bracket: BracketStyleOpts = BRACKET,
    output: OutputOpts = OUTPUT,
    ):
    """Boxplot with a labelled bracket over each compared pair.

    ``comparisons`` is a list of ``(left_label, right_label, text)``; the two
    labels must appear in ``lab_order``, which fixes the x positions. Brackets
    stack upward in the order given.

    Returns ``(fig, ax)``.
    """
    # Saving is suppressed on the inner call so the bracket lands on the
    # figure before it is written out; be_scan saves and closes in one step.
    fig, ax = boxplot_figure(
        df, x_col, y_col, palette, lab_order=lab_order,
        axis=axis, style=style,
        output=replace(output, save=False, show=False),
    )

    order = list(lab_order) if lab_order is not None else list(df[x_col].unique())
    bottom, top = ax.get_ylim()
    span = top - bottom

    for index, (left, right, text) in enumerate(comparisons):
        if left not in order or right not in order:
            continue
        x_left, x_right = order.index(left), order.index(right)
        y = bottom + span * (bracket.start + index * bracket.step)
        tick = span * bracket.tick
        ax.plot([x_left, x_left, x_right, x_right],
                [y - tick, y, y, y - tick],
                color=bracket.color, linewidth=bracket.linewidth,
                clip_on=False, solid_joinstyle='miter')
        ax.text((x_left + x_right) / 2, y + span * bracket.text_pad, text,
                ha='center', va='bottom', clip_on=False,
                fontproperties=arial_font6)

    finish_figure(fig, output, bbox_inches="tight")
    return fig, ax
