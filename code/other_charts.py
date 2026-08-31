"""The five figure types the non-screen manuscript figures are drawn with.

============================  ==========================================
``grid_plot``                 competitive growth line panels
``plot_heatmap``              MEK inhibitor fold-change heatmaps
``plot_mean_sd_curves``       CellTiter-Glo response, BaF3 growth curves
``plot_dose_response``        4-parameter logistic fits
``plot_site_boxplots``        published DMS scores at selected residues
============================  ==========================================

Between them they produce all 24 panels. ``plot_mean_sd_curves`` replaces four
near-duplicate functions spread across the original notebooks
(``plot_drugA_no_errorbars``, ``plot_drugA_with_sd_errorbars`` and two
versions of ``plot_manual_with_errorbars``); they differed only in styling,
which is exposed here as arguments. The 4PL toolkit had been copy-pasted into
four notebooks and is now fitted once, in ``other_figures``.

These are matplotlib and seaborn directly rather than be_scan primitives: the
figures are PDFs of assay data, not the SVG position-axis panels be_scan's
``figure_plot`` is built for. Each returns its figure so a caller can adjust
it further, and each saves only when given an ``output_path``.
"""

from __future__ import annotations

import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator

from . import config as cfg
from .other_figures import fit_dose_response, four_param_logistic


def apply_figure_style() -> None:
    """Set the shared 6 pt Arial style, with text left editable in the PDF.

    Matplotlib defaults are restored first. ``be_scan.figure_plot``, which
    ``code/__init__.py`` imports, sets its own tick sizes, pads and line
    widths on the global rcParams at import time. Those are right for the SVG
    position-axis panels it draws and wrong for these figures: inheriting them
    moves every axes box by several points. Resetting makes the result
    independent of what else has been imported.

    The backend is left alone - resetting it would detach the figures from
    the notebook's inline renderer.
    """
    mpl.rcParams.update({key: value for key, value in mpl.rcParamsDefault.items()
                         if key not in ("backend", "interactive")})
    mpl.rcParams.update(cfg.OTHER_FIG_RCPARAMS)


# ---------------------------------------------------------------------------
# COMPETITIVE GROWTH LINE PANELS
# ---------------------------------------------------------------------------


def grid_plot(
    df, groups, group_col="sgRNA_label", y_col="GFP+ Fraction", x_col="Day",
    condition_col="drug_label", palette=None, legend_order=None,
    nrows=1, ncols=None, panel_size=(1.1, 1.5),
    x_ticks=None, x_pad=0.0, y_lim=(0, 1), y_ticks=None,
    line_width=None, output_path=None,
):
    """Grid of GFP+ fraction vs. day line plots, one panel per group.

    Points are individual replicates; the line follows the mean and error bars
    show the standard deviation across replicates.

    palette      dict of colour per condition, or a list in `legend_order` order.
    legend_order condition order for colour assignment and the legend.
    x_ticks      explicit tick positions; also fixes the x limits.
    x_pad        extra space to the left of the first x tick.
    line_width   line/spine width; None keeps the matplotlib default.
    """
    conditions = list(legend_order) if legend_order is not None \
        else list(df[condition_col].dropna().unique())

    if palette is None:
        palette = dict(zip(conditions, sns.color_palette("pastel", len(conditions))))
    elif not isinstance(palette, dict):
        palette = dict(zip(conditions, list(palette)))
    missing = [c for c in conditions if c not in palette]
    if missing:
        raise ValueError(f"No colour supplied for: {missing}")

    ncols = ncols if ncols is not None else len(groups)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(panel_size[0] * ncols, panel_size[1] * nrows),
        sharey=True, sharex=x_ticks is None,
    )
    axes = np.atleast_1d(axes).flatten()
    line_kw = {} if line_width is None else {"linewidth": line_width}

    for ax, group in zip(axes, groups):
        panel = df[df[group_col] == group]
        shared = dict(data=panel, x=x_col, y=y_col, hue=condition_col,
                      hue_order=conditions, palette=palette, ax=ax, legend=False)
        sns.lineplot(**shared, errorbar="sd", err_style="bars", **line_kw)
        sns.scatterplot(**shared, alpha=0.5, s=10, edgecolor="none", linewidth=0)

        ax.set_ylim(*y_lim)
        ax.set_title(group)
        sns.despine(ax=ax)
        if y_ticks is not None:
            ax.set_yticks(y_ticks)
        if x_ticks is None:
            ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        else:
            ax.set_xticks(list(x_ticks))
            ax.set_xlim(min(x_ticks) - x_pad, max(x_ticks))
            ax.tick_params(labelbottom=True)
        if line_width is not None:
            for spine in ax.spines.values():
                spine.set_linewidth(line_width)
            ax.tick_params(width=line_width)

    for ax in axes[len(groups):]:
        ax.set_visible(False)

    # Legend drawn by hand so it shows plain markers rather than error bars.
    handles = [Line2D([0], [0], color=palette[c], marker="o", linestyle="-",
                      markersize=4, linewidth=line_width or 1) for c in conditions]
    fig.legend(handles=handles, labels=conditions, loc="center left",
               bbox_to_anchor=(0.99, 0.5), frameon=False,
               handletextpad=0.4, borderaxespad=0.2, labelspacing=0.3)

    plt.tight_layout(rect=[0, 0, 0.97, 1])
    if output_path is not None:
        plt.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.show()
    return fig


# ---------------------------------------------------------------------------
# FOLD-CHANGE HEATMAPS
# ---------------------------------------------------------------------------


def plot_heatmap(
    df, row_order, col_order, value_col="FC from d0 compared to NT",
    row_col="sgRNA_mut_names", col_col="drug_names",
    figsize=(3.7, 2.7), cmap="bwr", vmin=0.4, vmax=1.4, center=1.0,
    label_fontsize=6, tick_fontsize=6, output_path=None,
):
    """Heatmap of a value column, rows and columns in the order given.

    The colour scale is centred on `center` (1.0 = no change vs. the NT
    control), so depletion and enrichment are visually symmetric even with an
    asymmetric range.
    """
    subset = df[df[row_col].isin(row_order) & df[col_col].isin(col_order)]
    matrix = (subset.pivot_table(index=row_col, columns=col_col,
                                 values=value_col, aggfunc="mean")
              .reindex(index=row_order, columns=col_order))

    data = matrix.values
    vmin = np.nanmin(data) if vmin is None else vmin
    vmax = np.nanmax(data) if vmax is None else vmax

    plt.figure(figsize=figsize)
    im = plt.imshow(data, aspect="auto", cmap=cmap,
                    norm=TwoSlopeNorm(vmin=vmin, vcenter=center, vmax=vmax))
    plt.xticks(range(len(matrix.columns)), matrix.columns,
               rotation=45, ha="right", fontsize=tick_fontsize)
    plt.yticks(range(len(matrix.index)), matrix.index, fontsize=tick_fontsize)
    plt.xlabel(col_col, fontsize=tick_fontsize)
    plt.ylabel(row_col, fontsize=tick_fontsize)

    cbar = plt.colorbar(im)
    cbar.set_label(value_col, fontsize=label_fontsize)
    cbar.ax.tick_params(labelsize=tick_fontsize)

    plt.tight_layout()
    if output_path is not None:
        plt.savefig(output_path, bbox_inches="tight")
    plt.show()
    return matrix


# ---------------------------------------------------------------------------
# REPLICATE POINTS WITH A MEAN +- SD LINE
# ---------------------------------------------------------------------------


def _sci_notation(x, _pos):
    """Compact tick labels: 3E5 for large or tiny values, plain digits otherwise."""
    if x == 0:
        return "0"
    if abs(x) < 0.01 or abs(x) >= 10000:
        exponent = int(np.floor(np.log10(abs(x))))
        coefficient = x / (10 ** exponent)
        if abs(coefficient - round(coefficient)) < 1e-10:
            coefficient = int(round(coefficient))
        return f"{coefficient}E{exponent}"
    if abs(x - round(x)) < 1e-10:
        return str(int(round(x)))
    return f"{x:.2f}".rstrip("0").rstrip(".")


def plot_mean_sd_curves(
    df, x_col, y_col, group_col, group_order, colors, legend_labels=None,
    x_scale="linear", scale_factor=1.0, sci_y=False,
    figsize=(1.5, 1.5), title=None, x_label=None, y_label=None,
    x_lim=None, y_lim=None, x_ticks=None, y_ticks=None, x_tick_labels=None,
    point_alpha=0.55, point_size=10, marker_size=2.5,
    line_width=0.8, elinewidth=0.6, capsize=2, cap_thick=0.6,
    spine_width=0.5, tick_length=2,
    legend_loc="center left", legend_anchor=(1.02, 0.5),
    layout="tight", legend_right_margin=0.72, transparent=False,
    output_path=None,
):
    """Replicate scatter plus a mean line with SD error bars, one series per group.

    Covers both the CellTiter-Glo drug-response panels (x = drug
    concentration) and the BaF3 growth curves (x = day); the two differ only in
    styling.

    scale_factor  divide y by this before plotting (e.g. 1e6 for cells per ml).
    sci_y         use the compact E-notation y tick formatter.
    layout        "tight" runs tight_layout and saves with bbox_inches='tight';
                  "reserve" instead reserves a fixed right margin for the legend,
                  which keeps the plotting area identical across figures.
    """
    data = df.dropna(subset=[x_col, y_col, group_col]).copy()
    data[x_col] = pd.to_numeric(data[x_col], errors="coerce")
    data["_y"] = pd.to_numeric(data[y_col], errors="coerce") / scale_factor
    data = data.dropna(subset=[x_col, "_y"])
    if x_scale == "log":
        data = data[data[x_col] > 0]

    colors = dict(zip(group_order, colors)) if not isinstance(colors, dict) else colors
    if legend_labels is None:
        labels = {g: str(g) for g in group_order}
    elif isinstance(legend_labels, dict):
        labels = {g: str(legend_labels.get(g, g)) for g in group_order}
    else:
        labels = dict(zip(group_order, map(str, legend_labels)))

    summary = (data.groupby([group_col, x_col], observed=True)["_y"]
               .agg(mean="mean", sd="std").reset_index())

    rc = dict(cfg.OTHER_FIG_RCPARAMS,
              **{"axes.linewidth": spine_width, "lines.linewidth": line_width})
    with mpl.rc_context(rc):
        fig, ax = plt.subplots(figsize=figsize)
        drawn = []

        for group in group_order:
            points = data[data[group_col] == group].sort_values(x_col)
            means = summary[summary[group_col] == group].sort_values(x_col)
            if points.empty or means.empty:
                continue
            drawn.append(group)
            ax.scatter(points[x_col], points["_y"], color=colors[group],
                       s=point_size, alpha=point_alpha, linewidths=0,
                       edgecolors="none", zorder=2)
            ax.errorbar(means[x_col].to_numpy(float), means["mean"].to_numpy(float),
                        yerr=means["sd"].to_numpy(float), color=colors[group],
                        marker="o", markersize=marker_size, linestyle="-",
                        linewidth=line_width, elinewidth=elinewidth,
                        capsize=capsize, capthick=cap_thick, zorder=3)

        ax.set_xscale(x_scale)
        ax.set_xlabel(x_label if x_label is not None else x_col)
        ax.set_ylabel(y_label if y_label is not None else y_col)
        if title is not None:
            ax.set_title(title)
        if sci_y:
            ax.yaxis.set_major_formatter(FuncFormatter(_sci_notation))
        if x_lim is not None:
            ax.set_xlim(x_lim)
        if y_lim is not None:
            ax.set_ylim(y_lim)
        if x_ticks is not None:
            ax.set_xticks(list(x_ticks))
        elif x_scale == "linear":
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        if x_tick_labels is not None:
            ax.set_xticklabels(list(x_tick_labels))
        if y_ticks is not None:
            ax.set_yticks(list(y_ticks))

        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_linewidth(spine_width)
            ax.spines[side].set_color("black")
        ax.tick_params(axis="both", which="major", direction="out",
                       length=tick_length, width=spine_width, labelsize=6)

        handles = [Line2D([0], [0], color=colors[g], marker="o", linestyle="-",
                          markersize=3, linewidth=line_width) for g in drawn]
        ax.legend(handles, [labels[g] for g in drawn], loc=legend_loc,
                  bbox_to_anchor=legend_anchor, frameon=False,
                  borderaxespad=0, handlelength=1.5)

        if layout == "tight":
            fig.tight_layout()
            save_kw = {"bbox_inches": "tight"}
        else:
            fig.subplots_adjust(right=legend_right_margin)
            save_kw = {}

        if output_path is not None:
            fig.savefig(output_path, format="pdf", transparent=transparent, **save_kw)
        plt.show()
    return fig, ax


# ---------------------------------------------------------------------------
# 4-PARAMETER LOGISTIC DOSE RESPONSE
# ---------------------------------------------------------------------------


def plot_dose_response(
    df, conc_col, response_col, group_col, group_order, colors,
    label_map=None, x_label="Concentration", y_label="Response", title=None,
    x_tick_labels=None, y_tick_labels=None, figsize=(3, 2.5),
    point_alpha=0.9, point_size=8, capsize=2, output_path=None,
):
    """Replicate points, mean +- SEM, and a fitted 4PL curve per group, log x.

    Returns ``(fig, ax, fits)``; the fits are the same dict
    ``other_figures.fit_dose_response`` returns, so ``ecx_table`` can be run on
    them to read EC50s off the drawn curves.

    x_tick_labels / y_tick_labels take (position, label) pairs.
    """
    label_map = label_map or {}
    colors = dict(zip(group_order, colors)) if not isinstance(colors, dict) else colors

    data = df[[conc_col, response_col, group_col]].dropna()
    data = data[data[conc_col] > 0]
    fits = fit_dose_response(data, conc_col, response_col, group_col, group_order)

    fig, ax = plt.subplots(figsize=figsize)

    for group in group_order:
        sub = data[data[group_col] == group]
        if sub.empty:
            continue
        color = colors[group]

        ax.scatter(sub[conc_col], sub[response_col], color=color,
                   alpha=point_alpha, s=point_size, label="_nolegend_")

        summary = (sub.groupby(conc_col)[response_col]
                   .agg(["mean", "std", "count"]).reset_index().sort_values(conc_col))
        summary["sem"] = summary["std"] / np.sqrt(summary["count"])
        ax.errorbar(summary[conc_col], summary["mean"], yerr=summary["sem"],
                    fmt="o", color=color, capsize=capsize, markersize=3,
                    linewidth=1, label="_nolegend_")

        fit = fits.get(group)
        if fit is not None:
            x_fit = np.logspace(np.log10(sub[conc_col].min()),
                                np.log10(sub[conc_col].max()), 300)
            ax.plot(x_fit, four_param_logistic(x_fit, fit["bottom"], fit["top"],
                                               fit["logEC50"], fit["hill_slope"]),
                    color=color, linewidth=1.5, label=label_map.get(group, group))
        else:
            ax.scatter([], [], color=color, label=label_map.get(group, group))

    ax.set_xscale("log")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if title is not None:
        ax.set_title(title)
    if x_tick_labels is not None:
        positions, labels = zip(*x_tick_labels)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels)
    if y_tick_labels is not None:
        positions, labels = zip(*y_tick_labels)
        ax.set_yticks(positions)
        ax.set_yticklabels(labels)

    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

    fig.tight_layout()
    if output_path is not None:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()
    return fig, ax, fits


# ---------------------------------------------------------------------------
# DMS SCORES AT SELECTED RESIDUES
# ---------------------------------------------------------------------------


def _subset_site(df, site, id_col):
    """Rows whose mutation_id starts at `site` ('Q61' matches Q61R, not Q610)."""
    match = re.match(r"^([A-Za-z])(\d+)$", site.strip())
    if not match:
        raise ValueError("site must look like 'Q61' or 'Y4'")
    site_norm = f"{match.group(1).upper()}{match.group(2)}"
    pattern = re.compile(rf"^{re.escape(site_norm)}(?!\d)")
    ids = df[id_col].fillna("").astype(str).str.strip().str.upper()
    return df[ids.apply(lambda s: bool(pattern.match(s)))].copy(), site_norm


def plot_site_boxplots(
    df, sites, highlight_muts, id_col="mutation_id", label_map=None,
    figsize=(1, 2.7), title=None, highlight_color="#0197C1",
    seed=None, output_path=None,
):
    """Stacked boxplots of z-scored assay values, one row per residue.

    Each box is one published assay; each point is one amino-acid substitution
    at that residue, with ``highlight_muts[i]`` coloured. Point x-jitter is
    random, so the seed is fixed to keep the figure reproducible.
    """
    if len(sites) != len(highlight_muts):
        raise ValueError("sites and highlight_muts must be the same length")

    seed = cfg.OTHER_FIG_JITTER_SEED if seed is None else seed
    label_map = label_map or {}
    score_cols = [c for c in df.select_dtypes(include=[np.number]).columns if c != id_col]
    positions = np.arange(len(score_cols))
    rng = np.random.default_rng(seed)

    fig, axes = plt.subplots(nrows=len(sites), ncols=1, figsize=figsize, sharex=True)
    axes = np.atleast_1d(axes)

    for ax, site, highlight in zip(axes, sites, highlight_muts):
        site_df, site_norm = _subset_site(df, site, id_col)
        if site_df.empty:
            raise ValueError(f"No data found for site '{site_norm}'")

        ax.boxplot([site_df[c].dropna().values for c in score_cols],
                   positions=positions, widths=0.6, showfliers=False,
                   patch_artist=False, boxprops=dict(linewidth=1),
                   whiskerprops=dict(linewidth=1), capprops=dict(linewidth=1),
                   medianprops=dict(linewidth=1, color="black"))

        for j, col in enumerate(score_cols):
            col_data = site_df[[id_col, col]].dropna(subset=[col])
            ids = col_data[id_col].fillna("").astype(str).str.strip().str.upper()
            is_highlight = ids.eq(str(highlight).strip().upper()) if highlight \
                else pd.Series(False, index=col_data.index)
            x = rng.normal(loc=j, scale=0.03, size=len(col_data))

            ax.scatter(x[~is_highlight.values], col_data.loc[~is_highlight, col].values,
                       c="gray", s=8, alpha=0.9, linewidths=0, zorder=3)
            ax.scatter(x[is_highlight.values], col_data.loc[is_highlight, col].values,
                       c=highlight_color, s=12, linewidths=0, zorder=5)

        ax.axhline(0, color="black", linewidth=0.5, linestyle=(0, (5, 3)))
        ax.set_ylabel(site_norm)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(0.8)
        ax.spines["bottom"].set_linewidth(0.8)
        if ax is not axes[-1]:
            ax.tick_params(axis="x", labelbottom=False)

        handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor="gray",
                          markeredgecolor="None", markersize=4,
                          label="Other Mutations")]
        if highlight:
            handles.append(Line2D([0], [0], marker="o", color="none",
                                  markerfacecolor=highlight_color,
                                  markeredgecolor="None", markersize=4,
                                  label=str(highlight)))
        ax.legend(handles=handles, frameon=False, loc="center left",
                  bbox_to_anchor=(1.02, 0.5), borderaxespad=0)

    axes[-1].set_xticks(positions)
    axes[-1].set_xticklabels([label_map.get(c, c) for c in score_cols],
                             rotation=90, ha="right")
    axes[-1].tick_params(axis="x", labelbottom=True)
    if title:
        fig.suptitle(title, y=0.98)

    fig.patch.set_alpha(0)
    for ax in axes:
        ax.patch.set_alpha(0)

    if output_path is not None:
        plt.savefig(output_path, bbox_inches="tight", transparent=True,
                    facecolor="none", edgecolor="none")
    plt.show()
    return fig
