"""Gaussian mixture clustering with a minimum-spanning-tree trajectory.

Guides are clustered on their condition Z-score profiles with a Gaussian
mixture model. The component means are joined by a minimum spanning tree, and
breadth-first distance from a chosen root component gives each component a
pseudotime; every guide inherits a soft pseudotime from its posterior
probabilities. PCA is used for display only, never for fitting.

Components are relabelled by pseudotime rank, so "Cluster 1" is the root and
numbering follows the trajectory. Cluster colors follow the same rank, so the
colormap runs along the trajectory rather than along arbitrary GMM component
indices.
"""

from __future__ import annotations

from collections import Counter, deque
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba
from matplotlib.patches import Ellipse
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture

from be_scan.figure_plot import arial_font6

from . import config as cfg
from .figure_options import (
    AXIS, OUTPUT, TRAJECTORYSTYLE,
    AxisLabelOpts, OutputOpts, TrajectoryStyleOpts,
)

def _apply_spine_style(ax) -> None:
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', 2))
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(arial_font6)


def _save(fig, path: Path, dpi: int, out_type: str, **kwargs) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, format=out_type, dpi=dpi, transparent=True, **kwargs)


def _covariance_for(gmm, k: int, covariance_type: str, n_features: int) -> np.ndarray:
    if covariance_type == 'full':
        return gmm.covariances_[k]
    if covariance_type == 'tied':
        return gmm.covariances_
    if covariance_type == 'diag':
        return np.diag(gmm.covariances_[k])
    return np.eye(n_features) * gmm.covariances_[k]


def fit_trajectory(
    df: pd.DataFrame,
    n_components: int,
    trajectory_root: int | None = None,
    covariance_type: str = 'full',
    random_state: int = cfg.RANDOM_STATE,
    n_init: int = 5,
) -> dict:
    """Fit the mixture, build the trajectory, and return every derived array.

    Separated from plotting so the numbers can be inspected or re-plotted
    without refitting.
    """
    values = np.asarray(df)

    gmm = GaussianMixture(
        n_components=n_components, covariance_type=covariance_type,
        random_state=random_state, n_init=n_init, max_iter=300,
    )
    gmm.fit(values)

    hard_labels = gmm.predict(values)
    soft_probs = gmm.predict_proba(values)
    confidence = soft_probs.max(axis=1)

    pca = PCA(n_components=2, random_state=random_state)
    coords = pca.fit_transform(values)
    means_pca = pca.transform(gmm.means_)

    weights = cdist(gmm.means_, gmm.means_, metric='euclidean')
    mst = minimum_spanning_tree(weights).toarray()
    mst_sym = mst + mst.T

    if trajectory_root is None:
        intra_var = {
            k: (np.var(values[hard_labels == k], axis=0).mean()
                if np.sum(hard_labels == k) > 1 else 0.0)
            for k in range(n_components)
        }
        root_idx = int(max(intra_var, key=intra_var.get))
    else:
        root_idx = int(trajectory_root)

    pseudotime = np.full(n_components, np.inf)
    pseudotime[root_idx] = 0.0
    visited = [False] * n_components
    visited[root_idx] = True
    queue = deque([root_idx])
    while queue:
        node = queue.popleft()
        for neighbor in np.where(mst_sym[node] > 0)[0]:
            if not visited[neighbor]:
                visited[neighbor] = True
                pseudotime[neighbor] = pseudotime[node] + mst_sym[node, neighbor]
                queue.append(neighbor)

    finite = pseudotime[pseudotime < np.inf]
    pseudotime = (pseudotime - finite.min()) / (finite.max() - finite.min() + 1e-9)
    point_pseudotime = soft_probs @ pseudotime

    trajectory_order = np.argsort(pseudotime)
    comp_to_rank = np.empty(n_components, dtype=int)
    comp_to_rank[trajectory_order] = np.arange(n_components)

    return {
        'gmm': gmm, 'pca': pca, 'values': values,
        'columns': list(df.columns),
        'hard_labels': hard_labels, 'soft_probs': soft_probs,
        'confidence': confidence, 'coords': coords, 'means_pca': means_pca,
        'mst': mst_sym, 'root_idx': root_idx,
        'pseudotime_component': pseudotime, 'pseudotime_points': point_pseudotime,
        'trajectory_order': trajectory_order,
        'comp_to_rank': comp_to_rank, 'comp_to_label': comp_to_rank + 1,
        'explained_variance': pca.explained_variance_ratio_,
        'converged': gmm.converged_, 'bic': gmm.bic(values), 'aic': gmm.aic(values),
        'log_likelihood': gmm.score(values),
        'covariance_type': covariance_type, 'n_components': n_components,
    }


def trajectory_summary(fit: dict) -> pd.DataFrame:
    """One-row-per-cluster summary of the fitted trajectory."""
    rows = []
    for component in fit['trajectory_order']:
        members = fit['hard_labels'] == component
        rows.append({
            'trajectory_cluster': int(fit['comp_to_label'][component]),
            'gmm_component': int(component),
            'is_root': bool(component == fit['root_idx']),
            'n_guides': int(members.sum()),
            'pseudotime': float(fit['pseudotime_component'][component]),
            'mean_confidence': float(fit['confidence'][members].mean()) if members.any() else np.nan,
        })
    return pd.DataFrame(rows)


def cluster_gene_composition(fit: dict, genes_list) -> pd.DataFrame:
    """Gene makeup of each cluster, in trajectory order."""
    genes_arr = np.asarray(genes_list)
    rows = []
    for component in fit['trajectory_order']:
        members = fit['hard_labels'] == component
        counts = Counter(genes_arr[members])
        total = int(members.sum())
        for gene, count in counts.most_common():
            rows.append({
                'trajectory_cluster': int(fit['comp_to_label'][component]),
                'Gene': gene, 'n_guides': count, 'n_cluster': total,
                'fraction': count / total if total else np.nan,
            })
    return pd.DataFrame(rows)


def plot_trajectory(
    fit, genes_list,
    trajectory_clusters=None,
    annotations=None,
    per_cluster_panels=True,

    axis: AxisLabelOpts = AXIS,
    style: TrajectoryStyleOpts = TRAJECTORYSTYLE,
    output: OutputOpts = OUTPUT,
    ):
    """Draw every trajectory panel and heatmap; return the per-guide table.

    ``output.path`` is treated as a stem, as everywhere else: its parent is the
    destination directory and its name prefixes each panel's filename.
    """
    stem = Path(output.path)
    out_dir, prefix = stem.parent, stem.name
    out_dir.mkdir(parents=True, exist_ok=True)

    dpi, out_type = output.dpi, output.out_type
    xlim, ylim = axis.xlim or (-12, 18), axis.ylim or (-9, 9)
    xticks = axis.xticks or (-12, -6, 0, 6, 12, 18)
    yticks = axis.yticks or (-9, -6, -3, 0, 3, 6, 9)
    trajectory_cmap = style.trajectory_cmap

    values = fit['values']
    coords = fit['coords']
    means_pca = fit['means_pca']
    mst_sym = fit['mst']
    hard_labels = fit['hard_labels']
    pseudotime = fit['pseudotime_component']
    point_pseudotime = fit['pseudotime_points']
    comp_to_label = fit['comp_to_label']
    comp_to_rank = fit['comp_to_rank']
    n = fit['n_components']
    root_idx = fit['root_idx']
    ev = fit['explained_variance']
    col_names = fit['columns']

    cluster_cmap = plt.get_cmap(style.cluster_cmap)
    pt_cmap = plt.get_cmap(trajectory_cmap)
    ambiguous = fit['confidence'] < style.prob_threshold

    def cluster_color(k):
        return cluster_cmap(comp_to_rank[k] / max(n - 1, 1))

    def pca_covariance(k):
        components = fit['pca'].components_
        cov = _covariance_for(fit['gmm'], k, fit['covariance_type'], values.shape[1])
        return components @ cov @ components.T

    def _draw_ellipses(ax, alpha=None, color_by='pseudotime', subset=None):
        alpha = style.ellipse_alpha if alpha is None else alpha
        for k in (range(n) if subset is None else subset):
            vals, vecs = np.linalg.eigh(pca_covariance(k))
            order = vals.argsort()[::-1]
            vals, vecs = vals[order], vecs[:, order]
            angle = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
            width, height = 2 * 2 * np.sqrt(vals)
            color = pt_cmap(pseudotime[k]) if color_by == 'pseudotime' else cluster_color(k)
            ax.add_patch(Ellipse(
                xy=means_pca[k], width=width, height=height, angle=angle,
                facecolor=color, alpha=alpha, edgecolor=color,
                linewidth=0, zorder=2,
            ))

    def draw_trajectory(ax, subset=None):
        nodes = list(range(n)) if subset is None else list(subset)
        nodeset = set(nodes)
        for i in range(n):
            for j in range(i + 1, n):
                if mst_sym[i, j] > 0 and i in nodeset and j in nodeset:
                    ax.plot([means_pca[i, 0], means_pca[j, 0]],
                            [means_pca[i, 1], means_pca[j, 1]],
                            **style.edge_kws, zorder=4)
                    mid = (means_pca[i] + means_pca[j]) / 2
                    src, tgt = ((means_pca[i], means_pca[j])
                                if pseudotime[i] < pseudotime[j]
                                else (means_pca[j], means_pca[i]))
                    delta = (tgt - src) * 0.01
                    ax.annotate('', xy=mid + delta, xytext=mid - delta,
                                arrowprops=dict(arrowstyle='-|>', color='black',
                                                lw=0.5, mutation_scale=5),
                                zorder=5)
        ax.scatter(means_pca[nodes, 0], means_pca[nodes, 1],
                   c=[pt_cmap(pseudotime[k]) for k in nodes],
                   s=style.node_size, edgecolors='black', linewidths=0.25, zorder=6)
        if root_idx in nodeset:
            ax.scatter(means_pca[root_idx, 0], means_pca[root_idx, 1],
                       s=style.root_size, facecolors='none', edgecolors='gold',
                       linewidths=0.25, zorder=7)
        for k in nodes:
            ax.annotate(f"Cluster {comp_to_label[k]}", means_pca[k],
                        fontproperties=arial_font6, ha='center', va='center',
                        color='black', zorder=8)

    def new_axes(title):
        fig, ax = plt.subplots(1, 1, figsize=style.scatter_figsize)
        if ambiguous.any():
            ax.scatter(coords[ambiguous, 0], coords[ambiguous, 1],
                       c='lightgrey', s=5, alpha=0.3, zorder=1)
        ax.set_xlabel(f"PC-1 ({ev[0] * 100:.1f}% var)", fontproperties=arial_font6)
        ax.set_ylabel(f"PC-2 ({ev[1] * 100:.1f}% var)", fontproperties=arial_font6)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xticks(list(xticks))
        ax.set_yticks(list(yticks))
        ax.set_title(title, fontproperties=arial_font6)
        return fig, ax

    # ── Panel: all clusters, colored by cluster ─────────────────────────────
    fig, ax = new_axes('Trajectory - component ID')
    if style.draw_ellipses:
        _draw_ellipses(ax, color_by='cluster')
    handles = []
    for k in range(n):
        mask = (~ambiguous) & (hard_labels == k)
        ax.scatter(coords[mask, 0], coords[mask, 1], color=cluster_color(k),
                   s=style.point_size, alpha=style.point_alpha, zorder=3, rasterized=True)
        if comp_to_label[k] <= 20:
            handles.append(plt.Line2D([0], [0], marker='o', linestyle='',
                                      color=cluster_color(k), markersize=3,
                                      label=f"Cluster {comp_to_label[k]}"))
    draw_trajectory(ax)
    _apply_spine_style(ax)
    fig.tight_layout()
    _save(fig, out_dir / f"{prefix}-Scatter-ClusterColor.{out_type}", dpi, out_type)
    plt.show()
    plt.close(fig)

    # Cluster legend, written separately for figure assembly.
    handles.sort(key=lambda h: int(h.get_label().split()[1]))
    fig_leg, ax_leg = plt.subplots(figsize=(1, max(1.2, len(handles) * 0.09)))
    ax_leg.axis('off')
    legend = ax_leg.legend(
        handles=handles, labels=[h.get_label().split(' ')[1] for h in handles],
        loc='center', ncol=1, frameon=False, title='Clusters', title_fontsize=6,
        prop=arial_font6, markerscale=1.2, handlelength=1.0,
        handletextpad=0.4, borderpad=0.3,
    )
    legend.get_title().set_fontproperties(arial_font6)
    fig_leg.tight_layout()
    _save(fig_leg, out_dir / f"{prefix}-ClusterColor-Legend.{out_type}", dpi, out_type)
    plt.show()
    plt.close(fig_leg)

    # ── Panels: one per cluster ─────────────────────────────────────────────
    if per_cluster_panels:
        for k in range(n):
            fig, ax = new_axes(f"Trajectory - cluster {comp_to_label[k]}")
            if style.draw_ellipses:
                _draw_ellipses(ax, color_by='cluster')
            mask = (~ambiguous) & (hard_labels == k)
            ax.scatter(coords[mask, 0], coords[mask, 1], color=cluster_color(k),
                       s=style.point_size, alpha=style.point_alpha, zorder=3, rasterized=True)
            draw_trajectory(ax)
            _apply_spine_style(ax)
            fig.tight_layout()
            _save(fig, out_dir / f"{prefix}-Scatter-ClusterColor-Cluster{comp_to_label[k]}.{out_type}",
                  dpi, out_type)
            plt.show()
            plt.close(fig)

    # ── Panels: trajectory branches ─────────────────────────────────────────
    if trajectory_clusters:
        for group in trajectory_clusters:
            group = list(group)
            label = 'to'.join(str(c) for c in group)
            fig, ax = new_axes(f"Trajectory - clusters {label}")
            if style.draw_ellipses:
                _draw_ellipses(ax, color_by='cluster', subset=group)
            for k in group:
                mask = (~ambiguous) & (hard_labels == k)
                ax.scatter(coords[mask, 0], coords[mask, 1], color=cluster_color(k),
                           s=style.point_size, alpha=style.point_alpha, zorder=3, rasterized=True)
            draw_trajectory(ax, subset=group)
            _apply_spine_style(ax)
            fig.tight_layout()
            _save(fig, out_dir / f"{prefix}-Scatter-ClusterColor-Traj_{label}.{out_type}",
                  dpi, out_type)
            plt.show()
            plt.close(fig)

    # ── Panel: pseudotime coloring ──────────────────────────────────────────
    fig, ax = new_axes('Trajectory pseudotime')
    if style.draw_ellipses:
        _draw_ellipses(ax, color_by='pseudotime')
    keep = ~ambiguous
    ax.scatter(coords[keep, 0], coords[keep, 1], c=point_pseudotime[keep],
               cmap=trajectory_cmap, s=style.point_size, alpha=0.85, vmin=0, vmax=1,
               zorder=3, rasterized=True)
    draw_trajectory(ax)
    _apply_spine_style(ax)
    fig.tight_layout()
    _save(fig, out_dir / f"{prefix}-Scatter-PseudotimeColor.{out_type}", dpi, out_type)
    plt.show()
    plt.close(fig)

    fig_cb, ax_cb = plt.subplots(figsize=(0.15, 0.5))
    colorbar = mpl.colorbar.ColorbarBase(
        ax_cb, cmap=plt.get_cmap(trajectory_cmap),
        norm=mcolors.Normalize(vmin=0, vmax=1), orientation='vertical')
    colorbar.set_label('Pseudotime', fontproperties=arial_font6)
    for label in ax_cb.get_yticklabels():
        label.set_fontproperties(arial_font6)
    ax_cb.tick_params(width=0.5, length=2)
    fig_cb.tight_layout()
    _save(fig_cb, out_dir / f"{prefix}-PseudotimeColor-ColorBar.{out_type}", dpi, out_type)
    plt.show()
    plt.close(fig_cb)

    # ── Heatmaps ────────────────────────────────────────────────────────────
    vmax = np.abs(values).max()
    genes_ordered = list(cfg.COLOR_MAP)
    genes_arr = np.asarray(genes_list)

    def draw_heatmap(ordered_components, suffix=''):
        blocks = []
        for k in ordered_components:
            indices = np.where(hard_labels == k)[0]
            if len(indices) == 0:
                continue
            ordered = indices[np.argsort(point_pseudotime[indices])]
            blocks.append({'component': k, 'data': values[ordered],
                           'genes': genes_arr[ordered], 'n': len(ordered)})
        if not blocks:
            return

        total_cols = len(genes_ordered) + 1
        fig = plt.figure(figsize=style.heatmap_figsize if suffix == '' else style.small_heatmap_figsize)
        grid = fig.add_gridspec(
            nrows=len(blocks), ncols=total_cols,
            height_ratios=[max(b['n'] * style.row_height, 0.1) for b in blocks],
            width_ratios=[style.heatmap_width] + [style.col_width] * len(genes_ordered),
            wspace=0, hspace=style.hspace,
        )

        gene_bar_axes = []
        for row, block in enumerate(blocks):
            ax_hm = fig.add_subplot(grid[row, 0])
            ax_hm.imshow(block['data'], aspect='auto', cmap=style.heatmap_cmap,
                         vmin=-vmax, vmax=vmax, interpolation='nearest')
            for spine in ax_hm.spines.values():
                spine.set_linewidth(style.spine_linewidth)
            ax_hm.spines['top'].set_visible(False)
            ax_hm.spines['right'].set_visible(False)
            ax_hm.set_yticks([block['n'] / 2 - 0.5])
            ax_hm.set_yticklabels(
                [f"Cluster {comp_to_label[block['component']]}  (n={block['n']})"],
                fontproperties=arial_font6)
            ax_hm.tick_params(axis='y', length=0, pad=2)

            if row < len(blocks) - 1:
                ax_hm.set_xticks([])
            else:
                ax_hm.set_xticks(np.arange(len(col_names)))
                ax_hm.set_xticklabels(col_names, rotation=90, fontproperties=arial_font6)
                ax_hm.tick_params(axis='x', length=2, width=style.spine_linewidth, pad=1)

            ax_hm.add_patch(mpatches.Rectangle(
                (-0.5, -0.5), len(col_names), block['n'], fill=False,
                edgecolor='black', linewidth=style.border_linewidth, zorder=10, clip_on=False))
            for c in range(1, len(col_names)):
                ax_hm.axvline(c - 0.5, color='black', lw=0.2, alpha=0.3)

            for j, gene in enumerate(genes_ordered):
                on_color = style.gene_bar_gray if style.gray_gene_bars else to_rgba(cfg.COLOR_MAP[gene])
                strip = np.array([on_color if g == gene else style.off_color
                                  for g in block['genes']]).reshape(-1, 1, 4)
                ax_gene = fig.add_subplot(grid[row, 1 + j])
                ax_gene.imshow(strip, aspect='auto', interpolation='nearest')
                ax_gene.axis('off')
                if row == 0:
                    gene_bar_axes.append(ax_gene)

        for index, ax_gene in enumerate(gene_bar_axes):
            box = ax_gene.get_position()
            fig.text(box.x0 + box.width / 2, 0.01, genes_ordered[index],
                     ha='center', va='top', fontproperties=arial_font6,
                     color='k', rotation=90)

        _save(fig, out_dir / f"{prefix}-Heatmap{suffix}.{out_type}", dpi, out_type,
              bbox_inches='tight')
        plt.show()
        plt.close(fig)

    draw_heatmap(list(fit['trajectory_order']), suffix='')
    if trajectory_clusters:
        for group in trajectory_clusters:
            draw_heatmap(list(group),
                         suffix='-Traj_' + 'to'.join(str(c) for c in group))

    fig_cb, ax_cb = plt.subplots(figsize=(0.15, 1.0))
    colorbar = mpl.colorbar.ColorbarBase(
        ax_cb, cmap=plt.cm.bwr,
        norm=mcolors.Normalize(vmin=-vmax, vmax=vmax), orientation='vertical')
    colorbar.set_label('Z-Score', fontproperties=arial_font6)
    for label in ax_cb.get_yticklabels():
        label.set_fontproperties(arial_font6)
    ax_cb.tick_params(width=style.spine_linewidth, length=2)
    _save(fig_cb, out_dir / f"{prefix}-Heatmap-Colorbar.{out_type}", dpi, out_type,
          bbox_inches='tight')
    plt.show()
    plt.close(fig_cb)

    fig_leg, ax_leg = plt.subplots(figsize=(0.5, max(1.0, len(genes_ordered) * 0.09)))
    ax_leg.axis('off')
    legend = ax_leg.legend(
        handles=[mpatches.Patch(color=cfg.COLOR_MAP[g], label=g) for g in genes_ordered],
        loc='center', ncol=1, frameon=False, title='Gene',
        prop=arial_font6, handlelength=0.6, handleheight=0.6)
    legend.get_title().set_fontproperties(arial_font6)
    _save(fig_leg, out_dir / f"{prefix}-Heatmap-GeneLegend.{out_type}", dpi, out_type,
          bbox_inches='tight')
    plt.show()
    plt.close(fig_leg)

    # ── Per-guide table ─────────────────────────────────────────────────────
    table = pd.DataFrame(values, columns=col_names)
    table['gene'] = np.asarray(genes_list)
    table['gmm_component'] = hard_labels
    table['trajectory_cluster'] = comp_to_label[hard_labels]
    table['confidence'] = fit['confidence']
    table['PC1'] = coords[:, 0]
    table['PC2'] = coords[:, 1]
    table['pseudotime'] = point_pseudotime
    if annotations is not None:
        table = pd.concat([table, annotations.reset_index(drop=True)], axis=1)
    table.to_csv(out_dir / f"{prefix}-Data.csv", index=False)
    return table
