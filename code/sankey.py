"""Sankey diagram of validated resistance hits, inhibitor to gene.

Each link carries the number of guides that scored as a hit in the primary GOF
screen and reproduced in the validation screen. Links are colored by where the
hit gene sits relative to the inhibited node in the pathway: the inhibitor's
own target is gray, genes downstream are pink, genes upstream are blue.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

from . import config as cfg

# Rank in the pathway, used to call a hit upstream or downstream of its
# inhibitor. Both gene symbols and protein names are listed so either
# naming resolves.
PATHWAY_RANK = {
    'SHC1': 0, 'GRB2': 0,
    'SHP2i': 1, 'PTPN11': 1, 'SHP2': 1,
    'KRASi': 2, 'KRAS': 2, 'NRAS': 2, 'MRAS': 2, 'HRAS': 2,
    'PPP1CA': 3, 'PPP1CB': 3, 'PPP1CC': 3, 'SHOC2': 3,
    'RAFi': 4, 'panRAFi': 4, 'ARAF': 4, 'BRAF': 4, 'RAF1': 4, 'CRAF': 4,
    'KSR1': 5, 'KSR2': 5,
    'MEKi': 6, 'MAP2K1': 6, 'MAP2K2': 6, 'MEK1': 6, 'MEK2': 6,
    'DUSP4': 7, 'DUSP6': 7,
    'ERKi': 8, 'MAPK1': 8, 'MAPK3': 8, 'ERK1': 8, 'ERK2': 8,
}

# Each inhibitor's own target(s), whose links are drawn gray.
ON_TARGET = {
    'SHP2i': {'PTPN11', 'SHP2'},
    'KRASi': {'KRAS'},
    'RAFi': {'ARAF', 'BRAF', 'RAF1', 'CRAF'},
    'MEKi': {'MAP2K1', 'MAP2K2', 'MEK1', 'MEK2'},
    'ERKi': {'MAPK1', 'MAPK3', 'ERK1', 'ERK2'},
}

LINK_ON_TARGET = '#D6D6D6'
LINK_DOWNSTREAM = 'pink'
LINK_UPSTREAM = 'lightblue'


def build_hit_connections(
    gof_valid_abe: pd.DataFrame,
    gof_valid_cbe: pd.DataFrame,
    gof_columns: list[str],
    valid_columns: list[str],
    genes: list[str] | None = None,
    hit_cutoff: float = cfg.GOF_HIT_CUTOFF,
    valid_cutoff: float = cfg.SANKEY_VALID_CUTOFF,
) -> pd.DataFrame:
    """Count validated resistance hits for every inhibitor and gene pair.

    A guide counts when its primary GOF Z-score is at least ``hit_cutoff`` and
    its validation Z-score is at least ``valid_cutoff`` in the same condition.
    ABE and CBE counts are summed.
    """
    genes = genes or cfg.GENES
    rows = []
    for gof_col, valid_col in zip(gof_columns, valid_columns):
        for gene in genes:
            count = 0
            for frame in (gof_valid_abe, gof_valid_cbe):
                hits = frame[
                    (frame['Gene'] == gene)
                    & (frame[gof_col] >= hit_cutoff)
                    & (frame[valid_col] >= valid_cutoff)
                ]
                count += len(hits)
            rows.append({
                'inhibitor': cfg.INHIBITOR_LABELS[gof_col],
                'gene': gene,
                'hit_count': count,
            })
    return pd.DataFrame(rows)


def link_color(inhibitor: str, gene: str) -> str:
    """Gray for the inhibitor's own target, pink downstream, blue upstream."""
    if gene in ON_TARGET.get(inhibitor, set()):
        return LINK_ON_TARGET
    if PATHWAY_RANK.get(inhibitor, 0) < PATHWAY_RANK.get(gene, 0):
        return LINK_DOWNSTREAM
    return LINK_UPSTREAM


def directional_summary(connections: pd.DataFrame) -> pd.DataFrame:
    """Per-inhibitor totals of on-target, upstream and downstream hits."""
    coloured = connections.assign(
        color=[link_color(i, g) for i, g in zip(connections['inhibitor'], connections['gene'])]
    )
    lookup = {LINK_ON_TARGET: 'on_target',
              LINK_UPSTREAM: 'upstream',
              LINK_DOWNSTREAM: 'downstream'}
    coloured['direction'] = coloured['color'].map(lookup)
    summary = (coloured.pivot_table(index='inhibitor', columns='direction',
                                    values='hit_count', aggfunc='sum', fill_value=0)
               .reindex(columns=['on_target', 'upstream', 'downstream'], fill_value=0))
    summary['total'] = summary.sum(axis=1)
    return summary.reset_index()


def _spread(n: int, top_margin: float = 0.05, bottom_margin: float = 0.05) -> list[float]:
    """Evenly spaced y positions between the margins."""
    if n == 1:
        return [0.5]
    span = 1 - top_margin - bottom_margin
    return [top_margin + i * span / (n - 1) for i in range(n)]


def sankey_figure(
    connections: pd.DataFrame,
    min_hit_count: int = cfg.SANKEY_MIN_HIT_COUNT,
    use_protein_names: bool = True,
    fixed_layout: bool = True,
    width: int = 300, height: int = 300,
    font_size: int = 6,
):
    """Build the Sankey figure from a hit-connection table.

    Nodes are laid out in pathway order rather than left to Plotly, so the
    diagram reads top-to-bottom along the cascade.
    """
    subset = connections[connections['hit_count'] >= min_hit_count].copy()
    if use_protein_names:
        subset['gene'] = subset['gene'].map(lambda g: cfg.GENE_PROTEIN_MAP.get(g, g))

    inhibitors = [i for i in ['SHP2i', 'KRASi', 'RAFi', 'MEKi', 'ERKi']
                  if i in set(subset['inhibitor'])]
    ordered_genes = cfg.PROTEINS if use_protein_names else cfg.GENES
    targets = [g for g in ordered_genes if g in set(subset['gene'])]

    labels = inhibitors + targets
    index = {label: i for i, label in enumerate(labels)}

    inhibitor_node_colors = [
        cfg.PROTEIN_COLOR_MAP[t]
        for t in (['SHP2', 'KRAS', 'BRAF', 'MEK1', 'ERK2'] if use_protein_names
                  else ['PTPN11', 'KRAS', 'BRAF', 'MAP2K1', 'MAPK1'])
    ][:len(inhibitors)]
    node_colors = inhibitor_node_colors + [cfg.PROTEIN_COLOR_MAP[g] for g in targets]

    link_colors = [link_color(i, g) for i, g in zip(subset['inhibitor'], subset['gene'])]

    node = dict(label=labels, color=node_colors)
    if fixed_layout:
        left_y = _spread(len(inhibitors), top_margin=0.1, bottom_margin=0.1)
        right_y = _spread(len(targets), top_margin=0.02, bottom_margin=0.02)
        positions = {label: (0.01, y) for label, y in zip(inhibitors, left_y)}
        positions.update({label: (0.99, y) for label, y in zip(targets, right_y)})
        node['x'] = [positions[label][0] for label in labels]
        node['y'] = [positions[label][1] for label in labels]

    figure = go.Figure(go.Sankey(
        arrangement='fixed' if fixed_layout else 'snap',
        node=node,
        link=dict(
            arrowlen=15,
            source=[index[i] for i in subset['inhibitor']],
            target=[index[g] for g in subset['gene']],
            value=subset['hit_count'].tolist(),
            color=link_colors,
        ),
    ))
    figure.update_layout(
        width=width, height=height,
        paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
        font=dict(family='Arial', size=font_size, color='black'),
    )
    return figure


def save_sankey(figure, out_path: Path | str, write_svg: bool = True) -> None:
    """Write the interactive HTML, and an SVG when kaleido is available."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # plotly names the container div with a fresh uuid on every write, so
    # without div_id the same figure produces a different file each run and
    # the output tree cannot be diffed. See errors.md §23.
    figure.write_html(f"{out_path}.html", div_id=out_path.name)
    if not write_svg:
        return
    try:
        svg = figure.to_image(format='svg').decode('utf-8')
        # plotly suffixes the two <defs> ids with six random hex digits on
        # every export, which is the only thing that differs between two
        # renders of the same figure. Naming them keeps the SVG reproducible.
        svg = re.sub(r'(id="(?:top)?defs-)[0-9a-f]{6}"', r'\g<1>sankey"', svg)
        Path(f"{out_path}.svg").write_text(svg)
    except Exception as exc:                                  # kaleido missing
        print(f"  SVG export skipped ({type(exc).__name__}); HTML written")
