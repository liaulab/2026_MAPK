"""Analysis steps shared by the static and interactive notebooks.

Everything here computes; nothing draws. Both `MAPK_Master_Analysis.ipynb`
(print figures) and `MAPK_Interactive_Analysis.ipynb` (HTML figures) import
these, so the two notebooks are guaranteed to be plotting the same numbers.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from . import config as cfg
from . import data as dat
from . import structures as struct

# ---------------------------------------------------------------------------
# SCREEN LOADING
# ---------------------------------------------------------------------------


def load_screens() -> dict:
    """Every screen frame, with the standard renames applied.

    Returns a dict keyed by short name, plus a `table` summarising what was
    loaded.
    """
    screens = dat.load_all_screens()

    frames = {
        'gof_abe': screens[('gof', 'ABE')], 'gof_cbe': screens[('gof', 'CBE')],
        'lof_abe': screens[('lof', 'ABE')], 'lof_cbe': screens[('lof', 'CBE')],
        'meki_abe': screens[('meki', 'ABE')], 'meki_cbe': screens[('meki', 'CBE')],
        'valid_abe': dat.rename_validation_abe(screens[('valid', 'ABE')]),
        'valid_cbe': screens[('valid', 'CBE')],
        'lof_erk_abe': screens[('lof_erk', 'ABE')].rename(
            columns={'Day14-ABE-Dox-UT-Z': cfg.LOF_ERK_COLUMN}),
        'lof_erk_cbe': screens[('lof_erk', 'CBE')].rename(
            columns={'Day14-CBE-Dox-UT-Z': cfg.LOF_ERK_COLUMN}),
    }

    frames['table'] = pd.DataFrame(
        [{'screen': s, 'editor': e, 'guides': len(df), 'genes': df['Gene'].nunique(),
          'conditions': sum(c.endswith('-Z') for c in df.columns)}
         for (s, e), df in screens.items()])
    return frames


# ---------------------------------------------------------------------------
# SCATTERPLOT INPUTS
# ---------------------------------------------------------------------------


def gene_panels(abe, cbe, ranges, conditions, genes=None):
    """Yield (gene, gene_range, subset, present_conditions) for each panel."""
    genes = genes or cfg.GENES
    ranges = ranges if isinstance(ranges, dict) else dict(zip(cfg.GENES, ranges))

    for gene in genes:
        gene_range = ranges.get(gene)
        if not gene_range:
            continue
        subset = dat.combine_editors(abe, cbe, gene=gene)
        if subset.empty:
            continue
        present = [c for c in conditions if c in subset.columns]
        if not present:
            continue
        yield gene, gene_range, subset, present


def melt_conditions(df, conditions, x_col='pos'):
    """Long form: one row per guide and condition, score in a 'score' column."""
    frames = []
    for condition in conditions:
        slim = df[['Gene', 'Editor', 'muttype', x_col, condition]].copy()
        if 'sgRNA_ID' in df.columns:
            slim['sgRNA_ID'] = df['sgRNA_ID'].to_numpy()
        if 'mutations' in df.columns:
            slim['mutations'] = df['mutations'].to_numpy()
        slim = slim.rename(columns={condition: 'score'})
        slim['Condition'] = condition
        frames.append(slim)
    return pd.concat(frames, ignore_index=True)


def combined_panel_frame(frames, gene, random_state=cfg.RANDOM_STATE):
    """All-conditions long frame for one gene, shuffled for even overplotting.

    ``frames`` is a list of (abe, cbe, conditions) triples, so the LOF figure
    can draw its DMSO column from the GOF screen alongside the LOF score.
    """
    pieces = []
    for abe, cbe, conditions in frames:
        subset = dat.combine_editors(abe, cbe, gene=gene)
        present = [c for c in conditions if c in subset.columns]
        if subset.empty or not present:
            continue
        pieces.append(melt_conditions(subset, present))
    if not pieces:
        return None
    plot_df = pd.concat(pieces, ignore_index=True)
    return plot_df.sample(frac=1, random_state=random_state)


# ---------------------------------------------------------------------------
# GOF VERSUS VALIDATION
# ---------------------------------------------------------------------------


def gof_validation_frames(gof_abe, gof_cbe, valid_abe, valid_cbe):
    """Validation frames with the matching GOF scores joined on, per editor."""
    return (
        dat.merge_gof_validation(gof_abe, valid_abe, cfg.GOF_SCREEN_COLUMNS,
                                 cfg.VALID_CBE_SCREEN_COLUMNS),
        dat.merge_gof_validation(gof_cbe, valid_cbe, cfg.GOF_SCREEN_COLUMNS,
                                 cfg.VALID_CBE_SCREEN_COLUMNS),
    )


GOF_CONDITIONS = cfg.COMBINED_GOF_COLUMNS
VALID_CONDITIONS = [f'Valid-{c}' for c in GOF_CONDITIONS]


def hit_recovery(frame, editor):
    """Hit and validated-hit counts per condition, in both directions."""
    rows = []
    for gof_col, valid_col in zip(GOF_CONDITIONS, VALID_CONDITIONS):
        pos = frame[gof_col] > cfg.GOF_HIT_CUTOFF
        neg = frame[gof_col] < -cfg.GOF_HIT_CUTOFF
        pos_valid = pos & (frame[valid_col] > cfg.VALID_HIT_CUTOFF)
        neg_valid = neg & (frame[valid_col] < -cfg.VALID_HIT_CUTOFF)
        rows.append({
            'editor': editor, 'condition': gof_col,
            'pos_hits': int(pos.sum()), 'pos_validated': int(pos_valid.sum()),
            'neg_hits': int(neg.sum()), 'neg_validated': int(neg_valid.sum()),
            'pos_rate': pos_valid.sum() / pos.sum() if pos.sum() else np.nan,
            'neg_rate': neg_valid.sum() / neg.sum() if neg.sum() else np.nan,
        })
    return pd.DataFrame(rows)


def correlation_frame(frame, editor):
    """Hits and controls in long form, with the per-condition Pearson r."""
    pieces, stats = [], []
    for gof_col, valid_col in zip(GOF_CONDITIONS, VALID_CONDITIONS):
        keep = ((frame[gof_col] > cfg.GOF_HIT_CUTOFF)
                | (frame[gof_col] < -cfg.GOF_HIT_CUTOFF)
                | (frame['Gene'] == 'Control'))
        piece = frame.loc[keep, ['Gene', gof_col, valid_col]].rename(
            columns={gof_col: 'GOF Score', valid_col: 'Valid Score'})
        piece['Condition'] = gof_col
        pieces.append(piece)
        stats.append({
            'editor': editor, 'condition': gof_col, 'n': len(piece),
            'pearson_r': piece['GOF Score'].corr(piece['Valid Score'], method='pearson'),
        })
    return pd.concat(pieces, ignore_index=True), pd.DataFrame(stats)


# ---------------------------------------------------------------------------
# POCKET ANALYSIS
# ---------------------------------------------------------------------------


def node_hits(abe, cbe, gene, condition, cutoff=cfg.POCKET_HIT_CUTOFF):
    """Resistance hits at one node, both editors, with a real position."""
    subset = dat.combine_editors(abe, cbe, gene=gene, keep_muttypes=False)
    subset = subset[subset[condition] > cutoff]
    return dat.drop_unpositioned(subset)


def off_node_hits(abe, cbe, gene, condition, paralogs, cutoff=cfg.POCKET_HIT_CUTOFF):
    """Hits in the same condition that fall outside the inhibited node."""
    counts = {'paralog': 0, 'other': 0}
    for frame in (abe, cbe):
        hits = dat.drop_unpositioned(
            frame[(frame['Gene'] != gene) & (frame[condition] > cutoff)])
        counts['paralog'] += int(hits['Gene'].isin(paralogs).sum())
        counts['other'] += int((~hits['Gene'].isin(paralogs)).sum())
    return counts


def pocket_analysis(structure_specs, paralog_map, distance_tables, abe, cbe):
    """Classify every node's hits as in-pocket, outside-pocket or unresolved."""
    records = []
    for gene, structure, inhibitor, *_rest, condition in structure_specs:
        distances = distance_tables[structure]
        record = {
            'gene': gene, 'structure': structure, 'inhibitor': inhibitor,
            'condition': condition, 'label': f'{gene} {structure}',
            'residues_in_pocket': int((distances['Distance'] <= cfg.POCKET_DISTANCE_CUTOFF).sum()),
            'residues_outside': int((distances['Distance'] > cfg.POCKET_DISTANCE_CUTOFF).sum()),
            'distances': distances,
            'off_node': off_node_hits(abe, cbe, gene, condition, paralog_map[gene]),
        }
        for editor, pair in (('ABE', (abe, cbe.iloc[0:0])), ('CBE', (abe.iloc[0:0], cbe))):
            hits = node_hits(*pair, gene, condition)
            split = struct.split_in_out_pocket(hits, distances, pos_col='pos')
            for key, part in split.items():
                record[f'{editor}_{key}'] = part
                record[f'n_{editor}_{key}'] = len(part)
        records.append(record)
    return records


def pocket_summary(records):
    """Flat table of the pocket classification."""
    return pd.DataFrame([{
        'label': r['label'], 'gene': r['gene'], 'structure': r['structure'],
        'inhibitor': r['inhibitor'], 'condition': r['condition'],
        'residues_in_pocket': r['residues_in_pocket'],
        'residues_outside': r['residues_outside'],
        **{f'n_{e}_{k}': r[f'n_{e}_{k}']
           for e in ('ABE', 'CBE') for k in ('InPocket', 'OutPocket', 'Unresolved')},
        'n_paralog_hits': r['off_node']['paralog'],
        'n_other_gene_hits': r['off_node']['other'],
    } for r in records])


def rearrange(values):
    """Interleave the two halves of a list, pairing each structure's editors."""
    mid = len(values) // 2
    return [x for pair in zip(values[:mid], values[mid:]) for x in pair]


def stacked_series(records, keys):
    """One interleaved ABE/CBE series per pocket class."""
    return [rearrange([r[f'n_ABE_{k}'] for r in records]
                      + [r[f'n_CBE_{k}'] for r in records]) for k in keys]


def structure_labels(records):
    return rearrange([f"{r['label']} ABE" for r in records]
                     + [f"{r['label']} CBE" for r in records])


def outside_node_series(records):
    return rearrange([r['off_node']['paralog'] + r['off_node']['other']
                      for r in records] * 2)


def hit_distance_frame(record):
    """Ligand distance for each hit, by editor."""
    lookup = dict(zip(record['distances']['Position'], record['distances']['Distance']))
    rows = []
    for editor in ('ABE', 'CBE'):
        for key in ('InPocket', 'OutPocket', 'Unresolved'):
            for pos in record[f'{editor}_{key}']['pos']:
                if int(pos) in lookup:
                    rows.append({'distance': lookup[int(pos)], 'editor': editor,
                                 'pocket': key})
    return pd.DataFrame(rows, columns=['distance', 'editor', 'pocket'])


def hit_score_frame(record, keys):
    """Hit Z-scores labelled by editor and pocket class."""
    rows = []
    for editor in ('ABE', 'CBE'):
        for key in keys:
            part = record[f'{editor}_{key}']
            for value, gene in zip(part[record['condition']], part['Gene']):
                rows.append({'Editor': editor, 'pocket': key,
                             'Z-Score': value, 'Gene': gene})
    return pd.DataFrame(rows, columns=['Editor', 'pocket', 'Z-Score', 'Gene'])


PARALOG_FAMILIES = {
    'AllRAS': ['KRAS', 'NRAS', 'HRAS', 'MRAS'],
    'AllRAF': ['ARAF', 'BRAF', 'RAF1'],
    'AllMEK': ['MAP2K1', 'MAP2K2'],
    'AllERK': ['MAPK1', 'MAPK3'],
}


# ---------------------------------------------------------------------------
# TRAJECTORY INPUTS
# ---------------------------------------------------------------------------


def prepare_trajectory_input(abe, cbe, conditions, lower, upper,
                             clip_columns=None, two_sided=True):
    """Guides responding in at least one condition, clipped to the plot range.

    ``two_sided`` keeps guides moving in either direction. The MEKi trajectory
    is one-sided: it maps resistance only, so a depleted guide under MEK
    inhibition is not informative about the landscape being described.
    """
    clip_cols = clip_columns or conditions
    pieces = []
    for frame in (abe, cbe):
        scores = np.abs(frame[conditions]) if two_sided else frame[conditions]
        subset = frame[scores.gt(lower).any(axis=1)].copy()
        subset[clip_cols] = subset[clip_cols].clip(lower=-upper, upper=upper)
        pieces.append(subset)
    return pd.concat(pieces, ignore_index=True)


def validation_trajectory_input(valid_abe, valid_cbe, gof_abe, gof_cbe):
    """Guides entering the validation trajectory, ordered reproducibly.

    Day-0 representation and control filtering are applied first, annotations
    are taken from the GOF screen under the same filtering, and rows are
    ordered by editing-window position. That ordering matters: k-means
    initialization is order-dependent and this fit is multimodal.
    """
    abe = dat.drop_controls(dat.filter_day0(valid_abe, 'valid', 'ABE'))
    cbe = dat.drop_controls(dat.filter_day0(valid_cbe, 'valid', 'CBE'))
    gof_a = dat.drop_controls(dat.filter_day0(gof_abe, 'gof', 'ABE'))
    gof_c = dat.drop_controls(dat.filter_day0(gof_cbe, 'gof', 'CBE'))

    abe = dat.annotate_validation(abe, gof_a)
    cbe = dat.annotate_validation(cbe, gof_c)

    combined = prepare_trajectory_input(
        abe, cbe, cfg.TRAJECTORY_VALID_CONDITIONS,
        cfg.TRAJECTORY_VALID_LOWER_CUTOFF, cfg.TRAJECTORY_VALID_UPPER_CLIP)
    return dat.order_by_edit_site(combined)


def meki_trajectory_input(meki_abe, meki_cbe, gof_abe, gof_cbe):
    """Guides entering the MEKi trajectory, with GOF conditions joined on."""
    abe = dat.drop_controls(dat.filter_day0(meki_abe, 'meki', 'ABE'))
    cbe = dat.drop_controls(dat.filter_day0(meki_cbe, 'meki', 'CBE'))
    gof_a = dat.drop_controls(dat.filter_day0(gof_abe, 'gof', 'ABE'))
    gof_c = dat.drop_controls(dat.filter_day0(gof_cbe, 'gof', 'CBE'))

    abe = dat.merge_gof_into_meki(abe, gof_a)
    cbe = dat.merge_gof_into_meki(cbe, gof_c)

    # The response filter runs on the eight MEK inhibitors plus SHP2i.
    filter_columns = cfg.TRAJECTORY_MEKI_CONDITIONS[1:10]
    combined = prepare_trajectory_input(
        abe, cbe, filter_columns,
        cfg.TRAJECTORY_MEKI_LOWER_CUTOFF, cfg.TRAJECTORY_MEKI_UPPER_CLIP,
        clip_columns=cfg.TRAJECTORY_MEKI_CONDITIONS, two_sided=False)
    combined = combined.dropna(subset=cfg.TRAJECTORY_MEKI_CONDITIONS).reset_index(drop=True)
    return combined.rename(columns=dict(zip(cfg.TRAJECTORY_MEKI_CONDITIONS,
                                            cfg.TRAJECTORY_MEKI_LABELS)))


# ---------------------------------------------------------------------------
# PER-GENE BOXPLOTS
# ---------------------------------------------------------------------------


def boxplot_palette(genes):
    """Gene colors plus the two control categories, keyed by category.

    be_scan's boxplot figures look colors up by category name, so every
    category on the axis has to be named. Passing a bare list of gene colors,
    as the original notebook did, left the two control categories to fall off
    the end of the palette and be drawn in recycled gene colors.
    """
    return {**{gene: cfg.COLOR_MAP[gene] for gene in genes},
            **cfg.BOXPLOT_CONTROL_COLORS}


def boxplot_order(genes):
    """Axis order: the screen's genes in pathway order, then the controls."""
    return list(genes) + list(cfg.BOXPLOT_CONTROL_ORDER)


def boxplot_conditions(screen, frame):
    """Condition columns of a screen that are present in a frame."""
    columns = {
        'gof': cfg.GOF_SCREEN_COLUMNS,
        'meki': cfg.MEKI_SCREEN_COLUMNS,
        'lof': cfg.LOF_SCREEN_COLUMNS,
    }[screen]
    return [c for c in columns if c in frame.columns]


def control_boxplot_frame(abe, cbe, column):
    """Controls and essential-splice-site guides from both editors, stacked.

    Returns None when ``column`` is missing from either frame, so a screen
    that has not yet been re-exported with a new condition is skipped rather
    than raising.
    """
    if column not in abe.columns or column not in cbe.columns:
        return None

    pieces = []
    for frame, editor in ((abe, 'ABE'), (cbe, 'CBE')):
        subset = frame.loc[frame['Gene'].isin(cfg.CONTROL_GENES), ['Gene', column]].copy()
        subset['Gene'] = np.where(subset['Gene'] == 'Control',
                                  f'{editor} Control', f'{editor} Essential')
        pieces.append(subset)
    return pd.concat(pieces, ignore_index=True)


def boxplot_summary(frame, screen, editor, condition, genes):
    """Median and spread per category, so the figures have numbers behind them."""
    order = boxplot_order(genes)
    grouped = frame[frame['Gene'].isin(order)].groupby('Gene')[condition]
    table = grouped.agg(n='count', median='median', q1=lambda s: s.quantile(0.25),
                        q3=lambda s: s.quantile(0.75)).reindex(order)
    return table.reset_index().assign(screen=screen, editor=editor,
                                      condition=condition)


def pvalue_stars(p: float) -> str:
    """Significance label for a p-value, or 'ns'."""
    for cutoff, label in cfg.PVALUE_STARS:
        if p < cutoff:
            return label
    return cfg.PVALUE_NS


def control_test(frame, condition, control='Control',
                 essential='Essential splice site', gene_col='Gene'):
    """Non-targeting against essential-splice-site guides, in one condition.

    Two-sided Mann-Whitney U; see ``config.BOXPLOT_CONTROL_TEST`` for why the
    test is rank-based. Returns None when either group is empty, so a screen
    missing a control class is skipped rather than reported as a failed test.
    """
    a = frame.loc[frame[gene_col] == control, condition].dropna()
    b = frame.loc[frame[gene_col] == essential, condition].dropna()
    if a.empty or b.empty:
        return None

    statistic, p = mannwhitneyu(a, b, alternative='two-sided')
    return {'n_control': len(a), 'n_essential': len(b),
            'median_control': float(a.median()),
            'median_essential': float(b.median()),
            'U': float(statistic), 'p': float(p), 'stars': pvalue_stars(p)}


def control_pvalue_conditions(screen, frame):
    """Every column of a screen the control test should cover.

    The per-gene panels' conditions, plus the column the screen's Figure 1
    control panel is drawn from. Those are not the same list: GOF draws its
    control panel from ``DMSO-Day0-Z``, which is one of its screen columns
    anyway, but LOF draws its from ``UT-Day0-Z``, which is not in
    ``LOF_SCREEN_COLUMNS``. Without this the LOF control panel would be drawn
    with a bracket that had no row behind it in the summary table.
    """
    conditions = list(boxplot_conditions(screen, frame))
    for control_screen, column, *_ in cfg.BOXPLOT_CONTROL_FIGURES:
        if control_screen == screen and column in frame.columns \
                and column not in conditions:
            conditions.append(column)
    return conditions


def control_pvalues(box_frames):
    """The control test for every screen, editor and condition.

    ``box_frames`` is keyed by (screen, editor), as the boxplot cells build it.
    One row per panel, so the table lines up with the figures.
    """
    rows = []
    for (screen, editor), frame in box_frames.items():
        for condition in control_pvalue_conditions(screen, frame):
            result = control_test(frame, condition)
            if result is None:
                continue
            rows.append({'screen': cfg.BOXPLOT_SCREEN_LABELS[screen],
                         'editor': editor, 'condition': condition,
                         'test': cfg.BOXPLOT_CONTROL_TEST, **result})
    return pd.DataFrame(rows)


def control_frame_comparisons(frame, column, pairs=None, gene_col='Gene'):
    """Bracket labels for a control panel, one per editor.

    The control panel relabels the two control classes per editor
    ('ABE Control', 'ABE Essential', ...), so each pair is tested within the
    panel's own frame. Returns ``(comparisons, results)``: the first is what
    the figure draws, the second is the numbers behind it.
    """
    comparisons, results = [], []
    for control, essential in (pairs or cfg.BOXPLOT_CONTROL_PAIRS):
        result = control_test(frame, column, control=control,
                              essential=essential, gene_col=gene_col)
        if result is None:
            continue
        comparisons.append((control, essential, result['stars']))
        results.append({'comparison': f'{control} vs {essential}',
                        'condition': column, 'test': cfg.BOXPLOT_CONTROL_TEST,
                        **result})
    return comparisons, results
