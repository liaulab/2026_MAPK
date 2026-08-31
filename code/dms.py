"""Deep mutational scanning reference data, mapped onto screen guides.

Three published DMS datasets are used, one per gene:

PTPN11  Jiang et al., Nature Communications 2025 (position x amino-acid matrix)
KRAS    integrated KRAS dataset, HCC827 G12D log fold change
MAPK1   Brenan et al., Cell Reports 2016 (ERK2, ETP vs DOX log fold change)

All three score depletion negatively. A guide making several scorable edits is
summarized by ``config.DMS_AGGREGATE`` across them.

All three are read from ``Inputs-DMS.xlsx``, one sheet of Gene/Variant/Score
rows built by ``build_inputs`` from the published files under
``Inputs/DMS-Data``; each dataset arrives in a different shape and that module
is where the three shapes are understood.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from . import config as cfg


def load_dms_table(path=None) -> pd.DataFrame:
    """The stacked DMS sheet: Gene, Variant, Score, Dataset."""
    return pd.read_excel(path or cfg.DMS_XLSX, sheet_name=cfg.DMS_SHEET)


def load_all_dms(path=None) -> dict[str, dict[str, float]]:
    """DMS lookup per gene symbol.

    ``Inputs-DMS.xlsx`` holds the three published datasets stacked into one
    sheet, already reduced to the variant/score pairs the scoring uses; the
    per-dataset readers that produced it live in ``build_inputs``.
    """
    table = load_dms_table(path)
    return {gene: dict(zip(group['Variant'], group['Score']))
            for gene, group in table.groupby('Gene', sort=False)}


def dms_summary_stats(dms_dicts: dict[str, dict[str, float]]) -> pd.DataFrame:
    """Mean and standard deviation of each DMS dataset, for Z-scoring.

    Statistics come from the whole published dataset, not from the subset of
    variants the library happens to make, so a guide's Z-score is relative to
    the full mutational landscape.
    """
    rows = []
    for gene, scores in dms_dicts.items():
        values = [v for v in scores.values()
                  if v is not None and not (isinstance(v, float) and math.isnan(v))]
        rows.append({'Gene': gene, 'n_variants': len(values),
                     'mean': float(np.mean(values)), 'stdev': float(np.std(values))})
    return pd.DataFrame(rows).set_index('Gene')


AGGREGATORS = {
    'min': min,
    'mean': lambda values: float(np.mean(values)),
}


def score_guide(mutations: str, dms_dict: dict[str, float],
                aggregate: str | None = None) -> float | None:
    """DMS score for a guide, aggregated across the coding edits it makes.

    ``mutations`` is the semicolon/slash separated edit string from the
    library annotation. Silent edits (same start and end amino acid) carry no
    DMS score and are skipped; a guide with no scorable edit returns None
    rather than 0, so 'not measured' stays distinguishable from 'measured as
    zero'.

    ``aggregate`` defaults to ``config.DMS_AGGREGATE``. Roughly a third of
    scored guides make more than one scorable edit, so which aggregator is
    used changes the score materially; see errors.md §24.
    """
    if not isinstance(mutations, str) or not mutations:
        return None

    edits = {m for m in mutations.replace('/', ';').split(';') if m}
    scores = [dms_dict[m] for m in edits
              if m and m[0] != m[-1] and m in dms_dict]
    scores = [s for s in scores if s is not None and not (isinstance(s, float) and math.isnan(s))]
    if not scores:
        return None
    return AGGREGATORS[aggregate or cfg.DMS_AGGREGATE](scores)


def add_dms_scores(df: pd.DataFrame, dms_dicts: dict[str, dict[str, float]],
                   stats: pd.DataFrame, mutation_col: str = 'mutations',
                   gene_col: str = 'Gene', aggregate: str | None = None) -> pd.DataFrame:
    """Add per-gene and combined DMS score columns to a screen frame.

    Each guide is scored only against its own gene's dataset. The combined
    ``DMS-Score`` / ``DMS-Score-Z`` columns are NaN where no dataset applies,
    which is why they are built by coalescing rather than summing: summing
    turned an unscored guide into a real zero in the original code.
    """
    df = df.copy()

    raw_cols, z_cols = [], []
    for gene, dms_dict in dms_dicts.items():
        raw_col, z_col = f'{gene}-DMS-Score', f'{gene}-DMS-Score-Z'
        matches = df[gene_col] == gene
        df[raw_col] = np.where(
            matches,
            df[mutation_col].map(lambda m: score_guide(m, dms_dict, aggregate)),
            np.nan,
        )
        df[raw_col] = pd.to_numeric(df[raw_col], errors='coerce')
        df[z_col] = (df[raw_col] - stats.loc[gene, 'mean']) / stats.loc[gene, 'stdev']
        raw_cols.append(raw_col)
        z_cols.append(z_col)

    # A guide belongs to exactly one gene, so coalescing is unambiguous.
    df['DMS-Score'] = df[raw_cols].bfill(axis=1).iloc[:, 0]
    df['DMS-Score-Z'] = df[z_cols].bfill(axis=1).iloc[:, 0]
    return df


def prepare_dms_frame(screen_long: pd.DataFrame,
                      dms_dicts: dict[str, dict[str, float]],
                      stats: pd.DataFrame,
                      genes: list[str] | None = None) -> pd.DataFrame:
    """Screen frame restricted to the DMS genes and annotated with DMS scores."""
    genes = genes or cfg.DMS_GENES
    subset = screen_long[screen_long['Gene'].isin(genes)]
    subset = subset[subset['mutations'].notna()]
    return add_dms_scores(subset, dms_dicts, stats)
