"""Loading and reshaping of the MAPK screen data.

Everything reads from ``TableS1-ScreenData.xlsx``. The ten original notebooks
each repeated the same load-and-clean preamble; it lives here once.

Editor-agnostic columns
-----------------------
The ABE sheets annotate edits in ``AtoG_*`` columns and the CBE sheets in
``CtoT_*``. Rather than making every downstream figure branch on the editor,
:func:`load_screen` adds three normalized columns to every frame:

``Editor``      'ABE' or 'CBE'
``pos``         amino-acid position of the edit (``AtoG_pos`` / ``CtoT_pos``)
``muttype``     ``AtoG_muttype`` / ``CtoT_muttype``
``mutations``   ``AtoG_mutations`` / ``CtoT_mutations``

The original notebooks used a guide-level ``edit_site`` column for some x-axes.
It is not in TableS1, so ``pos`` is used throughout instead; guides with no
coding edit carry ``pos`` of -1 or NaN and are dropped by
:func:`drop_unpositioned` wherever a position axis is required.
"""

from __future__ import annotations

import pandas as pd

from . import config as cfg

# Columns Excel round-trips in that carry no data.
_JUNK_COLUMNS = ['Unnamed: 0', '-']

_EDITOR_COLUMNS = {
    'ABE': {'pos': 'AtoG_pos', 'muttype': 'AtoG_muttype',
            'mutations': 'AtoG_mutations', 'muttypes': 'AtoG_muttypes'},
    'CBE': {'pos': 'CtoT_pos', 'muttype': 'CtoT_muttype',
            'mutations': 'CtoT_mutations', 'muttypes': 'CtoT_muttypes'},
}

ANNOTATION_COLUMNS = ['pos', 'muttype', 'mutations', 'muttypes']


def _read_sheet(sheet_name: str) -> pd.DataFrame:
    df = pd.read_excel(cfg.SCREEN_DATA_XLSX, sheet_name=sheet_name)
    return df.drop(columns=[c for c in _JUNK_COLUMNS if c in df.columns])


def deduplicate(df: pd.DataFrame) -> pd.DataFrame:
    """Drop every row of any (Gene, sgRNA_seq) group seen more than once.

    A handful of guides appear twice per sheet with different annotations, so
    neither copy can be trusted; the original notebooks discarded both, and
    that behaviour is preserved.
    """
    return df[~df.duplicated(subset=['Gene', 'sgRNA_seq'], keep=False)].copy()


def add_editor_columns(df: pd.DataFrame, editor: str) -> pd.DataFrame:
    """Add Editor plus the editor-agnostic pos/muttype/mutations columns."""
    df = df.copy()
    df['Editor'] = editor
    for target, source in _EDITOR_COLUMNS[editor].items():
        df[target] = df[source] if source in df.columns else pd.NA
    return df


def load_screen(screen: str, editor: str, dedup: bool = True) -> pd.DataFrame:
    """One screen/editor combination, cleaned and annotated.

    screen is one of 'gof', 'lof', 'meki', 'valid', 'lof_erk'.
    """
    df = _read_sheet(cfg.SCREEN_SHEETS[(screen, editor)])
    if dedup:
        df = deduplicate(df)
    return add_editor_columns(df, editor)


def load_all_screens(dedup: bool = True) -> dict[tuple[str, str], pd.DataFrame]:
    """Every screen/editor frame, keyed by (screen, editor)."""
    return {(screen, editor): load_screen(screen, editor, dedup=dedup)
            for screen, editor in cfg.SCREEN_SHEETS}


def load_edit_site_key() -> pd.DataFrame:
    """Guide-level editing-window position, keyed by (Gene, sgRNA_seq).

    Used only to order rows before fitting the GMM trajectory. scikit-learn's
    k-means initialization depends on row order, and the 14-component fit on
    the validation screen is multimodal enough that a different ordering lands
    on a different local optimum. Ordering by this column reproduces the
    published clustering exactly.

    Read from the ``edit-site`` sheet of TableS1; it was previously a separate
    CSV under ``Inputs/``.
    """
    return _read_sheet(cfg.EDIT_SITE_SHEET)


def order_by_edit_site(df: pd.DataFrame) -> pd.DataFrame:
    """Sort by gene then editing-window position, for a reproducible GMM fit.

    Guides with no recorded editing-window position sort last, consistently.
    """
    key = load_edit_site_key()
    merged = df.merge(key, how='left', on=['Gene', 'sgRNA_seq'])
    return merged.sort_values(['Gene', 'edit_site'], na_position='last').reset_index(drop=True)


def drop_unpositioned(df: pd.DataFrame, pos_col: str = 'pos') -> pd.DataFrame:
    """Rows with a real amino-acid position.

    Guides whose editing window produces no coding change carry a position of
    -1 or NaN. Any figure with a position x-axis must exclude them.
    """
    return df[df[pos_col].notna() & (df[pos_col] > 0)].copy()


def combine_editors(abe: pd.DataFrame, cbe: pd.DataFrame,
                    gene: str | None = None,
                    keep_muttypes: bool = True) -> pd.DataFrame:
    """Stack the ABE and CBE frames, optionally for a single gene.

    Restricting to the mutation types in ``MUT_PAL`` drops UTR-annotated guides,
    matching the original scatterplot filtering.
    """
    if gene is not None:
        abe = abe[abe['Gene'] == gene]
        cbe = cbe[cbe['Gene'] == gene]
    combined = pd.concat([abe, cbe], ignore_index=True)
    if keep_muttypes:
        combined = combined[combined['muttype'].isin(cfg.MUT_PAL)]
    return combined.copy()


def filter_day0(df: pd.DataFrame, screen: str, editor: str,
                cutoff: int = cfg.DAY0_COUNT_CUTOFF) -> pd.DataFrame:
    """Guides represented above ``cutoff`` reads in the Day 0 sample."""
    column = cfg.DAY0_COUNT_COL[screen]
    if isinstance(column, dict):
        column = column[editor]
    return df[df[column] > cutoff].copy()


def drop_controls(df: pd.DataFrame) -> pd.DataFrame:
    """Remove non-targeting and essential-splice-site control guides."""
    return df[~df['Gene'].isin(cfg.CONTROL_GENES)].copy()


def splice_guides(df: pd.DataFrame) -> pd.DataFrame:
    """The splice-targeting validation guides, identified by sgRNA_ID."""
    return df[df['sgRNA_ID'].str.contains('splice_pool', na=False)].copy()


def to_protein_names(df: pd.DataFrame, gene_col: str = 'Gene') -> pd.DataFrame:
    """Relabel gene symbols with the protein names used in the figures."""
    df = df.copy()
    df[gene_col] = df[gene_col].map(lambda g: cfg.GENE_PROTEIN_MAP.get(g, g))
    return df


# ---------------------------------------------------------------------------
# DERIVED FRAMES
# ---------------------------------------------------------------------------


def rename_validation_abe(df: pd.DataFrame) -> pd.DataFrame:
    """Put Validation-ABE Day18 columns onto the shared condition names.

    Validation-ABE is a Day18/Day32 timecourse while Validation-CBE has a
    single timepoint. Renaming the Day18 columns lets both editors share one
    set of condition names. Day32 columns are left untouched.
    """
    return df.rename(columns=cfg.VALID_ABE_TO_COMMON)


def annotate_validation(valid: pd.DataFrame, gof: pd.DataFrame) -> pd.DataFrame:
    """Attach GOF mutation annotations to a validation frame.

    Must be called per editor: the GOF-ABE sheet carries only ``AtoG_*``
    columns and GOF-CBE only ``CtoT_*``, so a validation frame can only be
    annotated from the GOF sheet of its own editor.
    """
    editor = valid['Editor'].iloc[0]
    source_cols = ['Gene', 'sgRNA_seq'] + [
        c for c in _EDITOR_COLUMNS[editor].values() if c in gof.columns
    ]
    merged = valid.drop(columns=ANNOTATION_COLUMNS, errors='ignore').merge(
        gof[source_cols], how='left', on=['Gene', 'sgRNA_seq'],
    )
    return add_editor_columns(merged, editor)


def merge_gof_validation(gof: pd.DataFrame, valid: pd.DataFrame,
                         gof_columns: list[str],
                         valid_columns: list[str],
                         prefix: str = 'Valid-') -> pd.DataFrame:
    """Validation frame with the matching GOF scores joined on.

    Validation is the left frame so every validated guide is kept even when it
    has no GOF counterpart, which is how the hit-recovery barplots count.
    """
    renamed = dict(zip(valid_columns, [f'{prefix}{c}' for c in gof_columns]))
    valid_slim = valid[['sgRNA_seq', 'Gene'] + valid_columns].rename(columns=renamed)
    gof_slim = gof[['sgRNA_seq', 'Gene'] + gof_columns]
    return valid_slim.merge(gof_slim, how='left', on=['sgRNA_seq', 'Gene'])


def merge_lof(base: pd.DataFrame, lof: pd.DataFrame,
              column: str = 'Dox-UT-Z', fillna: float | None = None) -> pd.DataFrame:
    """Join the LOF hyperactivation Z-score onto another screen's frame."""
    merged = base.merge(
        lof[['Gene', 'sgRNA_seq', column]], how='left', on=['Gene', 'sgRNA_seq'],
    )
    if fillna is not None:
        merged[column] = merged[column].fillna(fillna)
    return merged


def merge_gof_into_meki(meki: pd.DataFrame, gof: pd.DataFrame) -> pd.DataFrame:
    """Join the GOF conditions onto the MEKi frame for the MEKi trajectory.

    DMSO-Day0 and Trametinib exist in both screens, so the GOF copies are
    prefixed rather than overwriting the MEKi values.
    """
    renames = {'DMSO-Day0-Z': 'GOF-DMSO-Day0-Z', 'Trametinib-DMSO-Z': 'GOF-Trametinib-DMSO-Z'}
    gof_renamed = gof.rename(columns=renames)
    columns = ['Gene', 'sgRNA_seq'] + list(renames.values()) + [
        'SHP2i-DMSO-Z', 'KRASi-DMSO-Z', 'RAFi-DMSO-Z', 'ERKi-DMSO-Z',
    ]
    return meki.merge(gof_renamed[columns], how='left', on=['Gene', 'sgRNA_seq'])


def build_lof_gof_long(gof_abe: pd.DataFrame, gof_cbe: pd.DataFrame,
                       lof_abe: pd.DataFrame, lof_cbe: pd.DataFrame) -> pd.DataFrame:
    """Long-format GOF frame carrying the LOF hyperactivation score.

    This reconstructs what the original notebooks read from the pre-merged
    ``Merged_LOF_data.csv``: each GOF sheet left-joined with its editor's LOF
    ``Dox-UT-Z``, then stacked with an ``Editor`` column. Guides absent from
    the LOF library (MAPK1, MAPK3) keep a NaN score.
    """
    frames = [
        merge_lof(gof_abe, lof_abe),
        merge_lof(gof_cbe, lof_cbe),
    ]
    return pd.concat(frames, ignore_index=True)
