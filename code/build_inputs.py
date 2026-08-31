"""Build the three consolidated Excel inputs from the raw files under ``Inputs/``.

The analysis reads three workbooks in the repository root rather than the
directory trees they were built from:

``Inputs-cBioPortal.xlsx``       one sheet per gene, plus a transcript summary
``Inputs-PocketDistances.xlsx``  one sheet per structure, plus a summary
``Inputs-DMS.xlsx``              one sheet, all three DMS datasets stacked

Consolidating them means a reviewer downloads three files instead of sixty, and
the expensive step — parsing eighteen mmCIF files with Biopython to measure
every residue's distance to its ligand — runs once here rather than on every
notebook execution. Nothing in the pipeline reads ``Inputs/`` any more.

This module is the only place the raw formats are understood. It is a build
step, not part of a run: execute it when a raw input changes.

    python -m code.build_inputs

The readers below are the ones the analysis modules used to carry, moved here
unchanged, so the workbooks hold exactly what the previous code loaded.
"""

from __future__ import annotations

import re

import pandas as pd

from . import config as cfg
from . import structures as struct

# xlsxwriter would otherwise turn a cell beginning with '=' into a formula and
# anything URL-shaped into a hyperlink; every value here is literal text.
_WRITER_KWARGS = dict(
    engine='xlsxwriter',
    engine_kwargs={'options': {'strings_to_formulas': False,
                               'strings_to_urls': False}},
)


# ---------------------------------------------------------------------------
# cBIOPORTAL
# ---------------------------------------------------------------------------


def transcript_ids(filename: str) -> list[str]:
    """Every ENST accession in a cBioPortal export's filename.

    Two exports were downloaded across a pair of transcripts and are named for
    both; the accessions are in the filename only, never in the table.
    """
    return re.findall(r'ENST\d+', filename)


def build_cbioportal(path=None) -> pd.DataFrame:
    """Write the per-gene workbook; return the transcript summary sheet.

    Sheet order follows ``config.CBIOPORTAL_FILES``, which is pathway order.
    Where an export covers two transcripts the first accession is the one
    recorded as canonical, matching how the files were requested.
    """
    path = path or cfg.CBIOPORTAL_XLSX

    tables, rows = {}, []
    for gene, filename in cfg.CBIOPORTAL_FILES.items():
        table = pd.read_csv(cfg.CBIOPORTAL_DIR / filename, sep='\t')
        tables[gene] = table
        accessions = transcript_ids(filename)
        rows.append({
            'Gene': gene,
            'Transcript': accessions[0],
            'All transcripts': ';'.join(accessions),
            'Source file': filename,
            'Mutations': len(table),
            'Samples': table['Sample ID'].nunique(),
        })
    summary = pd.DataFrame(rows)

    with pd.ExcelWriter(path, **_WRITER_KWARGS) as writer:
        summary.to_excel(writer, sheet_name=cfg.CBIOPORTAL_TRANSCRIPT_SHEET,
                         index=False)
        for gene, table in tables.items():
            table.to_excel(writer, sheet_name=gene, index=False)
    return summary


# ---------------------------------------------------------------------------
# POCKET DISTANCES
# ---------------------------------------------------------------------------

# Every structure in Inputs/Pocket-Structures, tagged by the role it plays.
# 9O0R is a second KRAS/adagrasib co-crystal that no figure uses; it is
# measured and stored anyway so the workbook covers the whole directory.
STRUCTURE_SOURCES = (
    ('Co-crystal', cfg.PDB_STRUCTURES),
    ('AlphaFold', cfg.AF_STRUCTURES),
    ('Unused', cfg.EXTRA_STRUCTURES),
)


def build_pocket_distances(path=None, csv_dir=None) -> pd.DataFrame:
    """Measure every structure, write the CSVs and the workbook.

    Each structure's table is stored exactly as
    :func:`structures.residue_ligand_distances` produced it — sorted by
    distance, index reset — so reading a sheet back is interchangeable with
    re-parsing the mmCIF.
    """
    path = path or cfg.POCKET_DISTANCES_XLSX
    csv_dir = csv_dir or cfg.POCKET_DISTANCES_DIR
    csv_dir.mkdir(parents=True, exist_ok=True)

    tables, rows = {}, []
    for source, specs in STRUCTURE_SOURCES:
        for gene, structure, inhibitor, lig_name, lig_chain, prot_chain, condition in specs:
            table = struct.residue_ligand_distances(
                cfg.STRUCTURES_DIR / f'{structure}.cif',
                lig_name, lig_chain, prot_chain,
            )
            table.to_csv(csv_dir / f'{structure}-distances.csv', index=False)
            tables[structure] = table
            rows.append({
                'Gene': gene,
                'Structure': structure,
                'Source': source,
                'Inhibitor': inhibitor,
                'Ligand': lig_name,
                'Ligand chain': lig_chain,
                'Protein chain': prot_chain,
                'Condition': condition,
                'Residues': len(table),
                'Min distance': table['Distance'].min(),
                'Max distance': table['Distance'].max(),
                'Sequence': table.attrs.get('sequence', ''),
            })
    summary = pd.DataFrame(rows)

    with pd.ExcelWriter(path, **_WRITER_KWARGS) as writer:
        summary.to_excel(writer, sheet_name=cfg.POCKET_STRUCTURE_SHEET, index=False)
        for structure, table in tables.items():
            table.to_excel(writer, sheet_name=structure, index=False)
    return summary


# ---------------------------------------------------------------------------
# DMS DATASETS
# ---------------------------------------------------------------------------


def read_erk2_dms() -> dict[str, float]:
    """ERK2 mutant -> log fold change."""
    key_col, score_col = cfg.DMS_SCORE_COLUMNS['MAPK1']
    df = pd.read_csv(cfg.DMS_FILES['MAPK1'])
    return dict(zip(df[key_col], df[score_col]))


def read_kras_dms() -> dict[str, float]:
    """KRAS variant -> HCC827 G12D log fold change."""
    key_col, score_col = cfg.DMS_SCORE_COLUMNS['KRAS']
    df = pd.read_csv(cfg.DMS_FILES['KRAS'])
    return dict(zip(df[key_col], df[score_col]))


def read_shp2_dms() -> dict[str, float]:
    """SHP2 variant -> score, flattened from a position x amino-acid matrix.

    The file is indexed by mutant amino acid with one column per residue
    position, and its first row holds the wild-type amino acid. Combining
    those gives the usual 'M1K' variant key.
    """
    df = pd.read_csv(cfg.DMS_FILES['PTPN11'], index_col=0)
    df.columns = df.iloc[0] + df.columns          # wild-type AA + position
    df = df.drop(index=df.index[0])               # drop the wild-type row

    scores = {}
    for column in df.columns:
        for mutant_aa, value in df[column].items():
            scores[f"{column}{mutant_aa}"] = float(value)
    return scores


RAW_DMS_READERS = {
    'MAPK1': read_erk2_dms,
    'KRAS': read_kras_dms,
    'PTPN11': read_shp2_dms,
}


def build_dms(path=None) -> pd.DataFrame:
    """Write the stacked DMS sheet; return it.

    One row per (gene, variant), in the order the source files list them, so
    that rebuilding the per-gene lookup preserves the original last-one-wins
    resolution of the handful of variants a dataset reports twice.
    """
    path = path or cfg.DMS_XLSX

    frames = []
    for gene, reader in RAW_DMS_READERS.items():
        scores = reader()
        frames.append(pd.DataFrame({
            'Gene': gene,
            'Variant': list(scores.keys()),
            'Score': list(scores.values()),
            'Dataset': cfg.DMS_DATASETS[gene],
        }))
    table = pd.concat(frames, ignore_index=True)

    with pd.ExcelWriter(path, **_WRITER_KWARGS) as writer:
        table.to_excel(writer, sheet_name=cfg.DMS_SHEET, index=False)
    return table


# ---------------------------------------------------------------------------


def build_all() -> dict[str, pd.DataFrame]:
    """Rebuild all three workbooks."""
    return {
        'cBioPortal': build_cbioportal(),
        'PocketDistances': build_pocket_distances(),
        'DMS': build_dms(),
    }


if __name__ == '__main__':
    for name, summary in build_all().items():
        print(f'{name}: {len(summary)} rows')
