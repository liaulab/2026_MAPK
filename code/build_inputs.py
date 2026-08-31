"""Build the consolidated Excel inputs from the raw files they came from.

Every notebook reads workbooks in the repository root rather than the directory
trees they were built from:

``TableS3-cBioPortal-Annotations.xlsx``       one sheet per gene, plus a summary
``TableS5-PocketDistances-Annotations.xlsx``  one sheet per structure, plus a summary
``TableS4-DMS-Annotations.xlsx``              one sheet, all three DMS datasets stacked
``TableS6-LibraryDesign.xlsx``                the sgRNA library design inputs

Consolidating them means a reviewer downloads a handful of files instead of
sixty, and the expensive steps — parsing eighteen mmCIF files with Biopython to
measure every residue's distance to its ligand, reading four large design and
score tables — run once here rather than on every notebook execution. Nothing
in the pipeline reads ``Inputs/`` or ``library_design_JW/inputs/`` any more.

Two root workbooks are not built here because they arrive already consolidated:
``TableS1-ScreenData.xlsx`` and ``TableS7-OtherFiguresData.xlsx``.

This module is the only place the raw formats are understood. It is a build
step, not part of a run: execute it when a raw input changes.

    python -m code.build_inputs

The readers below are the ones the analysis modules used to carry, moved here
unchanged, so the workbooks hold exactly what the previous code loaded.
"""

from __future__ import annotations

import re

import pandas as pd
from Bio import SeqIO

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
# LIBRARY DESIGN
# ---------------------------------------------------------------------------


def read_mane_transcripts() -> pd.DataFrame:
    """The target genes mapped to their MANE Select transcripts.

    Restricting the MANE v1.0 summary to ``MANE Select`` yields exactly one
    transcript per gene. The version suffix is dropped
    (ENST00000644974.5 -> ENST00000644974) because that is the form the design
    script expects.

    Only the 22 target genes are kept. The other 19,040 MANE Select rows are
    not part of this manuscript, and the pipeline used them only as a lookup.
    """
    mane = pd.read_table(cfg.LIB_MANE_SUMMARY)
    mane = mane[mane['MANE_status'] == 'MANE Select']

    gene_info = mane[['symbol', 'name', 'Ensembl_nuc', 'RefSeq_nuc']].copy()
    gene_info.columns = ['Gene', 'name', 'Ensembl_nuc', 'RefSeq_nuc']
    gene_info['Ensembl_nuc'] = [re.match(r'ENST\d+', x).group(0)
                                for x in gene_info['Ensembl_nuc']]

    targets = pd.read_csv(cfg.LIB_INITIAL_GENE_LIST)
    transcripts = targets.merge(gene_info, 'left', on='Gene')
    assert transcripts['Ensembl_nuc'].notna().all(), \
        'some genes did not map to a MANE transcript'
    return transcripts


def read_exon_records() -> pd.DataFrame:
    """The UCSC exon FASTA as one row per record.

    Record IDs look like ``hg38_ncbiRefSeq_NM_001654.5_3``; the RefSeq
    accession is parsed back out when BEhive looks a guide's context up, so it
    is stored here rather than the parsed form. Only the forward strand is
    stored - the reverse complement is derived at read time.
    """
    records = [{'Record': record.id, 'Sequence': str(record.seq)}
               for record in SeqIO.parse(str(cfg.LIB_EXON_FASTA), 'fasta')]
    return pd.DataFrame(records)


def build_library_design(path=None) -> pd.DataFrame:
    """Write the library-design workbook; return the base-editor summary sheet.

    Sheets, in order: the two base editors and their design parameters, the
    target-gene transcript mapping, one design table per base editor, the
    precomputed BEhive efficiency scores, the FlashFry specificity scores, and
    the exon sequences BEhive reads a guide's genomic context from.

    ClinVar's ``variant_summary.txt`` is deliberately absent: it is 3.9 GB, it
    is only read when step 2 regenerates designs, and the designs it annotated
    are stored here already.
    """
    path = path or cfg.LIBRARY_DESIGN_XLSX

    base_editors = pd.read_table(cfg.LIB_BE_PARAMETERS)
    transcripts = read_mane_transcripts()
    exons = read_exon_records()

    designs = {be: pd.read_table(cfg.LIB_SGRNA_DESIGNS_DIR /
                                 f'sgrna_designs_MAPK_{be}.txt')
               for be in (cfg.LIB_ABE_NAME, cfg.LIB_CBE_NAME)}

    # Written with an index and read back with index_col=0, so the CSV carries
    # a second, nameless index column. It holds nothing.
    behive = pd.read_csv(cfg.LIB_BEHIVE_PRECOMPUTED, index_col=0)
    behive = behive.drop(columns=[c for c in behive.columns
                                  if c.startswith('Unnamed:')])
    flashfry = pd.read_table(cfg.LIB_FLASHFRY_PRECOMPUTED)

    summary = base_editors.assign(
        **{'Designed guides': [len(designs[be]) for be in base_editors['BEs']]})

    with pd.ExcelWriter(path, **_WRITER_KWARGS) as writer:
        summary.to_excel(writer, sheet_name=cfg.LIB_BASE_EDITOR_SHEET, index=False)
        transcripts.to_excel(writer, sheet_name=cfg.LIB_TRANSCRIPT_SHEET, index=False)
        for be, table in designs.items():
            table.to_excel(writer, sheet_name=be, index=False)
        behive.to_excel(writer, sheet_name=cfg.LIB_BEHIVE_SHEET, index=False)
        flashfry.to_excel(writer, sheet_name=cfg.LIB_FLASHFRY_SHEET, index=False)
        exons.to_excel(writer, sheet_name=cfg.LIB_EXON_SHEET, index=False)
    return summary


# ---------------------------------------------------------------------------


def build_all() -> dict[str, pd.DataFrame]:
    """Rebuild every workbook that is built from raw files."""
    return {
        'cBioPortal': build_cbioportal(),
        'PocketDistances': build_pocket_distances(),
        'DMS': build_dms(),
        'LibraryDesign': build_library_design(),
    }


if __name__ == '__main__':
    for name, summary in build_all().items():
        print(f'{name}: {len(summary)} rows')
