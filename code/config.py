"""Shared configuration for the MAPK base-editing screen analysis.

Every constant that was duplicated across the ten original notebooks lives here
exactly once. Nothing in this module reads data or draws anything.

Path resolution
---------------
``REPO_ROOT`` is discovered at import time so the same code runs unchanged in
Google Colab (where the repository is cloned into ``/content``) and locally.
No absolute Google Drive paths appear anywhere.
"""

from pathlib import Path

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------


def _find_repo_root() -> Path:
    """Directory holding TableS1-ScreenData.xlsx, searched upward from here."""
    for candidate in [Path(__file__).resolve().parent, *Path(__file__).resolve().parents]:
        if (candidate / "TableS1-ScreenData.xlsx").exists():
            return candidate
    # Fall back to the parent of this file's directory (code/ -> repo root).
    return Path(__file__).resolve().parent.parent


REPO_ROOT = _find_repo_root()

SCREEN_DATA_XLSX = REPO_ROOT / "TableS1-ScreenData.xlsx"
INPUTS_DIR = REPO_ROOT / "Inputs"
DMS_DIR = INPUTS_DIR / "DMS-Data"
CBIOPORTAL_DIR = INPUTS_DIR / "cBioPortal"
STRUCTURES_DIR = INPUTS_DIR / "Pocket-Structures"
# Residue-to-ligand distances measured from the mmCIF files above, written by
# code/build_inputs.py.
POCKET_DISTANCES_DIR = INPUTS_DIR / "Pocket-Distances"

# CONSOLIDATED EXCEL INPUTS
#
# The analysis reads these three workbooks, not the directories above. Each is
# built by `python -m code.build_inputs` from the raw files, which stay in
# Inputs/ as the provenance of what the workbooks contain. Reading them keeps
# a run to four spreadsheets and takes the eighteen-structure Biopython parse
# out of the hot path.
CBIOPORTAL_XLSX = REPO_ROOT / "TableS3-cBioPortal-Annotations.xlsx"
POCKET_DISTANCES_XLSX = REPO_ROOT / "TableS5-PocketDistances-Annotations.xlsx"
DMS_XLSX = REPO_ROOT / "TableS4-DMS-Annotations.xlsx"

# Summary sheets, first in their workbook; every other sheet is data.
CBIOPORTAL_TRANSCRIPT_SHEET = "Transcripts"
POCKET_STRUCTURE_SHEET = "Structures"
DMS_SHEET = "DMS"

# Guide-level editing-window position, kept only as a deterministic row
# ordering key for the GMM trajectory. It is not used as a plot axis; every
# position axis uses AtoG_pos / CtoT_pos via the `pos` column. See errors.md.
#
# This is now a sheet of TableS1 rather than Inputs/guide-edit-site.csv, so the
# ordering key ships with the screen data it orders. The retired CSV stays in
# Inputs/ as provenance; the sheet holds the same 15,364 rows in the same
# order. See errors.md §1.
EDIT_SITE_SHEET = "edit-site"
OUTPUTS_DIR = REPO_ROOT / "Outputs"

# Interactive companion site, written by MAPK_Interactive_Analysis.ipynb.
#
# `docs/` is the site root because GitHub Pages can only publish from the
# repository root, `/docs`, or a gh-pages branch, and every file the pages
# link to has to sit underneath whichever of those is chosen. The plots
# therefore live in their own folder *inside* the site root rather than
# beside `Outputs/`: OUT_INTERACTIVE holds every interactive figure and is
# entirely separate from OUTPUTS_DIR, which holds the print SVGs.
SITE_DIR = REPO_ROOT / "docs"
OUT_INTERACTIVE = SITE_DIR / "figures"
SITE_MANIFEST = OUT_INTERACTIVE / "manifest.json"

# Backwards-compatible alias.
SITE_FIGURES = OUT_INTERACTIVE

# One subdirectory per figure family, mirroring the original notebooks.
OUT_SCATTER = OUTPUTS_DIR / "01-Scatterplots"
OUT_SCATTER_COMBINED = OUTPUTS_DIR / "02-Scatterplots-Combined"
OUT_GOF_VS_VALID = OUTPUTS_DIR / "03-GOF-vs-Validation"
OUT_POCKET_PDB = OUTPUTS_DIR / "04-Pocket-Hits-PDB"
OUT_POCKET_AF = OUTPUTS_DIR / "05-Pocket-Hits-AlphaFold"
OUT_SPLICE_HEATMAP = OUTPUTS_DIR / "06-Splice-Validation-Heatmaps"
OUT_SANKEY = OUTPUTS_DIR / "07-Sankey"
OUT_TRAJECTORY = OUTPUTS_DIR / "08-GMM-Trajectory"
OUT_DMS = OUTPUTS_DIR / "09-LOF-vs-DMS"
OUT_LOLLIPOP = OUTPUTS_DIR / "10-cBioPortal-Lollipops"
OUT_BOXPLOT = OUTPUTS_DIR / "11-Boxplots"
OUT_SUMMARIES = OUTPUTS_DIR / "00-Summaries"

ALL_OUTPUT_DIRS = [
    OUT_SUMMARIES, OUT_SCATTER, OUT_SCATTER_COMBINED, OUT_GOF_VS_VALID,
    OUT_POCKET_PDB, OUT_POCKET_AF, OUT_SPLICE_HEATMAP, OUT_SANKEY,
    OUT_TRAJECTORY, OUT_DMS, OUT_LOLLIPOP, OUT_BOXPLOT,
]


def make_output_dirs() -> None:
    """Create every output subdirectory. Safe to call repeatedly."""
    for directory in ALL_OUTPUT_DIRS:
        directory.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# GENES AND COLORS
# ---------------------------------------------------------------------------

# Pathway order, N-terminal input to terminal kinase. Every per-gene list in
# this module (y-ranges, cBioPortal axis limits) is indexed in this order.
GENES = [
    'SHC1', 'GRB2', 'PTPN11',
    'KRAS', 'NRAS', 'HRAS', 'MRAS',
    'PPP1CA', 'PPP1CB', 'PPP1CC',
    'SHOC2',
    'ARAF', 'BRAF', 'RAF1',
    'KSR1', 'KSR2',
    'MAP2K1', 'MAP2K2',
    'DUSP4', 'DUSP6',
    'MAPK1', 'MAPK3',
]

GENE_COLORS = [
    '#E099ED', '#F59DC5', '#7584E6',
    '#75CBE6', '#75A7E6', '#A6B1E9', '#B6C1E9',
    '#75E687', '#8DE675', '#BCE675',
    '#95E6DC',
    '#E69F75', '#E68C75', '#E67575',
    '#E6CF75', '#E6E475',
    '#75E6AE', '#BBE8B7',
    '#E8CFB0', '#CFB79D',
    '#E9B5B5', '#D9A9A9',
]

DARKGRAY = '#767676'
LIGHTGRAY = '#a9a9a9'

# Centimetres to inches. Every figure size in this module is written in cm,
# as the original notebooks were.
CM = 1 / 2.54

COLOR_MAP = dict(zip(GENES, GENE_COLORS))

# Gene symbol -> protein name, used by the splice heatmaps and the Sankey.
GENE_PROTEIN_MAP = {
    'PTPN11': 'SHP2',
    'MAP2K1': 'MEK1',
    'MAP2K2': 'MEK2',
    'MAPK1': 'ERK2',
    'MAPK3': 'ERK1',
    'RAF1': 'CRAF',
}

PROTEINS = [GENE_PROTEIN_MAP.get(g, g) for g in GENES]

# Colors keyed by protein name as well, so Sankey nodes resolve either way.
PROTEIN_COLOR_MAP = {**COLOR_MAP, **{p: COLOR_MAP[g] for g, p in GENE_PROTEIN_MAP.items()}}

# Mutation-type palette shared by every scatterplot.
MUT_PAL = {
    "Missense": "#FDB462",
    "No Mutation": "#80B1D3",
    "Silent": "#80B1D3",
    "Nonsense": "#C03221",
}

# Rows whose Gene is one of these are library controls, not pathway genes.
CONTROL_GENES = ['Control', 'Essential splice site']

# ---------------------------------------------------------------------------
# PROTEIN LENGTHS AND DOMAINS
# ---------------------------------------------------------------------------

XMAX_DICT = {
    'SHC1': 584, 'GRB2': 217,
    'PTPN11': 593, 'KRAS': 188, 'NRAS': 189, 'HRAS': 189,
    'MRAS': 208, 'PPP1CA': 330, 'PPP1CB': 327, 'PPP1CC': 323,
    'SHOC2': 582, 'ARAF': 606, 'BRAF': 766, 'RAF1': 648,
    'KSR1': 923, 'KSR2': 950, 'MAP2K1': 393, 'MAP2K2': 400,
    'DUSP4': 394, 'DUSP6': 381, 'MAPK1': 360, 'MAPK3': 379,
}

# The revised domain boundaries, marked "### Changed after edits" in the three
# most recently edited notebooks. The older boundaries that appeared in the
# other seven notebooks are deliberately not carried over; see errors.md for
# the list of figures whose domain shading this changes.
DOMAINS_LIST = {
    'SHC1': [
        {'description': 'PID', 'start': 156, 'end': 339, 'color': COLOR_MAP['SHC1']},
        {'description': 'SH2', 'start': 488, 'end': 579, 'color': COLOR_MAP['SHC1']},
    ],
    'GRB2': [
        {'description': 'SH3', 'start': 1, 'end': 58, 'color': COLOR_MAP['GRB2']},
        {'description': 'SH2', 'start': 60, 'end': 152, 'color': COLOR_MAP['GRB2']},
        {'description': 'SH3', 'start': 156, 'end': 215, 'color': COLOR_MAP['GRB2']},
    ],
    'PTPN11': [
        {'description': 'N-SH2', 'start': 6, 'end': 102, 'color': COLOR_MAP['PTPN11']},
        {'description': 'C-SH2', 'start': 112, 'end': 216, 'color': COLOR_MAP['PTPN11']},
        {'description': 'PTPase', 'start': 247, 'end': 517, 'color': COLOR_MAP['PTPN11']},
    ],
    'KRAS': [
        {'description': 'GTPase', 'start': 1, 'end': 189, 'color': COLOR_MAP['KRAS']},
    ],
    'NRAS': [
        {'description': 'GTPase', 'start': 1, 'end': 189, 'color': COLOR_MAP['NRAS']},
    ],
    'HRAS': [
        {'description': 'GTPase', 'start': 1, 'end': 188, 'color': COLOR_MAP['HRAS']},
    ],
    'MRAS': [
        {'description': 'GTPase', 'start': 9, 'end': 208, 'color': COLOR_MAP['MRAS']},
    ],
    'PPP1CA': [
        {'description': 'PP1/PPKL metallophosphatase', 'start': 12, 'end': 305, 'color': COLOR_MAP['PPP1CA']},
    ],
    'PPP1CB': [
        {'description': 'PP1/PPKL metallophosphatase', 'start': 23, 'end': 302, 'color': COLOR_MAP['PPP1CB']},
    ],
    'PPP1CC': [
        {'description': 'PP1/PPKL metallophosphatase', 'start': 8, 'end': 298, 'color': COLOR_MAP['PPP1CC']},
    ],
    'SHOC2': [
        {'description': 'LRR', 'start': 101, 'end': 122, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR', 'start': 124, 'end': 145, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR', 'start': 147, 'end': 169, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR-SDS22', 'start': 170, 'end': 191, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR', 'start': 193, 'end': 214, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR-SDS22', 'start': 216, 'end': 237, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR', 'start': 239, 'end': 260, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR', 'start': 262, 'end': 283, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR', 'start': 285, 'end': 307, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR', 'start': 308, 'end': 329, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR', 'start': 332, 'end': 353, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR', 'start': 356, 'end': 377, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR', 'start': 380, 'end': 400, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR', 'start': 403, 'end': 424, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR-SDS22', 'start': 426, 'end': 448, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR', 'start': 449, 'end': 470, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR', 'start': 472, 'end': 494, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR', 'start': 495, 'end': 516, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR', 'start': 518, 'end': 540, 'color': COLOR_MAP['SHOC2']},
        {'description': 'LRR', 'start': 542, 'end': 563, 'color': COLOR_MAP['SHOC2']},
    ],
    'ARAF': [
        {'description': 'RBD', 'start': 19, 'end': 91, 'color': COLOR_MAP['ARAF']},
        {'description': 'CRD', 'start': 98, 'end': 144, 'color': COLOR_MAP['ARAF']},
        {'description': 'Inhib pSite', 'start': 214, 'end': 214, 'color': COLOR_MAP['ARAF']},
        {'description': 'Kinase', 'start': 310, 'end': 570, 'color': COLOR_MAP['ARAF']},
    ],
    'BRAF': [
        {'description': 'RBD', 'start': 155, 'end': 227, 'color': COLOR_MAP['BRAF']},
        {'description': 'CRD', 'start': 234, 'end': 280, 'color': COLOR_MAP['BRAF']},
        {'description': 'Inhib pSite', 'start': 365, 'end': 365, 'color': COLOR_MAP['BRAF']},
        {'description': 'Kinase', 'start': 457, 'end': 717, 'color': COLOR_MAP['BRAF']},
    ],
    'RAF1': [
        {'description': 'RBD', 'start': 56, 'end': 131, 'color': COLOR_MAP['RAF1']},
        {'description': 'CRD', 'start': 138, 'end': 184, 'color': COLOR_MAP['RAF1']},
        {'description': 'Inhib pSite', 'start': 259, 'end': 259, 'color': COLOR_MAP['RAF1']},
        {'description': 'Kinase', 'start': 349, 'end': 609, 'color': COLOR_MAP['RAF1']},
    ],
    'KSR1': [
        {'description': 'CRD', 'start': 347, 'end': 391, 'color': COLOR_MAP['KSR1']},
        {'description': 'Kinase', 'start': 613, 'end': 883, 'color': COLOR_MAP['KSR1']},
    ],
    'KSR2': [
        {'description': 'CRD', 'start': 412, 'end': 456, 'color': COLOR_MAP['KSR2']},
        {'description': 'Kinase', 'start': 666, 'end': 931, 'color': COLOR_MAP['KSR2']},
    ],
    'MAP2K1': [
        {'description': 'Kinase', 'start': 68, 'end': 361, 'color': COLOR_MAP['MAP2K1']},
    ],
    'MAP2K2': [
        {'description': 'Kinase', 'start': 72, 'end': 369, 'color': COLOR_MAP['MAP2K2']},
    ],
    'DUSP4': [
        {'description': 'Rhodanese-like', 'start': 41, 'end': 159, 'color': COLOR_MAP['DUSP4']},
        {'description': 'Phosphatase', 'start': 195, 'end': 336, 'color': COLOR_MAP['DUSP4']},
    ],
    'DUSP6': [
        {'description': 'Rhodanese-like', 'start': 30, 'end': 148, 'color': COLOR_MAP['DUSP6']},
        {'description': 'Phosphatase', 'start': 206, 'end': 349, 'color': COLOR_MAP['DUSP6']},
    ],
    'MAPK1': [
        {'description': 'Kinase', 'start': 25, 'end': 313, 'color': COLOR_MAP['MAPK1']},
    ],
    'MAPK3': [
        {'description': 'Kinase', 'start': 44, 'end': 304, 'color': COLOR_MAP['MAPK3']},
    ],
}

# ---------------------------------------------------------------------------
# SCREEN SHEETS AND CONDITION COLUMNS
# ---------------------------------------------------------------------------

# Sheet name in TableS1-ScreenData.xlsx, keyed by (screen, editor).
SCREEN_SHEETS = {
    ('gof', 'ABE'): 'GOF-ABE',
    ('gof', 'CBE'): 'GOF-CBE',
    ('lof', 'ABE'): 'LOF-ABE',
    ('lof', 'CBE'): 'LOF-CBE',
    ('meki', 'ABE'): 'MEKi-ABE',
    ('meki', 'CBE'): 'MEKi-CBE',
    ('valid', 'ABE'): 'Validation-ABE',
    ('valid', 'CBE'): 'Validation-CBE',
    ('lof_erk', 'ABE'): 'LOF-ERK-ABE',
    ('lof_erk', 'CBE'): 'LOF-ERK-CBE',
}

# Read-count column used by the Day-0 representation filter, per screen.
DAY0_COUNT_COL = {
    'gof': 'Day 0',
    'lof': 'Day 0',
    'meki': 'Day 0',
    'valid': 'Day0',
    'lof_erk': {'ABE': 'Day0-ABE', 'CBE': 'Day0-CBE'},
}

DAY0_COUNT_CUTOFF = 100

GOF_SCREEN_COLUMNS = [
    'DMSO-Day0-Z',
    'SHP2i-DMSO-Z', 'KRASi-DMSO-Z', 'RAFi-DMSO-Z', 'Trametinib-DMSO-Z', 'ERKi-DMSO-Z',
    'SHP2i-Day0-Z', 'KRASi-Day0-Z', 'RAFi-Day0-Z', 'Trametinib-Day0-Z', 'ERKi-Day0-Z',
]

# The five drug-vs-DMSO comparisons, used wherever "the inhibitor conditions"
# are meant rather than all eleven GOF columns.
GOF_DMSO_COLUMNS = GOF_SCREEN_COLUMNS[1:6]

MEKI_SCREEN_COLUMNS = [
    'DMSO-Day0-Z',
    'Trametinib-DMSO-Z', 'Trametiglue-DMSO-Z', 'GDC0623-DMSO-Z', 'Avutametinib-DMSO-Z',
    'Selumetinib-DMSO-Z', 'Mirdametinib-DMSO-Z', 'MAP855-DMSO-Z', 'Cobimetinib-DMSO-Z',
    'APS_Adagrasib-DMSO-Z',
    'Trametinib-Day0-Z', 'Trametiglue-Day0-Z', 'GDC0623-Day0-Z', 'Avutametinib-Day0-Z',
    'Selumetinib-Day0-Z', 'Mirdametinib-Day0-Z', 'MAP855-Day0-Z', 'Cobimetinib-Day0-Z',
    'APS_Adagrasib-Day0-Z',
]

MEKI_DMSO_COLUMNS = MEKI_SCREEN_COLUMNS[1:10]

LOF_SCREEN_COLUMNS = ['Dox-UT-Z']

# Dox-treated versus Day 0, added to the LOF sheets after the rest of this
# analysis was written. Only the Figure 1 control boxplot uses it, and it
# reads the depletion of the essential-splice-site guides against the
# starting library rather than against the untreated arm.
LOF_DAY0_COLUMN = 'UT-Day0-Z'

LOF_ERK_COLUMN = 'Day14-Dox-UT-Z'

VALID_ABE_SCREEN_COLUMNS = [
    'DMSO-Day18-Day0-Z',
    'SHP2i-Day18-DMSO-Z', 'KRASi-Day18-DMSO-Z', 'RAFi-Day18-DMSO-Z',
    'Trametinib-Day18-DMSO-Z', 'ERKi-Day18-DMSO-Z',
    'SHP2i-Day18-Day0-Z', 'KRASi-Day18-Day0-Z', 'RAFi-Day18-Day0-Z',
    'Trametinib-Day18-Day0-Z', 'ERKi-Day18-Day0-Z',
]

VALID_CBE_SCREEN_COLUMNS = [
    'DMSO-Day0-Z',
    'SHP2i-DMSO-Z', 'KRASi-DMSO-Z', 'RAFi-DMSO-Z', 'Trametinib-DMSO-Z', 'ERKi-DMSO-Z',
    'SHP2i-Day0-Z', 'KRASi-Day0-Z', 'RAFi-Day0-Z', 'Trametinib-Day0-Z', 'ERKi-Day0-Z',
]

# Validation-ABE is a Day18/Day32 timecourse; renaming its Day18 columns onto
# the CBE naming lets both editors share one set of condition names.
VALID_ABE_TO_COMMON = dict(zip(VALID_ABE_SCREEN_COLUMNS, VALID_CBE_SCREEN_COLUMNS))

# Human-readable condition labels for the splice heatmaps.
HEATMAP_CONDITION_LABELS = ['DMSO', 'SHP2i', 'KRASi', 'RAFi', 'MEKi', 'ERKi']

# The splice heatmaps are drawn twice: once across the whole pathway, and
# once across the three RAF paralogs alone. The RAF panel keeps the full
# panel's width so the two stack in the figure, at half its height.
SPLICE_HEATMAP_FIGSIZE = (2.0, 1.6)
SPLICE_HEATMAP_RAF_FIGSIZE = (SPLICE_HEATMAP_FIGSIZE[0],
                             SPLICE_HEATMAP_FIGSIZE[1] / 2)

# Protein names, since the splice heatmaps relabel RAF1 as CRAF.
SPLICE_HEATMAP_RAF_PROTEINS = ['ARAF', 'BRAF', 'CRAF']

# Condition -> inhibitor label, used by the Sankey.
INHIBITOR_LABELS = {
    'SHP2i-DMSO-Z': 'SHP2i',
    'KRASi-DMSO-Z': 'KRASi',
    'RAFi-DMSO-Z': 'RAFi',
    'Trametinib-DMSO-Z': 'MEKi',
    'ERKi-DMSO-Z': 'ERKi',
}

# ---------------------------------------------------------------------------
# PER-GENE BOXPLOTS
# ---------------------------------------------------------------------------
# One boxplot per screen, editor and condition, showing the Z-score
# distribution of every guide targeting each gene, with the two control
# categories alongside. Unlike most other figures these keep the controls, so
# each gene reads against the library background it was scored on.

BOXPLOT_CONTROL_ORDER = ['Essential splice site', 'Control']

# Screen names as the figures label them; 'MEKi' is not upper-case.
BOXPLOT_SCREEN_LABELS = {'gof': 'GOF', 'meki': 'MEKi', 'lof': 'LOF'}

# The originals passed a 22-color list next to a 24-category order, so seaborn
# recycled the palette and drew the two control categories in the SHC1 and
# GRB2 colors. Naming every category removes the collision; see errors.md.
BOXPLOT_CONTROL_COLORS = {
    'Essential splice site': DARKGRAY,
    'Control': LIGHTGRAY,
}

# Genes present in each screen's library, in pathway order.
BOXPLOT_GENES = {
    'gof': GENES,
    'meki': GENES[11:18],        # ARAF through MAP2K2
    'lof': GENES[:20],           # everything but MAPK1 and MAPK3
}

# (ylim, n_ticks) per screen and editor, as drawn in the original notebook.
BOXPLOT_LIMITS = {
    ('gof', 'ABE'): ((-6, 6), 7),
    ('gof', 'CBE'): ((-6, 4), 6),
    ('meki', 'ABE'): ((-6, 6), 7),
    ('meki', 'CBE'): ((-6, 6), 5),
    ('lof', 'ABE'): ((-5, 5), 5),
    ('lof', 'CBE'): ((-5, 5), 5),
}

BOXPLOT_FIGSIZES = {
    'gof': (10 * CM, 1.0),
    'meki': (5 * CM, 1.3),
    'lof': (10 * CM, 1.0),
}

# Figure 1 controls: essential-splice-site guides against non-targeting
# controls, both editors on one axis, showing that the library separates a
# real depletion phenotype from background.
BOXPLOT_CONTROL_FIGURES = [
    # (screen, column, ylim, n_ticks)
    ('gof', 'DMSO-Day0-Z', (-6, 4), 6),
    ('lof', LOF_DAY0_COLUMN, (-8, 4), 7),
]

BOXPLOT_CONTROL_LABELS = ['ABE Control', 'ABE Essential', 'CBE Control', 'CBE Essential']
BOXPLOT_CONTROL_PALETTE = dict(zip(BOXPLOT_CONTROL_LABELS,
                                   ['gray', 'red', 'gray', 'red']))
BOXPLOT_CONTROL_FIGSIZE = (2 * CM, 1.0)

# CONTROL-VERSUS-ESSENTIAL TEST
#
# The control panel asserts that the library separates a real depletion
# phenotype from background, so each editor's non-targeting guides are tested
# against its essential-splice-site guides. Rank-based: the essential guides
# are strongly left-skewed in every screen and the two groups differ in size
# by an order of magnitude, neither of which a t-test tolerates.
BOXPLOT_CONTROL_TEST = 'Mann-Whitney U, two-sided'
BOXPLOT_CONTROL_PAIRS = [
    ('ABE Control', 'ABE Essential'),
    ('CBE Control', 'CBE Essential'),
]
# Significance labels, most stringent first.
PVALUE_STARS = [(1e-4, '****'), (1e-3, '***'), (1e-2, '**'), (0.05, '*')]
PVALUE_NS = 'ns'

# ---------------------------------------------------------------------------
# HIT-CALLING CUTOFFS
# ---------------------------------------------------------------------------

GOF_HIT_CUTOFF = 3.0
VALID_HIT_CUTOFF = 2.0
POCKET_HIT_CUTOFF = 3.0
POCKET_DISTANCE_CUTOFF = 4.0     # angstroms
GUIDE_POSITION_WINDOW = 1        # residues either side of a hit position
SANKEY_VALID_CUTOFF = 1.0
SANKEY_MIN_HIT_COUNT = 2

# ---------------------------------------------------------------------------
# CONDITION PALETTES (combined scatterplots)
# ---------------------------------------------------------------------------

COMBINED_GOF_COLUMNS = GOF_SCREEN_COLUMNS[:6]
COMBINED_MEKI_COLUMNS = MEKI_SCREEN_COLUMNS[:10]
COMBINED_LOF_COLUMNS = ['DMSO-Day0-Z', 'Dox-UT-Z']

GOF_CONDITION_PAL = dict(zip(
    COMBINED_GOF_COLUMNS,
    [LIGHTGRAY, '#7584E6', '#75CBE6', '#E69F75', '#75E6AE', '#E9B5B5'],
))

MEKI_CONDITION_PAL = dict(zip(
    COMBINED_MEKI_COLUMNS,
    [LIGHTGRAY, '#75E6AE', '#BBE8B7', '#66C897', '#A4CBA0',
     '#396F54', '#57AA80', '#8EAF8A', '#488E6B', '#779374'],
))

LOF_CONDITION_PAL = dict(zip(COMBINED_LOF_COLUMNS, [LIGHTGRAY, '#722F37']))

# GOF-vs-Validation correlation scatterplot, one color per GOF condition.
CORR_SCATTER_COLORS = [
    '#888888',
    COLOR_MAP['PTPN11'],
    COLOR_MAP['KRAS'],
    COLOR_MAP['BRAF'],
    COLOR_MAP['MAP2K1'],
    COLOR_MAP['MAPK1'],
]

# ---------------------------------------------------------------------------
# PER-GENE AXIS RANGES
# ---------------------------------------------------------------------------
# Each entry is (ymin, ymax, n_ticks) in GENES order. An empty tuple means the
# gene is absent from that screen and is skipped. These differ between the
# single-condition and combined-condition figures because they are different
# figures with different point densities; they are kept separate on purpose.

SCATTER_RANGES_GOF = [
    (-7, 7, 5), (-6, 6, 3), (-8, 24, 5),
    (-8, 24, 5), (-5, 15, 5), (-6, 6, 3), (-7, 7, 5),
    (-7, 7, 5), (-7, 7, 5), (-6, 6, 5),
    (-9, 6, 6), (-15, 20, 8), (-5, 10, 4), (-7, 14, 4),
    (-8, 8, 5), (-8, 8, 5), (-14, 21, 5), (-10, 10, 5),
    (-8, 8, 5), (-6, 6, 3), (-10, 15, 6), (-8, 12, 6),
]

SCATTER_RANGES_MEKI = [
    (), (), (),
    (), (), (), (),
    (), (), (),
    (), (-6, 6, 3), (-12, 12, 5), (-8, 8, 5),
    (-8, 8, 5), (-8, 8, 5), (-8, 16, 4), (-8, 8, 5),
    (), (), (), (),
]

SCATTER_RANGES_LOF = [
    (-6, 12, 4), (-6, 12, 4), (-6, 6, 3),
    (-6, 12, 4), (-6, 6, 3), (-6, 6, 3), (-6, 6, 3),
    (-8, 8, 5), (-6, 6, 3), (-6, 6, 3),
    (-6, 6, 3), (-6, 6, 3), (-6, 6, 3), (-6, 12, 4),
    (-6, 6, 3), (-6, 6, 3), (-6, 6, 3), (-6, 12, 4),
    (-8, 8, 5), (-6, 6, 3), (), (),
]

# LOF-ERK screen, MAPK1 only; the entry for every other gene is unused.
SCATTER_RANGES_LOF_ERK = {'MAPK1': (-4, 6, 6)}

COMBINED_RANGES_GOF = [
    (-8, 8, 5), (-6, 6, 5), (-8, 24, 5),
    (-8, 24, 5), (-5, 15, 5), (-6, 6, 5), (-8, 8, 5),
    (-8, 8, 5), (-8, 8, 5), (-6, 6, 5),
    (-9, 9, 7), (-10, 20, 4), (-6, 12, 4), (-7, 14, 4),
    (-8, 8, 5), (-8, 8, 5), (-14, 21, 6), (-10, 10, 5),
    (-8, 8, 5), (-6, 6, 3), (-10, 15, 6), (-8, 12, 6),
]

COMBINED_RANGES_MEKI = [
    (), (), (),
    (), (), (), (),
    (), (), (),
    (), (-6, 6, 3), (-12, 12, 5), (-9, 9, 7),
    (-8, 8, 5), (-8, 8, 5), (-8, 16, 4), (-9, 9, 7),
    (), (), (), (),
]

COMBINED_RANGES_LOF = [
    (-8, 8, 5), (-6, 8, 8), (-6, 6, 5),
    (-8, 12, 6), (-6, 6, 5), (-6, 6, 5), (-6, 6, 5),
    (-8, 8, 5), (-8, 8, 5), (-6, 6, 5),
    (-8, 8, 5), (-8, 8, 5), (-6, 6, 5), (-8, 8, 5),
    (-8, 8, 5), (-8, 8, 5), (-6, 6, 5), (-6, 9, 6),
    (-8, 8, 5), (-6, 6, 5), (-6, 6, 5), (-6, 6, 5),
]

# X-axis padding around the protein, in residues.
SCATTER_SPACING = 20
COMBINED_SPACING = 15

# Figure sizes, in inches (the notebooks wrote these as cm / 2.54).
FIGSIZE_FIVE_ACROSS = (4.3 * CM, 2.8 * CM)
FIGSIZE_FIVE_ACROSS_LARGE = (1.5 * 4.3 * CM, 1.5 * 2.8 * CM)

# ---------------------------------------------------------------------------
# STRUCTURES
# ---------------------------------------------------------------------------
# Co-crystal structures: an inhibitor bound to its target. Distances are
# measured from the ligand's heavy atoms to every residue of the protein chain.

PDB_STRUCTURES = [
    # (gene, structure, inhibitor, ligand_resname, ligand_chain, protein_chain, condition)
    ('PTPN11', '7JVM', 'Batoprotafib', 'VKS', 'A', 'A', 'SHP2i-DMSO-Z'),
    ('KRAS',   '6UT0', 'Adagrasib',    'M1X', 'A', 'A', 'KRASi-DMSO-Z'),
    ('BRAF',   '5C9C', 'LY3009120',    '4Z5', 'A', 'A', 'RAFi-DMSO-Z'),
    ('MAP2K1', '7JUR', 'Trametinib',   'QOM', 'C', 'C', 'Trametinib-DMSO-Z'),
    ('MAPK1',  '6RQ4', 'Temuterkib',   'KE8', 'A', 'A', 'ERKi-DMSO-Z'),
]

# Paralogs of each node, for the "hits outside this node" breakdown.
PDB_PARALOGS = {
    'PTPN11': [],
    'KRAS': ['NRAS', 'HRAS', 'MRAS'],
    'BRAF': ['ARAF', 'RAF1'],
    'MAP2K1': ['MAP2K2'],
    'MAPK1': ['MAPK3'],
}

# AlphaFold models with the corresponding inhibitor aligned into the pocket,
# so every paralog gets its own distance profile.
AF_STRUCTURES = [
    ('PTPN11', 'PTPN11-AF', 'Batoprotafib', 'VKS', 'A', 'A', 'SHP2i-DMSO-Z'),
    ('KRAS',   'KRAS-AF',   'Adagrasib',    'M1X', 'A', 'A', 'KRASi-DMSO-Z'),
    ('NRAS',   'NRAS-AF',   'Adagrasib',    'M1X', 'A', 'A', 'KRASi-DMSO-Z'),
    ('HRAS',   'HRAS-AF',   'Adagrasib',    'M1X', 'A', 'A', 'KRASi-DMSO-Z'),
    ('MRAS',   'MRAS-AF',   'Adagrasib',    'M1X', 'A', 'A', 'KRASi-DMSO-Z'),
    ('ARAF',   'ARAF-AF',   'LY3009120',    '4Z5', 'A', 'A', 'RAFi-DMSO-Z'),
    ('BRAF',   'BRAF-AF',   'LY3009120',    '4Z5', 'A', 'A', 'RAFi-DMSO-Z'),
    ('RAF1',   'RAF1-AF',   'LY3009120',    '4Z5', 'A', 'A', 'RAFi-DMSO-Z'),
    ('MAP2K1', 'MAP2K1-AF', 'Trametinib',   'QOM', 'C', 'A', 'Trametinib-DMSO-Z'),
    ('MAP2K2', 'MAP2K2-AF', 'Trametinib',   'QOM', 'C', 'A', 'Trametinib-DMSO-Z'),
    ('MAPK1',  'MAPK1-AF',  'Temuterkib',   'KE8', 'A', 'A', 'ERKi-DMSO-Z'),
    ('MAPK3',  'MAPK3-AF',  'Temuterkib',   'KE8', 'A', 'A', 'ERKi-DMSO-Z'),
]

AF_PARALOGS = {
    'PTPN11': [],
    'KRAS': ['NRAS', 'HRAS', 'MRAS'],
    'NRAS': ['KRAS', 'HRAS', 'MRAS'],
    'HRAS': ['NRAS', 'KRAS', 'MRAS'],
    'MRAS': ['NRAS', 'HRAS', 'KRAS'],
    'ARAF': ['BRAF', 'RAF1'],
    'BRAF': ['ARAF', 'RAF1'],
    'RAF1': ['ARAF', 'BRAF'],
    'MAP2K1': ['MAP2K2'],
    'MAP2K2': ['MAP2K1'],
    'MAPK1': ['MAPK3'],
    'MAPK3': ['MAPK1'],
}

# Structures present in Inputs/Pocket-Structures but used by no figure. 9O0R
# is a second KRAS/adagrasib co-crystal; 6UT0 is the one the pocket analysis
# uses. It is listed so build_inputs measures the whole directory.
EXTRA_STRUCTURES = [
    ('KRAS', '9O0R', 'Adagrasib', 'A1B7W', 'A', 'A', 'KRASi-DMSO-Z'),
]

# Pocket-figure palettes.
POCKET_PDB_COLORS = {'InPocket': '#f7c297', 'OutPocket': '#ffecb8', 'Unresolved': '#888888'}
POCKET_AF_COLORS = {'InPocket': '#C19EE0', 'OutPocket': '#7CCDF4'}
POCKET_EXTRA_COLORS = {'OutsideNode': '#D6D6D6', 'InParalogs': '#90E0EF', 'OtherGenes': '#FFCCCC'}
EDITOR_PALETTE = {'ABE': '#90E0EF', 'CBE': '#FFCCCC'}
EDITOR_MARKERS = {'ABE': 'o', 'CBE': 'D'}

# ---------------------------------------------------------------------------
# cBIOPORTAL
# ---------------------------------------------------------------------------

CBIOPORTAL_FILES = {
    'SHC1': 'SHC1_ENST00000448116.tsv',
    'GRB2': 'GRB2_ENST00000316804.tsv',
    'PTPN11': 'PTPN11_ENST00000351677.tsv',
    'KRAS': 'KRAS_ENST00000311936.tsv',
    'NRAS': 'NRAS_ENST00000369535.tsv',
    'HRAS': 'HRAS_ENST00000311189.tsv',
    'MRAS': 'MRAS_ENST00000423968.tsv',
    'PPP1CA': 'PPP1CA_ENST00000376745.tsv',
    'PPP1CB': 'PPP1CB_ENST00000395366.tsv',
    'PPP1CC': 'PPP1CC_ENST00000335007.tsv',
    'SHOC2': 'SHOC2_ENST00000369452.tsv',
    'ARAF': 'ARAF_ENST00000377045.tsv',
    'BRAF': 'BRAF_ENST00000646891_ENST00000288602.tsv',
    'RAF1': 'RAF1_ENST00000251849.tsv',
    'KSR1': 'KSR1_ENST00000644974_ENST00000319524.tsv',
    'KSR2': 'KSR2_ENST00000339824.tsv',
    'MAP2K1': 'MAP2K1_ENST00000307102.tsv',
    'MAP2K2': 'MAP2K2_ENST00000262948.tsv',
    'DUSP4': 'DUSP4_ENST00000240100.tsv',
    'DUSP6': 'DUSP6_ENST00000279488.tsv',
    'MAPK1': 'MAPK1_ENST00000215832.tsv',
    'MAPK3': 'MAPK3_ENST00000263025.tsv',
}

CBIOPORTAL_MUTATION_TYPES = [
    'Missense_Mutation', 'Nonsense_Mutation',
    'Frame_Shift_Del', 'In_Frame_Del', 'Frame_Shift_Ins', 'In_Frame_Ins',
]

# Lollipop axis limits, in GENES order. min_score is the floor of the
# (inverted) cBioPortal count panel; break is the y-range hidden by the axis
# break, or None for an unbroken axis.
LOLLIPOP_MIN_SCORES = dict(zip(GENES, [
    -5, -4, -8, -580, -220, -56, -4,
    -4, -4, -3, -5, -12, -600, -15,
    -13, -10, -20, -6, -5, -4, -20, -4,
]))

LOLLIPOP_MAX_SCORES = dict(zip(GENES, [
    6, 6, 28, 24, 18, 6, 8,
    8, 8, 6, 12, 30, 12, 18,
    8, 6, 24, 10, 12, 6, 18, 12,
]))

LOLLIPOP_BREAK_SCORES = dict(zip(GENES, [
    None, None, None, (-550, -100), (-200, -35), (-16, -8), None,
    None, None, None, None, (-9, -6), (-500, -25), (-12, -4),
    (-11, -10), (-8, -5), (-16, -8), None, None, None, (-18, -4), None,
]))

# Distance-panel ticks: deeper pockets need a taller axis.
LOLLIPOP_LINE_YTICKS_TALL = [0, 45, 90]
LOLLIPOP_LINE_YTICKS_SHORT = [0, 25, 50]
LOLLIPOP_TALL_GENES = ['PTPN11', 'ARAF', 'BRAF', 'RAF1']

LOLLIPOP_POS_COLORS = [
    COLOR_MAP['PTPN11'], COLOR_MAP['KRAS'], COLOR_MAP['BRAF'],
    COLOR_MAP['MAP2K1'], COLOR_MAP['MAPK1'],
]
LOLLIPOP_NEG_COLORS = ['#888888']
LOLLIPOP_POS_THRESHOLD = 2
LOLLIPOP_NEG_THRESHOLD = 0

# ---------------------------------------------------------------------------
# DMS DATASETS
# ---------------------------------------------------------------------------

DMS_FILES = {
    'MAPK1': DMS_DIR / 'Brenan-Johannessen-CellReports2016.csv',
    'KRAS': DMS_DIR / 'full_KRAS_DMS_integrated_dataset_from_Jason.csv',
    'PTPN11': DMS_DIR / 'Jiang-Shah-NatureComms2025.csv',
}

# Column carrying the fitness/abundance score in each DMS dataset. Negative
# means depleted in all three, so mutation scores are aggregated with min().
DMS_SCORE_COLUMNS = {
    'MAPK1': ('ERK2_Mutant', 'LFC (ETP vs. DOX)'),
    'KRAS': ('VAR', 'HCC827G12D_LFC'),
    # PTPN11 is a position x amino-acid matrix rather than a long table.
}

# Citation carried through to the consolidated DMS sheet, so a row states
# which published dataset scored it.
DMS_DATASETS = {
    'MAPK1': 'Brenan et al., Cell Reports 2016 (ERK2)',
    'KRAS': 'Integrated KRAS dataset, HCC827 G12D',
    'PTPN11': 'Jiang et al., Nature Communications 2025 (SHP2)',
}

DMS_GENES = ['PTPN11', 'KRAS', 'MAPK1']

# How a guide making several scorable edits is summarized into one score.
#
# 'min' takes the most depleting edit, on the assumption that the strongest
# effect dominates the phenotype. 'mean' averages them, which is the right
# model when the edits act additively. This is not a cosmetic choice: 585 of
# the 1,604 scored guides (36.5%) make more than one scorable edit, up to six.
DMS_AGGREGATE = 'mean'

# ---------------------------------------------------------------------------
# GMM TRAJECTORY
# ---------------------------------------------------------------------------

TRAJECTORY_VALID_CONDITIONS = [
    'SHP2i-DMSO-Z', 'KRASi-DMSO-Z', 'RAFi-DMSO-Z', 'Trametinib-DMSO-Z', 'ERKi-DMSO-Z',
]
TRAJECTORY_VALID_LOWER_CUTOFF = 2
TRAJECTORY_VALID_UPPER_CLIP = 10
TRAJECTORY_VALID_N_COMPONENTS = 14
TRAJECTORY_VALID_ROOT = 13
TRAJECTORY_VALID_SUBSETS = [
    [13, 2, 1, 7, 11],
    [13, 2, 1, 6],
    [13, 2, 1, 4],
    [4, 9, 5],
    [4, 8, 10, 0],
    [4, 8, 10, 3],
    [13, 2, 1, 6, 7, 11, 4, 8, 10, 9, 5, 0, 12, 3],
]

TRAJECTORY_MEKI_CONDITIONS = [
    'DMSO-Day0-Z',
    'Avutametinib-DMSO-Z', 'GDC0623-DMSO-Z', 'Trametiglue-DMSO-Z', 'Trametinib-DMSO-Z',
    'Selumetinib-DMSO-Z', 'Mirdametinib-DMSO-Z', 'Cobimetinib-DMSO-Z', 'MAP855-DMSO-Z',
    'SHP2i-DMSO-Z', 'KRASi-DMSO-Z',
]
TRAJECTORY_MEKI_LABELS = [
    'DMSO-Day0-Z',
    'Avutametinib', 'GDC0623', 'Trametiglue', 'Trametinib',
    'Selumetinib', 'Mirdametinib', 'Cobimetinib', 'MAP855',
    'SHP2i', 'KRASi',
]
TRAJECTORY_MEKI_LOWER_CUTOFF = 4
TRAJECTORY_MEKI_UPPER_CLIP = 10
TRAJECTORY_MEKI_N_COMPONENTS = 10
TRAJECTORY_MEKI_ROOT = 3

RANDOM_STATE = 0
