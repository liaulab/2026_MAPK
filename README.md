# 2026_MAPK
All code used for MAPK submission analysis and figures

Three notebooks reproduce every figure and the library they were screened with. 
Each reads Excel tables available in the manuscript's Supplemental Materials and writes to `Outputs/`.

Plotting primitives come from [`be_scan`](https://github.com/liaulab/be-scan). 
Everything specific to this manuscript lives in `code/`.

---

## Quick start

### Google Colab

Open any of the three notebooks and run section 0. 
It installs the pinned dependency set and `be_scan` into the runtime, so the environment is built from scratch on every run and nothing needs to persist between sessions.

Section 0 needs this repository present. Either set `REPO_URL` in that cell to clone it, or mount Drive and point `REPO_DIR` at the copy there.

Then run the notebook top to bottom.

### Locally

```bash
conda env create -f environment.yml
conda activate mapk-analysis
jupyter lab
```

Skip section 0 when running locally — the conda environment already provides everything.

---

## Notebooks

| Notebook | What it does | Reads | Writes |
|---|---|---|---|
| **`MAPK_Analysis.ipynb`** | Every screen figure. Scatterplots, GOF-vs-validation, inhibitor pocket hits, splice-sgRNA heatmaps, the gene-inhibitor Sankey, the GMM clustering and trajectory, hyperactivation vs DMS, cBioPortal plots | TableS1, S3, S4, S5 | `Outputs/00-Summaries/` to `Outputs/11-Boxplots/` |
| **`MAPK_LibraryDesign.ipynb`** | Designs the tiling ABE/CBE sgRNA library from transcript mapping through sgRNA design, filtering, BEhive efficiency and FlashFry specificity | TableS6 | `Outputs/12-Library-Design/` |
| **`MAPK_Other.ipynb`** | Non-screen figures: competitive growth, CellTiter-Glo viability, inhibitor and kinase dose–response, BaF3 IL3-independent growth, published KRAS DMS meta-analysis | TableS7 | `Outputs/13-Other-Figures/` |

For all 3 notebooks, section 0 sets up the environment, section 1 imports, section 2 holds the shared configuration, and each later section is one analysis. 
In `MAPK_Analysis.ipynb` and `MAPK_Other.ipynb` those later sections are independent and can be run selectively once sections 0–3 have run. 
`MAPK_LibraryDesign.ipynb` runs in order, each step feeding the next. 

---

## Input tables

| Table | Sheets | Contents |
|---|---|---|
| `TableS1-ScreenData.xlsx` | 11 | All screen data plus the `edit-site` guide position key |
| `TableS3-cBioPortal-Annotations.xlsx` | 23 | Patient mutation tables with one sheet per gene plus a `Transcripts` summary |
| `TableS4-DMS-Annotations.xlsx` | 1 | All three published deep mutational scanning datasets |
| `TableS5-PocketDistances-Annotations.xlsx` | 19 | Minimum heavy-atom distance from the inhibitor to every residue with one sheet per structure plus a `Structures` summary |
| `TableS6-LibraryDesign.xlsx` | 7 | Everything the library design needs: base-editor parameters, target-gene transcripts, one sgRNA design table per editor, BEhive efficiency scores, FlashFry specificity scores, exon sequences |
| `TableS7-OtherFiguresData.xlsx` | 16 | One tidy tab per non-screen experiment plus the competitive-growth label keys |

---

## Repository layout

```
MAPK_Analysis.ipynb               screen figures        -> Outputs/00-11
MAPK_LibraryDesign.ipynb          the sgRNA library     -> Outputs/12-Library-Design
MAPK_Other.ipynb                  non-screen figures    -> Outputs/13-Other-Figures

code/                        analysis code not in be_scan
  config.py                  every shared variable (see below)
  data.py                    loading and reshaping TableS1
  analysis.py                the screen analysis itself
  build_inputs.py            builds TableS3, S4, S5 and S6 from their raw files
  figure_options.py          option dataclasses, in be_scan's style
  structures.py              ligand-to-residue distances, pocket assignment
  charts.py                  bar / histogram / pie / heatmap figures
  boxplot.py                 control boxplot with its significance bracket
  lollipop.py                split lollipop with a pocket-distance panel
  trajectory.py              GMM clustering with an MST pseudotime trajectory
  dms.py                     deep mutational scanning reference data
  sankey.py                  inhibitor-to-gene validated-hit Sankey
  library_design.py          the seven-step sgRNA library design pipeline
  other_figures.py           loading and reshaping for the non-screen figures
  other_charts.py            the five figure types those figures are drawn with
  interactive.py             Plotly twin of every screen figure
  website.py                 static site generator for docs/

requirements.txt             Colab dependency set
```

---

## Shared variables

Everything used by more than one figure lives in **`code/config.py`**, once, and every notebook reads it from there. 

| Group | Holds |
|---|---|
| Paths | Repository root, the seven tables, `Inputs/`, and every `Outputs/` subdirectory |
| Genes and colors | Pathway gene order, the gene palette, mutation-type and editor palettes |
| Protein lengths and domains | Per-gene sequence length and domain boundaries for the position axes |
| Screen sheets and condition columns | Which sheet and which Z-score columns belong to each screen |
| Hit-calling cutoffs | GOF, validation and hyperactivation thresholds |
| Per-gene axis ranges | y-ranges and tick counts for every scatterplot and lollipop panel |
| Structures | Co-crystal and AlphaFold structures, their ligands, chains and paralog groupings |
| DMS datasets | Dataset citations, gene list, and how a guide with several edits is summarized |
| GMM trajectory | Conditions, component counts, roots and cluster subsets for both trajectories |
| Library design | Base editors, the Ensembl endpoint, guide filters and the published counts a run is checked against |
| Other figures | The shared figure style, palette and jitter seed |

---

## The interactive site

`docs/` is the finished site. Open `docs/index.html` in a browser to read it locally.
This is also published on its own [`Github Page`](https://github.com/liaulab/be-scan](https://liaulab.github.io/2026_MAPK/ ).
