"""Loading and reshaping for the manuscript figures outside the screens.

Competitive growth, CellTiter-Glo viability, inhibitor and kinase dose
response, BaF3 IL3-independent growth, and a meta-analysis of published KRAS
deep mutational scanning. Eleven separate analysis projects; one tidy input
workbook, ``Other_figures_JW/inputs/other_figures_raw_data.xlsx``, one tab per
experiment plus the key tabs the competitive growth tabs are labelled from.
The workbook's ``README`` tab records what each tab holds.

The reshaping that does not depend on the experiment was done when the
workbook was built: 384-well plates flattened well by well and joined to their
layout keys, Prism-style wide sheets melted to long format with a trailing
``*`` becoming an ``exclude`` flag, FCS Express exports parsed out of their
filenames. What remains - normalisation, gate correction, unit conversion,
z-scoring, curve fitting - is assay-specific and lives here.

``other_charts`` draws what this module computes. Nothing here prints or
plots; every function returns what the notebook reports.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from . import config as cfg


def read_tab(name: str, path=None) -> pd.DataFrame:
    """Read one tab of the input workbook."""
    return pd.read_excel(path or cfg.OTHER_FIGURES_XLSX, sheet_name=name)


# ---------------------------------------------------------------------------
# COMPETITIVE GROWTH
# ---------------------------------------------------------------------------
#
# GFP-labelled edited cells are mixed with unlabelled parental cells; the GFP+
# fraction over time reports the relative fitness of the edit under each drug.


def load_competitive_growth(tab: str) -> pd.DataFrame:
    """One competitive growth tab, with sgRNA and drug display labels attached.

    The key tabs are indexed by experiment without the ``CG_`` prefix, and
    both merges are validated many-to-one so a duplicated key row cannot
    silently multiply the readouts.
    """
    df = read_tab(tab)
    name = tab.removeprefix("CG_")
    sgrna_key = read_tab("CG_sgRNA_key").query("Experiment == @name")
    drug_key = read_tab("CG_drug_key").query("Experiment == @name")

    df["sgRNA"] = df["sgRNA"].astype(str)
    sgrna_key = sgrna_key.assign(sgRNA=sgrna_key["sgRNA"].astype(str))

    df = df.merge(sgrna_key[["BE", "sgRNA", "label"]]
                  .rename(columns={"label": "sgRNA_label"}),
                  on=["BE", "sgRNA"], how="left", validate="many_to_one")
    df = df.merge(drug_key[["Condition", "label", "order"]]
                  .rename(columns={"label": "drug_label", "order": "drug_order"}),
                  on="Condition", how="left", validate="many_to_one")
    df = df.rename(columns={"M4-2 % Parent": "GFP+ Fraction"})

    assert df["sgRNA_label"].isna().sum() == 0, f"unmatched sgRNAs in {name}"
    return df


def expand_baseline(df: pd.DataFrame, baseline_day: int = 0) -> pd.DataFrame:
    """Copy the shared day-0 read into every drug arm.

    Day 0 is measured once, before the culture is split across drugs, so those
    rows have a blank ``Condition`` and the same baseline value has to start
    each condition's growth curve. Copied rows are flagged ``pseudo_day0``.
    """
    baseline = df[df["Day"] == baseline_day]
    treated = df[df["Day"] != baseline_day]
    conditions = treated[["Condition", "drug_label", "drug_order"]].drop_duplicates()
    expanded = (baseline.drop(columns=["Condition", "drug_label", "drug_order"])
                .merge(conditions, how="cross").assign(pseudo_day0=True))
    return pd.concat([expanded, treated.assign(pseudo_day0=False)], ignore_index=True)


def fold_change_vs_nt(df: pd.DataFrame, baseline_day: int = 0,
                      endpoint_day: int = 12, nt_sgrna: str = "NT") -> pd.DataFrame:
    """Endpoint / baseline GFP+ fraction, normalised to the NT mean.

    Divided by the matched day-0 read, then by the mean of the non-targeting
    controls for the same base editor and drug, then averaged over replicates.
    Below 1 means the edited cells were depleted relative to NT; above 1 means
    enriched.
    """
    baseline = (df[df["Day"] == baseline_day][["BE", "sgRNA", "GFP+ Fraction"]]
                .rename(columns={"GFP+ Fraction": "GFP+ Fraction d0"}))
    endpoint = df[df["Day"] == endpoint_day].merge(baseline, on=["BE", "sgRNA"],
                                                   how="left")
    endpoint["FC from d0"] = endpoint["GFP+ Fraction"] / endpoint["GFP+ Fraction d0"]

    nt_mean = (endpoint[endpoint["sgRNA"] == nt_sgrna]
               .groupby(["BE", "Condition"], observed=True)["FC from d0"].mean()
               .rename("FC from d0 NT mean").reset_index())
    endpoint = endpoint.merge(nt_mean, on=["BE", "Condition"], how="left")
    endpoint["FC from d0 compared to NT"] = (endpoint["FC from d0"]
                                             / endpoint["FC from d0 NT mean"])

    return (endpoint.groupby(["BE", "Condition", "drug_label", "sgRNA", "sgRNA_label"],
                             observed=True)["FC from d0 compared to NT"]
            .mean().reset_index())


# ---------------------------------------------------------------------------
# CELLTITER-GLO VIABILITY
# ---------------------------------------------------------------------------
#
# 384-well viability readouts, normalised two ways. Dividing by the no-Drug-A
# wells within each series starts every curve at 1, so the shapes are directly
# comparable. Dividing by the wells with neither drug, per plate, keeps the
# vertical offset between series, which is the Drug B effect.


def normalize_to_zero_dose(df: pd.DataFrame, dose_col: str = "Drug A nM",
                           value_col: str = "Value",
                           group_cols=("category", "Plate"),
                           out_col: str = "Value_norm") -> pd.DataFrame:
    """Divide by the mean value at zero Drug A within each group."""
    df = df.copy()
    df[dose_col] = pd.to_numeric(df[dose_col], errors="coerce")
    baseline = (df.loc[df[dose_col] == 0].groupby(list(group_cols))[value_col]
                .mean().rename("_baseline"))
    df = df.merge(baseline, on=list(group_cols), how="left")
    df[out_col] = df[value_col] / df["_baseline"]
    return df.drop(columns="_baseline")


def normalize_to_zero_both(df: pd.DataFrame, dose_a: str = "Drug A nM",
                           dose_b: str = "Drug B nM", value_col: str = "Value",
                           group_cols=("Plate",),
                           out_col: str = "Value_norm_zeroAB") -> pd.DataFrame:
    """Divide by the mean of the untreated wells (both drugs zero) per plate."""
    df = df.copy()
    for col in (dose_a, dose_b, value_col):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    untreated = df.loc[df[dose_a].eq(0) & df[dose_b].eq(0)]
    if untreated.empty:
        raise ValueError("No wells with both drugs at zero to normalise against.")
    baseline = untreated.groupby(list(group_cols))[value_col].mean().rename("_baseline")
    df = df.merge(baseline, on=list(group_cols), how="left")
    df[out_col] = df[value_col] / df["_baseline"]
    return df.drop(columns="_baseline")


# ---------------------------------------------------------------------------
# DOSE RESPONSE
# ---------------------------------------------------------------------------


def load_dose_response(tab: str) -> pd.DataFrame:
    """A dose-response tab, dropping the values flagged for exclusion.

    Values carrying a trailing ``*`` in the original Prism sheets became
    ``exclude = yes`` when the workbook was built; they are left out of the
    fits, as they were in the published figures.
    """
    df = read_tab(tab)
    return df[df["exclude"] == "no"].copy()


def four_param_logistic(x, bottom, top, log_ec50, hill_slope):
    """Standard 4PL curve: response as a function of log concentration."""
    x = np.asarray(x, dtype=float)
    return bottom + (top - bottom) / (1 + 10 ** ((log_ec50 - np.log10(x)) * hill_slope))


def fit_dose_response(df: pd.DataFrame, conc_col: str, response_col: str,
                      group_col: str, group_order) -> dict:
    """Fit one 4PL curve per group. Returns ``{group: fitted parameters}``.

    Zero and negative concentrations are dropped: the curve is a function of
    log concentration. A group whose fit does not converge maps to None rather
    than aborting the panel.
    """
    data = df[[conc_col, response_col, group_col]].dropna()
    data = data[data[conc_col] > 0]

    fits = {}
    for group in group_order:
        sub = data[data[group_col] == group]
        if sub.empty:
            continue
        x = sub[conc_col].to_numpy(float)
        y = sub[response_col].to_numpy(float)
        p0 = [np.nanmin(y), np.nanmax(y), np.median(np.log10(x)), 1.0]
        try:
            popt, pcov = curve_fit(four_param_logistic, x, y, p0=p0, maxfev=10000)
            bottom, top, log_ec50, hill = popt
            fits[group] = {"bottom": bottom, "top": top, "logEC50": log_ec50,
                           "EC50": 10 ** log_ec50, "hill_slope": hill,
                           "covariance": pcov}
        except RuntimeError:
            fits[group] = None
    return fits


def ecx_table(fits: dict, ecx_values=(50, 90, 10)) -> pd.DataFrame:
    """ECx values from fitted parameters, x% of the way from bottom to top."""
    rows = []
    for group, fit in fits.items():
        row = {"group": group}
        if fit is None:
            row.update({f"EC{x}": np.nan for x in ecx_values})
        else:
            for x in ecx_values:
                fraction = x / 100
                log_ecx = (fit["logEC50"]
                           - (1 / fit["hill_slope"]) * np.log10(1 / fraction - 1))
                row[f"EC{x}"] = 10 ** log_ecx
            row.update({k: fit[k] for k in ("bottom", "top", "hill_slope")})
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# BAF3 IL3-INDEPENDENT GROWTH
# ---------------------------------------------------------------------------
#
# BaF3 cells need IL3 unless a transduced construct drives proliferation on
# its own, so cell count over time is a direct readout of how activating each
# variant is.

# Sampled volume in ml, and the dilution factor the sample was read at.
BAF3_SAMPLE_VOLUME_ML = 0.02
BAF3_DILUTION = 2
# Below this many events a gate is treated as having failed.
BAF3_MIN_EVENTS = 10


def load_baf3_jw311() -> pd.DataFrame:
    """The JW311 counts, gate-corrected and converted to cells per ml.

    Gate 2 is the live-cell gate. A handful of wells failed it but have a
    usable "shifted singles" count, which is substituted for those wells only.
    """
    baf3 = read_tab("BaF3_JW311")
    baf3["Cells"] = np.where(
        (baf3["Gate 2 events"] < BAF3_MIN_EVENTS)
        & (baf3["Shifted singles events"] > BAF3_MIN_EVENTS),
        baf3["Shifted singles events"], baf3["Gate 2 events"])
    baf3["Cells/ml"] = (baf3["Cells"] / BAF3_SAMPLE_VOLUME_ML) * BAF3_DILUTION
    return baf3


def baf3_subset(baf3: pd.DataFrame, labels, gene: str = "KRAS",
                promoter: str = "Strong",
                include_no_il3: bool = False) -> pd.DataFrame:
    """Rows for one construct set, optionally with the no-IL3 floor appended."""
    sub = baf3[(baf3["Gene"] == gene) & (baf3["Promoter"] == promoter)
               & (baf3["label"].isin(labels))]
    if include_no_il3:
        sub = pd.concat([baf3[baf3["label"] == "no IL3"], sub])
    return sub


# ---------------------------------------------------------------------------
# PUBLISHED DMS META-ANALYSIS
# ---------------------------------------------------------------------------


def zscore_columns(df: pd.DataFrame, id_col: str = "mutation_id") -> pd.DataFrame:
    """Column-wise z-score of every numeric column, so assays become comparable.

    The published screens are on very different scales; z-scoring each assay
    across all variants puts them on one axis. A column with no spread scores
    NaN rather than dividing by zero.
    """
    out = df.copy()
    for col in out.select_dtypes(include=[np.number]).columns:
        if col == id_col:
            continue
        std = out[col].std(skipna=True)
        out[col] = (np.nan if (pd.isna(std) or std == 0)
                    else (out[col] - out[col].mean(skipna=True)) / std)
    return out
