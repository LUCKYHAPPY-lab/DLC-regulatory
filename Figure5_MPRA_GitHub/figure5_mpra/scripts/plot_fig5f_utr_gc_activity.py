#!/usr/bin/env python3
"""GC content vs 3'UTR MPRA activity (Fig. 5f style).

Display rule — LINEAR RESCALE of the lower tail
------------------------------------------------
Values with log2FoldChange < THRESH are linearly rescaled into
[THRESH - RESCALE_RANGE, THRESH]:

    y_plot = THRESH - RESCALE_RANGE * (THRESH - y) / (THRESH - min_obs)

where min_obs is the smallest observed value among the rescaled points.
The mapping is monotone and fixes both ends:
    THRESH  -> THRESH                 (threshold points stay put)
    min_obs -> THRESH - RESCALE_RANGE (extreme points land on the bottom edge)
Rank order among the rescaled points is preserved.

This only compresses the empty region that the extreme tail would otherwise
open at the bottom of the panel. It is a pure display transform.

Statistics rule
---------------
Pearson correlation and the fitted regression line are computed on the RAW
(untransformed) log2FoldChange, so the reported r and the drawn line always
agree with each other and with the untransformed data.

Outputs (per cell) : <CELL>_GC_vs_log2FC_rescale.pdf
                     SourceData_<CELL>_GC_vs_log2FC_rescale.xlsx
                     SourceData_<CELL>_GC_vs_log2FC_rescale.tsv
                     figure_legend_<CELL>.txt
Plus comparison_<CELL>.png — rescale vs floor vs raw, for review only.
"""

from __future__ import annotations

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr

parser = argparse.ArgumentParser(description="Plot 3'UTR GC content versus MPRA activity.")
parser.add_argument("--pk15", required=True)
parser.add_argument("--psc", required=True)
parser.add_argument("--pearson-helper", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()
from pathlib import Path
import importlib.util

helper_path = Path(args.pearson_helper).resolve()
if not helper_path.is_file():
    raise FileNotFoundError(f"Pearson helper not found: {helper_path}")
spec = importlib.util.spec_from_file_location("pearson_p_exact", helper_path)
if spec is None or spec.loader is None:
    raise ImportError(f"Cannot load Pearson helper: {helper_path}")
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
pearson_p_value = helper.pearson_p_value

OUT_DIR = args.outdir
os.makedirs(OUT_DIR, exist_ok=True)
INPUTS = {"PK15": args.pk15, "PSC": args.psc}
THRESH = -8.0
RESCALE_RANGE = 0.5
CELLS = ["PK15", "PSC"]


def load(cell: str) -> pd.DataFrame:
    src = INPUTS[cell]
    df = pd.read_excel(src)
    n_raw = len(df)
    df = df.dropna(subset=["GC_percent", "log2FoldChange"])
    n_drop = len(df)
    df = df.groupby("ID", as_index=False).mean(numeric_only=True)
    if df["GC_percent"].max() <= 1.0:
        df["GC_percent"] = df["GC_percent"] * 100
    df["log2FoldChange"] = df["log2FoldChange"].astype(float)
    df["n_rows_input"] = n_raw
    df["n_rows_after_dropna"] = n_drop
    return df


def rescale_y(values, thresh=THRESH, rng=RESCALE_RANGE):
    """Linearly map values below thresh into [thresh - rng, thresh]."""
    v = values.astype(float).copy()
    m = v < thresh
    if m.any():
        mn = v[m].min()
        v[m] = thresh - rng * (thresh - v[m]) / (thresh - mn)
    return v


def floor_y(values, thresh=THRESH):
    """For comparison, place values below thresh at thresh."""
    return np.maximum(values.astype(float), thresh)


def fmt_p(p):
    return "p < 2.2e-16" if p < 2.2e-16 else f"p = {p:.2e}"


for cell in CELLS:
    df = load(cell)
    x = df["GC_percent"].to_numpy(float)
    y_raw = df["log2FoldChange"].to_numpy(float)

    y_plot = rescale_y(y_raw)
    y_floor = floor_y(y_raw)

    # Calculate statistics on raw values.
    r, p = pearsonr(x, y_raw)
    slope, intercept = np.polyfit(x, y_raw, 1)
    # Use the helper to retain the order of magnitude if the SciPy P value underflows.
    p_exact, log10p, p_text = pearson_p_value(r, len(df))

    mask = y_raw < THRESH
    n_below = int(mask.sum())
    min_obs = float(y_raw[mask].min()) if n_below else np.nan
    n_total = len(df)

    # ---------- Main plot: linear tail rescaling ----------
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(x, y_plot, s=5, alpha=0.15, color="black", rasterized=True, zorder=1)
    sns.kdeplot(x=x, y=y_plot, fill=True, levels=35, thresh=0.05,
                cmap="YlOrRd", alpha=0.85, ax=ax, zorder=2)
    xs = np.linspace(x.min(), x.max(), 200)
    ax.plot(xs, slope * xs + intercept, color="#0047AB", lw=2.5, zorder=3)

    ax.set_ylim(bottom=THRESH - RESCALE_RANGE)
    ax.set_xlabel("GC content (%)", fontsize=11, labelpad=8)
    ax.set_ylabel("log2FoldChange of luciferase activity", fontsize=11, labelpad=8)
    ax.set_title(f"Pearson correlation of GC content and luciferase activity ({cell})",
                 fontsize=12, pad=12)
    ax.text(0.05, 0.95, f"r = {r:.2f}\n{fmt_p(p)}",
            transform=ax.transAxes, fontsize=12, fontweight="semibold", va="top")
    sns.despine()
    plt.tight_layout()

    # Place the note below the x-axis to keep the points visible.
    fig.text(0.5, -0.035,
             f"Values with log₂FC < {THRESH:g} were linearly rescaled to "
             f"[{THRESH - RESCALE_RANGE:g}, {THRESH:g}] "
             f"(n = {n_below} of {n_total}). Statistics use untransformed values.",
             ha="center", va="top", fontsize=7.0, color="#444444")

    out_pdf = os.path.join(OUT_DIR, f"{cell}_GC_vs_log2FC_rescale.pdf")
    plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close()

    # ---------- source data ----------
    points = pd.DataFrame({
        "ID": df["ID"],
        "GC_percent": df["GC_percent"],
        "log2FoldChange_raw": y_raw,
        "log2FoldChange_plotted": y_plot,
        "rescaled": mask,
    }).sort_values("GC_percent").reset_index(drop=True)

    stats = pd.DataFrame({
        "Statistic": [
            "N_points_plotted", "N_rows_input", "N_rows_after_dropna",
            "Y_display_transform", "Rescale_threshold", "Rescale_range",
            "Rescale_min_observed", "Rescale_formula",
            "N_points_rescaled",
            "Pearson_r", "Pearson_p", "Pearson_log10P", "Pearson_on",
            "Regression_slope", "Regression_intercept", "Regression_on",
            "GC_percent_min", "GC_percent_max",
            "log2FoldChange_min", "log2FoldChange_max",
        ],
        "Value": [
            n_total, int(df["n_rows_input"].iloc[0]), int(df["n_rows_after_dropna"].iloc[0]),
            "linear rescale of lower tail",
            THRESH, RESCALE_RANGE, min_obs,
            "y_plotted = THRESH - RESCALE_RANGE*(THRESH-y)/(THRESH-min_obs) for y < THRESH; "
            "y_plotted = y otherwise",
            n_below,
            r, p_text, round(-log10p, 6), "raw log2FoldChange",
            slope, intercept, "raw log2FoldChange",
            df["GC_percent"].min(), df["GC_percent"].max(),
            y_raw.min(), y_raw.max(),
        ],
    })

    readme = pd.DataFrame({"Field": [
        "Figure", "Panel", "Test", "x", "y",
        "Deduplication", "Missing values",
        "y-axis display transform", "Statistics",
        "Density layer", "Regression line", "Input file",
    ], "Description": [
        f"{cell}_GC_vs_log2FC_rescale.pdf",
        "Single panel: scatter of per-sequence GC content vs luciferase activity, "
        "overlaid with a 2D KDE density surface and a linear regression line",
        "Pearson correlation on the RAW log2FoldChange (scipy.stats.pearsonr); "
        "Pearson_p is the exact two-sided p computed in log space (pearson_p_exact.py)"
        " — scipy double precision underflows to 0 below ~1e-308 — and"
        " Pearson_log10P = -log10(Pearson_p)",
        "GC_percent — GC content of the 101 bp oligo insert (%)",
        f"log2FoldChange — {cell} luciferase activity from MPRAmodel",
        "Averaged per ID (pandas groupby mean over numeric columns), so each sequence contributes exactly one point",
        "Rows with missing GC_percent or log2FoldChange dropped before averaging",
        f"DISPLAY ONLY: values with log2FoldChange < {THRESH:g} are linearly rescaled into "
        f"[{THRESH - RESCALE_RANGE:g}, {THRESH:g}] as "
        f"y' = {THRESH:g} - {RESCALE_RANGE:g}*({THRESH:g}-y)/({THRESH:g}-min_obs), "
        f"where min_obs = {min_obs} is the smallest observed value among those points "
        f"(n = {n_below} of {n_total}). The threshold {THRESH:g} maps to itself and "
        f"min_obs maps to {THRESH - RESCALE_RANGE:g}, so rank order is preserved within the tail. "
        f"Raw range is [{y_raw.min():.3f}, {y_raw.max():.3f}].",
        "Pearson r and the fitted regression are computed on the untransformed log2FoldChange; "
        "neither uses the rescaled values",
        "seaborn kdeplot, levels=35, thresh=0.05, cmap=YlOrRd — derived from the plotted points, no extra data",
        "OLS of raw log2FoldChange on GC_percent (numpy.polyfit); same fit as the OLS underlying the reported r",
        f"3UTR_all_info_annotion_encodeid_pass_with_GC_MFE_activity_merged_{cell}.xlsx",
    ]})

    xlsx = os.path.join(OUT_DIR, f"SourceData_{cell}_GC_vs_log2FC_rescale.xlsx")
    with pd.ExcelWriter(xlsx) as w:
        readme.to_excel(w, sheet_name="README", index=False)
        stats.to_excel(w, sheet_name="Fig_stats", index=False)
        points.to_excel(w, sheet_name="Fig_points", index=False)
    points.to_csv(os.path.join(OUT_DIR, f"SourceData_{cell}_GC_vs_log2FC_rescale.tsv"),
                  sep="\t", index=False, float_format="%.6g")

    # ---------- Figure legend text ----------
    legend = (
        "Correlation between oligonucleotide GC content and MPRA transcriptional activity "
        f"for all 3′UTR oligonucleotides in {cell} cells. The x-axis shows GC content percentage, "
        "and the y-axis shows MPRA transcriptional activity, measured as log2FC activity. "
        f"For display, log2FC values below {THRESH:g} were linearly rescaled to the interval "
        f"[{THRESH - RESCALE_RANGE:g}, {THRESH:g}] (n = {n_below} of {n_total:,} "
        "oligonucleotides); Pearson correlation and the blue linear regression line were "
        "computed on untransformed log2FC values. The blue line represents the linear regression fit."
    )
    with open(os.path.join(OUT_DIR, f"figure_legend_{cell}.txt"), "w", encoding="utf-8") as fh:
        fh.write(legend + "\n")

    # ---------- Compare rescaled, clipped and raw displays ----------
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), sharex=True)
    variants = [
        (y_plot, "linear rescale to [−8.5, −8]  (chosen)",
         f"n(<−8)={n_below} rescaled", THRESH - RESCALE_RANGE),
        (y_floor, "floor at −8", f"n(<−8)={n_below} set to −8", THRESH - 0.55),
        (y_raw, "raw  (no transform)", f"n(<−8)={n_below} shown as is", None),
    ]
    for ax, (yv, title, sub, ybot) in zip(axes, variants):
        ax.scatter(x, yv, s=5, alpha=0.15, color="black", rasterized=True, zorder=1)
        sns.kdeplot(x=x, y=yv, fill=True, levels=35, thresh=0.05,
                    cmap="YlOrRd", alpha=0.85, ax=ax, zorder=2)
        ax.plot(xs, slope * xs + intercept, color="#0047AB", lw=2.2, zorder=3)
        ax.set_title(title, fontsize=9.5)
        ax.set_xlabel("GC content (%)")
        if ybot is not None:
            ax.set_ylim(bottom=ybot)
        ax.text(0.05, 0.95, f"r = {r:.2f}\n{sub}", transform=ax.transAxes, fontsize=8, va="top")
        sns.despine(ax=ax)
    axes[0].set_ylabel("log2FoldChange")
    fig.suptitle(f"{cell} — display transform comparison (all fits on raw y)", fontsize=11, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, f"comparison_{cell}.png"), dpi=130, bbox_inches="tight")
    plt.close()

    print(f"[{cell}] points={n_total}  n_rescaled={n_below}  min_obs={min_obs:.4f}")
    print(f"       Pearson r={r:.6f}  slope={slope:.6f}  intercept={intercept:.6f}")
    print(f"       wrote {out_pdf}")
