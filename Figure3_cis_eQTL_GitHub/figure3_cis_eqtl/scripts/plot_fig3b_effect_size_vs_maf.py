#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.stats import spearmanr, gaussian_kde
from matplotlib.colors import LinearSegmentedColormap


# ============================================================
# 1. Input and output
# ============================================================

parser = argparse.ArgumentParser(description="Generate a DLC cis-eQTL figure panel.")
parser.add_argument("--eqtl", required=True, help="DLC cis-eQTL table")
parser.add_argument("--outdir", required=True)
args = parser.parse_args()
input_file = args.eqtl

out_dir = args.outdir
os.makedirs(out_dir, exist_ok=True)

output_pdf = os.path.join(
    out_dir,
    "eQTL_MAF005_maxMAF_ge005_vs_absBeta_density_red_xmax1p5_triangle.pdf"
)

output_png = os.path.join(
    out_dir,
    "eQTL_MAF005_maxMAF_ge005_vs_absBeta_density_red_xmax1p5_triangle.png"
)

source_excel = os.path.join(
    out_dir,
    "eQTL_MAF005_maxMAF_ge005_vs_absBeta_density_red_xmax1p5_triangle_source_data.xlsx"
)

csv_dir = os.path.join(
    out_dir,
    "eQTL_MAF005_maxMAF_ge005_vs_absBeta_density_red_xmax1p5_triangle_source_csv"
)

os.makedirs(csv_dir, exist_ok=True)


# ============================================================
# 2. Analysis and plotting parameters
# ============================================================

# Use the analysis range for correlation, fitting and KDE.
x_min, x_max = 0, 2.0

# The main plot axis extends to 1.5.
# Values above 1.5 are shown at the right edge.
x_display_max = 1.5
x_outlier_display_position = 1.53

# Apply a lower threshold of 0.05 to max_MAF.
maf_min = 0.05

kde_bw_method = 0.25
grid_n = 300

# Set the cutoff for the high-density region.
# A cutoff of 80 shows the top 20% density region.
# Lower cutoffs such as 70 or 75 display a larger region.
density_percentile_cutoff = 80

n_contour_levels = 12

# Cap the number of background points to limit PDF size.
max_scatter_points = 15000
random_seed = 1

fig_width = 4.2
fig_height = 4.2

panel_label = ""


# ============================================================
# 3. Read the data
# ============================================================

df = pd.read_csv(
    input_file,
    sep="\t",
    low_memory=False
)

required_cols = ["max_MAF", "beta"]

missing_cols = [
    c for c in required_cols
    if c not in df.columns
]

if missing_cols:
    raise ValueError(
        f"Missing required columns: {missing_cols}\n"
        f"Available columns: {list(df.columns)}"
    )

df["max_MAF"] = pd.to_numeric(
    df["max_MAF"],
    errors="coerce"
)

df["beta"] = pd.to_numeric(
    df["beta"],
    errors="coerce"
)

df = df.dropna(
    subset=[
        "max_MAF",
        "beta"
    ]
).copy()

df = df[
    np.isfinite(df["max_MAF"]) &
    np.isfinite(df["beta"])
].copy()

df["abs_beta"] = df["beta"].abs()

# Use records within the analysis limits for correlation and KDE.
plot_df_all = df[
    (df["abs_beta"] >= x_min) &
    (df["abs_beta"] <= x_max) &
    (df["max_MAF"] >= maf_min)
].copy()

x = plot_df_all["abs_beta"].values
y = plot_df_all["max_MAF"].values

# Set the upper y-axis limit from the observed values.
y_min = maf_min
y_max_data = float(np.nanmax(y))
y_max = y_max_data + max(0.01, y_max_data * 0.02)

if len(plot_df_all) < 3:
    raise ValueError("At least three valid data points are required to calculate correlation.")

print(f"Total clean points: {len(df):,}")
print(f"Points within plot range: {len(plot_df_all):,}")
print("X variable: |beta|")
print("Y variable: max_MAF")
print(f"MAF filter: max_MAF >= {maf_min}")


# ============================================================
# 4. Calculate Spearman correlation
# ============================================================

rho, pval = spearmanr(x, y)

if pval == 0:
    p_display = r"$P < 1.0 \times 10^{-308}$"
    p_display_plain = "P < 1.0e-308"
elif pval < 0.001:
    exponent = int(np.floor(np.log10(pval)))
    mantissa = pval / (10 ** exponent)
    p_display = rf"$P = {mantissa:.1f} \times 10^{{{exponent}}}$"
    p_display_plain = f"P = {mantissa:.1f}e{exponent}"
else:
    p_display = rf"$P = {pval:.3g}$"
    p_display_plain = f"P = {pval:.3g}"

print(f"Spearman rho = {rho:.4f}")
print(f"P-value = {p_display_plain}")


# ============================================================
# 5. Estimate kernel density
# ============================================================

xy = np.vstack(
    [
        x,
        y
    ]
)

kde = gaussian_kde(
    xy,
    bw_method=kde_bw_method
)

x_grid = np.linspace(
    x_min,
    x_max,
    grid_n
)

y_grid = np.linspace(
    y_min,
    y_max,
    grid_n
)

xx, yy = np.meshgrid(
    x_grid,
    y_grid
)

grid_coords = np.vstack(
    [
        xx.ravel(),
        yy.ravel()
    ]
)

zz = kde(
    grid_coords
).reshape(
    xx.shape
)

z_cut = np.percentile(
    zz,
    density_percentile_cutoff
)

zz_masked = np.ma.masked_less(
    zz,
    z_cut
)

levels = np.linspace(
    z_cut,
    zz.max(),
    n_contour_levels
)


# ============================================================
# 6. Sample points for plotting and Source Data
# ============================================================

if len(plot_df_all) > max_scatter_points:
    scatter_df = plot_df_all.sample(
        max_scatter_points,
        random_state=random_seed
    ).copy()
else:
    scatter_df = plot_df_all.copy()

scatter_df = scatter_df.reset_index(drop=True)


# Flag points beyond the main display range.
scatter_df["display_category"] = np.where(
    scatter_df["abs_beta"] > x_display_max,
    "abs_beta_gt_1.5_display_as_triangle",
    "abs_beta_le_1.5"
)

# ============================================================
# 6.1 Position points for display
# ------------------------------------------------------------
# For points with abs_beta > 1.5:
# Retain the original values in Source Data.
# Retain these points in the analysis.
# Move their display positions to the right of 1.5.
# ============================================================

scatter_df["abs_beta_plot"] = scatter_df["abs_beta"]

scatter_df.loc[
    scatter_df["abs_beta_plot"] > x_display_max,
    "abs_beta_plot"
] = x_outlier_display_position


# ============================================================
# 7. Fit the trend line
# ============================================================

coef = np.polyfit(
    x,
    y,
    1
)

poly = np.poly1d(
    coef
)

x_line = np.linspace(
    x_min,
    x_max,
    300
)

y_line = poly(
    x_line
)

trendline_df = pd.DataFrame({
    "abs_beta": x_line,
    "predicted_max_MAF": y_line
})


# ============================================================
# 8. Export Source Data
# ============================================================

# Compress the full point data as gzip CSV.
full_point_csv_gz = os.path.join(
    csv_dir,
    "Full_Point_SourceData.csv.gz"
)

plot_df_all.to_csv(
    full_point_csv_gz,
    index=False,
    compression="gzip"
)

scatter_csv = os.path.join(
    csv_dir,
    "Scatter_SourceData_used_for_plot.csv"
)

scatter_df.to_csv(
    scatter_csv,
    index=False
)

density_grid_df = pd.DataFrame({
    "abs_beta_grid": xx.ravel(),
    "max_MAF_grid": yy.ravel(),
    "KDE_density": zz.ravel(),
    "KDE_density_masked_top_region": np.where(
        zz.ravel() >= z_cut,
        zz.ravel(),
        np.nan
    )
})

density_grid_csv_gz = os.path.join(
    csv_dir,
    "KDE_Density_Grid_SourceData.csv.gz"
)

density_grid_df.to_csv(
    density_grid_csv_gz,
    index=False,
    compression="gzip"
)

contour_levels_df = pd.DataFrame({
    "Level_index": np.arange(1, len(levels) + 1),
    "KDE_density_level": levels
})

trendline_csv = os.path.join(
    csv_dir,
    "Trendline_SourceData.csv"
)

trendline_df.to_csv(
    trendline_csv,
    index=False
)

summary_df = pd.DataFrame({
    "Item": [
        "input_file",
        "n_total_clean_points",
        "n_points_within_plot_range",
        "n_scatter_points_used_for_plot",
        "x_variable",
        "y_variable",
        "Spearman_rho",
        "Spearman_pvalue_scipy",
        "Spearman_pvalue_display",
        "KDE_bw_method",
        "KDE_grid_n",
        "density_percentile_cutoff",
        "z_cut",
        "z_max",
        "n_contour_levels",
        "linear_fit_slope",
        "linear_fit_intercept",
        "x_min",
        "x_max",
        "maf_filter_rule",
        "y_min",
        "y_max_display",
        "output_pdf",
        "output_png",
        "source_excel",
        "full_point_csv_gz",
        "density_grid_csv_gz"
    ],
    "Value": [
        input_file,
        len(df),
        len(plot_df_all),
        len(scatter_df),
        "abs_beta = |beta|",
        "max_MAF",
        rho,
        pval,
        p_display_plain,
        kde_bw_method,
        grid_n,
        density_percentile_cutoff,
        z_cut,
        zz.max(),
        n_contour_levels,
        coef[0],
        coef[1],
        x_min,
        x_max,
        f"max_MAF >= {maf_min}",
        y_min,
        y_max,
        output_pdf,
        output_png,
        source_excel,
        full_point_csv_gz,
        density_grid_csv_gz
    ]
})

plot_config_df = pd.DataFrame({
    "Item": [
        "figure_width",
        "figure_height",
        "scatter_color",
        "scatter_alpha",
        "density_colormap",
        "density_color_low",
        "density_color_mid",
        "density_color_high",
        "contour_line_color",
        "trendline_color",
        "trendline_style",
        "x_label",
        "y_label",
        "y_source_column"
    ],
    "Value": [
        fig_width,
        fig_height,
        "#fdd0c2",
        0.22,
        "custom_light_to_red",
        "#fff5f0",
        "#fb6a4a",
        "#a50f15",
        "white",
        "#444444",
        "dashed",
        "|Effect Size|",
        "MAF",
        "max_MAF"
    ]
})

readme_df = pd.DataFrame({
    "Sheet_or_file": [
        "README",
        "Summary",
        "Scatter_SourceData",
        "KDE_Density_Grid_SourceData",
        "Contour_Levels",
        "Trendline_SourceData",
        "Plot_Config",
        "Full_Point_SourceData.csv.gz"
    ],
    "Description": [
        "Description of each Source Data file.",
        "Input and output paths, point counts, correlation, KDE settings and trend-line parameters.",
        "Sampled points used for the background scatterplot.",
        "Two-dimensional KDE grid (|beta| × max_MAF) for recreating the filled contours.",
        "Density levels used for filled contours and contour lines.",
        "Linear trend-line coordinates.",
        "Plot settings and colors.",
        "All points within the analysis limits, saved as gzip CSV."
    ]
})

# Keep the Excel file compact by omitting the full point and density tables.
# The full tables are available as compressed CSV files.
with pd.ExcelWriter(
    source_excel,
    engine="openpyxl"
) as writer:
    readme_df.to_excel(
        writer,
        sheet_name="README",
        index=False
    )

    summary_df.to_excel(
        writer,
        sheet_name="Summary",
        index=False
    )

    scatter_df.to_excel(
        writer,
        sheet_name="Scatter_SourceData",
        index=False
    )

    contour_levels_df.to_excel(
        writer,
        sheet_name="Contour_Levels",
        index=False
    )

    trendline_df.to_excel(
        writer,
        sheet_name="Trendline_SourceData",
        index=False
    )

    plot_config_df.to_excel(
        writer,
        sheet_name="Plot_Config",
        index=False
    )

print("Source data saved:")
print("Excel:", source_excel)
print("CSV dir:", csv_dir)
print("Full point data:", full_point_csv_gz)
print("KDE density grid:", density_grid_csv_gz)


# ============================================================
# 9. Draw the figure
# ============================================================

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.unicode_minus"] = False

# Red gradient from low to high density.
# Low density: near-white pink.
# Medium density: pink-orange to bright red.
# High density: dark red.
red_density_cmap = LinearSegmentedColormap.from_list(
    "soft_to_deep_red",
    [
        "#fff5f0",
        "#fee0d2",
        "#fcbba1",
        "#fc9272",
        "#fb6a4a",
        "#de2d26",
        "#a50f15"
    ]
)

fig, ax = plt.subplots(
    figsize=(fig_width, fig_height)
)

# 1. Draw translucent pale-pink background points.
# Draw the points within the main display range.
normal_scatter = scatter_df[
    scatter_df["abs_beta"] <= x_display_max
]

ax.scatter(
    normal_scatter["abs_beta_plot"],
    normal_scatter["max_MAF"],
    s=10,
    alpha=0.22,
    color="#fdd0c2",
    edgecolors="none",
    rasterized=True,
    zorder=1
)

# Handle points beyond the main display range.
# Retain their data and mark them with triangles at the right edge.
outlier_scatter = scatter_df[
    scatter_df["abs_beta"] > x_display_max
]

if len(outlier_scatter) > 0:
    ax.scatter(
        np.repeat(
            x_outlier_display_position,
            len(outlier_scatter)
        ),
        outlier_scatter["max_MAF"],
        s=28,
        marker=">",
        alpha=0.55,
        color="#fdd0c2",
        edgecolors="none",
        rasterized=True,
        zorder=2
    )

    ax.text(
        x_outlier_display_position,
        y_max * 0.96,
        ">|1.5",
        ha="center",
        va="top",
        fontsize=8,
        rotation=90
    )

# 2. Draw red filled contours.
cf = ax.contourf(
    xx,
    yy,
    zz_masked,
    levels=levels,
    cmap=red_density_cmap,
    antialiased=True,
    alpha=0.90,
    zorder=2
)

# 3. Outline density contours in white.
ax.contour(
    xx,
    yy,
    zz_masked,
    levels=levels,
    colors="white",
    linewidths=0.65,
    zorder=3
)

# 4. Draw the trend line.
ax.plot(
    x_line,
    y_line,
    color="#444444",
    linewidth=1.5,
    linestyle="--",
    zorder=4
)

# 5. Annotate the correlation.
ax.text(
    0.50,
    1.03,
    f"Spearman rho = {rho:.3f}, {p_display}\n$n$ = {len(plot_df_all):,}",
    transform=ax.transAxes,
    ha="center",
    va="bottom",
    fontsize=10,
    linespacing=1.5
)

# Optional panel label.
if panel_label:
    ax.text(
        -0.18,
        1.05,
        panel_label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        ha="left",
        va="bottom"
    )

# Configure the axes.
ax.set_xlabel(
    "|Effect Size|",
    fontsize=12
)

ax.set_ylabel(
    "MAF",
    fontsize=12
)

ax.set_xlim(
    x_min,
    x_outlier_display_position + 0.02
)

# Show major ticks through 1.5.
ax.set_xticks(
    [0, 0.5, 1.0, 1.5]
)

ax.set_ylim(
    y_min,
    y_max
)

# Configure the axis borders.
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

for spine in ax.spines.values():
    spine.set_linewidth(1)

ax.tick_params(
    axis="both",
    labelsize=10,
    width=1,
    length=4
)

ax.grid(False)

plt.tight_layout()

plt.savefig(
    output_pdf,
    format="pdf",
    bbox_inches="tight"
)

plt.savefig(
    output_png,
    dpi=600,
    bbox_inches="tight"
)

plt.close()

print("The x-axis shows 0 to 1.5. Larger values are displayed at the right edge.")
print("Completed.")
print(f"Number of points: {len(plot_df_all):,}")
print(f"Spearman rho: {rho:.4f}")
print(f"P-value: {p_display_plain}")
print("PDF:", output_pdf)
print("PNG:", output_png)
print("Source data Excel:", source_excel)
print("Source data CSV dir:", csv_dir)
