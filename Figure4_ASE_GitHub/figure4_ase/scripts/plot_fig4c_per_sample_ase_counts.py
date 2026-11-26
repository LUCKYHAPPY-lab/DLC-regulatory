#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.stats import pearsonr, t


# ============================================================
# 1. Input and output
# ============================================================

parser = argparse.ArgumentParser(description="Generate Figure 4c: per sample ase counts.")
parser.add_argument("--sample-counts", required=True, type=Path)
parser.add_argument("--outdir", required=True, type=Path)
args = parser.parse_args()
INPUT = args.sample_counts
OUTDIR = args.outdir

OUT_PDF = OUTDIR / "ASE_scatter_Pearson_95PI.pdf"
OUT_PNG = OUTDIR / "ASE_scatter_Pearson_95PI.png"
OUT_WHITE = OUTDIR / "ASE_scatter_Pearson_95PI_white.png"
OUT_XLSX = OUTDIR / "ASE_scatter_95PI_SourceData_single_sheet.xlsx"
OUT_POINTS = OUTDIR / "Plot_points.tsv"
OUT_FIT = OUTDIR / "Regression_95PI.tsv"
OUT_STATS = OUTDIR / "Statistics.tsv"

X_COL = "biallelic_count_MAF005"
Y_COL = "ASE_list_observed_count_MAF005"


# ============================================================
# 2. Read and validate per-sample counts
# ============================================================

if not INPUT.is_file():
    raise FileNotFoundError(f"Count table not found: {INPUT}")

OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(
    INPUT,
    sep="\t",
    dtype={"sample_id": "string"}
)

df.columns = df.columns.astype(str).str.strip()

required = {"sample_id", X_COL, Y_COL}
missing = required - set(df.columns)

if missing:
    raise ValueError(
        f"Missing columns: {sorted(missing)}\n"
        f"Available columns: {list(df.columns)}"
    )

df["sample_id"] = df["sample_id"].str.strip()

values = df[[X_COL, Y_COL]].apply(
    pd.to_numeric,
    errors="coerce"
)

valid_sample_id = (
    df["sample_id"].notna()
    & ~df["sample_id"].isin(["", "NA", "nan"])
)

finite_values = np.isfinite(
    values.to_numpy(dtype=float)
).all(axis=1)

valid = pd.Series(
    finite_values,
    index=df.index
) & valid_sample_id

excluded = df.loc[~valid].copy()

excluded.to_csv(
    OUTDIR / "Excluded_samples.tsv",
    sep="\t",
    index=False,
    na_rep="NA"
)

points = df.loc[valid, ["sample_id"]].copy()

points["Both_alleles_observed_SNPs"] = values.loc[valid, X_COL]
points["Observed_ASE_SNPs"] = values.loc[valid, Y_COL]

if points["sample_id"].duplicated().any():
    raise ValueError("Duplicate sample IDs found in the count table.")

x = points["Both_alleles_observed_SNPs"].to_numpy(dtype=float)
y = points["Observed_ASE_SNPs"].to_numpy(dtype=float)

n = len(points)

if n < 3:
    raise ValueError("At least three valid samples are required.")

if np.any(x < 0) or np.any(y < 0):
    raise ValueError("SNP counts must be nonnegative.")

if np.any(x != np.floor(x)) or np.any(y != np.floor(y)):
    raise ValueError("SNP counts must be integers.")

if np.any(y > x):
    raise ValueError(
        "Some samples have more ASE SNPs than RNA-covered biallelic SNPs. "
        "Check the upstream count table."
    )

if np.ptp(x) == 0 or np.ptp(y) == 0:
    raise ValueError("Pearson correlation requires variation in both variables.")


# ============================================================
# 3. Pearson correlation and ordinary least-squares regression
# ============================================================

r, p = pearsonr(x, y)

x_mean = x.mean()
y_mean = y.mean()

sxx = np.sum((x - x_mean) ** 2)

slope = (
    np.sum((x - x_mean) * (y - y_mean))
    / sxx
)

intercept = y_mean - slope * x_mean

fitted_y = intercept + slope * x
residuals = y - fitted_y

degrees_of_freedom = n - 2
sse = np.sum(residuals ** 2)

residual_se = np.sqrt(
    sse / degrees_of_freedom
)


# ============================================================
# 4. Pointwise 95% prediction interval
# ============================================================

# Classical linear regression prediction interval for a new observation at a given x:
#
# y_hat ± t(0.975, n−2) × residual_se
#         × sqrt[1 + 1/n + (x−x_mean)^2 / Sxx]
#
# The additional 1 under the square root accounts for the variability of a new observation.
# The interval assumes independent, homoscedastic and approximately normal errors.

x_line = np.linspace(
    x.min(),
    x.max(),
    300
)

y_line = intercept + slope * x_line

prediction_se = residual_se * np.sqrt(
    1.0
    + 1.0 / n
    + (x_line - x_mean) ** 2 / sxx
)

t_critical = t.ppf(
    0.975,
    df=degrees_of_freedom
)

pi_lower = y_line - t_critical * prediction_se
pi_upper = y_line + t_critical * prediction_se


# ============================================================
# 5. Plotting
# ============================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.unicode_minus": False,
})

fig, ax = plt.subplots(figsize=(4.4, 3.8))

# 95% prediction interval
ax.fill_between(
    x_line,
    pi_lower,
    pi_upper,
    color="#00AAE3",
    alpha=0.18,
    linewidth=0,
    zorder=1
)

# Include all valid samples, including outliers.
ax.scatter(
    x,
    y,
    s=15,
    color="#00AAE3",
    alpha=0.75,
    edgecolors="none",
    zorder=2
)

# Fitted regression line
ax.plot(
    x_line,
    y_line,
    color="#1F4E9F",
    linestyle="--",
    linewidth=1.4,
    zorder=3
)

# Display the two-sided P value in scientific notation.
if p > 0:
    mantissa, exponent = f"{p:.2e}".split("e")
    p_label = (
        rf"$P = {float(mantissa):.2f}"
        rf"\times 10^{{{int(exponent)}}}$"
    )
else:
    # Handle numerical underflow when formatting the P value.
    p_label = r"$P$: numerical underflow"

ax.set_title(
    rf"$r = {r:.2f}$    {p_label}    $n = {n}$",
    fontsize=10,
    pad=10
)

ax.set_xlabel("SNPs with both alleles observed")
ax.set_ylabel("Observed ASE SNPs")

formatter = FuncFormatter(
    lambda value, position: f"{value:,.0f}"
)

ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_formatter(formatter)

ax.tick_params(
    direction="out",
    length=3,
    width=0.8
)

for spine in ax.spines.values():
    spine.set_linewidth(0.8)

ax.grid(False)
ax.margins(x=0.05, y=0.08)

fig.tight_layout()

# Export a transparent PDF with TrueType text.
fig.savefig(
    OUT_PDF,
    bbox_inches="tight",
    transparent=True
)

# Export a transparent PNG.
fig.savefig(
    OUT_PNG,
    dpi=600,
    bbox_inches="tight",
    transparent=True
)

# Export a white-background PNG for previewing.
fig.savefig(
    OUT_WHITE,
    dpi=600,
    bbox_inches="tight",
    transparent=False,
    facecolor="white"
)

plt.close(fig)


# ============================================================
# 6. Assemble Source Data
# ============================================================

points["Fitted_y"] = fitted_y
points["Residual"] = residuals

fit_data = pd.DataFrame({
    "X": x_line,
    "Fitted_y": y_line,
    "Prediction_standard_error": prediction_se,
    "PI95_lower": pi_lower,
    "PI95_upper": pi_upper,
})

summary = pd.DataFrame([
    ("Input_file", str(INPUT)),
    ("X_input_column", X_COL),
    ("Y_input_column", Y_COL),
    ("Input_sample_rows", len(df)),
    ("Excluded_invalid_sample_rows", len(excluded)),
    ("N_samples", n),
    ("Pearson_r", r),
    ("Pearson_two_sided_P", p),
    ("P_numeric_underflow", bool(p == 0)),
    ("Regression_slope", slope),
    ("Regression_intercept", intercept),
    ("Regression_R_squared", r ** 2),
    ("Residual_standard_error", residual_se),
    ("Degrees_of_freedom", degrees_of_freedom),
    ("t_critical_0.975", t_critical),
    (
        "Band_definition",
        "Pointwise 95% prediction interval for a new observation"
    ),
    (
        "Band_formula",
        "y_hat +/- t_critical * s * sqrt(1 + 1/n + (x-x_mean)^2/Sxx)"
    ),
    ("Outlier_exclusion", "None"),
    (
        "X_definition",
        "RNA-covered SNPs with ref_count >= 1 and alt_count >= 1 "
        "after DLC site-max MAF filtering"
    ),
    (
        "Y_definition",
        "X sites present in the final ASE list; "
        "not individual-sample ASE significance"
    ),
], columns=["Statistic", "Value"])

points.to_csv(
    OUT_POINTS,
    sep="\t",
    index=False,
    na_rep="NA"
)

fit_data.to_csv(
    OUT_FIT,
    sep="\t",
    index=False,
    na_rep="NA"
)

summary.to_csv(
    OUT_STATS,
    sep="\t",
    index=False,
    na_rep="NA"
)

# Use Record_type to distinguish sample, fit and statistics records on one worksheet.
source_data = pd.concat([
    points.assign(Record_type="Sample"),
    fit_data.assign(Record_type="Prediction_interval"),
    summary.assign(Record_type="Statistic"),
], ignore_index=True, sort=False)

source_data = source_data[
    ["Record_type"]
    + [
        col for col in source_data.columns
        if col != "Record_type"
    ]
]

with pd.ExcelWriter(
    OUT_XLSX,
    engine="openpyxl"
) as writer:
    source_data.to_excel(
        writer,
        sheet_name="Source_Data",
        index=False,
        na_rep="NA"
    )


# ============================================================
# 7. Report results
# ============================================================

print(f"Valid samples: {n}")
print(f"Invalid sample rows excluded: {len(excluded)}")
print(f"Pearson r：{r:.10f}")
print(f"Two-sided P value: {p:.10e}")
print(f"Regression slope: {slope:.10f}")
print(f"Regression intercept: {intercept:.10f}")
print(f"Residual standard error: {residual_se:.10f}")
print("Shading shows the pointwise 95% prediction interval.")
print("Correlation and regression use all valid samples, including outliers.")

for path in [
    OUT_PDF,
    OUT_PNG,
    OUT_WHITE,
    OUT_XLSX,
    OUT_POINTS,
    OUT_FIT,
    OUT_STATS,
]:
    print(f"Saved: {path}")