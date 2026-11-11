#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================
# 1. Input file settings
# ============================================================

# Use an eQTL table with one SNP–gene pair per row.
parser = argparse.ArgumentParser(description="Generate a DLC cis-eQTL figure panel.")
parser.add_argument("--eqtl", required=True, help="DLC cis-eQTL table")
parser.add_argument("--outdir", required=True)
args = parser.parse_args()
input_file = args.eqtl

# The input contains SNP–gene pairs.
INPUT_MODE = "pairs"

# Set these names to match the input columns.
gene_col = "Gene.ID"
snp_col = "SNP_id"

# Used only with INPUT_MODE = "gene_count".
count_col = "n_unique_SNPs"

# Use False to analyze observed data.
USE_SIMULATED_DATA = False

# ============================================================
# 2. Output paths
# ============================================================

out_dir = args.outdir
os.makedirs(out_dir, exist_ok=True)

source_excel = os.path.join(
    out_dir,
    "eQTL_MAF005_unique_SNPs_per_gene_source_data.xlsx"
)

csv_dir = os.path.join(
    out_dir,
    "eQTL_MAF005_unique_SNPs_per_gene_source_csv"
)
os.makedirs(csv_dir, exist_ok=True)

out_pdf = os.path.join(
    out_dir,
    "eQTL_MAF005_unique_SNPs_per_gene_histogram.pdf"
)

out_png = os.path.join(
    out_dir,
    "eQTL_MAF005_unique_SNPs_per_gene_histogram.png"
)

# ============================================================
# 3. Plot settings
# ============================================================

# Log-spaced bins from 1 to 100,000.
bins = np.logspace(0, 5, 60)

bar_color = "#318195"

x_label = "Number of unique SNPs per gene"
y_label = "Number of genes"

x_ticks = [1, 10, 100, 1000, 10000, 100000]
x_tick_labels = ["1", "10", "100", "1K", "10K", "100K"]

# Figure dimensions.
fig_width = 3.2
fig_height = 3.0

# ============================================================
# 4. Utility functions
# ============================================================

def read_table_auto(path):
    """
    Read XLSX, CSV, TSV, TXT or an extensionless tab-separated file.
    Use tab separation for files without a recognized extension.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    lower = path.lower()

    if lower.endswith(".xlsx") or lower.endswith(".xls"):
        return pd.read_excel(path)

    elif lower.endswith(".csv"):
        return pd.read_csv(path, low_memory=False)

    else:
        return pd.read_csv(path, sep="\t", low_memory=False)


def write_df_to_excel_safely(writer, df, sheet_name, index=False):
    """
    An Excel worksheet supports up to 1,048,576 rows.
    Split larger tables across worksheets.
    """
    max_excel_rows = 1048576
    max_data_rows = max_excel_rows - 1

    sheet_name = sheet_name[:31]

    if len(df) <= max_data_rows:
        df.to_excel(writer, sheet_name=sheet_name, index=index)
    else:
        n_parts = int(np.ceil(len(df) / max_data_rows))

        for i in range(n_parts):
            start = i * max_data_rows
            end = min((i + 1) * max_data_rows, len(df))

            part_sheet = f"{sheet_name[:27]}_{i + 1}"

            df.iloc[start:end].to_excel(
                writer,
                sheet_name=part_sheet,
                index=index
            )


# ============================================================
# 5. Read eQTL data and count unique SNPs per gene
# ============================================================

if USE_SIMULATED_DATA:

    print(">>> Using simulated data for testing. Set USE_SIMULATED_DATA = False for study Source Data.")

    np.random.seed(42)
    data_eqtl = np.random.lognormal(5, 1.2, 2000)
    data_eqtl = data_eqtl[data_eqtl >= 1]

    gene_source_data = pd.DataFrame({
        "gene_id": [f"gene_{i+1}" for i in range(len(data_eqtl))],
        "n_unique_SNPs": data_eqtl
    })

    raw_pair_data = pd.DataFrame()

else:

    print(">>> Reading eQTL input file:")
    print(input_file)

    raw_df = read_table_auto(input_file)

    print(">>> Input columns:")
    print(list(raw_df.columns))

    if INPUT_MODE == "pairs":

        required_cols = [gene_col, snp_col]

        missing_cols = [
            c for c in required_cols
            if c not in raw_df.columns
        ]

        if missing_cols:
            raise ValueError(
                f"INPUT_MODE='pairs' requires these columns: {required_cols}\n"
                f"Missing columns: {missing_cols}\n"
                f"Available columns: {list(raw_df.columns)}"
            )

        raw_pair_data = raw_df[[gene_col, snp_col]].copy()

        raw_pair_data[gene_col] = raw_pair_data[gene_col].astype(str).str.strip()
        raw_pair_data[snp_col] = raw_pair_data[snp_col].astype(str).str.strip()

        raw_pair_data = raw_pair_data[
            raw_pair_data[gene_col].notna() &
            raw_pair_data[snp_col].notna() &
            (raw_pair_data[gene_col] != "") &
            (raw_pair_data[snp_col] != "") &
            (raw_pair_data[gene_col].str.lower() != "nan") &
            (raw_pair_data[snp_col].str.lower() != "nan")
        ].copy()

        raw_pair_data = raw_pair_data.rename(
            columns={
                gene_col: "gene_id",
                snp_col: "SNP_id"
            }
        )

        # Remove duplicate SNP–gene pairs.
        unique_pair_data = raw_pair_data.drop_duplicates(
            subset=["gene_id", "SNP_id"]
        ).copy()

        # Count unique SNPs per gene.
        gene_source_data = (
            unique_pair_data
            .groupby("gene_id", as_index=False)["SNP_id"]
            .nunique()
            .rename(columns={
                "SNP_id": "n_unique_SNPs"
            })
        )

    elif INPUT_MODE == "gene_count":

        required_cols = [gene_col, count_col]

        missing_cols = [
            c for c in required_cols
            if c not in raw_df.columns
        ]

        if missing_cols:
            raise ValueError(
                f"INPUT_MODE='gene_count' requires these columns: {required_cols}\n"
                f"Missing columns: {missing_cols}\n"
                f"Available columns: {list(raw_df.columns)}"
            )

        raw_pair_data = pd.DataFrame()

        gene_source_data = raw_df[[gene_col, count_col]].copy()

        gene_source_data = gene_source_data.rename(
            columns={
                gene_col: "gene_id",
                count_col: "n_unique_SNPs"
            }
        )

    else:
        raise ValueError("INPUT_MODE must be 'pairs' or 'gene_count'.")

# Clean numeric values.
gene_source_data["n_unique_SNPs"] = pd.to_numeric(
    gene_source_data["n_unique_SNPs"],
    errors="coerce"
)

gene_source_data = gene_source_data.dropna(
    subset=["gene_id", "n_unique_SNPs"]
).copy()

gene_source_data = gene_source_data[
    np.isfinite(gene_source_data["n_unique_SNPs"]) &
    (gene_source_data["n_unique_SNPs"] >= 1)
].copy()

gene_source_data = gene_source_data.sort_values(
    "n_unique_SNPs",
    ascending=False
).reset_index(drop=True)

gene_source_data.insert(
    0,
    "Gene_order",
    np.arange(1, len(gene_source_data) + 1)
)

data_eqtl = gene_source_data["n_unique_SNPs"].values

print(f">>> Number of genes used for plotting: {len(data_eqtl):,}")
print(f">>> Min unique SNPs per gene: {np.min(data_eqtl):.0f}")
print(f">>> Median unique SNPs per gene: {np.median(data_eqtl):.0f}")
print(f">>> Max unique SNPs per gene: {np.max(data_eqtl):.0f}")

# ============================================================
# 6. Calculate histogram Source Data
# ============================================================

counts_eqtl, bin_edges = np.histogram(
    data_eqtl,
    bins=bins
)

bin_centers = (bin_edges[:-1] * bin_edges[1:]) ** 0.5

hist_source_data = pd.DataFrame({
    "Bin_ID": np.arange(1, len(counts_eqtl) + 1),
    "Bin_left": bin_edges[:-1],
    "Bin_right": bin_edges[1:],
    "Bin_center_geometric": bin_centers,
    "Bin_width": np.diff(bin_edges),
    "Number_of_genes": counts_eqtl,
    "Bar_height": counts_eqtl,
    "Bar_color": bar_color
})

hist_nonzero_source_data = hist_source_data[
    hist_source_data["Number_of_genes"] > 0
].copy()

# Assign each gene to its histogram bin.
gene_source_data["Histogram_bin_ID"] = pd.cut(
    gene_source_data["n_unique_SNPs"],
    bins=bin_edges,
    labels=np.arange(1, len(bin_edges)),
    include_lowest=True
)

# ============================================================
# 7. Summary statistics
# ============================================================

summary_data = pd.DataFrame({
    "Item": [
        "input_file",
        "input_mode",
        "use_simulated_data",
        "number_of_genes",
        "number_of_unique_gene_snp_pairs",
        "min_unique_SNPs_per_gene",
        "q25_unique_SNPs_per_gene",
        "median_unique_SNPs_per_gene",
        "mean_unique_SNPs_per_gene",
        "q75_unique_SNPs_per_gene",
        "max_unique_SNPs_per_gene",
        "number_of_bins",
        "x_axis_min",
        "x_axis_max"
    ],
    "Value": [
        input_file,
        INPUT_MODE,
        USE_SIMULATED_DATA,
        len(gene_source_data),
        len(raw_pair_data.drop_duplicates()) if not raw_pair_data.empty else "NA",
        np.min(data_eqtl),
        np.percentile(data_eqtl, 25),
        np.median(data_eqtl),
        np.mean(data_eqtl),
        np.percentile(data_eqtl, 75),
        np.max(data_eqtl),
        len(counts_eqtl),
        bins[0],
        bins[-1]
    ]
})

plot_config = pd.DataFrame({
    "Item": [
        "source_excel",
        "pdf_output",
        "png_output",
        "x_axis",
        "y_axis",
        "x_scale",
        "bins",
        "bar_height",
        "bar_color",
        "x_ticks",
        "x_tick_labels",
        "removed_finemapping_result",
        "figure_width",
        "figure_height"
    ],
    "Value": [
        source_excel,
        out_pdf,
        out_png,
        "Number of unique SNPs per gene",
        "Number of genes",
        "log10",
        "np.logspace(0, 5, 60)",
        "Number_of_genes",
        bar_color,
        ",".join(map(str, x_ticks)),
        ",".join(x_tick_labels),
        "True",
        fig_width,
        fig_height
    ]
})

readme = pd.DataFrame({
    "Sheet": [
        "README",
        "Gene_SourceData",
        "Histogram_SourceData",
        "Histogram_Nonzero_Bins",
        "Raw_Unique_Pair_Data",
        "Summary",
        "Plot_Config"
    ],
    "Description": [
        "Description of each worksheet.",
        "Number of unique eQTL SNPs per gene.",
        "Full histogram table. Bin_center_geometric gives the x-coordinate and Number_of_genes/Bar_height gives the bar height.",
        "Histogram bins with Number_of_genes > 0.",
        "Deduplicated gene–SNP pairs for INPUT_MODE='pairs'.",
        "Summary of the input file, gene count and unique SNP counts.",
        "Plot parameters, output paths, axes and colors."
    ]
})

# ============================================================
# 8. Export CSV copies
# ============================================================

csv_outputs = {
    "README": os.path.join(csv_dir, "README.csv"),
    "Gene_SourceData": os.path.join(csv_dir, "Gene_SourceData.csv"),
    "Histogram_SourceData": os.path.join(csv_dir, "Histogram_SourceData.csv"),
    "Histogram_Nonzero_Bins": os.path.join(csv_dir, "Histogram_Nonzero_Bins.csv"),
    "Raw_Unique_Pair_Data": os.path.join(csv_dir, "Raw_Unique_Pair_Data.csv"),
    "Summary": os.path.join(csv_dir, "Summary.csv"),
    "Plot_Config": os.path.join(csv_dir, "Plot_Config.csv")
}

readme.to_csv(csv_outputs["README"], index=False)
gene_source_data.to_csv(csv_outputs["Gene_SourceData"], index=False)
hist_source_data.to_csv(csv_outputs["Histogram_SourceData"], index=False)
hist_nonzero_source_data.to_csv(csv_outputs["Histogram_Nonzero_Bins"], index=False)

if not USE_SIMULATED_DATA and INPUT_MODE == "pairs":
    unique_pair_data.to_csv(csv_outputs["Raw_Unique_Pair_Data"], index=False)
else:
    pd.DataFrame().to_csv(csv_outputs["Raw_Unique_Pair_Data"], index=False)

summary_data.to_csv(csv_outputs["Summary"], index=False)
plot_config.to_csv(csv_outputs["Plot_Config"], index=False)

print(f">>> CSV Source Data exported to: {csv_dir}")

# ============================================================
# 9. Export Excel Source Data
# ============================================================

try:
    with pd.ExcelWriter(source_excel, engine="openpyxl") as writer:
        write_df_to_excel_safely(writer, readme, "README", index=False)
        write_df_to_excel_safely(writer, gene_source_data, "Gene_SourceData", index=False)
        write_df_to_excel_safely(writer, hist_source_data, "Histogram_SourceData", index=False)
        write_df_to_excel_safely(writer, hist_nonzero_source_data, "Histogram_Nonzero_Bins", index=False)

        if not USE_SIMULATED_DATA and INPUT_MODE == "pairs":
            write_df_to_excel_safely(writer, unique_pair_data, "Raw_Unique_Pair_Data", index=False)

        write_df_to_excel_safely(writer, summary_data, "Summary", index=False)
        write_df_to_excel_safely(writer, plot_config, "Plot_Config", index=False)

    print(f">>> Excel Source Data exported to: {source_excel}")

except Exception as e:
    print(">>> Excel export failed:")
    print(e)
    print(f">>> CSV Source Data remain available at: {csv_dir}")
    print(">>> Install openpyxl with: pip install openpyxl")

# ============================================================
# 10. Plot the eQTL histogram
# ============================================================

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.unicode_minus"] = False

fig, ax = plt.subplots(
    figsize=(fig_width, fig_height),
    dpi=300
)

ax.bar(
    hist_source_data["Bin_center_geometric"],
    hist_source_data["Bar_height"],
    width=hist_source_data["Bin_width"],
    color=bar_color,
    alpha=0.85,
    edgecolor=bar_color,
    linewidth=0.2,
    align="center",
    label="eQTL"
)

ax.set_xscale("log")
ax.set_xlim(1, 100000)

ax.set_xticks(x_ticks)
ax.set_xticklabels(
    x_tick_labels,
    fontsize=9
)

ax.set_xlabel(
    x_label,
    fontsize=10
)

ax.set_ylabel(
    y_label,
    fontsize=10
)

ax.tick_params(
    axis="y",
    labelsize=9
)

ax.legend(
    frameon=False,
    fontsize=9,
    loc="upper right"
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.grid(False)

plt.tight_layout()

fig.savefig(
    out_pdf,
    format="pdf",
    bbox_inches="tight"
)

fig.savefig(
    out_png,
    dpi=600,
    bbox_inches="tight"
)

plt.close(fig)

print(">>> Completed.")
print(f">>> PDF: {out_pdf}")
print(f">>> PNG: {out_png}")
print(f">>> Source data Excel: {source_excel}")
print(f">>> CSV source data dir: {csv_dir}")

print("\n>>> Worksheets for recreating the figure:")
print("1. Gene_SourceData")
print("2. Histogram_SourceData")
print("3. Plot_Config")
