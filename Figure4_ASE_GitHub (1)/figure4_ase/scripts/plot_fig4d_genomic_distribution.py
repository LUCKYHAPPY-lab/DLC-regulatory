#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

parser = argparse.ArgumentParser(description="Generate Figure 4d: genomic distribution.")
parser.add_argument("--annotation-statistics", required=True, type=Path)
parser.add_argument("--outdir", required=True, type=Path)
args = parser.parse_args()
INPUT = args.annotation_statistics
OUT = args.outdir

OUT.mkdir(parents=True, exist_ok=False)

categories = ["3_prime_UTR", "CDS", "5_prime_UTR", "Noncoding_exon"]
labels = {
    "3_prime_UTR": "3′UTR",
    "CDS": "CDS",
    "5_prime_UTR": "5′UTR",
    "Noncoding_exon": "Noncoding exon",
}
colors = {
    "3_prime_UTR": "#159D86",
    "CDS": "#45AEC0",
    "5_prime_UTR": "#E96A38",
    "Noncoding_exon": "#B79AC8",
}

full = pd.read_csv(INPUT, sep="\t")
required = {"Category", "Background_SNP_count", "ASE_SNP_count"}
if not required.issubset(full.columns):
    raise ValueError(f"Missing columns: {required - set(full.columns)}")
if full["Category"].duplicated().any():
    raise ValueError("Duplicate Category values in the input.")
if not set(categories).issubset(set(full["Category"])):
    raise ValueError("The input must include all four plotted categories.")

for col in ["Background_SNP_count", "ASE_SNP_count"]:
    full[col] = pd.to_numeric(full[col], errors="raise")
    if full[col].isna().any() or (full[col] < 0).any():
        raise ValueError(f"{col} contains missing or negative values.")
    if (full[col] % 1 != 0).any():
        raise ValueError(f"{col} contains noninteger values.")
    full[col] = full[col].astype(int)

if (full["ASE_SNP_count"] > full["Background_SNP_count"]).any():
    raise ValueError("ASE counts exceed background counts in one or more categories.")

total_ase = int(full["ASE_SNP_count"].sum())
total_background = int(full["Background_SNP_count"].sum())
if total_ase == 0:
    raise ValueError("The total ASE count is zero.")

shown = full.set_index("Category").loc[categories].reset_index()
shown_background = int(shown["Background_SNP_count"].sum())
outside_background = total_background - shown_background
outside_ase = total_ase - int(shown["ASE_SNP_count"].sum())

# Calculate bar percentages using all final ASE SNPs, including sites in undisplayed categories.
shown["Percentage_of_ASE_SNPs"] = (
    shown["ASE_SNP_count"] / total_ase * 100
)
shown["ASE_over_background_percent"] = (
    shown["ASE_SNP_count"]
    / shown["Background_SNP_count"].replace(0, float("nan"))
    * 100
)

# Sort by decreasing bar height and retain the original order for ties.
shown = shown.sort_values(
    "Percentage_of_ASE_SNPs", ascending=False, kind="stable"
).reset_index(drop=True)
shown.insert(0, "Plot_order", range(1, len(shown) + 1))
shown["Display_label"] = shown["Category"].map(labels)
shown["Total_final_ASE_SNPs"] = total_ase
shown["Total_MAF_filtered_background_SNPs"] = total_background
shown["Four_category_background_SNPs"] = shown_background
shown["Background_SNPs_outside_four_categories"] = outside_background
shown["ASE_SNPs_outside_four_categories"] = outside_ase

# Export one Source Data worksheet with one row per plotted category.
shown.to_csv(
    OUT / "ASE_four_categories_SourceData.tsv",
    sep="\t", index=False, na_rep="NA",
)
shown.to_excel(
    OUT / "ASE_four_categories_SourceData.xlsx",
    sheet_name="Source_Data", index=False, na_rep="NA",
    engine="openpyxl",
)

# Retain counts for undisplayed categories for verification.
full.to_csv(
    OUT / "All_categories_statistics.tsv",
    sep="\t", index=False, na_rep="NA",
)

readme = f"""Input statistics table:
{INPUT}

Total MAF-filtered background SNPs: {total_background}
Background SNPs in the four plotted categories: {shown_background}
Background SNPs in unplotted categories: {outside_background}
Total final ASE SNPs: {total_ase}
ASE SNPs in unplotted categories: {outside_ase}

Bars show the percentage of all final ASE SNPs in each category.
The line shows the percentage of background SNPs classified as ASE within each category.
The calculation includes unplotted categories in the total ASE count.
The Source Data workbook has one worksheet with one row per plotted category.
The background uses the original MAF-only filter.
"""
(OUT / "README.txt").write_text(readme, encoding="utf-8")

print(f"Reading: {INPUT}", flush=True)
print(shown[[
    "Display_label", "Background_SNP_count", "ASE_SNP_count",
    "Percentage_of_ASE_SNPs", "ASE_over_background_percent"
]].to_string(index=False), flush=True)
print(f"\nSource Data exported to: {OUT}", flush=True)

# Save the statistics before plotting so plotting errors leave the Source Data available.
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })

    x = list(range(len(shown)))
    fig, ax1 = plt.subplots(figsize=(4.8, 3.7))

    bars = ax1.bar(
        x,
        shown["Percentage_of_ASE_SNPs"],
        color=shown["Category"].map(colors).tolist(),
        width=0.62,
        edgecolor="black",
        linewidth=0.5,
    )
    ax1.set_ylabel("Percentage of ASE SNPs (%)")
    ax1.set_xticks(x)
    ax1.set_xticklabels(shown["Display_label"], rotation=30, ha="right")
    ax1.set_ylim(0, max(5, shown["Percentage_of_ASE_SNPs"].max() * 1.18))
    ax1.tick_params(axis="x", length=0)

    for bar, value in zip(bars, shown["Percentage_of_ASE_SNPs"]):
        ax1.annotate(
            f"{value:.1f}%",
            (bar.get_x() + bar.get_width() / 2, bar.get_height()),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center", va="bottom", fontsize=8,
        )

    ax2 = ax1.twinx()
    ax2.plot(
        x, shown["ASE_over_background_percent"],
        color="#F05A4F", marker="o",
        linewidth=1.2, markersize=4,
    )
    ax2.set_ylabel("ASE / background SNPs (%)")
    ymax = shown["ASE_over_background_percent"].max()
    ax2.set_ylim(0, max(1, ymax * 1.18))

    ax1.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax1.grid(False)
    ax2.grid(False)
    fig.tight_layout()

    for ext in ["pdf", "svg", "png"]:
        fig.savefig(
            OUT / f"ASE_four_categories.{ext}",
            dpi=600, bbox_inches="tight", facecolor="white",
        )
    plt.close(fig)
    print("\nPDF, SVG and PNG files exported.", flush=True)

except Exception as exc:
    (OUT / "Plot_error.txt").write_text(repr(exc), encoding="utf-8")
    print(f"\nPlotting failed. Source Data was saved. Error: {exc}", flush=True)
    raise

print(f"\nCompleted. Output directory: {OUT}", flush=True)
