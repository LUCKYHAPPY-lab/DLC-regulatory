#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================
# 1. Input and output
# ============================================================

parser = argparse.ArgumentParser(description="Compare DLC and PigGTEx muscle cis-eQTLs.")
parser.add_argument("--dlc", type=Path, required=True)
parser.add_argument("--piggtex", type=Path, required=True)
parser.add_argument("--outdir", type=Path, required=True)
args = parser.parse_args()
DLC_FILE = args.dlc
PIGGTEX_FILE = args.piggtex
OUT_DIR = args.outdir
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_PDF = OUT_DIR / "Fig_DLC_vs_PigGTEx_reported_new_stacked_bar.pdf"
OUT_PNG = OUT_DIR / "Fig_DLC_vs_PigGTEx_reported_new_stacked_bar.png"
OUT_SVG = OUT_DIR / "Fig_DLC_vs_PigGTEx_reported_new_stacked_bar.svg"

SUMMARY_TSV = OUT_DIR / "01_DLC_vs_PigGTEx_reported_new_summary.tsv"
SOURCE_XLSX = OUT_DIR / "02_DLC_vs_PigGTEx_reported_new_sourcedata.xlsx"

SITE_ANNOTATED = OUT_DIR / "03_DLC_eQTL_sites_annotated_reported_new.tsv.gz"
PAIR_ANNOTATED = OUT_DIR / "04_DLC_eQTL_gene_pairs_annotated_reported_new.tsv.gz"
REPORTED_SITE = OUT_DIR / "05_reported_eQTL_sites.tsv.gz"
NEW_SITE = OUT_DIR / "06_DLC_new_eQTL_sites.tsv.gz"
REPORTED_PAIR = OUT_DIR / "07_reported_eQTL_gene_pairs.tsv.gz"
NEW_PAIR = OUT_DIR / "08_DLC_new_eQTL_gene_pairs.tsv.gz"

# ============================================================
# 2. Plot settings
# ============================================================

REPORTED_COLOR = "#38A9CC"
NEW_COLOR = "#C9483F"

# ============================================================
# 3. Utility functions
# ============================================================

def clean_chrom(x):
    s = str(x).strip()
    s = re.sub(r"^chr", "", s, flags=re.IGNORECASE)
    if s.endswith(".0"):
        s = s[:-2]
    return s


def strip_gene_version(x):
    if pd.isna(x):
        return np.nan
    s = str(x).strip()
    if s == "":
        return np.nan
    return re.sub(r"\.\d+$", "", s)


def parse_dlc_snp(x):
    """
    DLC example:
        14:48776571_A_G
    """
    if pd.isna(x):
        return None

    s = str(x).strip()

    m = re.match(
        r"^(?:chr)?([^:]+):(\d+)_([^_]+)_([^_]+)$",
        s,
        flags=re.IGNORECASE
    )

    if not m:
        return None

    chrom = clean_chrom(m.group(1))
    pos = int(m.group(2))
    a1 = m.group(3).upper()
    a2 = m.group(4).upper()

    alleles = sorted([a1, a2])

    return f"{chrom}:{pos}:{alleles[0]}:{alleles[1]}"


def parse_pig_variant(x):
    """
    PigGTEx examples:
        1_69021_G_A
        chr1_69021_G_A
        1:69021_G_A
    """
    if pd.isna(x):
        return None

    s = str(x).strip()

    patterns = [
        r"^(?:chr)?([^_]+)_(\d+)_([^_]+)_([^_]+)$",
        r"^(?:chr)?([^:]+):(\d+)_([^_]+)_([^_]+)$",
        r"^(?:chr)?([^:]+):(\d+):([^:]+):([^:]+)$",
    ]

    for pattern in patterns:
        m = re.match(pattern, s, flags=re.IGNORECASE)
        if m:
            chrom = clean_chrom(m.group(1))
            pos = int(m.group(2))
            a1 = m.group(3).upper()
            a2 = m.group(4).upper()

            alleles = sorted([a1, a2])

            return f"{chrom}:{pos}:{alleles[0]}:{alleles[1]}"

    return None


def require_file(path, label):
    if not path.exists():
        raise FileNotFoundError(f"{label} input not found:\n{path}")


# ============================================================
# 4. Validate inputs
# ============================================================

require_file(DLC_FILE, "DLC eQTL file")
require_file(PIGGTEX_FILE, "PigGTEx muscle eQTL file")

# ============================================================
# 5. Read DLC data
# ============================================================

print("============================================================")
print("DLC vs PigGTEx reported/new eQTL analysis")
print("============================================================")

print("\n[1] Reading DLC eQTL...")

dlc = pd.read_csv(
    DLC_FILE,
    sep="\t",
    low_memory=False
)

required_dlc = ["SNP_id", "Gene.ID"]
missing = [c for c in required_dlc if c not in dlc.columns]

if missing:
    raise ValueError(
        f"DLC input is missing columns: {missing}\n"
        f"Available columns: {list(dlc.columns)}"
    )

dlc["variant_key"] = dlc["SNP_id"].map(parse_dlc_snp)
dlc["gene_key"] = dlc["Gene.ID"].map(strip_gene_version)

dlc = dlc[
    dlc["variant_key"].notna()
    &
    dlc["gene_key"].notna()
].copy()

dlc["pair_key"] = (
    dlc["variant_key"].astype(str)
    + "|"
    + dlc["gene_key"].astype(str)
)

# site-level: unique variant
dlc_sites = (
    dlc
    .sort_values(["variant_key", "gene_key"])
    .drop_duplicates(subset=["variant_key"], keep="first")
    .copy()
)

# pair-level: unique variant-gene
dlc_pairs = (
    dlc
    .sort_values("pair_key")
    .drop_duplicates(subset=["pair_key"], keep="first")
    .copy()
)

print(f">>> DLC unique eQTL sites: {len(dlc_sites):,}")
print(f">>> DLC unique eQTL-gene pairs: {len(dlc_pairs):,}")

# ============================================================
# 6. Read PigGTEx data
# ============================================================

print("\n[2] Reading PigGTEx muscle significant eQTL...")

pig = pd.read_csv(
    PIGGTEX_FILE,
    sep="\t",
    compression="gzip",
    low_memory=False
)

if (
    "variant_id" not in pig.columns
    or
    "phenotype_id" not in pig.columns
):
    pig = pd.read_csv(
        PIGGTEX_FILE,
        sep=r"\s+",
        compression="gzip",
        low_memory=False
    )

required_pig = ["variant_id", "phenotype_id"]
missing = [c for c in required_pig if c not in pig.columns]

if missing:
    raise ValueError(
        f"PigGTEx input is missing columns: {missing}\n"
        f"Available columns: {list(pig.columns)}"
    )

pig["variant_key"] = pig["variant_id"].map(parse_pig_variant)
pig["gene_key"] = pig["phenotype_id"].map(strip_gene_version)

pig = pig[
    pig["variant_key"].notna()
    &
    pig["gene_key"].notna()
].copy()

pig["pair_key"] = (
    pig["variant_key"].astype(str)
    + "|"
    + pig["gene_key"].astype(str)
)

pig_site_set = set(pig["variant_key"].unique())
pig_pair_set = set(pig["pair_key"].unique())

print(f">>> PigGTEx unique eQTL sites: {len(pig_site_set):,}")
print(f">>> PigGTEx unique eQTL-gene pairs: {len(pig_pair_set):,}")

# ============================================================
# 7. Reported / DLC-new
# ============================================================

print("\n[3] Classifying Reported vs DLC-new...")

dlc_sites["Status"] = np.where(
    dlc_sites["variant_key"].isin(pig_site_set),
    "Reported",
    "DLC-new"
)

dlc_pairs["Status"] = np.where(
    dlc_pairs["pair_key"].isin(pig_pair_set),
    "Reported",
    "DLC-new"
)

# ============================================================
# 8. Summarize overlap
# ============================================================

def summarize(df, level):
    total = len(df)
    reported = int((df["Status"] == "Reported").sum())
    new = int((df["Status"] == "DLC-new").sum())

    return {
        "Level": level,
        "Total_DLC": total,
        "Reported_count": reported,
        "Reported_percentage": reported / total * 100 if total else np.nan,
        "DLC_new_count": new,
        "DLC_new_percentage": new / total * 100 if total else np.nan,
    }


summary = pd.DataFrame(
    [
        summarize(dlc_sites, "eQTL sites"),
        summarize(dlc_pairs, "eQTL-gene pairs")
    ]
)

summary.to_csv(
    SUMMARY_TSV,
    sep="\t",
    index=False
)

print("\nSummary:")
print(summary.to_string(index=False))

# ============================================================
# 9. Export detailed results
# ============================================================

dlc_sites.to_csv(
    SITE_ANNOTATED,
    sep="\t",
    index=False,
    compression="gzip"
)

dlc_pairs.to_csv(
    PAIR_ANNOTATED,
    sep="\t",
    index=False,
    compression="gzip"
)

dlc_sites[
    dlc_sites["Status"] == "Reported"
].to_csv(
    REPORTED_SITE,
    sep="\t",
    index=False,
    compression="gzip"
)

dlc_sites[
    dlc_sites["Status"] == "DLC-new"
].to_csv(
    NEW_SITE,
    sep="\t",
    index=False,
    compression="gzip"
)

dlc_pairs[
    dlc_pairs["Status"] == "Reported"
].to_csv(
    REPORTED_PAIR,
    sep="\t",
    index=False,
    compression="gzip"
)

dlc_pairs[
    dlc_pairs["Status"] == "DLC-new"
].to_csv(
    NEW_PAIR,
    sep="\t",
    index=False,
    compression="gzip"
)

# ============================================================
# 10. Source data
# ============================================================

plot_source = pd.DataFrame(
    {
        "Category": [
            "eQTL sites",
            "eQTL sites",
            "eQTL-gene pairs",
            "eQTL-gene pairs"
        ],
        "Status": [
            "Reported",
            "DLC-new",
            "Reported",
            "DLC-new"
        ],
        "Count": [
            int(summary.loc[0, "Reported_count"]),
            int(summary.loc[0, "DLC_new_count"]),
            int(summary.loc[1, "Reported_count"]),
            int(summary.loc[1, "DLC_new_count"])
        ],
        "Percentage": [
            float(summary.loc[0, "Reported_percentage"]),
            float(summary.loc[0, "DLC_new_percentage"]),
            float(summary.loc[1, "Reported_percentage"]),
            float(summary.loc[1, "DLC_new_percentage"])
        ]
    }
)

qc = pd.DataFrame(
    {
        "Item": [
            "DLC input",
            "PigGTEx input",
            "DLC unique eQTL sites",
            "DLC unique eQTL-gene pairs",
            "PigGTEx unique eQTL sites",
            "PigGTEx unique eQTL-gene pairs",
            "Variant harmonization",
            "Gene harmonization",
            "Pair definition"
        ],
        "Value": [
            str(DLC_FILE),
            str(PIGGTEX_FILE),
            len(dlc_sites),
            len(dlc_pairs),
            len(pig_site_set),
            len(pig_pair_set),
            "chromosome + position + unordered allele pair",
            "Ensembl gene ID with version removed",
            "harmonized variant_key + harmonized gene_key"
        ]
    }
)

with pd.ExcelWriter(
    SOURCE_XLSX,
    engine="openpyxl"
) as writer:

    plot_source.to_excel(
        writer,
        sheet_name="Plot_SourceData",
        index=False
    )

    summary.to_excel(
        writer,
        sheet_name="Summary",
        index=False
    )

    qc.to_excel(
        writer,
        sheet_name="QC_and_Method",
        index=False
    )

    dlc_sites.head(100000).to_excel(
        writer,
        sheet_name="Sites_top100k",
        index=False
    )

    dlc_pairs.head(100000).to_excel(
        writer,
        sheet_name="Pairs_top100k",
        index=False
    )

# ============================================================
# 11. Draw the figure
# ============================================================

plt.rcParams.update(
    {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.family": "Arial",
        "font.size": 9,
        "axes.linewidth": 0.8,
        "axes.labelsize": 9,
        "xtick.labelsize": 9,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
    }
)

categories = ["eQTL sites", "eQTL-gene pairs"]

reported_pct = summary["Reported_percentage"].to_numpy(dtype=float)
new_pct = summary["DLC_new_percentage"].to_numpy(dtype=float)

reported_count = summary["Reported_count"].to_numpy(dtype=int)
new_count = summary["DLC_new_count"].to_numpy(dtype=int)

x = np.arange(len(categories))

fig, ax = plt.subplots(
    figsize=(3.15, 3.15)
)

ax.bar(
    x,
    reported_pct,
    width=0.56,
    color=REPORTED_COLOR,
    edgecolor="none",
    label="Reported"
)

ax.bar(
    x,
    new_pct,
    width=0.56,
    bottom=reported_pct,
    color=NEW_COLOR,
    edgecolor="none",
    label="DLC-new"
)

for i in range(len(categories)):

    if reported_pct[i] > 4:
        ax.text(
            x[i],
            reported_pct[i] / 2,
            f"{reported_count[i]:,}",
            ha="center",
            va="center",
            rotation=90,
            fontsize=8
        )

    if new_pct[i] > 4:
        ax.text(
            x[i],
            reported_pct[i] + new_pct[i] / 2,
            f"{new_count[i]:,}",
            ha="center",
            va="center",
            rotation=90,
            fontsize=8
        )

ax.set_ylim(0, 100)
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_ylabel("Percentage (%)")

ax.set_xticks(x)
ax.set_xticklabels(categories)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.legend(
    loc="upper center",
    bbox_to_anchor=(0.5, 1.12),
    ncol=2,
    frameon=False,
    handlelength=1.8,
    handletextpad=0.4,
    columnspacing=1.8
)

fig.tight_layout(pad=0.8)

fig.savefig(
    OUT_PDF,
    format="pdf",
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    OUT_PNG,
    format="png",
    dpi=600,
    bbox_inches="tight"
)

fig.savefig(
    OUT_SVG,
    format="svg",
    bbox_inches="tight"
)

plt.close(fig)

# ============================================================
# 12. Report results
# ============================================================

print("\n============================================================")
print("COMPLETED")
print("============================================================")
print(summary.to_string(index=False))

print(f"\nOutput directory:\n{OUT_DIR}")
print(f"\nFigure PDF:\n{OUT_PDF}")
print(f"\nSource data Excel:\n{SOURCE_XLSX}")
