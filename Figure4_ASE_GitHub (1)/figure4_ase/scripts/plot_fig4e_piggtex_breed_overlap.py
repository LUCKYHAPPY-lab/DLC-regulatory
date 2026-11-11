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
from matplotlib.ticker import MaxNLocator


# ============================================================
# 1. Input and output
# ============================================================

parser = argparse.ArgumentParser(description="Generate Figure 4e: piggtex breed overlap.")
parser.add_argument("--dlc-ase", required=True, type=Path)
parser.add_argument("--piggtex-ase", required=True, type=Path)
parser.add_argument("--outdir", required=True, type=Path)
args = parser.parse_args()
DLC_FILE = args.dlc_ase
PIGGTEX_FILE = args.piggtex_ase
OUTDIR = args.outdir
DLC_SHEET = "Final_ASE_sites"

PIGGTEX_FDR = 0.05
MIN_TOTAL_COUNT = 10

COLOR_DUROC = "#B8322A"
COLOR_OTHER = "#3C91B1"


# ============================================================
# 2. Normalize variant identifiers
# ============================================================

def clean_chr(value):
    value = re.sub(
        r"^chr", "", str(value).strip(), flags=re.I
    )
    if value.endswith(".0"):
        value = value[:-2]
    return {"23": "X", "M": "MT"}.get(value, value)


def parse_snp(value):
    if pd.isna(value):
        raise ValueError("The DLC catalog contains a missing SNP_id.")

    match = re.fullmatch(
        r"(?:chr)?([^:_]+)[:_](\d+)[:_]([^:_]+)[:_]([^:_]+)",
        str(value).strip(),
        flags=re.I
    )
    if not match:
        raise ValueError(f"Cannot parse SNP_id: {value}")

    chrom, pos, ref, alt = match.groups()

    return (
        clean_chr(chrom),
        int(pos),
        ref.upper(),
        alt.upper()
    )


def allele_key(chrom, pos, ref, alt):
    a1, a2 = sorted((ref, alt))
    return f"{chrom}:{int(pos)}:{a1}:{a2}"


def save_tsv(table, name):
    table.to_csv(
        OUTDIR / name,
        sep="\t",
        index=False,
        na_rep="NA"
    )


# ============================================================
# 3. Main workflow
# ============================================================

def main():
    for path in [DLC_FILE, PIGGTEX_FILE]:
        if not path.is_file():
            raise FileNotFoundError(f"File not found: {path}")

    OUTDIR.mkdir(parents=True, exist_ok=True)

    if (OUTDIR / "Summary.tsv").exists():
        raise FileExistsError(
            "The output directory already contains completed results. "
            "Use a different --outdir to preserve them."
        )

    # --------------------------------------------------------
    # 3.1 Load the final DLC ASE catalog
    # --------------------------------------------------------

    print("1. Reading the final DLC ASE catalog...", flush=True)

    dlc = pd.read_excel(
        DLC_FILE, sheet_name=DLC_SHEET
    )
    dlc.columns = dlc.columns.astype(str).str.strip()

    if "SNP_id" not in dlc.columns:
        raise ValueError("The DLC file is missing the SNP_id column.")

    parsed = dlc["SNP_id"].map(parse_snp)

    dlc["allele_key"] = [
        allele_key(*value) for value in parsed
    ]

    # Count each unique chromosome, position and unordered allele pair once.
    dlc = dlc.drop_duplicates("allele_key").copy()
    dlc_set = set(dlc["allele_key"])

    save_tsv(dlc, "00_DLC_final_ASE_input.tsv")

    print(f"Unique DLC ASE sites: {len(dlc_set):,}")

    # --------------------------------------------------------
    # 3.2 Read PigGTEx records in chunks
    # --------------------------------------------------------

    print("2. Reading PigGTEx muscle ASE records in chunks...", flush=True)

    # Detect the delimiter to preserve empty fields in tab-delimited input.
    import gzip

    with gzip.open(
        PIGGTEX_FILE, "rt", encoding="utf-8-sig"
    ) as handle:
        header = handle.readline()

    sep = "\t" if "\t" in header else r"\s+"

    required = {
        "chr", "position", "refAllele", "altAllele",
        "totalCount", "FDR", "SampleID", "Breed", "Tissue"
    }

    parts = []
    input_rows = 0
    invalid_identity_rows = 0

    for chunk in pd.read_csv(
        PIGGTEX_FILE,
        sep=sep,
        compression="gzip",
        dtype=str,
        chunksize=500_000
    ):
        chunk.columns = chunk.columns.astype(str).str.strip()

        missing = required - set(chunk.columns)
        if missing:
            raise ValueError(
                f"Missing PigGTEx columns: {sorted(missing)}"
            )

        input_rows += len(chunk)

        fdr = pd.to_numeric(
            chunk["FDR"], errors="coerce"
        )
        total = pd.to_numeric(
            chunk["totalCount"], errors="coerce"
        )

        breed = chunk["Breed"].astype("string").str.strip()
        tissue = chunk["Tissue"].astype("string").str.strip()

        keep = (
            fdr.between(0, 1)
            & (fdr < PIGGTEX_FDR)
            & np.isfinite(total)
            & (total >= MIN_TOTAL_COUNT)
            & breed.notna()
            & breed.ne("")
            & ~breed.str.casefold().eq("cross")
            & tissue.str.contains(
                "muscle", case=False, na=False, regex=False
            )
        ).fillna(False)

        sub = chunk.loc[keep].copy()
        if sub.empty:
            continue

        sub["Breed"] = breed.loc[keep]
        sub["Tissue"] = tissue.loc[keep]
        sub["FDR"] = fdr.loc[keep]
        sub["totalCount"] = total.loc[keep]

        pos = pd.to_numeric(
            sub["position"], errors="coerce"
        )

        valid = (
            pos.notna()
            & np.isfinite(pos)
            & pos.mod(1).eq(0)
            & (pos >= 1)
        )

        for col in ["chr", "refAllele", "altAllele"]:
            text = sub[col].astype("string").str.strip()
            valid &= (
                text.notna()
                & ~text.isin(["", ".", "NA", "nan"])
            )

        invalid_identity_rows += int((~valid).sum())
        sub = sub.loc[valid].copy()

        if sub.empty:
            continue

        sub["position"] = pos.loc[valid].astype("int64")
        sub["chr"] = sub["chr"].map(clean_chr)

        sub["refAllele"] = (
            sub["refAllele"].str.strip().str.upper()
        )
        sub["altAllele"] = (
            sub["altAllele"].str.strip().str.upper()
        )

        sub["allele_key"] = [
            allele_key(chrom, position, ref, alt)
            for chrom, position, ref, alt in zip(
                sub["chr"], sub["position"],
                sub["refAllele"], sub["altAllele"]
            )
        ]

        parts.append(sub)

    if not parts:
        raise ValueError("No PigGTEx muscle ASE records passed the filters.")

    pig = pd.concat(parts, ignore_index=True)
    pig_set = set(pig["allele_key"])

    # Save all filtered records in a compressed file for verification.
    pig.to_csv(
        OUTDIR / "01_PigGTEx_nonCross_muscle_filtered.tsv.gz",
        sep="\t", index=False,
        compression="gzip", na_rep="NA"
    )

    print(f"Records passing filters: {len(pig):,}")
    print(f"Breeds with muscle ASE records: {pig['Breed'].nunique():,}")
    print(f"Unique muscle ASE sites: {len(pig_set):,}")

    # --------------------------------------------------------
    # 3.3 Calculate overlap counts by breed
    # --------------------------------------------------------

    rows = []

    for breed, group in pig.groupby("Breed"):
        breed_set = set(group["allele_key"])
        shared = dlc_set & breed_set

        rows.append({
            "Breed": breed,
            "PigGTEx_sample_count": group["SampleID"].nunique(),
            "PigGTEx_ASE_record_count": len(group),
            "PigGTEx_unique_ASE_sites": len(breed_set),
            "DLC_unique_ASE_sites": len(dlc_set),
            "Overlap_ASE_sites": len(shared),
            "DLC_sites_not_in_this_breed": len(dlc_set - breed_set),
            "Overlap_percent_of_DLC": (
                100 * len(shared) / len(dlc_set)
            ),
        })

    breed_table = pd.DataFrame(rows).sort_values(
        ["Overlap_ASE_sites", "Breed"],
        ascending=[False, True]
    ).reset_index(drop=True)

    save_tsv(
        breed_table,
        "02_breed_overlap_counts_nonCross.tsv"
    )

    # Keep one row per breed and overlapping site.
    overlap_pairs = (
        pig.loc[
            pig["allele_key"].isin(dlc_set),
            ["Breed", "allele_key"]
        ]
        .drop_duplicates()
        .merge(
            dlc[["SNP_id", "allele_key"]],
            on="allele_key",
            how="left",
            validate="many_to_one"
        )
        .sort_values(["Breed", "allele_key"])
    )
    save_tsv(
        overlap_pairs,
        "03_breed_overlap_unique_sites.tsv"
    )

    overlap_records = pig.loc[
        pig["allele_key"].isin(dlc_set)
    ]
    overlap_records.to_csv(
        OUTDIR / "04_overlap_sample_records.tsv.gz",
        sep="\t", index=False,
        compression="gzip", na_rep="NA"
    )

    overall_overlap = dlc_set & pig_set

    # Flag DLC sites reported in muscle for any PigGTEx breed other than Cross.
    dlc_status = dlc.copy()
    dlc_status["Reported_in_PigGTEx_nonCross_muscle"] = (
        dlc_status["allele_key"].isin(pig_set)
    )

    breed_map = (
        overlap_pairs.groupby("allele_key")["Breed"]
        .agg(lambda values: ";".join(sorted(set(values))))
    )
    dlc_status["Overlapping_breeds"] = (
        dlc_status["allele_key"].map(breed_map)
    )
    save_tsv(
        dlc_status,
        "05_DLC_sites_overlap_status.tsv"
    )

    # --------------------------------------------------------
    # 3.4 Plot the number of overlapping sites by breed
    # --------------------------------------------------------

    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.linewidth": 0.8,
    })

    plot = breed_table.copy()
    plot["Plot_rank"] = np.arange(1, len(plot) + 1)

    colors = [
        COLOR_DUROC if breed.casefold() == "duroc"
        else COLOR_OTHER
        for breed in plot["Breed"]
    ]

    fig_height = max(3.2, 0.22 * len(plot) + 0.7)
    fig, ax = plt.subplots(
        figsize=(4.2, fig_height)
    )

    y = np.arange(len(plot))
    x = plot["Overlap_ASE_sites"].to_numpy()

    ax.barh(
        y, x,
        color=colors,
        height=0.75,
        edgecolor="none"
    )

    ax.set_yticks(y)
    ax.set_yticklabels(plot["Breed"])
    ax.invert_yaxis()

    max_x = max(1, int(x.max()))

    for position, count in zip(y, x):
        ax.text(
            count + max_x * 0.015,
            position,
            f"{count:,}",
            ha="left", va="center",
            fontsize=7
        )

    ax.set_xlim(0, max_x * 1.20)
    ax.set_xlabel("Number of overlapping ASE SNPs")
    ax.xaxis.set_major_locator(
        MaxNLocator(nbins=5, integer=True)
    )
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", direction="out")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)

    fig.tight_layout()

    fig.savefig(
        OUTDIR / "Breed_ASE_overlap_counts.pdf",
        bbox_inches="tight"
    )
    fig.savefig(
        OUTDIR / "Breed_ASE_overlap_counts.png",
        dpi=600, bbox_inches="tight",
        facecolor="white"
    )
    plt.close(fig)

    # --------------------------------------------------------
    # 3.5 Record parameters and export single-sheet Source Data
    # --------------------------------------------------------

    summary = pd.DataFrame([
        ("DLC_input", str(DLC_FILE)),
        ("DLC_sheet", DLC_SHEET),
        ("PigGTEx_input", str(PIGGTEX_FILE)),
        ("PigGTEx_FDR", "< 0.05"),
        ("PigGTEx_totalCount", ">= 10"),
        ("Excluded_breed", "Cross; case insensitive"),
        ("Tissue_filter", "Contains muscle; case insensitive"),
        ("Matching", "CHROM + POS + unordered allele pair"),
        ("PigGTEx_input_records", input_rows),
        ("Filtered_records_with_invalid_identity",
         invalid_identity_rows),
        ("DLC_unique_ASE_sites", len(dlc_set)),
        ("PigGTEx_nonCross_muscle_records", len(pig)),
        ("PigGTEx_nonCross_muscle_unique_sites", len(pig_set)),
        ("PigGTEx_nonCross_muscle_breeds",
         pig["Breed"].nunique()),
        ("DLC_overlap_with_any_breed", len(overall_overlap)),
        ("DLC_not_reported_in_filtered_catalog",
         len(dlc_set - pig_set)),
        ("DLC_overlap_percent",
         100 * len(overall_overlap) / len(dlc_set)),
        ("Figure_measure", "Unique overlapping ASE SNP count"),
    ], columns=["Statistic", "Value"])

    source = pd.concat([
        plot.assign(Record_type="Breed_bar"),
        summary.assign(Record_type="Summary"),
    ], ignore_index=True, sort=False)

    source = source[
        ["Record_type"]
        + [col for col in source.columns if col != "Record_type"]
    ]

    with pd.ExcelWriter(
        OUTDIR / "PigGTEx_ASE_overlap_SourceData_single_sheet.xlsx",
        engine="openpyxl"
    ) as writer:
        source.to_excel(
            writer,
            sheet_name="Source_Data",
            index=False,
            na_rep="NA"
        )

    save_tsv(summary, "Summary.tsv")

    print("\nOverall summary:")
    print(summary.to_string(index=False))

    print("\nOverlap by breed:")
    print(
        breed_table[
            ["Breed", "PigGTEx_sample_count", "Overlap_ASE_sites"]
        ].to_string(index=False)
    )

    print(f"\nResults directory: {OUTDIR}")


if __name__ == "__main__":
    main()