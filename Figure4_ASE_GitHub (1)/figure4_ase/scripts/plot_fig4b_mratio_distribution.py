#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
from collections import defaultdict
from bisect import bisect_right
import re

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ============================================================
# 1. Paths and parameters
# ============================================================

parser = argparse.ArgumentParser(description="Generate Figure 4b: mratio distribution.")
parser.add_argument("--ase-table", required=True, type=Path)
parser.add_argument("--gtf", required=True, type=Path)
parser.add_argument("--dlc-frequency", required=True, type=Path)
parser.add_argument("--outdir", required=True, type=Path)
args = parser.parse_args()
INPUT = args.ase_table
GTF = args.gtf
FREQUENCY = args.dlc_frequency
OUTDIR = args.outdir

FDR_THRESHOLD = 0.05
MAF_THRESHOLD = 0.05
MRATIO_THRESHOLD = 0.60

BIN_WIDTH = 0.0125
HIST_MIN = 0.45
HIST_MAX = 1.00

COLOR_LOW = "#4DBBD5"
COLOR_ASE = "#E64B35"


# ============================================================
# 2. Utility functions
# ============================================================

def normalize_chrom(value):
    value = re.sub(r"^chr", "", str(value).strip(), flags=re.I)
    if value.endswith(".0"):
        value = value[:-2]
    return {"23": "X", "M": "MT"}.get(value, value)


def normalize_gene(value):
    if pd.isna(value):
        return ""
    value = str(value).strip()
    if value in {"", ".", "-", "NA", "nan"}:
        return ""
    return re.sub(r"\.\d+$", "", value)


def parse_snp(value):
    value = str(value).strip()
    match = re.fullmatch(r"([^:]+):(\d+)_([^_]+)_([^_]+)", value)
    if not match:
        raise ValueError(f"Cannot parse SNP_id: {value}")
    chrom, pos, ref, alt = match.groups()
    return normalize_chrom(chrom), int(pos), ref.upper(), alt.upper()


def allele_pair(ref, alt):
    return tuple(sorted((str(ref).upper(), str(alt).upper())))


def save_tsv(table, name):
    table.to_csv(
        OUTDIR / name, sep="\t", index=False, na_rep="NA"
    )


def read_input(path):
    if path.suffix.lower() in {".xlsx", ".xlsm"}:
        return pd.read_excel(path, sheet_name=0)

    # Preserve empty fields in tab-delimited input.
    with path.open("r", encoding="utf-8-sig") as handle:
        header = handle.readline()

    if "\t" in header:
        return pd.read_csv(
            path, sep="\t", encoding="utf-8-sig",
            low_memory=False
        )

    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)

    print("No tabs found in the header. Reading whitespace-delimited input.")
    return pd.read_csv(
        path, sep=r"\s+", engine="python",
        encoding="utf-8-sig"
    )


def make_interval_index(records):
    records = sorted(records)
    starts = []
    prefix_max = []
    maximum = 0

    for start, end in records:
        starts.append(start)
        maximum = max(maximum, end)
        prefix_max.append(maximum)

    return records, starts, prefix_max


def overlaps_exon(index, pos):
    if index is None:
        return False

    records, starts, prefix_max = index
    i = bisect_right(starts, pos) - 1

    while i >= 0 and prefix_max[i] >= pos:
        if records[i][0] <= pos <= records[i][1]:
            return True
        i -= 1

    return False


# ============================================================
# 3. Main workflow
# ============================================================

def main():
    for path in [INPUT, GTF, FREQUENCY]:
        if not path.is_file():
            raise FileNotFoundError(f"File not found: {path}")

    OUTDIR.mkdir(parents=True, exist_ok=True)

    # Protect the results of a completed run from being overwritten.
    if (OUTDIR / "Summary.tsv").exists():
        raise FileExistsError(
            f"{OUTDIR / 'Summary.tsv'} already exists. "
            "Use a different --outdir to preserve the previous results."
        )

    summary = []

    def record(name, value):
        summary.append((name, value))

    print("1. Reading the full input table...", flush=True)
    df = read_input(INPUT)
    df.columns = [
        str(c).strip().removesuffix("\\") for c in df.columns
    ]

    if df.columns.duplicated().any():
        raise ValueError("The input contains duplicate column names.")

    required = {
        "SNP_id", "Gene", "FDR_binomtest", "median_Mratio"
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing columns: {sorted(missing)}\nAvailable columns: {list(df.columns)}"
        )

    df["Source_row"] = np.arange(2, len(df) + 2)
    df["SNP_id"] = df["SNP_id"].astype("string").str.strip()
    df["Gene_normalized"] = df["Gene"].map(normalize_gene)

    fdr = pd.to_numeric(df["FDR_binomtest"], errors="coerce")
    mratio = pd.to_numeric(df["median_Mratio"], errors="coerce")

    invalid_fdr = ~fdr.between(0, 1)
    invalid_mratio = ~mratio.between(0.5, 1)
    missing_snp = df["SNP_id"].isna() | df["SNP_id"].eq("")

    invalid = invalid_fdr | invalid_mratio | missing_snp

    bad = df.loc[invalid].copy()
    bad["Invalid_FDR"] = invalid_fdr.loc[invalid]
    bad["Invalid_Mratio"] = invalid_mratio.loc[invalid]
    bad["Missing_SNP_id"] = missing_snp.loc[invalid]
    save_tsv(bad, "01_invalid_filter_values.tsv")

    record("Input_annotation_rows", len(df))
    record("Input_unique_SNPs", df["SNP_id"].nunique())
    record("Invalid_FDR_rows", int(invalid_fdr.sum()))
    record("Invalid_Mratio_rows", int(invalid_mratio.sum()))
    record("Invalid_any_rows", int(invalid.sum()))

    eligible = df.loc[
        ~invalid & (fdr < FDR_THRESHOLD)
    ].copy()

    eligible["FDR_binomtest"] = fdr.loc[eligible.index]
    eligible["median_Mratio"] = mratio.loc[eligible.index]

    save_tsv(eligible, "02_FDR005_annotations.tsv")

    record("FDR005_annotation_rows", len(eligible))
    record("FDR005_unique_SNPs", eligible["SNP_id"].nunique())

    if eligible.empty:
        raise ValueError("No records passed the FDR filter.")

    # Check that repeated records for each SNP agree on FDR and Mratio before deduplication.
    conflicts = []
    for col in ["FDR_binomtest", "median_Mratio"]:
        counts = eligible.groupby("SNP_id")[col].nunique()
        for snp in counts[counts > 1].index:
            conflicts.append({"SNP_id": snp, "Column": col})

    if conflicts:
        save_tsv(
            pd.DataFrame(conflicts),
            "ERROR_conflicting_SNP_statistics.tsv"
        )
        raise ValueError("Conflicting FDR or Mratio values for the same SNP. Conflicting records were exported.")

    coordinates = {
        snp: parse_snp(snp)
        for snp in eligible["SNP_id"].drop_duplicates()
    }

    # Report the sample-size distribution without applying an additional filter.
    if "num_all" in eligible.columns:
        sample_info = (
            eligible[["SNP_id", "num_all"]]
            .drop_duplicates()
        )
        save_tsv(sample_info, "02b_SNP_num_all.tsv")

        sample_counts = pd.to_numeric(
            sample_info["num_all"], errors="coerce"
        )
        record(
            "FDR005_SNPs_with_num_all_below19",
            sample_info.loc[
                sample_counts < 19, "SNP_id"
            ].nunique()
        )

    # ========================================================
    # 4. Check GTF exon overlap for each SNP–gene pair
    # ========================================================

    print("2. Reading the GTF and checking exon overlaps...", flush=True)

    wanted_genes = set(eligible["Gene_normalized"]) - {""}
    exon_records = defaultdict(set)
    gene_names = {}

    with GTF.open("r") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue

            attrs = dict(re.findall(
                r'(\S+)\s+"([^"]*)"', fields[8]
            ))
            gene = normalize_gene(attrs.get("gene_id"))

            if gene not in wanted_genes:
                continue

            chrom = normalize_chrom(fields[0])
            exon_records[(chrom, gene)].add(
                (int(fields[3]), int(fields[4]))
            )
            gene_names[gene] = attrs.get("gene_name", gene)

    indexes = {
        key: make_interval_index(records)
        for key, records in exon_records.items()
    }

    pair_check = eligible[
        ["SNP_id", "Gene_normalized"]
    ].drop_duplicates().copy()

    pair_check["GTF_exonic"] = [
        overlaps_exon(
            indexes.get((coordinates[snp][0], gene)),
            coordinates[snp][1]
        )
        for snp, gene in zip(
            pair_check["SNP_id"],
            pair_check["Gene_normalized"]
        )
    ]

    save_tsv(pair_check, "03_all_SNP_gene_pairs_GTF_check.tsv")

    exonic_pairs = pair_check.loc[
        pair_check["GTF_exonic"]
    ].copy()

    exonic_annotations = eligible.merge(
        exonic_pairs,
        on=["SNP_id", "Gene_normalized"],
        how="inner", validate="many_to_one"
    )
    save_tsv(
        exonic_annotations,
        "04_GTF_exonic_annotations.tsv"
    )

    no_exon = eligible.loc[
        ~eligible["SNP_id"].isin(exonic_pairs["SNP_id"])
    ]
    save_tsv(
        no_exon,
        "04b_SNPs_without_supported_exonic_gene.tsv"
    )

    if exonic_annotations.empty:
        raise ValueError("No GTF exon overlaps found. Check gene IDs and chromosome names.")

    record(
        "FDR005_GTF_exonic_unique_SNPs",
        exonic_annotations["SNP_id"].nunique()
    )

    # Combine all associated genes into one row per site.
    plot_sites = (
        exonic_annotations[
            ["SNP_id", "FDR_binomtest", "median_Mratio"]
        ]
        .drop_duplicates("SNP_id")
        .copy()
    )

    gene_summary = (
        exonic_pairs.groupby("SNP_id")["Gene_normalized"]
        .agg(lambda values: ";".join(sorted(set(values))))
    )
    gene_count = (
        exonic_pairs.groupby("SNP_id")["Gene_normalized"]
        .nunique()
    )

    plot_sites["ASE_gene_ids"] = plot_sites["SNP_id"].map(gene_summary)
    plot_sites["Candidate_gene_count"] = plot_sites["SNP_id"].map(
        gene_count
    )

    # ========================================================
    # 5. Calculate the maximum DLC MAF for each chromosome, position and reference allele
    # ========================================================

    print("3. Reading the DLC frequency file in chunks...", flush=True)

    wanted_positions = {
        coordinates[snp][:2] for snp in plot_sites["SNP_id"]
    }
    parts = []

    for chunk in pd.read_csv(
        FREQUENCY, sep="\t", dtype=str,
        chunksize=500_000, keep_default_na=False
    ):
        chunk.columns = [
            str(c).strip().removesuffix("\\")
            for c in chunk.columns
        ]

        required_freq = {
            "POPULATION", "CHROM", "POS", "REF", "ALT", "MAF"
        }
        if not required_freq.issubset(chunk.columns):
            raise ValueError("The frequency file is missing required columns. Check the delimiter.")

        chunk = chunk.loc[
            chunk["POPULATION"].str.strip().eq("DLC")
        ].copy()
        chunk["CHROM"] = chunk["CHROM"].map(normalize_chrom)

        pos = pd.to_numeric(chunk["POS"], errors="coerce")
        valid_pos = pos.notna() & pos.mod(1).eq(0)
        chunk = chunk.loc[valid_pos].copy()
        chunk["POS"] = pos.loc[valid_pos].astype("int64")

        mask = [
            (chrom, position) in wanted_positions
            for chrom, position in zip(
                chunk["CHROM"], chunk["POS"]
            )
        ]

        subset = chunk.loc[
            mask, ["CHROM", "POS", "REF", "ALT", "MAF"]
        ].copy()

        if not subset.empty:
            parts.append(subset)

    if not parts:
        raise ValueError("No matching DLC frequency records were found.")

    freq = pd.concat(parts, ignore_index=True)

    for col in ["REF", "ALT"]:
        freq[col] = freq[col].str.strip().str.upper()

    freq["MAF"] = pd.to_numeric(freq["MAF"], errors="coerce")

    if (~freq["MAF"].between(0, 1)).any():
        save_tsv(freq, "ERROR_frequency_records.tsv")
        raise ValueError("Relevant frequency records contain missing or invalid MAF values.")

    freq_keys = ["CHROM", "POS", "REF", "ALT"]
    conflict = freq.groupby(freq_keys)["MAF"].nunique()
    if (conflict > 1).any():
        raise ValueError("Frequency records contain conflicting MAF values for the same allele.")

    freq = freq.drop_duplicates(freq_keys)
    save_tsv(freq, "05_matched_coordinate_frequency_records.tsv")

    max_maf = (
        freq.groupby(["CHROM", "POS", "REF"])["MAF"]
        .max().to_dict()
    )

    lookup = defaultdict(list)
    for row in freq.itertuples(index=False):
        key = (
            row.CHROM, row.POS,
            allele_pair(row.REF, row.ALT)
        )
        lookup[key].append((row.REF, row.ALT, row.MAF))

    audit = []
    for snp in plot_sites["SNP_id"]:
        chrom, pos, ref, alt = coordinates[snp]
        matches = lookup.get(
            (chrom, pos, allele_pair(ref, alt)), []
        )
        forward = [
            match for match in matches
            if match[0] == ref and match[1] == alt
        ]
        matches = forward if forward else matches

        if len(matches) > 1:
            raise ValueError(f"Multiple frequency matches for SNP: {snp}")

        if not matches:
            audit.append({
                "SNP_id": snp,
                "DLC_allele_MAF": np.nan,
                "DLC_site_max_MAF": np.nan,
                "Frequency_status": "Unmatched",
            })
        else:
            matched_ref, matched_alt, maf = matches[0]
            maximum = max_maf[(chrom, pos, matched_ref)]

            audit.append({
                "SNP_id": snp,
                "DLC_allele_MAF": maf,
                "DLC_site_max_MAF": maximum,
                "Frequency_status": (
                    "Retained" if maximum >= MAF_THRESHOLD
                    else "Below_MAF005"
                ),
            })

    plot_sites = plot_sites.merge(
        pd.DataFrame(audit),
        on="SNP_id", validate="one_to_one"
    )
    save_tsv(plot_sites, "06_all_sites_frequency_audit.tsv")

    excluded = plot_sites.loc[
        plot_sites["Frequency_status"] != "Retained"
    ]
    save_tsv(excluded, "06b_frequency_excluded_sites.tsv")

    plot_sites = plot_sites.loc[
        plot_sites["Frequency_status"] == "Retained"
    ].copy()

    if plot_sites.empty:
        raise ValueError("No sites passed the MAF filter.")

    plot_sites["ASE_status"] = np.where(
        plot_sites["median_Mratio"] >= MRATIO_THRESHOLD,
        "ASE", "Below_Mratio_threshold"
    )

    save_tsv(plot_sites, "07_plot_unique_SNPs.tsv")

    ase_sites = plot_sites.loc[
        plot_sites["ASE_status"] == "ASE"
    ].copy()
    save_tsv(ase_sites, "08_final_ASE_unique_SNPs.tsv")

    final_pairs = exonic_pairs.loc[
        exonic_pairs["SNP_id"].isin(ase_sites["SNP_id"])
    ].rename(columns={"Gene_normalized": "Gene"})

    final_pairs["Gene_name"] = final_pairs["Gene"].map(gene_names)
    save_tsv(final_pairs, "09_final_ASE_SNP_gene_pairs.tsv")

    # ========================================================
    # 6. Binning and summary statistics
    # ========================================================

    print("4. Generating the plot and Source Data...", flush=True)

    bins = np.round(
        np.arange(HIST_MIN, HIST_MAX + BIN_WIDTH / 2, BIN_WIDTH),
        10
    )
    counts, edges = np.histogram(
        plot_sites["median_Mratio"], bins=bins
    )

    if counts.sum() != len(plot_sites):
        raise ValueError("Some Mratio values fall outside the plotted range.")

    centers = (edges[:-1] + edges[1:]) / 2
    cumulative = np.cumsum(counts) / counts.sum() * 100

    colors = np.where(
        edges[:-1] >= MRATIO_THRESHOLD,
        COLOR_ASE, COLOR_LOW
    )

    histogram = pd.DataFrame({
        "Bin_left": edges[:-1],
        "Bin_right": edges[1:],
        "Bin_center": centers,
        "SNP_count": counts,
        "Cumulative_percent": cumulative,
    })
    save_tsv(histogram, "10_histogram_bins.tsv")

    n_total = len(plot_sites)
    n_ase = len(ase_sites)
    percentage = 100 * n_ase / n_total

    record("FDR_threshold", "< 0.05")
    record("DLC_site_max_MAF_threshold", ">= 0.05")
    record("Mratio_ASE_threshold", ">= 0.60")
    record("Additional_num_all_filter", "Not applied")
    record("Plot_unique_SNPs", n_total)
    record("Final_ASE_unique_SNPs", n_ase)
    record("Below_Mratio_threshold_SNPs", n_total - n_ase)
    record("ASE_percentage_within_plot_set", percentage)
    record("Final_ASE_SNP_gene_pairs", len(final_pairs))
    record("Final_ASE_candidate_genes", final_pairs["Gene"].nunique())
    record("Unmatched_frequency_sites",
           int((excluded["Frequency_status"] == "Unmatched").sum()))
    record("Histogram_bin_width", BIN_WIDTH)

    summary_df = pd.DataFrame(
        summary, columns=["Statistic", "Value"]
    )

    # Use Record_type to distinguish sections in the single-sheet Source Data workbook.
    source = pd.concat([
        plot_sites.assign(Record_type="SNP"),
        histogram.assign(Record_type="Histogram_bin"),
        summary_df.assign(Record_type="Summary"),
    ], ignore_index=True, sort=False)

    source = source[
        ["Record_type"]
        + [col for col in source.columns if col != "Record_type"]
    ]

    with pd.ExcelWriter(
        OUTDIR / "Mratio_SourceData_single_sheet.xlsx",
        engine="openpyxl"
    ) as writer:
        source.to_excel(
            writer, sheet_name="Source_Data",
            index=False, na_rep="NA"
        )

    # Export the final ASE results in a separate workbook.
    with pd.ExcelWriter(
        OUTDIR / "Final_ASE_FDR005_Mratio060_MAF005.xlsx",
        engine="openpyxl"
    ) as writer:
        ase_sites.to_excel(
            writer, sheet_name="Final_ASE_sites",
            index=False, na_rep="NA"
        )
        final_pairs.to_excel(
            writer, sheet_name="SNP_gene_pairs",
            index=False, na_rep="NA"
        )
        summary_df.to_excel(
            writer, sheet_name="Summary", index=False
        )

    # ========================================================
    # 7. Scale the left axis to keep the full height of every bar visible
    # ========================================================

    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "DejaVu Sans",
        "axes.unicode_minus": False,
    })

    fig, ax1 = plt.subplots(figsize=(8, 5))

    ax1.bar(
        edges[:-1], counts,
        width=np.diff(edges), align="edge",
        color=colors, edgecolor="black",
        linewidth=0.6, alpha=0.8
    )
    ax1.axvline(
        MRATIO_THRESHOLD, color="red",
        linestyle="--", linewidth=1.3
    )

    ax1.set_xlabel("Mratio")
    ax1.set_ylabel("Number of SNPs")
    ax1.set_xlim(HIST_MIN, HIST_MAX)

    # Use six evenly spaced ticks on each y-axis.
    maximum = max(5, float(counts.max()) * 1.2)
    magnitude = 10 ** np.floor(np.log10(maximum / 5))
    step = np.ceil((maximum / 5) / magnitude) * magnitude
    ymax = step * 5

    ax1.set_ylim(0, ymax)
    ax1.set_yticks(np.linspace(0, ymax, 6))

    ax2 = ax1.twinx()

    # Plot each cumulative percentage at the right edge of its bin.
    ax2.plot(
        edges,
        np.r_[0, cumulative],
        color="blue", linewidth=1.8
    )
    ax2.set_ylabel("Cumulative percentage of SNPs", color="blue")
    ax2.set_ylim(0, 100)
    ax2.set_yticks(np.arange(0, 101, 20))
    ax2.tick_params(axis="y", colors="blue")

    ax1.text(
        0.97, 0.93,
        f"ASE SNPs: {n_ase:,} ({percentage:.1f}%)",
        transform=ax1.transAxes,
        ha="right", va="top", color=COLOR_ASE
    )

    ax1.set_title(
        "Mratio distribution\n"
        "FDR < 0.05; DLC site max MAF ≥ 0.05; exonic sites"
    )

    fig.tight_layout()
    fig.savefig(
        OUTDIR / "Mratio_distribution.pdf",
        bbox_inches="tight"
    )
    fig.savefig(
        OUTDIR / "Mratio_distribution.png",
        dpi=600, bbox_inches="tight"
    )
    plt.close(fig)

    # Write the completion marker and summary last.
    save_tsv(summary_df, "Summary.tsv")

    print("\n" + summary_df.to_string(index=False))
    print(f"\nAll results saved to: {OUTDIR}")


if __name__ == "__main__":
    main()