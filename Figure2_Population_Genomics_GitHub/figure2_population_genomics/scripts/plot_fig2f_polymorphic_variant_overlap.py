#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Use venn3 for compatibility across matplotlib-venn versions.
# Some versions of venn3_unweighted pass both normalize_to
# and layout_algorithm, which causes an error.
from matplotlib_venn import venn3 as draw_venn3

UNWEIGHTED = False


# ============================================================
# 1. Input paths and filtering schemes
# ============================================================

parser = argparse.ArgumentParser(description="Plot overlap among three population position sets.")
parser.add_argument("--input-root", type=Path, required=True)
parser.add_argument("--outdir", type=Path, required=True)
parser.add_argument("--scheme", choices=("all", "1", "2", "3", "4", "5", "6"), default="all")
args = parser.parse_args()
BASE_DIR = args.input_root
OUTPUT_ROOT = args.outdir

SCHEMES = [
    {
        "id": "scheme1_NS_ge30",
        "short": "NS ≥ 30",
        "condition": "MAF ≥ 0.05, NS ≥ 30",
        "directory": "NS_ge30_no_genotype_count_requirement",
        "files": {
            "Duroc": "Duroc_MAF005_NS_ge30_pos.txt",
            "Luchuan": "Luchuan_MAF005_NS_ge30_pos.txt",
            "DLC": "DLC_MAF005_NS_ge30_pos.txt",
        },
    },
    {
        "id": "scheme2_NS_ge30_group_ge3",
        "short": "NS ≥ 30 + genotype groups ≥ 3",
        "condition": (
            "MAF ≥ 0.05, NS ≥ 30, each observed genotype group ≥ 3"
        ),
        "directory": "NS_ge30_observed_genotype_group_ge3",
        "files": {
            "Duroc": (
                "Duroc_MAF005_NS_ge30_"
                "observed_genotype_group_ge3_pos.txt"
            ),
            "Luchuan": (
                "Luchuan_MAF005_NS_ge30_"
                "observed_genotype_group_ge3_pos.txt"
            ),
            "DLC": (
                "DLC_MAF005_NS_ge30_"
                "observed_genotype_group_ge3_pos.txt"
            ),
        },
    },
    {
        "id": "scheme3_NS_ge30_ALTcarrier_ge3",
        "short": "NS ≥ 30 + ALT carriers ≥ 3",
        "condition": "MAF ≥ 0.05, NS ≥ 30, ALT carriers ≥ 3",
        "directory": "NS_ge30_ALTcarrier_ge3",
        "files": {
            "Duroc": "Duroc_MAF005_NS_ge30_ALTcarrier_ge3_pos.txt",
            "Luchuan": "Luchuan_MAF005_NS_ge30_ALTcarrier_ge3_pos.txt",
            "DLC": "DLC_MAF005_NS_ge30_ALTcarrier_ge3_pos.txt",
        },
    },
    {
        "id": "scheme4_NS_ge50",
        "short": "NS ≥ 50",
        "condition": "MAF ≥ 0.05, NS ≥ 50",
        "directory": "NS_ge50_no_genotype_count_requirement",
        "files": {
            "Duroc": "Duroc_MAF005_NS_ge50_pos.txt",
            "Luchuan": "Luchuan_MAF005_NS_ge50_pos.txt",
            "DLC": "DLC_MAF005_NS_ge50_pos.txt",
        },
    },
    {
        "id": "scheme5_NS_ge50_group_ge3",
        "short": "NS ≥ 50 + genotype groups ≥ 3",
        "condition": (
            "MAF ≥ 0.05, NS ≥ 50, each observed genotype group ≥ 3"
        ),
        "directory": "NS_ge50_observed_genotype_group_ge3",
        "files": {
            "Duroc": (
                "Duroc_MAF005_NS_ge50_"
                "observed_genotype_group_ge3_pos.txt"
            ),
            "Luchuan": (
                "Luchuan_MAF005_NS_ge50_"
                "observed_genotype_group_ge3_pos.txt"
            ),
            "DLC": (
                "DLC_MAF005_NS_ge50_"
                "observed_genotype_group_ge3_pos.txt"
            ),
        },
    },
    {
        "id": "scheme6_NS_ge50_ALTcarrier_ge3",
        "short": "NS ≥ 50 + ALT carriers ≥ 3",
        "condition": "MAF ≥ 0.05, NS ≥ 50, ALT carriers ≥ 3",
        "directory": "NS_ge50_ALTcarrier_ge3",
        "files": {
            "Duroc": "Duroc_MAF005_NS_ge50_ALTcarrier_ge3_pos.txt",
            "Luchuan": "Luchuan_MAF005_NS_ge50_ALTcarrier_ge3_pos.txt",
            "DLC": "DLC_MAF005_NS_ge50_ALTcarrier_ge3_pos.txt",
        },
    },
]

# Optionally export the full position lists for the seven disjoint regions.
# The lists are disabled by default to limit disk use.
EXPORT_REGION_LISTS = False

SET_ORDER = ("Duroc", "Luchuan", "DLC")

REGION_NAMES = {
    1: "Duroc_only",
    2: "Luchuan_only",
    3: "Duroc_Luchuan_only",
    4: "DLC_only",
    5: "Duroc_DLC_only",
    6: "Luchuan_DLC_only",
    7: "All_three",
}

REGION_FILENAMES = {
    code: f"{name}_positions.txt"
    for code, name in REGION_NAMES.items()
}


# ============================================================
# 2. Stream and classify positions
# ============================================================

def read_next_site(handle):
    while True:
        line = handle.readline()
        if line == "":
            return None
        site = line.strip()
        if site:
            return site


def stream_venn_classification(input_files, output_dir):
    """Input files must contain unique positions sorted with LC_ALL=C."""

    counts = {code: 0 for code in range(1, 8)}
    region_handles = {}

    if EXPORT_REGION_LISTS:
        region_dir = output_dir / "region_position_lists"
        region_dir.mkdir(parents=True, exist_ok=True)
        for code, filename in REGION_FILENAMES.items():
            region_handles[code] = open(
                region_dir / filename, "w", encoding="utf-8"
            )

    handles = []
    try:
        for population in SET_ORDER:
            handles.append(
                open(input_files[population], "r", encoding="utf-8")
            )

        current = [read_next_site(handle) for handle in handles]
        processed = 0

        while any(site is not None for site in current):
            smallest = min(site for site in current if site is not None)
            code = 0

            for index, site in enumerate(current):
                if site == smallest:
                    code |= 1 << index

            counts[code] += 1
            processed += 1

            if EXPORT_REGION_LISTS:
                region_handles[code].write(smallest + "\n")

            for index, site in enumerate(current):
                if site == smallest:
                    current[index] = read_next_site(handles[index])

            if processed % 5_000_000 == 0:
                print(f"    Union positions processed: {processed:,}", flush=True)

    finally:
        for handle in handles:
            handle.close()
        for handle in region_handles.values():
            handle.close()

    return counts


# ============================================================
# 3. Count and validate regions
# ============================================================

def calculate_totals(counts):
    return {
        "Duroc": counts[1] + counts[3] + counts[5] + counts[7],
        "Luchuan": counts[2] + counts[3] + counts[6] + counts[7],
        "DLC": counts[4] + counts[5] + counts[6] + counts[7],
        "Union": sum(counts.values()),
    }


def count_nonempty_lines(path):
    with open(path, "r", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip())


def validate_totals(totals, input_files):
    for population in SET_ORDER:
        expected = count_nonempty_lines(input_files[population])
        observed = totals[population]
        if observed != expected:
            raise RuntimeError(
                f"{population} total failed validation: "
                f"input={expected:,}, reconstructed from Venn regions={observed:,}"
            )


def save_count_table(scheme, counts, totals, output_file):
    with open(output_file, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Scheme", "Condition", "Region", "Count"])

        for code in range(1, 8):
            writer.writerow(
                [
                    scheme["id"],
                    scheme["condition"],
                    REGION_NAMES[code],
                    counts[code],
                ]
            )

        for population in ("Duroc", "Luchuan", "DLC", "Union"):
            writer.writerow(
                [
                    scheme["id"],
                    scheme["condition"],
                    f"{population}_total",
                    totals[population],
                ]
            )


# ============================================================
# 4. Draw the Venn diagrams
# ============================================================

def draw_venn(scheme, counts, totals, output_dir):
    subsets = (
        counts[1], counts[2], counts[3], counts[4],
        counts[5], counts[6], counts[7],
    )

    fig, ax = plt.subplots(figsize=(8.6, 8.2))

    venn = draw_venn3(
        subsets=subsets,
        set_labels=(
            f"Duroc\n(n = {totals['Duroc']:,})",
            f"Luchuan\n(n = {totals['Luchuan']:,})",
            f"DLC\n(n = {totals['DLC']:,})",
        ),
        set_colors=("#F5A623", "#58A5D8", "#DE2D26"),
        alpha=0.84,
        ax=ax,
    )

    for label in venn.set_labels or []:
        if label is not None:
            label.set_fontsize(13)

    # Set labels by region ID to include all seven counts.
    region_ids = ("100", "010", "110", "001", "101", "011", "111")
    values = subsets

    # Use fixed coordinates when the layout does not supply a region label.
    # The All-three region can crowd out the Duroc–Luchuan-only region (110).
    fallback_positions = {
        "100": (-0.48, 0.18),
        "010": (0.48, 0.18),
        "110": (0.00, 0.42),
        "001": (0.00, -0.50),
        "101": (-0.27, -0.13),
        "011": (0.27, -0.13),
        "111": (0.00, 0.02),
    }

    for region_id, value in zip(region_ids, values):
        label = venn.get_label_by_id(region_id)
        if label is not None:
            label.set_text(f"{value:,}")
            label.set_fontsize(10.5)
            label.set_color("black")
        else:
            x, y = fallback_positions[region_id]
            ax.text(
                x,
                y,
                f"{value:,}",
                ha="center",
                va="center",
                fontsize=10.5,
                color="black",
                zorder=10,
            )

    ax.set_title(
        "Distribution of variant positions among Duroc, Luchuan and DLC\n"
        f"({scheme['condition']})",
        fontsize=14,
        pad=22,
    )

    if UNWEIGHTED:
        fig.text(
            0.5,
            0.025,
            "Circle areas are not proportional to set sizes.",
            ha="center",
            fontsize=9,
            color="#555555",
        )

    png_file = output_dir / "Venn_3population_variant_positions.png"
    pdf_file = output_dir / "Venn_3population_variant_positions.pdf"

    fig.savefig(
        png_file, dpi=300, bbox_inches="tight", facecolor="white"
    )
    fig.savefig(
        pdf_file, bbox_inches="tight", facecolor="white"
    )
    plt.close(fig)

    return png_file, pdf_file


# ============================================================
# 5. Main workflow
# ============================================================

def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    combined_rows = []

    schemes = SCHEMES if args.scheme == "all" else [SCHEMES[int(args.scheme) - 1]]
    for scheme in schemes:
        print("=" * 80, flush=True)
        print(f"Processing: {scheme['id']}", flush=True)
        print(f"Condition: {scheme['condition']}", flush=True)

        input_dir = BASE_DIR / scheme["directory"] / "position_sets"
        output_dir = OUTPUT_ROOT / scheme["id"]
        output_dir.mkdir(parents=True, exist_ok=True)

        input_files = {
            population: input_dir / scheme["files"][population]
            for population in SET_ORDER
        }

        for population, path in input_files.items():
            if not path.is_file():
                raise FileNotFoundError(
                    f"Input file for {population} not found: {path}"
                )

        counts = stream_venn_classification(input_files, output_dir)
        totals = calculate_totals(counts)
        validate_totals(totals, input_files)

        table_file = output_dir / "Venn_region_counts.tsv"
        save_count_table(scheme, counts, totals, table_file)
        png_file, pdf_file = draw_venn(
            scheme, counts, totals, output_dir
        )

        print("    Seven disjoint regions:", flush=True)
        for code in range(1, 8):
            print(
                f"      {REGION_NAMES[code]}: {counts[code]:,}",
                flush=True,
            )

        print(f"    Union: {totals['Union']:,}", flush=True)
        print(f"    Count table: {table_file}", flush=True)
        print(f"    PNG: {png_file}", flush=True)
        print(f"    PDF: {pdf_file}", flush=True)

        for code in range(1, 8):
            combined_rows.append(
                [
                    scheme["id"], scheme["short"], scheme["condition"],
                    REGION_NAMES[code], counts[code],
                ]
            )

        for population in ("Duroc", "Luchuan", "DLC", "Union"):
            combined_rows.append(
                [
                    scheme["id"], scheme["short"], scheme["condition"],
                    f"{population}_total", totals[population],
                ]
            )

    summary_name = (
        "Venn_six_schemes_region_counts.tsv" if args.scheme == "all"
        else f"Venn_scheme{args.scheme}_region_counts.tsv"
    )
    combined_file = OUTPUT_ROOT / summary_name
    with open(combined_file, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["Scheme", "Short_label", "Condition", "Region", "Count"]
        )
        writer.writerows(combined_rows)

    print("=" * 80, flush=True)
    print("Completed the Venn diagrams for all filtering schemes.", flush=True)
    print(f"Combined summary: {combined_file}", flush=True)
    print(f"Output directory: {OUTPUT_ROOT}", flush=True)


if __name__ == "__main__":
    main()
