#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import importlib.util
import sys
import json



# ============================================================
# 1. Input and output paths
# ============================================================

parser = argparse.ArgumentParser(description="Generate Figure 4f: novel ase maf heatmap.")
parser.add_argument("--function-file", required=True, type=Path)
parser.add_argument("--ase-status", required=True, type=Path)
parser.add_argument("--outdir", required=True, type=Path)
parser.add_argument("--du-list", required=True, type=Path)
parser.add_argument("--du-frequency", required=True, type=Path)
args = parser.parse_args()
FUNCTION_FILE = args.function_file
ASE_STATUS_FILE = args.ase_status
OUT = args.outdir
DU_LIST = args.du_list
DU_FREQ = args.du_frequency

MAX_PLOT_ROWS = 10000


# ============================================================
# 2. Validate and load the original heatmap helper
# ============================================================

if FUNCTION_FILE.resolve() == Path(__file__).resolve():
    raise ValueError(
        "The helper path points to this script, which would cause a recursive import.\n"
        "Use the original helper that defines get_keys, scan_du and draw."
    )

if not FUNCTION_FILE.is_file():
    raise FileNotFoundError(
        f"Heatmap helper not found:\n{FUNCTION_FILE}\n"
        "Supply the original plot_all_DLC_MAF005.py with --function-file."
    )

spec = importlib.util.spec_from_file_location(
    "maf_functions",
    FUNCTION_FILE
)

if spec is None or spec.loader is None:
    raise ImportError(f"Cannot load the heatmap helper: {FUNCTION_FILE}")

mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)

required_attributes = [
    "pd", "np",
    "get_keys", "scan_du", "draw",
    "LABELS", "DU_LIST", "DU_FREQ", "log"
]

for name in required_attributes:
    if not hasattr(mod, name):
        raise AttributeError(
            f"The heatmap helper is missing: {name}\n"
            f"Helper file: {FUNCTION_FILE}"
        )

pd = mod.pd
np = mod.np

if len(mod.LABELS) != 2:
    raise ValueError(
        f"The helper's LABELS must contain two categories. Found: {mod.LABELS}"
    )

OUT.mkdir(parents=True, exist_ok=True)
mod.OUT = OUT
mod.DU_LIST = str(DU_LIST) if isinstance(mod.DU_LIST, str) else DU_LIST
mod.DU_FREQ = str(DU_FREQ) if isinstance(mod.DU_FREQ, str) else DU_FREQ


# ============================================================
# 3. Load ASE sites unreported in PigGTEx
# ============================================================

def load_new_ase():
    df = pd.read_csv(
        ASE_STATUS_FILE,
        sep="\t",
        low_memory=False
    )

    df.columns = df.columns.astype(str).str.strip()

    status_col = "Reported_in_PigGTEx_nonCross_muscle"

    required = {
        "SNP_id",
        status_col,
        "DLC_site_max_MAF",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing columns: {sorted(missing)}\n"
            f"Available columns: {list(df.columns)}"
        )

    if df["SNP_id"].isna().any():
        raise ValueError("The input contains missing SNP_id values.")

    df["SNP_id"] = df["SNP_id"].astype(str).str.strip()

    if df["SNP_id"].eq("").any():
        raise ValueError("The input contains empty SNP_id values.")

    if df["SNP_id"].duplicated().any():
        duplicated = df.loc[
            df["SNP_id"].duplicated(keep=False)
        ]

        duplicated.to_csv(
            OUT / "ASE_duplicate_SNP_ids.tsv",
            sep="\t",
            index=False
        )

        raise ValueError(
            "Expected one row per SNP. Duplicate SNP_id values were found. "
            "See ASE_duplicate_SNP_ids.tsv."
        )

    # Require an explicit unreported status. Flag unknown or missing values.
    status = (
        df[status_col]
        .astype("string")
        .str.strip()
        .str.lower()
    )

    known = status.isin(["true", "false", "1", "0"])

    if not known.all():
        df.loc[~known].to_csv(
            OUT / "ASE_invalid_PigGTEx_status.tsv",
            sep="\t",
            index=False
        )

        raise ValueError(
            "Missing or unrecognized PigGTEx status values found. "
            "See ASE_invalid_PigGTEx_status.tsv."
        )

    reported = status.isin(["true", "1"])

    input_total = len(df)
    input_reported = int(reported.sum())

    new = df.loc[~reported].copy()

    if new.empty:
        raise ValueError(
            "The input contains no ASE sites unreported in PigGTEx."
        )

    if "allele_key" in new.columns:
        new["PigGTEx_comparison_allele_key"] = new["allele_key"]

    new["variant_key"] = new["SNP_id"]

    # Convert SNP IDs to the position keys expected by the helper.
    new["site_key"] = mod.get_keys(
        new[["SNP_id"]].rename(
            columns={"SNP_id": "variant_id"}
        )
    )

    if new["site_key"].isna().any():
        new.loc[new["site_key"].isna()].to_csv(
            OUT / "ASE_invalid_site_keys.tsv",
            sep="\t",
            index=False
        )
        raise ValueError("Could not generate position keys for some SNPs. Affected records were exported.")

    new["DLC_MAF"] = pd.to_numeric(
        new["DLC_site_max_MAF"],
        errors="coerce"
    )

    invalid_maf = (
        new["DLC_MAF"].isna()
        | ~np.isfinite(new["DLC_MAF"])
        | new["DLC_MAF"].lt(0.05)
        | new["DLC_MAF"].gt(0.5)
    )

    if invalid_maf.any():
        new.loc[invalid_maf].to_csv(
            OUT / "ASE_invalid_DLC_MAF.tsv",
            sep="\t",
            index=False
        )

        raise ValueError(
            "The unreported subset contains invalid DLC MAF values or values below the threshold."
        )

    # Count heatmap sites by chromosome and position.
    duplicate_position = new["site_key"].duplicated(
        keep=False
    )

    if duplicate_position.any():
        new.loc[duplicate_position].to_csv(
            OUT / "ASE_multiple_variants_same_position.tsv",
            sep="\t",
            index=False
        )

        raise ValueError(
            "Multiple unreported variants share a position. Affected records were exported. "
            "Verify the position-level aggregation before proceeding."
        )

    new.to_csv(
        OUT / "00_new_ASE_input_sites.tsv",
        sep="\t",
        index=False,
        na_rep="NA"
    )

    mod.log(f"Final ASE sites: {input_total:,}")
    mod.log(f"Reported in PigGTEx: {input_reported:,}")
    mod.log(f"Unreported in PigGTEx: {len(new):,}")

    return new, input_total, input_reported


# ============================================================
# 4. Analyze the ASE subset
# ============================================================

def main():
    for path in [
        ASE_STATUS_FILE,
        Path(mod.DU_LIST),
        Path(mod.DU_FREQ),
    ]:
        if not path.is_file():
            raise FileNotFoundError(
                f"Input file not found: {path}"
            )

    mod.log(f"Heatmap helper: {FUNCTION_FILE}")
    mod.log(f"Output directory: {OUT}")

    df, input_total, input_reported = load_new_ase()

    targets = set(df["site_key"])

    mod.log(
        "\nExtracting Duroc frequencies for ASE sites unreported in PigGTEx..."
    )

    # Reuse the original Duroc frequency scan.
    # common contains positions classified as common in Duroc by the helper.
    # maxima contains the maximum MAF across Duroc allele records at each position.
    common, maxima = mod.scan_du(targets)

    df["DU_MAF"] = pd.to_numeric(
        df["site_key"].map(maxima),
        errors="coerce"
    )

    df["DU_common_MAF005"] = (
        df["site_key"].isin(common)
    )

    invalid_du = (
        df["DU_MAF"].notna()
        & (
            ~np.isfinite(df["DU_MAF"])
            | df["DU_MAF"].lt(0)
            | df["DU_MAF"].gt(0.5)
        )
    )

    if invalid_du.any():
        df.loc[invalid_du].to_csv(
            OUT / "ASE_invalid_DU_MAF.tsv",
            sep="\t",
            index=False
        )
        raise ValueError("Invalid Duroc frequency values found. Affected records were exported.")

    df["Pattern"] = np.where(
        df["DU_common_MAF005"],
        mod.LABELS[0],
        mod.LABELS[1]
    )

    df["DU_frequency_status"] = np.select(
        [
            df["DU_common_MAF005"]
            & df["DU_MAF"].notna(),

            df["DU_common_MAF005"],

            df["DU_MAF"].notna(),
        ],
        [
            "DU common; frequency available",
            "DU common-list support; frequency unavailable",
            "DU low; frequency available",
        ],
        default=(
            "DU frequency unavailable; "
            "no common-list support"
        )
    )

    df["PigGTEx_status"] = "DLC-new"

    counts = (
        df["Pattern"]
        .value_counts()
        .reindex(mod.LABELS, fill_value=0)
    )

    summary = pd.DataFrame({
        "Analysis": "ASE",
        "Subset": "PigGTEx-unreported",
        "Category": list(mod.LABELS),
        "Count": counts.to_numpy(),
        "Total_positions": len(df),
        "Percentage": counts.to_numpy() / len(df) * 100,
    })

    # --------------------------------------------------------
    # 4.1 Save the full site table and category-specific results
    # --------------------------------------------------------

    df.to_csv(
        OUT / "ASE_new_annotated_sites.tsv",
        sep="\t",
        index=False,
        na_rep="NA"
    )

    summary.to_csv(
        OUT / "ASE_new_pattern_summary.tsv",
        sep="\t",
        index=False
    )

    for i, label in enumerate(mod.LABELS, 1):
        df.loc[df["Pattern"].eq(label)].to_csv(
            OUT / f"ASE_new_category_{i}_full_sites.tsv",
            sep="\t",
            index=False,
            na_rep="NA"
        )

    status_counts = (
        df.groupby("DU_frequency_status")
        .size()
        .rename("Count")
        .reset_index()
    )

    status_counts.to_csv(
        OUT / "ASE_new_DU_frequency_status_counts.tsv",
        sep="\t",
        index=False
    )

    # --------------------------------------------------------
    # 4.2 Draw the ASE heatmap
    # --------------------------------------------------------

    mod.draw(
        df,
        summary,
        "PigGTEx-unreported ASE",
        OUT,
        MAX_PLOT_ROWS
    )

    # --------------------------------------------------------
    # 4.3 Export single-sheet Source Data
    # --------------------------------------------------------

    source = pd.concat([
        summary.assign(
            Record_type="Pattern_summary"
        ),
        status_counts.assign(
            Record_type="DU_frequency_status"
        ),
    ], ignore_index=True, sort=False)

    source = source[
        ["Record_type"]
        + [
            column for column in source.columns
            if column != "Record_type"
        ]
    ]

    with pd.ExcelWriter(
        OUT / "ASE_new_SourceData.xlsx",
        engine="openpyxl"
    ) as writer:
        source.to_excel(
            writer,
            sheet_name="Source Data",
            index=False,
            na_rep="NA"
        )

    # --------------------------------------------------------
    # 4.4 Record analysis settings
    # --------------------------------------------------------

    information = {
        "Analysis": "ASE only",
        "ASE_status_input": str(ASE_STATUS_FILE),
        "Function_script": str(FUNCTION_FILE),
        "Final_ASE_total": input_total,
        "PigGTEx_reported_ASE": input_reported,
        "PigGTEx_unreported_ASE": len(df),
        "PigGTEx_subset_definition": (
            "Not reported in the filtered "
            "non-Cross muscle ASE catalog"
        ),
        "DU_list": str(mod.DU_LIST),
        "DU_frequency": str(mod.DU_FREQ),
        "Analysis_unit": "CHROM:POS",
        "DLC_MAF_column": "DLC_site_max_MAF",
        "DU_frequency_rule": (
            "Reuse original scan_du(); maximum MAF "
            "across allele records at each position"
        ),
        "Additional_NS_filter": False,
        "Additional_QC_filter": False,
        "Maximum_plot_rows": MAX_PLOT_ROWS,
        "Gray_cells": "Missing DU numerical frequency",
        "eQTL_reanalysis": False,
    }

    (OUT / "analysis_information.json").write_text(
        json.dumps(
            information,
            ensure_ascii=False,
            indent=2
        ),
        encoding="utf-8"
    )

    mod.log("\n" + summary.to_string(index=False))
    mod.log("\n" + status_counts.to_string(index=False))
    mod.log(
        f"\nASE analysis complete. Results directory: {OUT}"
    )


if __name__ == "__main__":
    main()