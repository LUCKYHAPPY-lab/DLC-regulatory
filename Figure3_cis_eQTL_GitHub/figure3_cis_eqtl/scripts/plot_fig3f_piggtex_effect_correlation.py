#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Harmonize shared DLC and PigGTEx muscle cis-eQTL effect estimates.

Match chromosome, position, gene ID and the unordered allele pair.
DLC SNP IDs use CHROM:POS_ALLELE1_ALLELE2. Specify which allele
the DLC beta describes using --dlc-effect-allele-position. PigGTEx
variant IDs use CHROM_POS_REF_ALT, and slope describes ALT dosage.
Orient the DLC beta to PigGTEx ALT and retain the PigGTEx slope.
Matching uses the supplied alleles without strand complementation.

The DLC input is read by column position. Keep the original order
of the first 14 fields in eQTL_MAF005_kept.tsv.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, linregress, pearsonr, spearmanr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Match shared DLC/PigGTEx muscle cis-eQTL pairs, harmonize beta "
            "to the PigGTEx ALT allele, and calculate effect correlations."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--dlc", required=True, help="DLC cis-eQTL result")
    parser.add_argument(
        "--pig",
        required=True,
        help="PigGTEx muscle significant cis-eQTL file",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="New output directory; existing result files are not overwritten",
    )
    parser.add_argument(
        "--dlc-effect-allele-position",
        required=True,
        choices=("first", "second"),
        help=(
            "Which allele in DLC SNP_id is counted by the genotype matrix and "
            "therefore corresponds to DLC beta"
        ),
    )
    parser.add_argument(
        "--pig-chunksize",
        type=int,
        default=500_000,
        help="Rows per PigGTEx input chunk",
    )
    parser.add_argument(
        "--plot-quantile",
        type=float,
        default=1.0,
        help=(
            "Symmetric absolute quantile used only for plot limits. "
            "The default 1.0 displays all finite observations; "
            "statistics always use all finite observations."
        ),
    )
    parser.add_argument(
        "--expected-shared-pairs",
        type=int,
        default=171_641,
        help="Previously reported site+gene overlap, used only as a QC reference",
    )
    return parser.parse_args()


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def clean_chr(series: pd.Series) -> pd.Series:
    out = series.astype("string").str.strip()
    out = out.str.replace(r"^chr", "", regex=True, case=False)
    out = out.str.replace(r"\.0$", "", regex=True)
    return out.str.upper()


def clean_gene(series: pd.Series) -> pd.Series:
    return (
        series.astype("string")
        .str.strip()
        .str.replace(r"\.\d+$", "", regex=True)
    )


def canonical_allele_key(a1: pd.Series, a2: pd.Series) -> pd.Series:
    a1 = a1.astype("string").str.upper().str.strip()
    a2 = a2.astype("string").str.upper().str.strip()
    left = a1.where(a1 <= a2, a2)
    right = a2.where(a1 <= a2, a1)
    return left + "|" + right


def inspect_dlc_field_counts(path: Path) -> Tuple[int, int]:
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        first_data = next(reader)
    return len(header), len(first_data)


def load_dlc(path: Path, effect_position: str) -> Tuple[pd.DataFrame, Dict[str, int]]:
    header_nf, first_nf = inspect_dlc_field_counts(path)
    n_fields = max(header_nf, first_nf, 14)
    generic_names = [f"c{i:02d}" for i in range(n_fields)]

    # Only the first 14 positional fields are needed. This is robust to the
    # known malformed/missing header near Reference/Alternative/pCADD.
    usecols = [0, 1, 2, 3, 4, 9, 12, 13]
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        skiprows=1,
        names=generic_names,
        usecols=usecols,
        dtype="string",
        low_memory=False,
        on_bad_lines="warn",
    )
    df = df.rename(
        columns={
            "c00": "DLC_SNP_id",
            "c01": "DLC_CHROM_column",
            "c02": "DLC_POS_column",
            "c03": "Gene.ID",
            "c04": "Gene.Name",
            "c09": "DLC_beta",
            "c12": "DLC_pvalue",
            "c13": "DLC_BBFDR",
        }
    )

    qc: Dict[str, int] = {
        "DLC_header_field_count": header_nf,
        "DLC_first_data_field_count": first_nf,
        "DLC_input_rows": len(df),
    }

    parsed = df["DLC_SNP_id"].str.extract(
        r"^(?:chr)?([^:]+):([0-9]+)_([^_]+)_([^_]+)$",
        expand=True,
    )
    parsed.columns = [
        "DLC_CHROM_from_SNP_id",
        "DLC_POS_from_SNP_id",
        "DLC_allele1",
        "DLC_allele2",
    ]
    df = pd.concat([df, parsed], axis=1)

    for col in ("DLC_allele1", "DLC_allele2"):
        df[col] = df[col].astype("string").str.upper().str.strip()

    df["CHROM"] = clean_chr(df["DLC_CHROM_column"])
    df["POS"] = pd.to_numeric(df["DLC_POS_column"], errors="coerce").astype("Int64")
    df["Gene.ID"] = clean_gene(df["Gene.ID"])
    df["DLC_beta"] = pd.to_numeric(df["DLC_beta"], errors="coerce")
    df["DLC_pvalue"] = pd.to_numeric(df["DLC_pvalue"], errors="coerce")
    df["DLC_BBFDR"] = pd.to_numeric(df["DLC_BBFDR"], errors="coerce")

    parsed_chr = clean_chr(df["DLC_CHROM_from_SNP_id"])
    parsed_pos = pd.to_numeric(
        df["DLC_POS_from_SNP_id"], errors="coerce"
    ).astype("Int64")
    df["DLC_SNP_id_coordinate_matches_columns"] = (
        (parsed_chr == df["CHROM"]) & (parsed_pos == df["POS"])
    )

    qc["DLC_unparseable_SNP_id_rows"] = int(
        df[["DLC_CHROM_from_SNP_id", "DLC_POS_from_SNP_id",
            "DLC_allele1", "DLC_allele2"]]
        .isna()
        .any(axis=1)
        .sum()
    )
    qc["DLC_SNP_id_coordinate_mismatch_rows"] = int(
        (~df["DLC_SNP_id_coordinate_matches_columns"].fillna(False)).sum()
    )

    if effect_position == "first":
        df["DLC_effect_allele"] = df["DLC_allele1"]
        df["DLC_other_allele"] = df["DLC_allele2"]
    else:
        df["DLC_effect_allele"] = df["DLC_allele2"]
        df["DLC_other_allele"] = df["DLC_allele1"]

    df["allele_key"] = canonical_allele_key(
        df["DLC_allele1"], df["DLC_allele2"]
    )

    required = [
        "CHROM",
        "POS",
        "Gene.ID",
        "DLC_allele1",
        "DLC_allele2",
        "DLC_effect_allele",
        "DLC_other_allele",
        "DLC_beta",
        "allele_key",
    ]
    before = len(df)
    df = df.dropna(subset=required).copy()
    df = df[np.isfinite(df["DLC_beta"])].copy()
    df = df[
        df["DLC_SNP_id_coordinate_matches_columns"].fillna(False)
    ].copy()
    df = df[df["DLC_allele1"] != df["DLC_allele2"]].copy()
    qc["DLC_rows_removed_before_matching"] = before - len(df)

    join_cols = ["Gene.ID", "CHROM", "POS", "allele_key"]
    df = df.sort_values(
        ["DLC_BBFDR", "DLC_pvalue"],
        ascending=[True, True],
        na_position="last",
        kind="mergesort",
    )
    duplicate_n = int(df.duplicated(join_cols, keep="first").sum())
    qc["DLC_duplicate_join_keys_removed"] = duplicate_n
    df = df.drop_duplicates(join_cols, keep="first").copy()
    qc["DLC_rows_available_for_matching"] = len(df)

    keep = [
        "Gene.ID",
        "Gene.Name",
        "CHROM",
        "POS",
        "allele_key",
        "DLC_SNP_id",
        "DLC_allele1",
        "DLC_allele2",
        "DLC_effect_allele",
        "DLC_other_allele",
        "DLC_beta",
        "DLC_pvalue",
        "DLC_BBFDR",
    ]
    return df[keep].copy(), qc


def iter_pig_chunks(path: Path, chunksize: int) -> Iterable[pd.DataFrame]:
    required = [
        "phenotype_id",
        "variant_id",
        "af",
        "pval_nominal",
        "slope",
        "slope_se",
    ]
    yield from pd.read_csv(
        path,
        sep="\t",
        usecols=required,
        dtype={
            "phenotype_id": "string",
            "variant_id": "string",
            "af": "float64",
            "pval_nominal": "float64",
            "slope": "float64",
            "slope_se": "float64",
        },
        chunksize=chunksize,
        low_memory=False,
    )


def match_shared_pairs(
    dlc: pd.DataFrame,
    pig_path: Path,
    chunksize: int,
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    join_cols = ["Gene.ID", "CHROM", "POS", "allele_key"]
    matches: List[pd.DataFrame] = []

    pig_rows = 0
    pig_unparseable = 0
    chunk_count = 0

    for chunk in iter_pig_chunks(pig_path, chunksize):
        chunk_count += 1
        pig_rows += len(chunk)
        chunk["Gene.ID"] = clean_gene(chunk["phenotype_id"])

        parsed = chunk["variant_id"].str.extract(
            r"^(?:chr)?([^_]+)_([0-9]+)_([^_]+)_([^_]+)$",
            expand=True,
        )
        parsed.columns = [
            "PigGTEx_CHROM_from_variant_id",
            "PigGTEx_POS_from_variant_id",
            "PigGTEx_REF",
            "PigGTEx_ALT",
        ]
        chunk = pd.concat([chunk, parsed], axis=1)
        invalid = chunk[
            [
                "PigGTEx_CHROM_from_variant_id",
                "PigGTEx_POS_from_variant_id",
                "PigGTEx_REF",
                "PigGTEx_ALT",
            ]
        ].isna().any(axis=1)
        pig_unparseable += int(invalid.sum())
        chunk = chunk[~invalid].copy()

        chunk["CHROM"] = clean_chr(chunk["PigGTEx_CHROM_from_variant_id"])
        chunk["POS"] = pd.to_numeric(
            chunk["PigGTEx_POS_from_variant_id"], errors="coerce"
        ).astype("Int64")
        chunk["PigGTEx_REF"] = (
            chunk["PigGTEx_REF"].astype("string").str.upper().str.strip()
        )
        chunk["PigGTEx_ALT"] = (
            chunk["PigGTEx_ALT"].astype("string").str.upper().str.strip()
        )
        chunk["allele_key"] = canonical_allele_key(
            chunk["PigGTEx_REF"], chunk["PigGTEx_ALT"]
        )
        chunk = chunk.dropna(
            subset=["Gene.ID", "CHROM", "POS", "allele_key", "slope"]
        ).copy()
        chunk = chunk[np.isfinite(chunk["slope"])].copy()

        merged = chunk.merge(
            dlc,
            on=join_cols,
            how="inner",
            validate="many_to_one",
        )
        if len(merged):
            matches.append(merged)

    if not matches:
        raise RuntimeError(
            "No shared pairs remained after gene, coordinate, and allele matching."
        )

    shared = pd.concat(matches, ignore_index=True)
    shared = shared.rename(
        columns={
            "variant_id": "PigGTEx_variant_id",
            "af": "PigGTEx_ALT_AF",
            "pval_nominal": "PigGTEx_pvalue",
            "slope": "PigGTEx_slope",
            "slope_se": "PigGTEx_slope_se",
        }
    )
    shared["PigGTEx_effect_allele"] = shared["PigGTEx_ALT"]
    shared["PigGTEx_other_allele"] = shared["PigGTEx_REF"]

    same_effect = (
        (shared["DLC_effect_allele"] == shared["PigGTEx_effect_allele"])
        & (shared["DLC_other_allele"] == shared["PigGTEx_other_allele"])
    )
    opposite_effect = (
        (shared["DLC_effect_allele"] == shared["PigGTEx_other_allele"])
        & (shared["DLC_other_allele"] == shared["PigGTEx_effect_allele"])
    )

    shared["allele_alignment"] = np.select(
        [same_effect, opposite_effect],
        ["same_effect_allele", "opposite_effect_allele_flipped"],
        default="unresolved",
    )
    shared["DLC_beta_aligned_to_PigGTEx_ALT"] = np.select(
        [same_effect, opposite_effect],
        [shared["DLC_beta"], -shared["DLC_beta"]],
        default=np.nan,
    )

    unresolved = int((shared["allele_alignment"] == "unresolved").sum())
    shared = shared[shared["allele_alignment"] != "unresolved"].copy()

    dedup_cols = ["Gene.ID", "CHROM", "POS", "allele_key"]
    shared = shared.sort_values(
        ["PigGTEx_pvalue", "DLC_BBFDR", "DLC_pvalue"],
        ascending=[True, True, True],
        na_position="last",
        kind="mergesort",
    )
    duplicated = int(shared.duplicated(dedup_cols, keep="first").sum())
    shared = shared.drop_duplicates(dedup_cols, keep="first").copy()

    shared["pair_key"] = (
        shared["Gene.ID"].astype(str)
        + "|"
        + shared["CHROM"].astype(str)
        + ":"
        + shared["POS"].astype(str)
        + "_"
        + shared["PigGTEx_REF"].astype(str)
        + "_"
        + shared["PigGTEx_ALT"].astype(str)
    )
    shared["effect_direction_concordant"] = (
        np.sign(shared["DLC_beta_aligned_to_PigGTEx_ALT"])
        == np.sign(shared["PigGTEx_slope"])
    )

    qc = {
        "PigGTEx_input_rows": pig_rows,
        "PigGTEx_chunks_read": chunk_count,
        "PigGTEx_unparseable_variant_id_rows": pig_unparseable,
        "shared_rows_before_final_deduplication": len(shared) + duplicated,
        "shared_unresolved_alignment_rows_removed": unresolved,
        "shared_duplicate_pair_rows_removed": duplicated,
        "shared_harmonized_pairs": len(shared),
        "shared_same_effect_allele_pairs": int(
            (shared["allele_alignment"] == "same_effect_allele").sum()
        ),
        "shared_flipped_effect_allele_pairs": int(
            (
                shared["allele_alignment"]
                == "opposite_effect_allele_flipped"
            ).sum()
        ),
        "shared_unique_sites": int(
            shared[["CHROM", "POS", "allele_key"]].drop_duplicates().shape[0]
        ),
        "shared_unique_genes": int(shared["Gene.ID"].nunique()),
    }

    output_cols = [
        "pair_key",
        "Gene.ID",
        "Gene.Name",
        "CHROM",
        "POS",
        "PigGTEx_REF",
        "PigGTEx_ALT",
        "DLC_SNP_id",
        "DLC_allele1",
        "DLC_allele2",
        "DLC_effect_allele",
        "DLC_other_allele",
        "DLC_beta",
        "DLC_beta_aligned_to_PigGTEx_ALT",
        "DLC_pvalue",
        "DLC_BBFDR",
        "PigGTEx_variant_id",
        "PigGTEx_effect_allele",
        "PigGTEx_other_allele",
        "PigGTEx_slope",
        "PigGTEx_slope_se",
        "PigGTEx_pvalue",
        "PigGTEx_ALT_AF",
        "allele_alignment",
        "effect_direction_concordant",
    ]
    return shared[output_cols].copy(), qc


def calculate_stats(df: pd.DataFrame, analysis: str) -> Dict[str, object]:
    x = pd.to_numeric(
        df["DLC_beta_aligned_to_PigGTEx_ALT"], errors="coerce"
    ).to_numpy(dtype=float)
    y = pd.to_numeric(df["PigGTEx_slope"], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]

    if len(x) < 3:
        raise RuntimeError(f"{analysis}: fewer than 3 finite pairs.")

    pearson_r, pearson_p = pearsonr(x, y)
    spearman_rho, spearman_p = spearmanr(x, y)
    regression = linregress(x, y)

    nonzero = (np.sign(x) != 0) & (np.sign(y) != 0)
    concordant = int((np.sign(x[nonzero]) == np.sign(y[nonzero])).sum())
    sign_n = int(nonzero.sum())
    sign_rate = concordant / sign_n if sign_n else math.nan
    sign_p = (
        binomtest(concordant, sign_n, p=0.5, alternative="greater").pvalue
        if sign_n
        else math.nan
    )

    return {
        "analysis": analysis,
        "N_pairs": len(x),
        "N_unique_sites": int(
            df.loc[valid, ["CHROM", "POS"]].drop_duplicates().shape[0]
        ),
        "N_unique_genes": int(df.loc[valid, "Gene.ID"].nunique()),
        "Pearson_r": pearson_r,
        "Pearson_P": pearson_p,
        "Spearman_rho": spearman_rho,
        "Spearman_P": spearman_p,
        "sign_concordant_pairs": concordant,
        "sign_test_pairs": sign_n,
        "sign_concordance_rate": sign_rate,
        "sign_concordance_binomial_P": sign_p,
        "regression_slope": regression.slope,
        "regression_intercept": regression.intercept,
        "regression_R_squared": regression.rvalue**2,
        "regression_P": regression.pvalue,
    }


def p_text(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p == 0 or p < 1e-300:
        return "<1×10⁻³⁰⁰"
    return f"{p:.2e}"


def plot_correlation(
    df: pd.DataFrame,
    stats: Dict[str, object],
    out_prefix: Path,
    plot_quantile: float,
) -> None:
    x = pd.to_numeric(
        df["DLC_beta_aligned_to_PigGTEx_ALT"], errors="coerce"
    ).to_numpy(dtype=float)
    y = pd.to_numeric(df["PigGTEx_slope"], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]

    # Plot limits only affect visualization; correlation statistics above
    # always use all finite SNP-gene pairs.
    #
    # plot_quantile = 1.0:
    #   plot every finite pair and add a small margin so edge points are not
    #   visually clipped by the axes.
    #
    # plot_quantile < 1.0:
    #   retain the original quantile-based display behavior.
    if plot_quantile >= 1.0:
        x_plot = x
        y_plot = y

        x_max = float(np.max(np.abs(x_plot)))
        y_max = float(np.max(np.abs(y_plot)))

        if not np.isfinite(x_max) or x_max <= 0:
            x_max = 1.0
        if not np.isfinite(y_max) or y_max <= 0:
            y_max = 1.0

        # Three-percent padding prevents hexagons at the extrema from being
        # cut by the plotting frame.
        xlim = x_max * 1.03
        ylim = y_max * 1.03

    else:
        xlim_raw = float(np.quantile(np.abs(x), plot_quantile))
        ylim_raw = float(np.quantile(np.abs(y), plot_quantile))

        if not np.isfinite(xlim_raw) or xlim_raw <= 0:
            xlim_raw = float(np.max(np.abs(x)))
        if not np.isfinite(ylim_raw) or ylim_raw <= 0:
            ylim_raw = float(np.max(np.abs(y)))

        in_view = (
            (np.abs(x) <= xlim_raw)
            & (np.abs(y) <= ylim_raw)
        )
        x_plot = x[in_view]
        y_plot = y[in_view]

        xlim = xlim_raw * 1.03
        ylim = ylim_raw * 1.03

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 1.1,
        }
    )
    fig, ax = plt.subplots(figsize=(7, 7))
    hb = ax.hexbin(
        x_plot,
        y_plot,
        gridsize=75,
        bins="log",
        mincnt=1,
        cmap="YlGnBu_r",
        linewidths=0,
        rasterized=True,
    )
    colorbar = fig.colorbar(hb, ax=ax, fraction=0.046, pad=0.03)

    # Colors represent the number of SNP-gene pairs per hexagonal bin.
    # bins="log" applies a logarithmic color scale while the ticks remain counts.
    colorbar.set_label("")
    colorbar.ax.set_title("Pair count", fontsize=10, pad=5)
    colorbar.ax.tick_params(labelsize=9)

    reg_x = np.linspace(-xlim, xlim, 200)
    reg_y = (
        float(stats["regression_intercept"])
        + float(stats["regression_slope"]) * reg_x
    )
    ax.plot(reg_x, reg_y, color="#D7301F", lw=1.8, zorder=4)
    ax.axhline(0, color="#636363", lw=0.8, ls="--", alpha=0.65)
    ax.axvline(0, color="#636363", lw=0.8, ls="--", alpha=0.65)

    annotation = (
        f"N = {int(stats['N_pairs']):,}\n"
        f"Pearson r = {float(stats['Pearson_r']):.2f}\n"
        f"P = {p_text(float(stats['Pearson_P']))}\n"
        f"Spearman ρ = {float(stats['Spearman_rho']):.2f}\n"
        f"Sign concordance = "
        f"{float(stats['sign_concordance_rate']) * 100:.1f}%"
    )
    ax.text(
        0.04,
        0.96,
        annotation,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11,
    )

    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(-ylim, ylim)
    ax.set_xlabel("DLC cis-eQTL effect aligned to PigGTEx ALT", fontsize=12)
    ax.set_ylabel("PigGTEx muscle cis-eQTL slope (ALT effect)", fontsize=12)
    ax.tick_params(axis="both", labelsize=10, direction="out")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".png"), dpi=400, bbox_inches="tight")
    plt.close(fig)


def write_qc(
    path: Path,
    dlc_qc: Dict[str, int],
    match_qc: Dict[str, int],
    effect_position: str,
    expected_shared_pairs: int,
) -> None:
    records = []
    for section, values in (("DLC_loading", dlc_qc), ("matching", match_qc)):
        for metric, value in values.items():
            records.append(
                {"section": section, "metric": metric, "value": value}
            )
    observed = match_qc["shared_harmonized_pairs"]
    records.extend(
        [
            {
                "section": "settings",
                "metric": "DLC_effect_allele_position_in_SNP_id",
                "value": effect_position,
            },
            {
                "section": "reference",
                "metric": "previously_reported_shared_site_gene_pairs",
                "value": expected_shared_pairs,
            },
            {
                "section": "reference",
                "metric": "observed_minus_previously_reported_pairs",
                "value": observed - expected_shared_pairs,
            },
        ]
    )
    pd.DataFrame(records).to_csv(path, sep="\t", index=False)


def ensure_new_outputs(outdir: Path, output_names: Iterable[str]) -> None:
    existing = [outdir / name for name in output_names if (outdir / name).exists()]
    if existing:
        lines = "\n".join(str(path) for path in existing)
        raise FileExistsError(
            "Refusing to overwrite existing result files. "
            "Choose a new --outdir.\n" + lines
        )


def main() -> None:
    args = parse_args()
    dlc_path = Path(args.dlc).expanduser().resolve()
    pig_path = Path(args.pig).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()

    for label, path in (("DLC", dlc_path), ("PigGTEx", pig_path)):
        if not path.is_file():
            raise FileNotFoundError(f"{label} input not found: {path}")
    if not 0.90 <= args.plot_quantile <= 1.0:
        raise ValueError("--plot-quantile must be between 0.90 and 1.0.")

    outdir.mkdir(parents=True, exist_ok=True)
    names = [
        "DLC_PigGTEx_shared_pairs_beta_harmonized.tsv.gz",
        "DLC_PigGTEx_one_lead_pair_per_gene.tsv.gz",
        "DLC_PigGTEx_beta_correlation_statistics.tsv",
        "DLC_PigGTEx_beta_harmonization_QC.tsv",
        "DLC_PigGTEx_run_parameters.json",
        "DLC_PigGTEx_beta_correlation_all_pairs.pdf",
        "DLC_PigGTEx_beta_correlation_all_pairs.png",
        "DLC_PigGTEx_beta_correlation_one_pair_per_gene.pdf",
        "DLC_PigGTEx_beta_correlation_one_pair_per_gene.png",
    ]
    ensure_new_outputs(outdir, names)

    print("[1/6] Loading DLC cis-eQTL results...", flush=True)
    dlc, dlc_qc = load_dlc(
        dlc_path, args.dlc_effect_allele_position
    )
    print(
        f"      DLC rows available for matching: {len(dlc):,}",
        flush=True,
    )

    print("[2/6] Matching PigGTEx shared SNP-gene pairs...", flush=True)
    shared, match_qc = match_shared_pairs(
        dlc, pig_path, args.pig_chunksize
    )
    print(
        f"      Harmonized shared pairs: {len(shared):,}",
        flush=True,
    )

    shared_path = outdir / "DLC_PigGTEx_shared_pairs_beta_harmonized.tsv.gz"
    shared.to_csv(
        shared_path,
        sep="\t",
        index=False,
        compression="gzip",
    )

    print("[3/6] Selecting one lead shared pair per gene...", flush=True)
    one_per_gene = (
        shared.sort_values(
            ["DLC_pvalue", "DLC_BBFDR", "PigGTEx_pvalue"],
            ascending=[True, True, True],
            na_position="last",
            kind="mergesort",
        )
        .drop_duplicates("Gene.ID", keep="first")
        .copy()
    )
    one_per_gene.to_csv(
        outdir / "DLC_PigGTEx_one_lead_pair_per_gene.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    print("[4/6] Calculating correlation statistics...", flush=True)
    all_stats = calculate_stats(shared, "all_harmonized_shared_pairs")
    lead_stats = calculate_stats(one_per_gene, "one_lead_shared_pair_per_gene")
    stats_df = pd.DataFrame([all_stats, lead_stats])
    stats_df.to_csv(
        outdir / "DLC_PigGTEx_beta_correlation_statistics.tsv",
        sep="\t",
        index=False,
    )

    print("[5/6] Drawing PDF and PNG figures...", flush=True)
    plot_correlation(
        shared,
        all_stats,
        outdir / "DLC_PigGTEx_beta_correlation_all_pairs",
        args.plot_quantile,
    )
    plot_correlation(
        one_per_gene,
        lead_stats,
        outdir / "DLC_PigGTEx_beta_correlation_one_pair_per_gene",
        args.plot_quantile,
    )

    print("[6/6] Writing QC report...", flush=True)
    write_qc(
        outdir / "DLC_PigGTEx_beta_harmonization_QC.tsv",
        dlc_qc,
        match_qc,
        args.dlc_effect_allele_position,
        args.expected_shared_pairs,
    )
    parameters = {
        "DLC_input": str(dlc_path),
        "PigGTEx_input": str(pig_path),
        "output_directory": str(outdir),
        "DLC_effect_allele_position_in_SNP_id": (
            args.dlc_effect_allele_position
        ),
        "PigGTEx_effect_allele_definition": "ALT in variant_id",
        "matching_key": "Gene.ID + CHROM + POS + unordered allele set",
        "aligned_effect_allele": "PigGTEx ALT",
        "pig_chunksize": args.pig_chunksize,
        "plot_quantile": args.plot_quantile,
        "expected_shared_pairs_reference_only": args.expected_shared_pairs,
    }
    with (outdir / "DLC_PigGTEx_run_parameters.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(parameters, handle, indent=2, ensure_ascii=False)

    observed = match_qc["shared_harmonized_pairs"]
    expected = args.expected_shared_pairs
    print("\nAnalysis completed.", flush=True)
    print(f"Output directory: {outdir}", flush=True)
    print(
        f"Shared harmonized pairs: {observed:,} "
        f"(previous reference count: {expected:,})",
        flush=True,
    )
    print(
        f"Pearson r (all pairs): {all_stats['Pearson_r']:.4f}",
        flush=True,
    )
    print(
        f"Pearson r (one pair/gene): {lead_stats['Pearson_r']:.4f}",
        flush=True,
    )
    if abs(observed - expected) > max(100, int(expected * 0.01)):
        print(
            "WARNING: Harmonized pair count differs from the previous "
            "site+gene overlap by >1%. Inspect the QC table before "
            "interpreting the correlation.",
            file=sys.stderr,
            flush=True,
        )


if __name__ == "__main__":
    main()
