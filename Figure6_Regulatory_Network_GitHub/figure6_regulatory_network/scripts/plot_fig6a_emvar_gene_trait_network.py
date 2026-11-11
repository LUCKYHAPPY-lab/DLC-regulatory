#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Plot the filtered emVar–gene–trait network.

Exclude Exterior and Adaptation. Label the smallest raw gene-based GWAS
P value in each remaining category. Break ties by the number of external
evidence types, then by gene, trait and emVar. Run input checks by default.
Use --run to generate plots and Source Data.
"""

from pathlib import Path
import argparse
import re

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, PathPatch
from matplotlib.path import Path as MplPath

import numpy as np
import pandas as pd


# ============================================================
# Plot settings
# ============================================================

EXCLUDED_CATEGORIES = {"Exterior", "Adaptation"}

CATEGORY_ORDER = [
    "Meat and carcass",
    "Production",
    "Reproduction",
    "Health",
    "Unmapped",
]

CATEGORY_COLORS = {
    "Adaptation": "#D51317",
    "Exterior": "#F39200",
    "Health": "#007B3D",
    "Meat and carcass": "#0094CD",
    "Production": "#9632B8",
    "Reproduction": "#FC97AF",
    "Unmapped": "#8C8C8C",
}

CELL_COLORS = {
    "PSC only": "#E64B35",
    "PK15 only": "#3C78D8",
    "PSC & PK15": "#F2C94C",
    "Unresolved": "#BDBDBD",
}

GENE_SOURCE_COLORS = {
    "DLC ASE": "#7B2CBF",
    "DLC eQTL": "#00A6A6",
    "DLC ASE & eQTL": "#F28E2B",
    "Unresolved": "#666666",
}

LIBRARY_MARKERS = {
    "CRE": "s",
    "3UTR": "^",
    "CRE & 3UTR": "D",
    "Unresolved": "o",
}

EXTERNAL_EVIDENCE_COLUMNS = [
    "gene_based",
    "TWAS_Muscle",
    "SMR_Muscle",
    "enloc_Muscle",
    "SNP_GWAS_P_lt_5e8",
]

REQUIRED_PATHWAY_COLUMNS = {
    "emVar",
    "gene",
    "trait",
    "ASE",
    "cis_eQTL",
    "gene_based",
    "TWAS_Muscle",
    "SMR_Muscle",
    "enloc_Muscle",
    "SNP_GWAS_P_lt_5e8",
    "gene_based_min_p",
}


# ============================================================
# Utilities
# ============================================================

def as_bool(series):
    return (
        series.astype("string")
        .str.strip()
        .str.lower()
        .isin(["true", "1", "1.0", "yes", "y"])
        .fillna(False)
    )


def clean_name(value):
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def find_column(df, candidates, required=True):
    lookup = {clean_name(c): c for c in df.columns}

    for candidate in candidates:
        key = clean_name(candidate)
        if key in lookup:
            return lookup[key]

    if required:
        raise KeyError(
            f"Cannot find {candidates}; columns={list(df.columns)}"
        )

    return None


def natural_key(value):
    return tuple(
        (0, int(x)) if x.isdigit() else (1, x.lower())
        for x in re.split(r"(\d+)", str(value))
        if x
    )


def collapse_tokens(series):
    tokens = set()

    for value in series.dropna().astype(str):
        tokens.update(
            x.strip()
            for x in value.split(";")
            if x.strip()
        )

    return ";".join(sorted(tokens))


def token_group(series, kind):
    tokens = {
        x.upper()
        for x in collapse_tokens(series).split(";")
        if x
    }

    if kind == "cell":
        a = "PSC" in tokens
        b = "PK15" in tokens

        if a and b:
            return "PSC & PK15"
        if a:
            return "PSC only"
        if b:
            return "PK15 only"
        return "Unresolved"

    a = "CRE" in tokens
    b = "3UTR" in tokens

    if a and b:
        return "CRE & 3UTR"
    if a:
        return "CRE"
    if b:
        return "3UTR"
    return "Unresolved"


def source_label(eqtl, ase):
    if eqtl and ase:
        return "DLC ASE & eQTL"
    if eqtl:
        return "DLC eQTL"
    if ase:
        return "DLC ASE"
    return "Unresolved"


def canonical_emvar(pos, alleles):
    p = (
        str(pos)
        .replace("chr", "")
        .replace("CHR", "")
        .strip()
    )

    a = (
        str(alleles)
        .replace("|", "/")
        .replace("_", "/")
        .strip()
    )

    parts = [
        x.strip().upper()
        for x in a.split("/")
        if x.strip()
    ]

    if len(parts) != 2:
        return f"{p}:{a}"

    return f"{p}:{'/'.join(sorted(parts))}"


def first_gene_name(name, gene):
    if pd.isna(name):
        return str(gene)

    tokens = [
        x.strip()
        for x in str(name).split(";")
        if x.strip() and x.strip().lower() not in {"nan", "none", "na"}
    ]

    return tokens[0] if tokens else str(gene)


def load_trait_map(path):
    x = pd.read_excel(path)

    category = find_column(
        x,
        ["Main_Category", "Main Category"],
    )

    keys = [
        find_column(x, ["Trait"], False),
        find_column(
            x,
            ["Study Symbol Raw", "study_symbol"],
            False,
        ),
    ]

    pieces = []

    for key in dict.fromkeys(k for k in keys if k):
        z = x[[key, category]].dropna().copy()
        z.columns = ["trait", "Main_Category"]

        z["trait"] = z["trait"].astype(str).str.strip()
        z["Main_Category"] = (
            z["Main_Category"].astype(str).str.strip()
        )

        pieces.append(z)

    if not pieces:
        raise ValueError(
            "Trait-category file has no usable trait key"
        )

    m = pd.concat(pieces, ignore_index=True)

    m["_key"] = (
        m["trait"].astype(str).str.strip().str.casefold()
    )

    conflicts = m.groupby("_key")["Main_Category"].nunique()

    if conflicts.gt(1).any():
        bad_keys = conflicts[conflicts.gt(1)].index.tolist()
        raise ValueError(
            f"Conflicting trait categories detected: {bad_keys[:20]}"
        )

    return m.drop_duplicates("_key")[
        ["_key", "Main_Category"]
    ]


def curved_link(
    ax,
    x0,
    y0,
    x1,
    y1,
    color,
    alpha,
    linewidth,
    zorder=1,
):
    middle = (y0 + y1) / 2

    path = MplPath(
        [
            (x0, y0),
            (x0, middle),
            (x1, middle),
            (x1, y1),
        ],
        [
            MplPath.MOVETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
        ],
    )

    ax.add_patch(
        PathPatch(
            path,
            facecolor="none",
            edgecolor=color,
            alpha=alpha,
            linewidth=linewidth,
            capstyle="round",
            zorder=zorder,
        )
    )


# ============================================================
# Load inputs
# ============================================================

def load_data(pathways, v19, categories, cre_mpra, utr_mpra):
    p = pd.read_csv(
        pathways,
        sep="\t",
        low_memory=False,
    )

    missing = sorted(REQUIRED_PATHWAY_COLUMNS - set(p.columns))

    if missing:
        raise ValueError(
            f"Strict table missing columns: {missing}"
        )

    if p.duplicated(["emVar","gene","trait"]).any():
        raise ValueError("Duplicate input variant–gene–trait paths")
    print(f"Input network: {len(p):,} paths; {p.emVar.nunique()} variants; {p.gene.nunique()} genes; {p.trait.nunique()} study–traits", flush=True)

    for c in ["ASE", "cis_eQTL"] + EXTERNAL_EVIDENCE_COLUMNS:
        p[c] = as_bool(p[c])

    p["gene_based_pvalue"] = pd.to_numeric(
        p["gene_based_min_p"],
        errors="coerce",
    )

    # Keep paths supported by gene-based GWAS for the P-value panel.
    p = p.loc[p["gene_based"]].copy()

    bad = (
        p["gene_based_pvalue"].isna()
        | ~p["gene_based_pvalue"].between(
            0,
            1,
            inclusive="right",
        )
    )

    if bad.any():
        raise ValueError(
            f"{int(bad.sum())} gene-based pathways "
            "lack valid raw P values"
        )

    p = (
        p.sort_values("gene_based_pvalue", kind="mergesort")
        .drop_duplicates(["emVar", "gene", "trait"])
        .copy()
    )

    old = pd.read_csv(
        v19,
        sep="\t",
        low_memory=False,
    )

    needed = {"gene_id_key"}

    missing = sorted(needed - set(old.columns))

    if missing:
        raise ValueError(
            f"v19 metadata missing columns: {missing}"
        )

    old = old.rename(columns={"gene_id_key":"gene"})

    # Gene names describe genes; they must not depend on an old variant-trait path.
    if "Gene.Name" in old:
        names = old.groupby("gene")["Gene.Name"].agg(collapse_tokens)
        p["Gene.Name"] = p["gene"].map(names).fillna(p["gene"])
    else:
        p["Gene.Name"] = p["gene"]

    # Cell/library support describes an MPRA variant, not a gene-trait pair.
    frames = []
    for library, path in [
        ("CRE", cre_mpra),
        ("3UTR", utr_mpra),
    ]:
        x = pd.read_csv(path, sep="\t", low_memory=False)
        required = {"chr", "variant_pos", "ref_allele", "alt_allele", "PSC.emVar_act_qc", "PK15.emVar_act_qc"}
        if not required.issubset(x):
            raise ValueError(f"MPRA columns missing in {path}: {sorted(required-set(x))}")
        # Only qualifying emVar records supply cell/library support.
        psc = as_bool(x["PSC.emVar_act_qc"])
        pk = as_bool(x["PK15.emVar_act_qc"])
        x = x.loc[psc | pk].copy()
        pos = pd.to_numeric(x.variant_pos, errors="raise")
        if pos.isna().any() or pos.mod(1).ne(0).any():
            raise ValueError(f"Invalid MPRA positions: {path}")
        x["emVar"] = [canonical_emvar(f"{c}:{int(v)}", f"{a}/{b}") for c,v,a,b in zip(x.chr,pos,x.ref_allele,x.alt_allele)]
        x["MPRA_libraries"] = library
        x["MPRA_cell_support_all"] = ["PSC;PK15" if a and b else ("PSC" if a else "PK15") for a,b in zip(as_bool(x["PSC.emVar_act_qc"]),as_bool(x["PK15.emVar_act_qc"]))]
        frames.append(x[["emVar","MPRA_libraries","MPRA_cell_support_all"]])
    meta = pd.concat(frames,ignore_index=True).groupby("emVar",as_index=False).agg(
        MPRA_libraries=("MPRA_libraries",collapse_tokens),
        MPRA_cell_support_all=("MPRA_cell_support_all",collapse_tokens))
    p=p.drop(columns=["MPRA_libraries","MPRA_cell_support_all"],errors="ignore").merge(meta,on="emVar",how="left",validate="many_to_one")
    missing_meta=p.loc[p.MPRA_libraries.isna(),"emVar"].unique()
    if len(missing_meta):
        raise ValueError(f"Network emVars absent from act_qc MPRA metadata ({len(missing_meta)}): {list(missing_meta[:20])}")

    em = p.groupby("emVar", as_index=False).agg(
        cell_group=(
            "MPRA_cell_support_all",
            lambda s: token_group(s, "cell"),
        ),
        library_group=(
            "MPRA_libraries",
            lambda s: token_group(s, "library"),
        ),
    )

    p = p.drop(
        columns=["cell_group", "library_group"],
        errors="ignore",
    )

    p = p.merge(
        em,
        on="emVar",
        how="left",
        validate="many_to_one",
    )

    m = load_trait_map(categories)

    p["_key"] = (
        p["trait"].astype(str).str.strip().str.casefold()
    )

    p = p.drop(columns=["Main_Category"], errors="ignore")

    p = (
        p.merge(
            m,
            on="_key",
            how="left",
            validate="many_to_one",
        )
        .drop(columns="_key")
    )

    p["Main_Category"] = (
        p["Main_Category"]
        .fillna("Unmapped")
        .astype(str)
        .str.strip()
    )

    p["gene_source"] = [
        source_label(e, a)
        for e, a in zip(p["cis_eQTL"], p["ASE"])
    ]

    p["n_external_evidence"] = (
        p[EXTERNAL_EVIDENCE_COLUMNS]
        .sum(axis=1)
        .astype(int)
    )

    # Exclude Exterior and Adaptation.
    excluded = p["Main_Category"].isin(EXCLUDED_CATEGORIES)

    print(
        f"Gene-based-supported pathways before filtering: {len(p):,}",
        flush=True,
    )

    print(
        f"Excluded Exterior / Adaptation: "
        f"{int(excluded.sum()):,} pathways",
        flush=True,
    )

    p = p.loc[~excluded].copy()

    if p.empty:
        raise ValueError(
            "No paths remain after excluding Exterior and Adaptation."
        )

    print(
        f"Retained pathways: {len(p):,}",
        flush=True,
    )

    return p


# ============================================================
# Select one top point per category
# ============================================================

def choose_category_top_labels(pathways, points, cats):
    """
    Select the path with the smallest gene-based GWAS P value per category.

    A gene–trait pair may have multiple emVars.
    Prefer paths with more external evidence types, then sort by ID.
    """

    ranked = pathways.sort_values(
        [
            "gene_based_pvalue",
            "n_external_evidence",
            "gene",
            "trait",
            "emVar",
        ],
        ascending=[True, False, True, True, True],
        kind="mergesort",
    )

    # Keep one representative emVar per gene–trait point.
    representatives = ranked.drop_duplicates(
        ["gene", "trait"],
        keep="first",
    )

    # Keep one top point per category.
    top = representatives.drop_duplicates(
        "Main_Category",
        keep="first",
    ).copy()

    cat_rank = {cat: i for i, cat in enumerate(cats)}

    top["_category_rank"] = top["Main_Category"].map(cat_rank)

    top = (
        top.sort_values(
            ["_category_rank", "gene_based_pvalue"],
            kind="mergesort",
        )
        .drop(columns="_category_rank")
    )

    top["label_reason"] = "top1_gene_based_GWAS_within_category"

    point_xy = points[
        [
            "gene",
            "trait",
            "trait_x",
            "plot_logp",
            "minus_log10_p",
        ]
    ]

    top = top.merge(
        point_xy,
        on=["gene", "trait"],
        how="left",
        validate="one_to_one",
    )

    top["gene_label"] = [
        first_gene_name(name, gene)
        for name, gene in zip(top["Gene.Name"], top["gene"])
    ]

    # Use "/" for sorted allele pairs, which have no REF-to-ALT direction.
    top["SNP_ID_label"] = "chr" + top["emVar"].astype(str)

    top["label_rank"] = np.arange(1, len(top) + 1)

    return top


# ============================================================
# Plot and export
# ============================================================

def plot_network(p, output, width, height, dpi):
    pathways = p.drop_duplicates(
        ["emVar", "gene", "trait"]
    ).copy()

    # Apply the same path ranking used to select labels.
    ranked = pathways.sort_values(
        [
            "gene_based_pvalue",
            "n_external_evidence",
            "gene",
            "trait",
            "emVar",
        ],
        ascending=[True, False, True, True, True],
        kind="mergesort",
    )

    points = ranked.drop_duplicates(
        ["gene", "trait"]
    ).copy()

    edges = pathways.drop_duplicates(
        ["emVar", "gene"]
    ).copy()

    points["minus_log10_p"] = -np.log10(
        points["gene_based_pvalue"].clip(
            lower=np.finfo(float).tiny
        )
    )

    present = set(points["Main_Category"])

    cats = [
        x for x in CATEGORY_ORDER
        if x in present
    ]

    cats += sorted(
        present - set(CATEGORY_ORDER),
        key=natural_key,
    )

    traits = points.groupby("trait", as_index=False).agg(
        Main_Category=("Main_Category", "first"),
        min_p=("gene_based_pvalue", "min"),
    )

    traits["_cat"] = pd.Categorical(
        traits["Main_Category"],
        categories=cats,
        ordered=True,
    )

    traits = (
        traits.sort_values(
            ["_cat", "min_p", "trait"]
        )
        .reset_index(drop=True)
    )

    traits["x"] = np.arange(len(traits), dtype=float)

    trait_x = dict(zip(traits["trait"], traits["x"]))

    points["trait_x"] = points["trait"].map(trait_x)

    gene_order = (
        points.groupby("gene")["trait_x"]
        .median()
        .sort_values()
    )

    xmax = max(len(traits) - 1, 1)

    gene_x = dict(
        zip(
            gene_order.index,
            np.linspace(0, xmax, len(gene_order)),
        )
    )

    edges["gene_x"] = edges["gene"].map(gene_x)

    em_order = (
        edges.groupby("emVar")["gene_x"]
        .median()
        .sort_values()
    )

    em_x = dict(
        zip(
            em_order.index,
            np.linspace(0, xmax, len(em_order)),
        )
    )

    raw_min = float(points["minus_log10_p"].min())
    raw_max = float(points["minus_log10_p"].max())

    axis_start = max(
        0.0,
        np.floor(raw_min * 2) / 2 - 0.5,
    )

    baseline_gap = 0.55

    points["plot_logp"] = (
        points["minus_log10_p"]
        - axis_start
        + baseline_gap
    )

    plot_max = max(
        float(points["plot_logp"].max()),
        1.0,
    )

    gene_y = -0.11 * plot_max
    em_y = -0.22 * plot_max

    labels = choose_category_top_labels(
        pathways,
        points,
        cats,
    )

    fig, ax = plt.subplots(figsize=(width, height))

    # --------------------------------------------------------
    # Category backgrounds
    # --------------------------------------------------------

    for cat in cats:
        xs = traits.loc[
            traits["Main_Category"].eq(cat),
            "x",
        ]

        if xs.empty:
            continue

        left = float(xs.min()) - 0.5
        right = float(xs.max()) + 0.5
        color = CATEGORY_COLORS.get(cat, "#8C8C8C")

        ax.axvspan(
            left,
            right,
            color=color,
            alpha=0.07,
            lw=0,
        )

        if left > -0.5:
            ax.axvline(
                left,
                color=color,
                alpha=0.32,
                lw=0.7,
            )

    # --------------------------------------------------------
    # Network edges
    # --------------------------------------------------------

    for r in edges.itertuples(index=False):
        curved_link(
            ax,
            em_x[r.emVar],
            em_y,
            gene_x[r.gene],
            gene_y,
            CELL_COLORS.get(r.cell_group, "#9E9E9E"),
            0.46,
            0.72,
        )

    for r in points.itertuples(index=False):
        curved_link(
            ax,
            gene_x[r.gene],
            gene_y,
            r.trait_x,
            0,
            CATEGORY_COLORS.get(r.Main_Category, "#707070"),
            0.38,
            0.68,
        )

    # --------------------------------------------------------
    # emVar nodes
    # --------------------------------------------------------

    em_meta = (
        edges.groupby("emVar", as_index=False)
        .agg(
            cell_group=("cell_group", "first"),
            library_group=("library_group", "first"),
        )
        .set_index("emVar")
    )

    for library, marker in LIBRARY_MARKERS.items():
        ids = [
            e for e in em_x
            if em_meta.loc[e, "library_group"] == library
        ]

        if ids:
            ax.scatter(
                [em_x[e] for e in ids],
                [em_y] * len(ids),
                s=19,
                marker=marker,
                c=[
                    CELL_COLORS.get(
                        em_meta.loc[e, "cell_group"],
                        "#BDBDBD",
                    )
                    for e in ids
                ],
                edgecolors="white",
                linewidths=0.18,
                zorder=4,
            )

    # --------------------------------------------------------
    # Gene nodes
    # --------------------------------------------------------

    gene_meta = pathways.groupby("gene", as_index=False).agg(
        ASE=("ASE", "max"),
        cis_eQTL=("cis_eQTL", "max"),
    )

    gene_meta["source"] = [
        source_label(e, a)
        for e, a in zip(
            gene_meta["cis_eQTL"],
            gene_meta["ASE"],
        )
    ]

    gene_meta = gene_meta.set_index("gene")

    genes = list(gene_x)

    ax.scatter(
        [gene_x[g] for g in genes],
        [gene_y] * len(genes),
        s=22,
        c=[
            GENE_SOURCE_COLORS[gene_meta.loc[g, "source"]]
            for g in genes
        ],
        edgecolors="white",
        linewidths=0.22,
        zorder=4,
    )

    # --------------------------------------------------------
    # Trait nodes
    # --------------------------------------------------------

    ax.scatter(
        traits["x"],
        np.zeros(len(traits)),
        s=14,
        marker="s",
        c=traits["Main_Category"].map(
            lambda c: CATEGORY_COLORS.get(c, "#8C8C8C")
        ),
        edgecolors="white",
        linewidths=0.18,
        zorder=4,
    )

    # --------------------------------------------------------
    # Draw all gene-based GWAS points in the upper panel.
    # --------------------------------------------------------

    multi = points["n_external_evidence"].ge(2)

    ax.scatter(
        points.loc[~multi, "trait_x"],
        points.loc[~multi, "plot_logp"],
        c=points.loc[~multi, "Main_Category"].map(
            lambda c: CATEGORY_COLORS.get(c, "#8C8C8C")
        ),
        s=26,
        edgecolors="none",
        alpha=0.82,
        zorder=3,
    )

    ax.scatter(
        points.loc[multi, "trait_x"],
        points.loc[multi, "plot_logp"],
        c=points.loc[multi, "Main_Category"].map(
            lambda c: CATEGORY_COLORS.get(c, "#8C8C8C")
        ),
        s=37,
        edgecolors="black",
        linewidths=0.75,
        alpha=0.94,
        zorder=4,
    )

    # --------------------------------------------------------
    # Label one top point per category.
    # Center the label above its category and connect it to the point.
    # --------------------------------------------------------

    if not labels.empty:
        category_centers = (
            traits.groupby("Main_Category")["x"]
            .mean()
            .to_dict()
        )

        for row in labels.itertuples(index=False):
            text = (
                f"{row.SNP_ID_label}\n"
                "$\\downarrow$\n"
                f"{row.gene_label}\n"
                "$\\downarrow$\n"
                f"{row.trait}"
            )

            color = CATEGORY_COLORS.get(
                row.Main_Category,
                "#555555",
            )

            label_x = category_centers[row.Main_Category]

            # Stagger neighboring labels slightly.
            label_y = plot_max * (
                1.09 + 0.10 * ((int(row.label_rank) - 1) % 2)
            )

            ax.annotate(
                text,
                xy=(row.trait_x, row.plot_logp),
                xycoords="data",
                xytext=(label_x, label_y),
                textcoords="data",
                ha="center",
                va="bottom",
                fontsize=6.5,
                fontweight="bold",
                color="#202020",
                bbox=dict(
                    boxstyle="round,pad=0.28",
                    facecolor="white",
                    edgecolor=color,
                    linewidth=0.85,
                    alpha=0.96,
                ),
                arrowprops=dict(
                    arrowstyle="->",
                    color=color,
                    lw=0.85,
                    shrinkA=2,
                    shrinkB=3,
                ),
                annotation_clip=False,
                zorder=9,
            )

            ax.scatter(
                [row.trait_x],
                [row.plot_logp],
                s=70,
                facecolors="none",
                edgecolors=color,
                linewidths=1.2,
                zorder=8,
            )

    # --------------------------------------------------------
    # Axes
    # --------------------------------------------------------

    ax.axhline(
        0,
        color="#333333",
        lw=0.8,
    )

    ax.set_xlim(
        -0.5,
        len(traits) - 0.5,
    )

    ax.set_xticks([])

    ax.set_ylabel(
        r"Gene-based GWAS $-\log_{10}(P)$",
        fontsize=11,
    )

    lo = int(np.ceil(axis_start))
    hi = int(np.ceil(raw_max))

    step = max(
        1,
        int(np.ceil((hi - lo + 1) / 6)),
    )

    ticks = np.arange(lo, hi + 1, step)
    positions = ticks - axis_start + baseline_gap

    keep = positions >= 0

    ax.set_yticks(positions[keep])

    ax.set_yticklabels(
        [str(int(x)) for x in ticks[keep]]
    )

    # Leave room above the axis for category labels.
    ax.set_ylim(
        em_y - 0.07 * plot_max,
        plot_max * 1.48,
    )

    ax.grid(
        axis="y",
        color="#DDDDDD",
        lw=0.6,
        alpha=0.65,
    )

    ax.spines[["top", "right"]].set_visible(False)

    ax.text(
        -0.015,
        0,
        "Trait",
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="center",
        fontsize=9,
        fontweight="bold",
    )

    ax.text(
        -0.015,
        gene_y,
        "DLC gene",
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="center",
        fontsize=9,
        fontweight="bold",
    )

    ax.text(
        -0.015,
        em_y,
        "MPRA frVar",
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="center",
        fontsize=9,
        fontweight="bold",
    )

    # --------------------------------------------------------
    # Legend
    # --------------------------------------------------------

    cat_handles = [
        Patch(
            facecolor=CATEGORY_COLORS.get(c, "#8C8C8C"),
            label=c,
        )
        for c in cats
    ]

    ev_handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markersize=6,
            markerfacecolor="#888",
            markeredgecolor="none",
            label="1 evidence type",
        ),
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markersize=7,
            markerfacecolor="#888",
            markeredgecolor="black",
            label="≥2 evidence types",
        ),
    ]

    gene_handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markersize=6,
            markerfacecolor=GENE_SOURCE_COLORS[x],
            markeredgecolor="white",
            label=x,
        )
        for x in [
            "DLC ASE",
            "DLC eQTL",
            "DLC ASE & eQTL",
        ]
    ]

    cell_handles = [
        Line2D(
            [],
            [],
            marker="s",
            linestyle="none",
            markersize=6,
            markerfacecolor=CELL_COLORS[x],
            markeredgecolor="white",
            label=x,
        )
        for x in [
            "PSC only",
            "PK15 only",
            "PSC & PK15",
        ]
    ]

    lib_handles = [
        Line2D(
            [],
            [],
            marker=LIBRARY_MARKERS[x],
            linestyle="none",
            markersize=6,
            markerfacecolor="#777",
            markeredgecolor="white",
            label=x,
        )
        for x in [
            "CRE",
            "3UTR",
            "CRE & 3UTR",
        ]
    ]

    legend_specs = [
        (cat_handles, "Trait category", 1.00),
        (ev_handles, "External evidence", 0.68),
        (gene_handles, "DLC regulatory source", 0.51),
        (cell_handles, "MPRA cell support", 0.32),
        (lib_handles, "MPRA library", 0.13),
    ]

    for i, (handles, title, y) in enumerate(legend_specs):
        leg = ax.legend(
            handles=handles,
            title=title,
            bbox_to_anchor=(1.01, y),
            loc="upper left",
            frameon=False,
            fontsize=7,
            title_fontsize=8,
            borderaxespad=0,
        )

        if i < len(legend_specs) - 1:
            ax.add_artist(leg)

    fig.subplots_adjust(
        left=0.075,
        right=0.83,
        bottom=0.08,
        top=0.97,
    )

    # --------------------------------------------------------
    # Save figures
    # --------------------------------------------------------

    pdf = (
        output
        / "strict_emVar_gene_trait_network_category_top1_new_ASE.pdf"
    )

    png = (
        output
        / "strict_emVar_gene_trait_network_category_top1_new_ASE.png"
    )

    fig.savefig(
        pdf,
        bbox_inches="tight",
        transparent=True,
    )

    fig.savefig(
        png,
        dpi=dpi,
        bbox_inches="tight",
        transparent=True,
    )

    plt.close(fig)

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    summary = pd.DataFrame(
        {
            "metric": [
                "plotted_pathways",
                "emVars",
                "genes",
                "traits",
                "gene_trait_points",
                "emVar_gene_edges",
                "labeled_category_top_points",
            ],
            "value": [
                len(pathways),
                pathways["emVar"].nunique(),
                pathways["gene"].nunique(),
                pathways["trait"].nunique(),
                len(points),
                len(edges),
                len(labels),
            ],
        }
    )

    category_summary = (
        pathways.groupby("Main_Category", as_index=False)
        .agg(
            pathways=("emVar", "size"),
            emVars=("emVar", "nunique"),
            genes=("gene", "nunique"),
            traits=("trait", "nunique"),
            min_gene_based_pvalue=("gene_based_pvalue", "min"),
        )
    )

    category_rank = {
        cat: i for i, cat in enumerate(cats)
    }

    category_summary["_rank"] = (
        category_summary["Main_Category"].map(category_rank)
    )

    category_summary = (
        category_summary.sort_values("_rank")
        .drop(columns="_rank")
    )

    trait_export = traits.drop(columns="_cat")

    # --------------------------------------------------------
    # Export Source Data
    # --------------------------------------------------------

    pathways.to_csv(
        output / "source_data_plotted_pathways.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    points.to_csv(
        output / "source_data_gene_trait_points.tsv",
        sep="\t",
        index=False,
    )

    edges.to_csv(
        output / "source_data_emVar_gene_edges.tsv",
        sep="\t",
        index=False,
    )

    trait_export.to_csv(
        output / "source_data_trait_nodes.tsv",
        sep="\t",
        index=False,
    )

    labels.to_csv(
        output / "source_data_category_top1_labels.tsv",
        sep="\t",
        index=False,
    )

    summary.to_csv(
        output / "network_summary_v34.tsv",
        sep="\t",
        index=False,
    )

    category_summary.to_csv(
        output / "category_summary_v34.tsv",
        sep="\t",
        index=False,
    )

    excel = output / "figure_source_data_v34.xlsx"

    blocks = []
    for kind, frame in [
        ("Pathway",pathways),("Gene_trait_point",points),("Variant_gene_edge",edges),
        ("Trait_node",trait_export),("Category_top_label",labels),
        ("Summary",summary),("Category_summary",category_summary),
    ]:
        z=frame.copy();z.insert(0,"Record_type",kind);blocks.append(z)
    source=pd.concat(blocks,ignore_index=True,sort=False)
    source.to_excel(excel,sheet_name="Source_Data",index=False,na_rep="NA")

    return summary, labels, pdf, png, excel


# ============================================================
# Entry point
# ============================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__)

    ap.add_argument(
        "--pathways",
        type=Path,
        required=True,
    )

    ap.add_argument(
        "--v19-metadata",
        type=Path,
        required=True,
    )

    ap.add_argument(
        "--categories",
        type=Path,
        required=True,
    )

    ap.add_argument(
        "--output",
        type=Path,
        required=True,
    )

    ap.add_argument(
        "--run",
        action="store_true",
    )

    ap.add_argument(
        "--force",
        action="store_true",
    )

    ap.add_argument(
        "--width",
        type=float,
        default=26,
    )

    ap.add_argument(
        "--height",
        type=float,
        default=9.5,
    )

    ap.add_argument(
        "--dpi",
        type=int,
        default=600,
    )

    ap.add_argument("--cre-mpra", type=Path, required=True)
    ap.add_argument("--utr-mpra", type=Path, required=True)
    args = ap.parse_args()

    if args.width <= 0 or args.height <= 0 or args.dpi <= 0:
        raise ValueError(
            "--width, --height and --dpi must be positive."
        )

    for path in [
        args.pathways,
        args.v19_metadata,
        args.categories,
        args.cre_mpra,
        args.utr_mpra,
    ]:
        if not path.is_file():
            raise FileNotFoundError(path)

    header = pd.read_csv(
        args.pathways,
        sep="\t",
        nrows=0,
    )

    missing = REQUIRED_PATHWAY_COLUMNS - set(header.columns)

    if missing:
        raise ValueError(
            f"Missing strict columns: {sorted(missing)}"
        )

    print(
        "PRECHECK PASSED: input files exist and strict-table "
        "required columns are present.",
        flush=True,
    )

    print(
        "Plot settings: exclude Exterior / Adaptation; "
        "label one strongest gene-based GWAS point per remaining category.",
        flush=True,
    )

    if not args.run:
        print(
            "Precheck only. Add --run to generate outputs.",
            flush=True,
        )
        return

    if args.output.exists():
        if not args.output.is_dir():
            raise NotADirectoryError(args.output)

        if any(args.output.iterdir()) and not args.force:
            raise FileExistsError(
                f"Output exists: {args.output}; "
                "use --force to replace same-named outputs"
            )

    print(
        "[1/3] Loading and filtering pathways...",
        flush=True,
    )

    data = load_data(
        args.pathways,
        args.v19_metadata,
        args.categories,
        args.cre_mpra,
        args.utr_mpra,
    )

    args.output.mkdir(
        parents=True,
        exist_ok=True,
    )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    print(
        "[2/3] Drawing network with one top label per category...",
        flush=True,
    )

    summary, labels, pdf, png, excel = plot_network(
        data,
        args.output,
        args.width,
        args.height,
        args.dpi,
    )

    print(
        "[3/3] Complete.\n"
        + summary.to_string(index=False),
        flush=True,
    )

    print(
        "\nTop labeled point in each category:",
        flush=True,
    )

    print(
        labels[
            [
                "Main_Category",
                "emVar",
                "gene_label",
                "trait",
                "gene_based_pvalue",
                "n_external_evidence",
            ]
        ].to_string(index=False),
        flush=True,
    )

    print("\nPDF:", pdf)
    print("PNG:", png)
    print("Excel:", excel)
    print("COMPLETED:", args.output)


if __name__ == "__main__":
    main()