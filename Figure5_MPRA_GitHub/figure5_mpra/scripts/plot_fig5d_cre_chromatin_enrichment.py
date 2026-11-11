#!/usr/bin/env python3
"""
Plot chromatin-feature enrichment in active versus inactive oligonucleotides.
Two-sided Fisher's exact test to detect both positive and negative enrichment.

Use activity BED intersections generated after cell-specific filtering.

BED file format:
- Column1-3: chr, start, end
- Column4: annotation (Active/Stable/Repressed)
- Column5: SEQ ID
- Column6-8: histone modification region
- Column9: feature type (ATAC, CTCF, H3K27ac, H3K4me3, H3K27me3, H3K4me1)
"""

import argparse
import os

import pandas as pd
import scipy.stats as stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from statsmodels.stats.multitest import multipletests

# =========================
# 1. Input and output paths
parser = argparse.ArgumentParser(description='Plot CRE chromatin-feature enrichment.')
parser.add_argument("--pk15-bed", required=True)
parser.add_argument("--psc-bed", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()
input_configs = [
    {"label": "PK15", "bed": args.pk15_bed},
    {"label": "PSC", "bed": args.psc_bed},
]
OUT_DIR = args.outdir
os.makedirs(OUT_DIR, exist_ok=True)

pdf_output = f"{OUT_DIR}/Active_vs_Stable_Histone_Enrichment_TwoSided.pdf"
xlsx_output = f"{OUT_DIR}/Active_vs_Stable_Histone_Enrichment_TwoSided_Statistics.xlsx"

# Export a Source Data workbook and a TSV of bar values.
source_xlsx = f"{OUT_DIR}/SourceData_Active_vs_Stable_Histone_Enrichment_TwoSided.xlsx"
source_tsv = f"{OUT_DIR}/SourceData_Active_vs_Stable_Histone_Enrichment_TwoSided.tsv"


# =========================
# 2. Combine activity groups
# =========================

def merge_groups(x):
    s = str(x).strip()
    if s == "Active":
        return "Active"
    else:
        return "Stable"


# =========================
# 2b. Compute small P values in log space
#     Double-precision Fisher P values can underflow to zero.
#     The CTCF tests can yield P values below the floating-point range.
# =========================

def log10_fisher_two_sided(a, b, c, d):
    """Compute log10 of the two-sided Fisher exact P value in log space."""
    from math import lgamma
    r1, c1 = a + b, a + c
    n = a + b + c + d

    def logpmf(x):
        return (lgamma(c1 + 1) - lgamma(x + 1) - lgamma(c1 - x + 1)
                + lgamma(n - c1 + 1) - lgamma(r1 - x + 1)
                - lgamma((n - c1) - (r1 - x) + 1)
                - (lgamma(n + 1) - lgamma(r1 + 1) - lgamma(n - r1 + 1)))

    lo, hi = max(0, r1 - (n - c1)), min(r1, c1)
    lp_obs = logpmf(a)
    logs = [logpmf(x) for x in range(lo, hi + 1) if logpmf(x) <= lp_obs + 1e-7]
    mx = max(logs)
    return (mx + np.log(sum(np.exp(t - mx) for t in logs))) / np.log(10)


def fmt10(log10p, sig=7):
    """Format a log10 P value in scientific notation without underflow."""
    if log10p is None or not np.isfinite(log10p):
        return ""
    e = int(np.floor(log10p))
    mant = 10 ** (log10p - e)
    if mant >= 9.9999995:          # Renormalize when rounding brings the mantissa to 10.
        mant /= 10
        e += 1
    return f"{mant:.{sig - 1}f}e{e}"


# =========================
# 3. Main analysis
# =========================

def run_analysis():

    all_plot_data = []
    all_stats_data = []
    all_seq_level = []      # One row per sequence with group and overlap flags.

    for config in input_configs:

        label = config["label"]
        bed_path = config["bed"]

        print(f"\n--- Processing {label} ---")

        # =========================
        # A. Read the BED intersection
        # =========================

        df = pd.read_csv(bed_path, sep='\t', header=None)
        df.columns = ['chr', 'start', 'end', 'annotation', 'seq_id',
                       'feat_chr', 'feat_start', 'feat_end', 'feature']

        print(f"Total rows: {len(df)}")
        print(f"Annotations: {df['annotation'].value_counts().to_dict()}")

        # =========================
        # B. Build the background set
        # =========================

        # Activity annotation for each sequence ID.
        seq_annotations = df[['seq_id', 'annotation']].drop_duplicates()
        seq_annotations['group_combined'] = seq_annotations['annotation'].apply(merge_groups)

        denoms = seq_annotations['group_combined'].value_counts().to_dict()
        n_active = denoms.get('Active', 0)
        n_stable = denoms.get('Stable', 0)

        print(f"Active total: {n_active}")
        print(f"Stable total: {n_stable}")

        # =========================
        # C. Deduplicate sequence–feature pairs
        # =========================

        df_hits = df[['seq_id', 'annotation', 'feature']].drop_duplicates()
        df_hits['group_combined'] = df_hits['annotation'].apply(merge_groups)

        features = df_hits['feature'].unique()
        print(f"Features: {features}")

        # =========================
        # C2. Source Data: binary sequence-by-feature overlap matrix
        #     Use the same overlap records as the enrichment test.
        # =========================

        overlap = (df_hits.assign(hit=1)
                          .pivot_table(index='seq_id', columns='feature',
                                       values='hit', aggfunc='max', fill_value=0))
        seq_tab = (seq_annotations.set_index('seq_id')
                                  .join(overlap.reindex(seq_annotations['seq_id']),
                                        how='left')
                                  .fillna(0))
        seq_tab[list(features)] = seq_tab[list(features)].astype(int)
        seq_tab.insert(0, 'CellType', label)
        all_seq_level.append(seq_tab.reset_index())

        # =========================
        # D. Calculate enrichment
        # =========================

        for feat in features:

            feat_df = df_hits[df_hits['feature'] == feat]

            hit_counts = feat_df['group_combined'].value_counts().to_dict()

            a = hit_counts.get('Active', 0)
            b = n_active - a
            c = hit_counts.get('Stable', 0)
            d = n_stable - c

            # enrichment ratio
            epsilon = 1e-6
            rate_active = (a + epsilon) / (n_active + epsilon)
            rate_stable = (c + epsilon) / (n_stable + epsilon)
            fold = rate_active / rate_stable
            log2fold = np.log2(fold)

            # Fisher exact test - two-sided to detect both enrichment and depletion
            # Compute P values in double precision and in log space.
            if (n_active == 0) or (n_stable == 0):
                p = np.nan
                log10p = np.nan
            else:
                _, p = stats.fisher_exact(
                    [[a, b], [c, d]],
                    alternative='two-sided'  # Two-sided test
                )
                log10p = log10_fisher_two_sided(a, b, c, d)

            all_plot_data.append({
                'Feature': feat,
                'Source': label,
                'Fold': fold,
                'log2Fold': log2fold,
                'P_value': p
            })

            all_stats_data.append({
                'CellType': label,
                'Feature': feat,
                'a_Active_hit': a,
                'b_Active_nonhit': b,
                'c_Stable_hit': c,
                'd_Stable_nonhit': d,
                'Active_total': n_active,
                'Stable_total': n_stable,
                'Fold': fold,
                'log2Fold': log2fold,
                'P_value': p,
                'log10P': -log10p if log10p == log10p else np.nan
            })

    # =========================
    # Assemble results in a DataFrame.
    # =========================

    res_plot_df = pd.DataFrame(all_plot_data)
    stats_df = pd.DataFrame(all_stats_data)

    # =========================
    # Adjust for multiple testing.
    # =========================

    if not stats_df.empty:
        pvals = stats_df['P_value'].fillna(1)
        _, fdr, _, _ = multipletests(pvals, method='fdr_bh')

        # Recompute BH adjustment in log space to retain values below floating-point range.
        # The log10P column stores -log10(P). Convert to log10(P) for adjustment.
        lp_true = -stats_df['log10P'].to_numpy(dtype=float)     # log10(P)
        m = len(lp_true)
        order = np.argsort(lp_true)          # Sort by increasing P value.
        run = np.inf
        log10q = np.empty(m)
        for r in range(m, 0, -1):            # Enforce BH monotonicity using a running minimum in reverse order.
            j = order[r - 1]
            run = min(run, np.log10(m) - np.log10(r) + lp_true[j])
            log10q[j] = run                  # log10(FDR)

        # Store representable P and FDR values as numbers and underflowed values as text.
        # Text retains values outside Excel numeric precision.
        stats_df['P_value'] = [
            p if (p == p and p > 0) else fmt10(l)
            for p, l in zip(pvals, lp_true)]
        stats_df['FDR'] = [
            q if (q == q and q > 0) else fmt10(l)
            for q, l in zip(fdr, log10q)]
        stats_df['log10FDR'] = -log10q       # Store -log10(FDR), following the log10P convention.

        stats_df['Star_FDR'] = pd.Series(fdr).apply(
            lambda x: '***' if x < 0.001 else
                      '**' if x < 0.01 else
                      '*' if x < 0.05 else ''
        ).values

    # =========================
    # Export statistics.
    # =========================

    stats_df.to_excel(xlsx_output, index=False)
    print(f"\nStatistics saved to: {xlsx_output}")

    # =========================
    # Export Source Data.
    # =========================

    seq_level_df = pd.concat(all_seq_level, ignore_index=True)

    # One row per bar, with bar height, test results and all four contingency counts.
    bars = stats_df[['CellType', 'Feature',
                     'a_Active_hit', 'b_Active_nonhit',
                     'c_Stable_hit', 'd_Stable_nonhit',
                     'Active_total', 'Stable_total',
                     'Fold', 'log2Fold',
                     'P_value', 'log10P', 'FDR', 'log10FDR', 'Star_FDR']].copy()
    bars['Active_hit_pct'] = bars.a_Active_hit / bars.Active_total * 100
    bars['Stable_hit_pct'] = bars.c_Stable_hit / bars.Stable_total * 100
    bars = bars.sort_values(['CellType', 'log2Fold'], ascending=[True, False])

    readme = pd.DataFrame({'Field': [
        'Figure', 'Panel', 'Test', 'Multiple testing', 'Bar height',
        'Grouping', 'Overlap definition', 'Denominator', 'Deduplication',
        'Star_FDR', 'P/FDR display',
        'Input TSV', 'Per-cell variant selection',
        'Oligo filter', 'Input BED (PK15)', 'Input BED (PSC)', 'Histone peaks',
    ], 'Description': [
        'Active_vs_Stable_Histone_Enrichment_TwoSided.pdf',
        'Single panel: grouped bars, x = feature, y = log2(enrichment ratio),'
        ' hue = cell type, features ordered by PK15 log2Fold descending',
        "Fisher's exact test, two-sided, on the 2x2 table"
        ' [[a_Active_hit, b_Active_nonhit], [c_Stable_hit, d_Stable_nonhit]]',
        'Benjamini-Hochberg across all cell type x feature tests in this figure',
        'log2Fold = log2[ (a/Active_total) / (c/Stable_total) ],'
        ' with epsilon = 1e-6 added to numerator and denominator',
        "MPRA annotation merged into two groups: 'Active' kept as is;"
        " every other label (Stable, Repressed) collapsed to 'Stable'",
        'bedtools intersect -wa -wb between the MPRA oligo interval and the'
        ' histone/ATAC/CTCF peak interval; any overlap of >=1 bp counts as a hit',
        'Sequences present in the intersect output, i.e. oligos overlapping at'
        ' least one feature of any kind (NOT the full candidate set)',
        'Counted once per (SEQ ID, feature) pair, so an oligo hit by several'
        ' peaks of the same mark contributes 1',
        '*** FDR < 0.001, ** FDR < 0.01, * FDR < 0.05, blank = not significant',
        'P_value / FDR computed exactly in log space; values within double range'
        ' stay numeric, but double precision underflows to 0 below ~1e-308'
        ' (PK15 CTCF p=3.16e-475, PSC CTCF p=4.70e-361), so those cells are text'
        ' with the exact value; log10P / log10FDR are numeric -log10() columns.'
        ' NOTE: pandas.read_excel coerces numeric-looking text back to float and'
        ' underflows it to 0 — read with dtype=str to see the exact values',
        'Processed CRE activity annotation table used to generate the BED inputs',
        'Each cell line uses only rows with <cell>.cell_filtering == "pass",'
        ' i.e. one row per (chr, variant_pos, strand, ref_allele, alt_allele)'
        ' group per cell (see VariantSelection_<cell>.tsv)',
        'Pass rows additionally require <cell>.count_ok == True;'
        ' ref_id gets act_ref_*, alt_id gets act_alt_*;'
        ' Active log2FC>1 & padj<0.01, Repressed log2FC<-1 & padj<0.01, else Stable',
        input_configs[0]['bed'],
        input_configs[1]['bed'],
        'Merged muscle chromatin-feature peak BED used for intersection',
    ]})

    with pd.ExcelWriter(source_xlsx) as w:
        readme.to_excel(w, sheet_name='README', index=False)
        bars.to_excel(w, sheet_name='Fig_bars', index=False)
        seq_level_df.to_excel(w, sheet_name='Per_sequence_overlap', index=False)

    bars.to_csv(source_tsv, sep='\t', index=False, float_format='%.6g')

    print(f"Source Data workbook saved to: {source_xlsx}")
    print(f"   (README / Fig_bars {len(bars)} rows / "
          f"Per_sequence_overlap {len(seq_level_df)} rows)")
    print(f"Source Data TSV saved to: {source_tsv}")

    # =========================
    # Sort by decreasing PK15 log2 enrichment.
    # =========================

    pk15_ref = res_plot_df[res_plot_df['Source'] == 'PK15'].sort_values(
        by='log2Fold', ascending=False
    )
    pk15_order = pk15_ref['Feature'].tolist()

    # =========================
    # Draw the plot.
    # =========================

    sns.set_theme(style="white")
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.linewidth'] = 1.2

    fig, ax = plt.subplots(figsize=(10, 5))

    palette = {
        "PK15": "#4C72B0",
        "PSC": "#DD8452"
    }

    h_order = ["PK15", "PSC"]

    # barplot
    sns.barplot(
        data=res_plot_df,
        x='Feature',
        y='log2Fold',
        hue='Source',
        palette=palette,
        order=pk15_order,
        hue_order=h_order,
        errorbar=None,
        linewidth=0.8,
        ax=ax
    )

    # baseline
    ax.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.7)

    # Add significance markers.
    num_sources = 2
    width = 0.8 / num_sources

    feature_to_idx = {feat: i for i, feat in enumerate(pk15_order)}

    for j, label in enumerate(h_order):
        subset = stats_df[stats_df['CellType'] == label]
        offset = (j - (num_sources - 1) / 2) * width

        for _, row in subset.iterrows():
            feat = row['Feature']
            if feat not in feature_to_idx:
                continue

            star = row['Star_FDR']
            if star == '':
                continue

            i = feature_to_idx[feat]
            y = row['log2Fold']

            # Position labels according to the sign of enrichment.
            if y >= 0:
                ypos = y + 0.06
                va = 'bottom'
            else:
                ypos = y - 0.06
                va = 'top'

            ax.text(
                i + offset, ypos, star,
                ha='center', va=va,
                fontsize=10, fontweight='bold'
            )

    # Format the axes.
    ax.set_xlabel("")
    ax.set_ylabel("log2(Enrichment Ratio)", fontsize=13)
    ax.set_title("Histone Modification Enrichment: Active vs Stable (Two-sided)", fontsize=14, pad=15)

    plt.xticks(rotation=45, ha='right')
    sns.despine()
    ax.legend(title="Cell type", frameon=False)
    plt.tight_layout()

    # Save the PDF.
    with PdfPages(pdf_output) as pdf:
        pdf.savefig(fig, bbox_inches='tight')

    print(f"\nFigure saved to: {pdf_output}")


# =========================
# Script entry point.
# =========================

if __name__ == "__main__":
    run_analysis()
