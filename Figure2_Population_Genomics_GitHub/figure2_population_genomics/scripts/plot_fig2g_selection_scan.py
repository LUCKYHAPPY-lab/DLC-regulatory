#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# =========================================================
# 0. Global plot settings
# =========================================================
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# Embed TrueType fonts to retain editable PDF text.
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


# =========================================================
# 1. Input and output paths
# =========================================================
parser = argparse.ArgumentParser(description="Plot the joint FST and nucleotide diversity selection scan.")
parser.add_argument("--fst", required=True, help="Window-level FST workbook")
parser.add_argument("--pi-ratio", required=True, help="Window-level DU/LC nucleotide diversity ratio workbook")
parser.add_argument("--outdir", required=True)
args = parser.parse_args()
fst_file = args.fst
pi_file = args.pi_ratio
output_dir = args.outdir

os.makedirs(output_dir, exist_ok=True)

out_luchuan = os.path.join(
    output_dir,
    'Luchuan_selected_regions_FST10_PI10.xlsx'
)

out_duroc = os.path.join(
    output_dir,
    'Duroc_selected_regions_FST10_PI10.xlsx'
)

out_all_windows = os.path.join(
    output_dir,
    'All_windows_with_FST10_PI10_category.xlsx'
)

out_thresholds = os.path.join(
    output_dir,
    'FST10_PI10_thresholds.txt'
)

output_pdf = os.path.join(
    output_dir,
    'Selection_Joint_Plot_FST10_PI10.pdf'
)

output_png = os.path.join(
    output_dir,
    'Selection_Joint_Plot_FST10_PI10.png'
)


# =========================================================
# 2. Read the data
# =========================================================
print("Reading inputs...")

df_fst = pd.read_excel(fst_file)
df_pi = pd.read_excel(pi_file)

# Trim whitespace from column names.
df_fst.columns = df_fst.columns.astype(str).str.strip()
df_pi.columns = df_pi.columns.astype(str).str.strip()

print(f"FST windows: {len(df_fst):,}")
print(f"Nucleotide diversity ratio windows: {len(df_pi):,}")

# Check for the Block column.
if 'Block' not in df_fst.columns:
    raise KeyError(
        "FST input is missing the Block column.\n"
        f"Available columns: {df_fst.columns.tolist()}"
    )

if 'Block' not in df_pi.columns:
    raise KeyError(
        "Nucleotide diversity ratio input is missing the Block column.\n"
        f"Available columns: {df_pi.columns.tolist()}"
    )

# Convert block IDs to strings and trim whitespace.
df_fst['Block'] = df_fst['Block'].astype(str).str.strip()
df_pi['Block'] = df_pi['Block'].astype(str).str.strip()

# Check for duplicate block IDs.
fst_duplicate_n = df_fst['Block'].duplicated().sum()
pi_duplicate_n = df_pi['Block'].duplicated().sum()

if fst_duplicate_n > 0:
    print(
        f"Warning: FST input contains "
        f"{fst_duplicate_n:,} duplicate Block entries."
    )

if pi_duplicate_n > 0:
    print(
        f"Warning: nucleotide diversity ratio input contains "
        f"{pi_duplicate_n:,} duplicate Block entries."
    )

# Join the inputs by block ID.
df = pd.merge(
    df_fst,
    df_pi,
    on='Block',
    how='inner',
    suffixes=('_fst', '_pi')
)

print(f"Matched windows: {len(df):,}")


# =========================================================
# 3. Identify the FST and nucleotide diversity ratio columns
# =========================================================
if 'WEIGHTED_FST' in df.columns:
    fst_col = 'WEIGHTED_FST'
elif 'MEAN_FST' in df.columns:
    fst_col = 'MEAN_FST'
else:
    raise KeyError(
        "Merged data lack both WEIGHTED_FST and MEAN_FST.\n"
        f"Available columns: {df.columns.tolist()}"
    )

if 'pi_ratio' in df.columns:
    pi_ratio_col = 'pi_ratio'
elif 'PI_ratio' in df.columns:
    pi_ratio_col = 'PI_ratio'
elif 'Pi_ratio' in df.columns:
    pi_ratio_col = 'Pi_ratio'
else:
    raise KeyError(
        "Merged data lack a pi_ratio column.\n"
        f"Available columns: {df.columns.tolist()}"
    )

print(f"FST column: {fst_col}")
print(f"Nucleotide diversity ratio column: {pi_ratio_col}")


# =========================================================
# 4. Clean and transform the data
# =========================================================
df[fst_col] = pd.to_numeric(
    df[fst_col],
    errors='coerce'
)

df[pi_ratio_col] = pd.to_numeric(
    df[pi_ratio_col],
    errors='coerce'
)

n_before_filter = len(df)

# Remove windows with missing FST or nucleotide diversity ratios.
df = df.dropna(
    subset=[fst_col, pi_ratio_col]
).copy()

# Require positive nucleotide diversity ratios for log10.
df = df[df[pi_ratio_col] > 0].copy()

# Remove nonfinite values.
df = df[
    np.isfinite(df[fst_col]) &
    np.isfinite(df[pi_ratio_col])
].copy()

n_after_filter = len(df)

print(
    f"Windows retained after filtering: "
    f"{n_after_filter:,}/{n_before_filter:,}"
)

if df.empty:
    raise ValueError(
        "No windows remain after filtering. Check the input files."
    )

# Define the x-axis as the log10 ratio.
# log10(theta_pi,DU / theta_pi,LC)
df['log10_pi_ratio'] = np.log10(
    df[pi_ratio_col]
)

x = df['log10_pi_ratio'].to_numpy()
y = df[fst_col].to_numpy()


# =========================================================
# 5. Calculate the upper FST decile and both nucleotide diversity ratio deciles
# =========================================================
# Lower decile of the nucleotide diversity ratio.
x_low_thresh = np.quantile(x, 0.10)

# Upper decile of the nucleotide diversity ratio.
x_high_thresh = np.quantile(x, 0.90)

# Use the 90th percentile of FST as the cutoff.
y_thresh = np.quantile(y, 0.90)

print("\n========================================")
print("Joint thresholds for FST and the nucleotide diversity ratio")
print("========================================")

print(
    f"Upper FST decile threshold: "
    f"{fst_col} >= {y_thresh:.6f}"
)

print(
    "Lower decile of log10(pi_DU/pi_LC): "
    f"<= {x_low_thresh:.6f}"
)

print(
    "Upper decile of log10(pi_DU/pi_LC): "
    f">= {x_high_thresh:.6f}"
)


# =========================================================
# 6. Apply joint filters and assign the selection direction
# =========================================================
# The x-axis represents:
# log10(theta_pi,DU / theta_pi,LC)
#
# Left tail:
# A low DU/LC nucleotide diversity ratio
# indicates relatively low diversity in Duroc
# and defines Duroc candidate windows when FST also passes.
#
# Right tail:
# A high DU/LC nucleotide diversity ratio
# indicates relatively low diversity in Luchuan
# and defines Luchuan candidate windows when FST also passes.

mask_du = (
    (df[fst_col] >= y_thresh) &
    (df['log10_pi_ratio'] <= x_low_thresh)
)

mask_lc = (
    (df[fst_col] >= y_thresh) &
    (df['log10_pi_ratio'] >= x_high_thresh)
)

# Initialize all windows as genome-wide background.
df['Category'] = 'Whole_genome'

df.loc[
    mask_du,
    'Category'
] = 'Selected_Duroc'

df.loc[
    mask_lc,
    'Category'
] = 'Selected_Luchuan'

# Record individual filter flags for verification.
df['FST_top10'] = (
    df[fst_col] >= y_thresh
)

df['Pi_left10_Duroc'] = (
    df['log10_pi_ratio'] <= x_low_thresh
)

df['Pi_right10_Luchuan'] = (
    df['log10_pi_ratio'] >= x_high_thresh
)

df['Joint_selected_Duroc'] = mask_du
df['Joint_selected_Luchuan'] = mask_lc


# =========================================================
# 7. Export results
# =========================================================
# Duroc candidate windows.
df.loc[mask_du].to_excel(
    out_duroc,
    index=False
)

# Luchuan candidate windows.
df.loc[mask_lc].to_excel(
    out_luchuan,
    index=False
)

# All windows with category labels.
df.to_excel(
    out_all_windows,
    index=False
)

# Save thresholds and window counts.
with open(
    out_thresholds,
    'w',
    encoding='utf-8'
) as f:

    f.write(
        "FST top 10% and Pi-ratio one-tailed 10%\n"
    )
    f.write(
        "========================================\n\n"
    )

    f.write(
        "Pi ratio definition: "
        "theta_pi_DU / theta_pi_LC\n"
    )

    f.write(
        "Transformation: "
        "log10(theta_pi_DU / theta_pi_LC)\n\n"
    )

    f.write(
        f"FST column: {fst_col}\n"
    )

    f.write(
        f"FST top 10% threshold: "
        f"{y_thresh:.10f}\n"
    )

    f.write(
        f"Pi ratio left 10% threshold: "
        f"{x_low_thresh:.10f}\n"
    )

    f.write(
        f"Pi ratio right 10% threshold: "
        f"{x_high_thresh:.10f}\n\n"
    )

    f.write(
        f"Total valid windows: {len(df)}\n"
    )

    f.write(
        f"Duroc selected windows: "
        f"{mask_du.sum()}\n"
    )

    f.write(
        f"Luchuan selected windows: "
        f"{mask_lc.sum()}\n"
    )

print("\n========================================")
print("Joint filtering results")
print("========================================")

print(
    f"Duroc candidate selection windows: "
    f"{mask_du.sum():,}"
)

print(
    f"Luchuan candidate selection windows: "
    f"{mask_lc.sum():,}"
)

print("\nOutput files:")
print(f"Duroc: {out_duroc}")
print(f"Luchuan: {out_luchuan}")
print(f"All classified windows: {out_all_windows}")
print(f"Thresholds: {out_thresholds}")


# =========================================================
# 8. Plot layout
# =========================================================
print("\nPlotting the joint selection scan...")

fig = plt.figure(
    figsize=(7, 7),
    dpi=300
)

gs = gridspec.GridSpec(
    2,
    2,
    width_ratios=[4, 1.2],
    height_ratios=[1.2, 4],
    wspace=0.05,
    hspace=0.05
)

ax_main = plt.subplot(
    gs[1, 0]
)

ax_histx = plt.subplot(
    gs[0, 0],
    sharex=ax_main
)

ax_histy = plt.subplot(
    gs[1, 1],
    sharey=ax_main
)


# =========================================================
# 9. Colors
# =========================================================
# Left tail: Duroc candidate windows.
color_du = '#2ca02c'

# Right tail: Luchuan candidate windows.
color_lc = '#1f77b4'

# Histogram for the upper FST decile.
color_fst = '#d62728'

# Genome-wide background.
color_genome = '#999999'

ax_main.set_rasterization_zorder(10)


# =========================================================
# 10. Main scatterplot
# =========================================================
mask_normal = (
    df['Category'] == 'Whole_genome'
)

# Genome-wide background.
ax_main.scatter(
    df.loc[
        mask_normal,
        'log10_pi_ratio'
    ],
    df.loc[
        mask_normal,
        fst_col
    ],
    c=color_genome,
    s=4,
    alpha=0.15,
    linewidths=0,
    label='Whole genome',
    rasterized=True,
    zorder=1
)

# Duroc candidate windows.
ax_main.scatter(
    df.loc[
        mask_du,
        'log10_pi_ratio'
    ],
    df.loc[
        mask_du,
        fst_col
    ],
    c=color_du,
    s=12,
    alpha=0.80,
    linewidths=0,
    label='Selected (Duroc)',
    rasterized=True,
    zorder=2
)

# Luchuan candidate windows.
ax_main.scatter(
    df.loc[
        mask_lc,
        'log10_pi_ratio'
    ],
    df.loc[
        mask_lc,
        fst_col
    ],
    c=color_lc,
    s=12,
    alpha=0.80,
    linewidths=0,
    label='Selected (Luchuan)',
    rasterized=True,
    zorder=2
)


# =========================================================
# 11. Label MC1R
# =========================================================
mc1r_block = '6:180001_220000'

mc1r_row = df[
    df['Block'].astype(str).str.strip() ==
    mc1r_block
]

if not mc1r_row.empty:

    mc1r_x = (
        mc1r_row['log10_pi_ratio'].iloc[0]
    )

    mc1r_y = (
        mc1r_row[fst_col].iloc[0]
    )

    ax_main.annotate(
        r'$MC1R$',
        xy=(
            mc1r_x,
            mc1r_y
        ),
        xytext=(
            mc1r_x + 0.30,
            mc1r_y - 0.05
        ),
        fontsize=10,
        fontweight='bold',
        fontstyle='italic',
        arrowprops=dict(
            arrowstyle='->',
            connectionstyle='arc3,rad=0.2',
            color='black',
            lw=0.8
        ),
        zorder=21
    )

    print("\nMC1R window:")

    print(
        mc1r_row[
            [
                'Block',
                fst_col,
                pi_ratio_col,
                'log10_pi_ratio',
                'Category'
            ]
        ].to_string(index=False)
    )

else:
    print(
        f"\nWarning: MC1R window not found: "
        f"{mc1r_block}"
    )


# =========================================================
# 12. Dashed threshold lines
# =========================================================
# Lower decile of the nucleotide diversity ratio.
ax_main.axvline(
    x_low_thresh,
    color='gray',
    linestyle='--',
    linewidth=0.8,
    zorder=15
)

# Upper decile of the nucleotide diversity ratio.
ax_main.axvline(
    x_high_thresh,
    color='gray',
    linestyle='--',
    linewidth=0.8,
    zorder=15
)

# Upper FST decile.
ax_main.axhline(
    y_thresh,
    color='gray',
    linestyle='--',
    linewidth=0.8,
    zorder=15
)


# =========================================================
# 13. Main plot axes
# =========================================================
ax_main.set_xlabel(
    r'$\log_{10}(\theta_{\pi,\mathrm{DU}}/'
    r'\theta_{\pi,\mathrm{LC}})$',
    fontsize=10
)

ax_main.set_ylabel(
    r'$F_{\mathrm{ST}}$',
    fontsize=10
)

ax_main.legend(
    loc='lower left',
    frameon=False,
    fontsize=8,
    markerscale=1.2
)

# Limit the display range so extreme x values do not compress the central distribution.
x_min_display = np.percentile(
    x,
    0.1
)

x_max_display = np.percentile(
    x,
    99.9
)

if x_min_display == x_max_display:
    x_min_display = x.min()
    x_max_display = x.max()

ax_main.set_xlim(
    x_min_display,
    x_max_display
)

# Display the full y range.
y_min_display = min(
    0,
    y.min()
)

y_max_display = y.max()

if y_max_display <= y_min_display:
    y_max_display = y_min_display + 1

y_margin = (
    y_max_display - y_min_display
) * 0.03

ax_main.set_ylim(
    y_min_display - y_margin,
    y_max_display + y_margin
)

ax_main.tick_params(
    axis='both',
    labelsize=8
)


# =========================================================
# 14. Top histogram of nucleotide diversity ratios
# =========================================================
bins_x = np.linspace(
    x_min_display,
    x_max_display,
    100
)

counts_x, edges_x = np.histogram(
    x,
    bins=bins_x
)

freq_x = (
    counts_x / len(x)
) * 100

for i in range(len(counts_x)):

    bin_left = edges_x[i]
    bin_right = edges_x[i + 1]

    if bin_right <= x_low_thresh:
        # Lower decile: Duroc.
        current_color = color_du

    elif bin_left >= x_high_thresh:
        # Upper decile: Luchuan.
        current_color = color_lc

    else:
        current_color = 'silver'

    ax_histx.bar(
        bin_left,
        freq_x[i],
        width=bin_right - bin_left,
        align='edge',
        color=current_color,
        edgecolor='none'
    )

# Top cumulative distribution.
ax_histx_cum = ax_histx.twinx()

sorted_x = np.sort(x)

cum_y_x = (
    np.arange(
        1,
        len(sorted_x) + 1
    ) /
    len(sorted_x)
) * 100

ax_histx_cum.plot(
    sorted_x,
    cum_y_x,
    color='black',
    linewidth=1
)

ax_histx_cum.set_ylim(
    0,
    100
)

ax_histx.set_ylabel(
    'Frequency (%)',
    fontsize=8
)

ax_histx_cum.set_ylabel(
    'Cumulative (%)',
    fontsize=8
)

ax_histx.axvline(
    x_low_thresh,
    color='gray',
    linestyle='--',
    linewidth=0.8
)

ax_histx.axvline(
    x_high_thresh,
    color='gray',
    linestyle='--',
    linewidth=0.8
)

ax_histx.tick_params(
    axis='y',
    labelsize=7
)

ax_histx_cum.tick_params(
    axis='y',
    labelsize=7
)

plt.setp(
    ax_histx.get_xticklabels(),
    visible=False
)


# =========================================================
# 15. Right-hand FST histogram
# =========================================================
bins_y = np.linspace(
    y_min_display,
    y_max_display,
    100
)

counts_y, edges_y = np.histogram(
    y,
    bins=bins_y
)

freq_y = (
    counts_y / len(y)
) * 100

for i in range(len(counts_y)):

    bin_low = edges_y[i]
    bin_high = edges_y[i + 1]

    # Highlight the upper FST decile in red.
    if bin_high >= y_thresh:
        current_color = color_fst
    else:
        current_color = 'silver'

    ax_histy.barh(
        bin_low,
        freq_y[i],
        height=bin_high - bin_low,
        align='edge',
        color=current_color,
        edgecolor='none'
    )

# Right-hand cumulative distribution.
ax_histy_cum = ax_histy.twiny()

sorted_y = np.sort(y)

cum_x_y = (
    np.arange(
        1,
        len(sorted_y) + 1
    ) /
    len(sorted_y)
) * 100

ax_histy_cum.plot(
    cum_x_y,
    sorted_y,
    color='black',
    linewidth=1
)

ax_histy_cum.set_xlim(
    0,
    100
)

ax_histy.set_xlabel(
    'Frequency (%)',
    fontsize=8
)

ax_histy_cum.set_xlabel(
    'Cumulative (%)',
    fontsize=8
)

ax_histy.axhline(
    y_thresh,
    color='gray',
    linestyle='--',
    linewidth=0.8
)

ax_histy.tick_params(
    axis='x',
    labelsize=7
)

ax_histy_cum.tick_params(
    axis='x',
    labelsize=7
)

plt.setp(
    ax_histy.get_yticklabels(),
    visible=False
)

plt.setp(
    ax_histy.get_xticklabels(),
    visible=False
)


# =========================================================
# 16. Plot borders
# =========================================================
for axis in [
    ax_main,
    ax_histx,
    ax_histy
]:
    for spine in axis.spines.values():
        spine.set_linewidth(0.8)


# =========================================================
# 17. Export PDF and PNG
# =========================================================
plt.savefig(
    output_pdf,
    format='pdf',
    dpi=300,
    bbox_inches='tight'
)

plt.savefig(
    output_png,
    format='png',
    dpi=600,
    bbox_inches='tight',
    transparent=False
)

plt.close(fig)

print("\n========================================")
print("Analysis and plotting complete.")
print("========================================")
print(f"PDF: {output_pdf}")
print(f"PNG: {output_png}")
print(f"Output directory: {output_dir}")
