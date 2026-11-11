#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot activity and allelic skew separately for the CRE and 3'UTR libraries.

Retain rows passing each cell's filter and classify them using the
precomputed count QC and emVar flags. Include untested rows in Source
Data. Plot tested rows using the signed activity with the larger absolute
allele-specific log2 fold change.
"""
import argparse
import os

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ==========================================
# 1. Input and output paths
parser = argparse.ArgumentParser(description='Plot MPRA activity and allelic skew.')
parser.add_argument("--cre", required=True)
parser.add_argument("--utr", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()
INPUT_CRE = args.cre
INPUT_UTR = args.utr
CELLS = ['PK15', 'PSC']
GROUP_KEY_COLS = ['chr', 'variant_pos', 'strand', 'ref_allele', 'alt_allele']


def process_cell(data, CELL):
    """Filter and classify one cell line, then return plot records and counts."""
    col_filter = f'{CELL}.cell_filtering'
    col_count_ok = f'{CELL}.count_ok'
    col_skew = f'{CELL}.Log2Skew'
    col_skew_p = f'{CELL}.skew_p'
    col_emvar_qc = f'{CELL}.emVar_all_qc'
    col_act = f'{CELL}.act'
    col_act_ref = f'{CELL}.act_ref_log2FC'
    col_act_alt = f'{CELL}.act_alt_log2FC'
    col_ref_id = f'{CELL}.ref_id'
    col_alt_id = f'{CELL}.alt_id'

    # Keep rows passing the filter for this cell line.
    df = data[data[col_filter] == 'pass'].copy()

    # Classify emVars, tested non-emVars and untested records.
    df['status_strict'] = 'non_emVar'
    df.loc[df[col_emvar_qc] == True, 'status_strict'] = 'emVar'
    df.loc[df[col_count_ok] != True, 'status_strict'] = 'no_tested'

    # Count groups within the cell line.
    tested = df[df['status_strict'] != 'no_tested']
    stats = pd.DataFrame([{
        'Cell': CELL,
        'n_groups': len(df),
        'n_emVar_groups': int((df['status_strict'] == 'emVar').sum()),
        'n_non_emVar_groups': int((df['status_strict'] == 'non_emVar').sum()),
        'n_no_tested_groups': int((df['status_strict'] == 'no_tested').sum()),
        'n_active_groups': int((tested[col_act] == True).sum()),
        'n_emVar_and_active': int(((df['status_strict'] == 'emVar') & (df[col_act] == True)).sum()),
    }])

    # Prepare plot columns.
    df['Log2Skew'] = df[col_skew]
    abs_ref = df[col_act_ref].abs()
    abs_alt = df[col_act_alt].abs()
    df['Max_Abs_Activity'] = np.where(abs_ref >= abs_alt, df[col_act_ref], df[col_act_alt])
    df['ref_id'] = df[col_ref_id]
    df['alt_id'] = df[col_alt_id]
    df['Cell'] = CELL
    return df, stats


def plot_joint(df, title, out_prefix):
    g = sns.JointGrid(data=df, x='Max_Abs_Activity', y='Log2Skew', hue='status_strict',
                      palette=my_pal, hue_order=['non_emVar', 'emVar'])
    g.plot_joint(sns.scatterplot, s=8, alpha=0.8, edgecolor='none', marker='o')
    g.plot_marginals(sns.kdeplot, fill=True, common_norm=False, alpha=0.6, bw_adjust=0.8)
    g.ax_joint.axhline(0, color='grey', lw=0.8, alpha=0.5)
    g.ax_joint.axvline(1.0, color='grey', lw=0.8, ls='--', alpha=0.4)
    g.ax_joint.axvline(-1.0, color='grey', lw=0.8, ls='--', alpha=0.4)
    g.set_axis_labels('Activity (log2FC)', 'Allelic effect (Δ log2FC)', fontsize=12)
    g.fig.suptitle(title, fontsize=11, y=1.02)
    plt.savefig(out_prefix + '.pdf', bbox_inches='tight', dpi=300)
    plt.savefig(out_prefix + '.png', bbox_inches='tight', dpi=300)
    plt.close()


# ==========================================
# 2. Read both inputs and analyze each library separately
# ==========================================
print('Reading input tables...')
frames = []
for name, path in [('cre', INPUT_CRE), ('utr', INPUT_UTR)]:
    d = pd.read_csv(path, sep='\t', low_memory=False)
    d['Dataset'] = name
    frames.append(d)
    print(f'  {name}: {len(d)} rows')
data = pd.concat(frames, ignore_index=True, sort=False)

sns.set_theme(style='ticks')
plt.rcParams['font.sans-serif'] = ['Arial']
my_pal = {'emVar': '#D62839', 'non_emVar': '#255C99'}

LIBS = [('cre', 'CRE'), ('utr', "3'UTR")]

for lib, lib_label in LIBS:
    outdir = os.path.join(args.outdir, lib)
    os.makedirs(outdir, exist_ok=True)
    sub = data[data['Dataset'] == lib]
    print(f'\n{"=" * 50}\nLibrary: {lib_label} ({len(sub)} rows). Output: {outdir}\n{"=" * 50}')

    all_dfs, all_stats = [], []
    for CELL in CELLS:
        df, stats = process_cell(sub, CELL)
        all_dfs.append(df)
        all_stats.append(stats)
        print(stats.to_string(index=False))

        tested = df[df['status_strict'] != 'no_tested']
        plot_joint(tested, f'{CELL} - Strict (emVar_all_qc) - {lib_label}',
                   os.path.join(outdir, f'MPRA_emVAR_Strict_{CELL}'))

    # Save per-cell-line counts.
    stats_all = pd.concat(all_stats, ignore_index=True)
    stats_all.to_csv(os.path.join(outdir, 'Stats_Groups_PerCell.tsv'), sep='\t', index=False)

    # Combine PK15 and PSC within this library.
    combined = pd.concat(all_dfs, ignore_index=True)
    combined_tested = combined[combined['status_strict'] != 'no_tested']
    print(f'  Combined tested points: {len(combined_tested)}'
          f' ({len(combined) - len(combined_tested)} untested rows retained only in Source Data)')
    plot_joint(combined_tested, f'PK15+PSC Combined - Strict (emVar_all_qc) - {lib_label}',
               os.path.join(outdir, 'MPRA_emVAR_Strict_Combined'))

    # source data
    ID_COLS = ['Dataset', 'key', 'chr', 'variant_pos', 'strand',
               'ref_allele', 'alt_allele', 'variant_type', 'INFO', 'project']
    cols = [c for c in ID_COLS if c in combined.columns] + [
        'Cell', 'ref_id', 'alt_id',
        'Max_Abs_Activity', 'Log2Skew', 'status_strict',
    ]
    points = combined[[c for c in cols if c in combined.columns]].copy()
    points = points.rename(columns={'Max_Abs_Activity': 'x_Activity_log2FC',
                                    'Log2Skew': 'y_Log2Skew'})
    points = points.sort_values(['Cell', 'status_strict', 'y_Log2Skew'],
                                ascending=[True, False, False])

    summary = (points.groupby(['Cell', 'status_strict']).size()
                     .rename('n_points').reset_index())
    tot = points.groupby('Cell').size().rename('n_total')
    summary = summary.merge(tot, on='Cell')
    summary['pct_of_cell'] = summary.n_points / summary.n_total * 100

    readme = pd.DataFrame({'Field': [
        'Figure', 'Library', 'Point color', 'Rows plotted', 'Filtering',
        'Grouping', 'Input file', 'Plotting script',
    ], 'Description': [
        'MPRA_emVAR_Strict_Combined.pdf / .png',
        f'{lib_label} panel only ({lib}); the CRE and 3UTR libraries are'
        ' plotted separately and never pooled',
        "Red = 'emVar', blue = 'non_emVar' (column status_strict); rows labeled"
        " 'no_tested' (count_ok != True, unreliable counts) are kept in Fig_points"
        ' but not drawn on the figure',
        'One row per (variant-group, cell line) after per-cell filtering within'
        ' this library; a variant tested in both cell lines contributes two points',
        "Per cell line, only rows with that cell's own cell_filtering == 'pass'"
        ' (PK15 uses PK15.cell_filtering, PSC uses PSC.cell_filtering); rows with'
        " count_ok != True are labeled 'no_tested' and excluded from the figure",
        'Group = (chr, variant_pos, strand, ref_allele, alt_allele), strand-aware;'
        ' one pass row per group within each panel',
        INPUT_CRE if lib == 'cre' else INPUT_UTR,
        os.path.abspath(__file__),
    ]})

    src_xlsx = os.path.join(outdir, 'SourceData_MPRA_emVAR_Strict_Combined.xlsx')
    src_tsv = os.path.join(outdir, 'SourceData_MPRA_emVAR_Strict_Combined.tsv')
    with pd.ExcelWriter(src_xlsx) as w:
        readme.to_excel(w, sheet_name='README', index=False)
        summary.to_excel(w, sheet_name='Summary_counts', index=False)
        points.to_excel(w, sheet_name='Fig_points', index=False)
    points.to_csv(src_tsv, sep='\t', index=False, float_format='%.6g')
    print(f'  ✅ {src_xlsx}')
    print(f'  ✅ {src_tsv}')
    print(summary.to_string(index=False))

print('\nCompleted.')
