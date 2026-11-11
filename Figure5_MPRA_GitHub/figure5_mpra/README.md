# Figure 5 — MPRA activity and reproducibility

Scripts for panels b, c, d and f of Figure 5 from *Multi-omic profiling and MPRA of hybrid pigs uncover extensive functional cis-regulatory variation undetectable in purebreds*. These scripts use processed barcode counts, MPRA result tables or feature intersections.

## Setup

Use Python 3.10 or newer and run these commands from this directory. Replace the example paths with the study inputs.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Panel b uses R, Bioconductor DESeq2, ggplot2, pheatmap and scales. The script can use data.table to speed up input. Supply the MPRASuite `MPRAmodel.R` used in the original analysis and install its dependencies.

## Fig. 5b — Replicate correlation

Provide four tab-separated barcode count tables in the command order shown below. Each table requires `Barcode` and `Oligo`, with plasmid replicate columns starting with `plasmid_` and RNA replicate columns starting with the matching cell-line prefix (`PK15_` or `PSC_`, case insensitive). The script aggregates counts per oligonucleotide and applies MPRASuite normalization with summit shift. The manuscript panel uses the CRE plasmid/PK15/PSC Pearson correlation heatmap from the log2(normalized count + 1) outputs. The script also writes alternate transformations, normalized counts and correlation matrices.

```bash
Rscript scripts/plot_fig5b_mpra_replicate_correlations.R \
  scripts/MPRAmodel.R data/cre_pk15.counts data/cre_psc.counts \
  data/utr_pk15.counts data/utr_psc.counts results/fig5b
```

## Fig. 5c — Activity and allelic skew

Provide the processed CRE and 3′UTR `filtered2` TSV files. Include the group columns `chr`, `variant_pos`, `strand`, `ref_allele`, `alt_allele` and the original per-cell columns. Required suffixes for both `PK15.` and `PSC.` are `cell_filtering`, `count_ok`, `emVar_all_qc`, `act`, `Log2Skew`, `act_ref_log2FC`, `act_alt_log2FC`, `ref_id` and `alt_id`. Keep the original result and annotation columns for complete Source Data.

```bash
python scripts/plot_fig5c_mpra_activity_and_skew.py \
  --cre data/RESULT_CROSSCELL_M3_cre_A0_F_eqtlase_filtered2.tsv \
  --utr data/RESULT_CROSSCELL_M3_utr_A0_F_eqtlase_filtered2.tsv \
  --outdir results/fig5c
```

The script uses the precomputed QC and emVar flags after per-cell filtering. It exports cell-specific and combined figures, group counts and Source Data under `cre/` and `utr/`. The combined CRE figure corresponds to the manuscript panel.

## Fig. 5d — CRE chromatin-feature enrichment

Provide the two headerless activity BED intersections. Their nine columns are oligonucleotide chromosome, start, end, activity label, sequence ID, peak chromosome, start, end and feature type. The analysis background comprises sequences present in the intersection files. The script groups labels other than Active into the inactive group, recorded as Stable in its tables. It counts overlaps once per sequence and feature, then applies two-sided Fisher exact tests with Benjamini–Hochberg adjustment.

```bash
python scripts/plot_fig5d_cre_chromatin_enrichment.py \
  --pk15-bed data/cre_pk15_activity_chromatin_intersect.bed \
  --psc-bed data/cre_psc_activity_chromatin_intersect.bed \
  --outdir results/fig5d
```

Outputs include the bar plot, enrichment statistics, bar-level Source Data and sequence-level overlap flags.

## Fig. 5f — GC content and 3′UTR activity

Provide one XLSX workbook per cell line with `ID`, `GC_percent` and `log2FoldChange`. The script averages numeric values per oligonucleotide ID, then calculates Pearson correlation and fits a linear regression on raw activity. Values below −8 undergo a display-only rescaling. Supply the original `pearson_p_exact.py` helper exposing `pearson_p_value(r, n)`, which returns a P value, its log10 and formatted text.

```bash
python scripts/plot_fig5f_utr_gc_activity.py \
  --pk15 data/utr_gc_activity_PK15.xlsx --psc data/utr_gc_activity_PSC.xlsx \
  --pearson-helper scripts/pearson_p_exact.py --outdir results/fig5f
```

The script writes plots, Source Data and legends for both cell lines. The manuscript panel uses the PK15 plot.
