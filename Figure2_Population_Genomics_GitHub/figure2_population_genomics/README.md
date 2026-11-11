# Figure 2 — Population genomics and selection scans

Scripts for panels e–g of Figure 2 from *Multi-omic profiling and MPRA of hybrid pigs uncover extensive functional cis-regulatory variation undetectable in purebreds*. DU and LC denote Duroc and Luchuan, respectively. DLC denotes their F1 cross. See the manuscript's Data Availability statement for the sequencing inputs.

## Setup

Use Python 3.10 or newer. Panel e also requires Bash at `/bin/bash` and bcftools with the `fill-tags` plugin.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Run the examples from this directory and replace the example paths with the study inputs.

## Fig. 2e — Per-site expected heterozygosity

Provide the combined, unimputed autosomal VCF and a directory containing `DU.txt`, `LC.txt` and `DLC.txt`. These lists must partition the VCF sample IDs into 100 DU, 100 LC and 189 DLC pigs.

```bash
python scripts/plot_fig2e_site_heterozygosity.py \
  --vcf data/combine_unimputed_Parents200_DLC189_chr1-18.vcf.gz \
  --sample-dir data/sample_lists --out-root results/fig2e
```

The script recalculates AC and AN within each population and computes Hp as 1 − sum(pᵢ²). It retains fixed records at Hp = 0 and excludes records with missing or invalid counts. Each run writes a timestamped directory with per-site TSV.gz files, summary statistics, run metadata and boxplots. The Kruskal–Wallis and Dunn tests and significance annotations described in the manuscript require a separate analysis script.

## Fig. 2f — Overlap of polymorphic variant positions

The manuscript uses scheme 4 (MAF ≥ 0.05 and at least 50 called individuals per population). Provide three prefiltered, headerless, deduplicated position lists, with one matching position key per line. Sort each list lexicographically with `LC_ALL=C sort -u`.

```text
data/position_sets/
  NS_ge50_no_genotype_count_requirement/
    position_sets/
      Duroc_MAF005_NS_ge50_pos.txt
      Luchuan_MAF005_NS_ge50_pos.txt
      DLC_MAF005_NS_ge50_pos.txt
```

```bash
python scripts/plot_fig2f_polymorphic_variant_overlap.py \
  --input-root data/position_sets --scheme 4 --outdir results/fig2f
```

Outputs include the seven disjoint region counts and Venn diagrams in the scheme subdirectory. The root summary for this command is `Venn_scheme4_region_counts.tsv`. The default `--scheme all` reproduces the original six-scheme workflow and requires the position lists for each scheme. Their directory and file names are defined in `SCHEMES` in the script.

## Fig. 2g — Joint FST and nucleotide diversity scan

Provide two window-level Excel workbooks with matching `Block` identifiers. The FST workbook requires `WEIGHTED_FST` or `MEAN_FST`. The diversity workbook requires the positive DU/LC ratio in `pi_ratio`, `PI_ratio` or `Pi_ratio`. The script joins the windows by `Block`, takes log10 of the diversity ratio, and recalculates quantiles from the valid matched windows.

```bash
python scripts/plot_fig2g_selection_scan.py \
  --fst data/100_LC_100_DU_fst_result.xlsx \
  --pi-ratio data/Pi_ratio_block_with_state_and_PI_top_split.xlsx \
  --outdir results/fig2g
```

Duroc candidate windows have FST at or above its 90th percentile and log10(πDU/πLC) at or below its 10th percentile. Luchuan candidate windows use the same FST cutoff and the upper diversity-ratio decile. Outputs include the joint plot, threshold report, all classified windows and per-breed candidate window tables. The adjacent-window merging step for the manuscript's sweep-region counts is upstream of the final region tables and is outside this plotting script.
