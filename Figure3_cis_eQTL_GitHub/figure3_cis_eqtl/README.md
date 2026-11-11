# Figure 3 — Cis-eQTLs in DLC pigs

Scripts accompanying Figure 3 of *Multi-omic profiling and MPRA of hybrid pigs uncover extensive functional cis-regulatory variation undetectable in purebreds*. DLC denotes the Duroc × Luchuan F1 population and DU denotes Duroc. The scripts use processed inputs. See the manuscript's Data Availability statement for study data and the [PigGTEx portal](https://piggtex.farmgtex.org/) for external muscle cis-eQTL results.

## Setup

Use Python 3.10 or newer for panels b–f, base R for panel a, and GNU `sort` for panel e preprocessing. Arial is used by the original plot settings in some panels and may be substituted by the plotting backend when unavailable.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Run the commands from this directory. Paths below are examples. Use separate output directories for each panel.

| Panel | Processed input | Main output |
| --- | --- | --- |
| 3a | Lead cis-eQTL scatterplot TSV with `Gene.ID`, `pos_to_TSS` (bp) and `BBFDR`. Use one lead variant per gene. | Scatterplot for distances within ±1 Mb and a kernel density TSV. |
| 3b | DLC cis-eQTL TSV with `max_MAF` and `beta`. | Effect size versus MAF plot and Source Data. The analysis uses 0 ≤ abs(beta) ≤ 2 and max_MAF ≥ 0.05. Up to 15,000 points are sampled for the scatter layer. |
| 3c | DLC cis-eQTL TSV with `SNP_id` and `Gene.ID`. | Log-binned histogram of unique associated SNP counts per eGene, with Excel and CSV Source Data. |
| 3d | DLC table with `SNP_id` and `Gene.ID`, and PigGTEx muscle significant pairs with `variant_id` and `phenotype_id`. | Stacked overlap bars, summary and annotated site/pair tables. |
| 3e (preprocessing) | DU allele-frequency TSV with `CHROM`, `POS` and `MAF`. | `DU_MAF005_only_pos.txt` and summary. |
| 3f | Original DLC cis-eQTL TSV and PigGTEx muscle significant pairs. | Harmonized shared pairs, statistics, QC and correlation plots. |

## Commands for panels a–e

```bash
Rscript scripts/plot_fig3a_lead_eqtl_tss_distance.R \
  data/eQTL_MAF005_distance_BBFDR_1Mb_scatter_plot_data.tsv results/fig3a

python scripts/plot_fig3b_effect_size_vs_maf.py \
  --eqtl data/eQTL_MAF005_kept.tsv --outdir results/fig3b

python scripts/plot_fig3c_unique_snps_per_gene.py \
  --eqtl data/eQTL_MAF005_kept.tsv --outdir results/fig3c

python scripts/plot_fig3d_piggtex_overlap.py \
  --dlc data/eQTL_MAF005_kept.tsv \
  --piggtex data/Muscle.cis_qtl_pairs.significant.txt.gz \
  --outdir results/fig3d

python scripts/fig3e_prepare_du_maf_positions.py \
  --du-frequency data/DU_allele_level_frequency.tsv --outdir results/fig3e
```

Panel e preprocessing retains a position if at least one DU allele record has MAF ≥ 0.05. The heatmap plotting script for panel e is an outstanding dependency.

## Panel f effect harmonization

Keep the original first 14 columns in the DLC TSV. The script reads one-based columns 1–5, 10, 13 and 14 as SNP ID, chromosome, position, gene ID, gene name, beta, P value and BBFDR. The PigGTEx table must contain `phenotype_id`, `variant_id`, `af`, `pval_nominal`, `slope` and `slope_se`.

Set `DLC_EFFECT_ALLELE_POSITION` to `first` or `second` according to the original DLC genotype coding. This tells the script which allele in `CHROM:POS_ALLELE1_ALLELE2` the DLC beta describes. The script matches by gene, chromosome, position and unordered allele pair, then orients DLC beta to the PigGTEx ALT allele.

```bash
python scripts/plot_fig3f_piggtex_effect_correlation.py \
  --dlc data/eQTL_MAF005_kept.tsv \
  --pig data/Muscle.cis_qtl_pairs.significant.txt.gz \
  --dlc-effect-allele-position "${DLC_EFFECT_ALLELE_POSITION:?Set to first or second from the original DLC genotype coding}" \
  --expected-shared-pairs 171544 \
  --outdir results/fig3f
```

The command uses the manuscript count (171,544) for QC. The original script defaults to 171,641 as a previous overlap reference. This option records a comparison and does not select or filter pairs.
