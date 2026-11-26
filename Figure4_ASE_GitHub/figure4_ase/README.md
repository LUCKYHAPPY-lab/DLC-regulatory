# Figure 4 — Allele-specific expression in DLC pigs

Scripts for panels a–f of Figure 4 from *Multi-omic profiling and MPRA of hybrid pigs uncover extensive functional cis-regulatory variation undetectable in purebreds*. DLC denotes the Duroc × Luchuan F1 population and DU denotes Duroc.

The scripts start from processed ASE and genotype summaries. The manuscript's Methods describe the upstream ASE workflow. For study data access, see the manuscript's Data Availability statement and [pigRegDb](https://gonglab.hzau.edu.cn/pigRegDb). The external comparison uses the [PigGTEx portal](https://piggtex.farmgtex.org/). Processed inputs and figure outputs are stored locally in the example paths below.

## Setup

Use Python 3.10 or newer for panels b–f. Panel a uses R with `ggplot2` and Cairo PDF support.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
Rscript -e 'install.packages("ggplot2", repos="https://cloud.r-project.org")'
```

Run the following commands from this directory. Replace the example input paths with the actual files. Use a fresh output directory for a new run, especially for panels b and d.

## Fig. 4a — Reference-allele ratio distribution

The sample map is comma-separated and contains an `Ear_ID` column. The count directory contains tab-separated per-sample files named `<Ear_ID>*ASEReadCounter_binom_test.txt`, each with a `ref_ratio` column. The script pools valid ratios across the listed individuals and writes `histogram_density_of_all_sample_balanced_bin.pdf`.

```bash
Rscript scripts/plot_fig4a_reference_allele_ratio.R \
  data/sample_map.csv data/per_sample_ase results/fig4a
```

## Fig. 4b — Mratio distribution and cumulative percentage

The input ASE table contains `SNP_id`, `Gene`, `FDR_binomtest` and `median_Mratio`. Use the Ensembl release 109 GTF for Sscrofa11.1 and an allele-level frequency TSV with `POPULATION`, `CHROM`, `POS`, `REF`, `ALT` and `MAF`. The plot uses FDR < 0.05, exonic annotation and DLC site maximum MAF ≥ 0.05. Mratio ≥ 0.60 defines the final ASE subset in the exported workbook. Outputs include `Mratio_distribution.pdf`, a single-sheet Source Data workbook and `Final_ASE_FDR005_Mratio060_MAF005.xlsx`.

```bash
python scripts/plot_fig4b_mratio_distribution.py \
  --ase-table data/ase_meta_analysis.tsv \
  --gtf data/Sus_scrofa.Sscrofa11.1.109.gtf \
  --dlc-frequency data/DLC_allele_level_frequency.tsv \
  --outdir results/fig4b
```

## Fig. 4c — Per-individual ASE counts

Use the processed `per_sample_snp_counts_MAF005.txt` table with `sample_id`, `biallelic_count_MAF005` and `ASE_list_observed_count_MAF005`. The script calculates Pearson correlation and fits a linear regression with a 95% prediction interval. It writes `ASE_scatter_Pearson_95PI.pdf` and the plotted points, regression grid and statistics as Source Data.

```bash
python scripts/plot_fig4c_per_sample_ase_counts.py \
  --sample-counts data/per_sample_snp_counts_MAF005.txt \
  --outdir results/fig4c
```

## Fig. 4d — Genomic distribution and density

Provide the full `Annotation_statistics.tsv` from the upstream ASE annotation step. It must contain `Category`, `Background_SNP_count` and `ASE_SNP_count`, including rows named `3_prime_UTR`, `CDS`, `5_prime_UTR` and `Noncoding_exon`. Bar percentages use the total ASE count across the input categories. The line uses ASE count divided by background SNP count within each displayed category. The script exports the four-category plot and its Source Data.

```bash
python scripts/plot_fig4d_genomic_distribution.py \
  --annotation-statistics data/Annotation_statistics.tsv \
  --outdir results/fig4d
```

## Fig. 4e — Muscle ASE overlap with PigGTEx

Use the study's final DLC ASE workbook, with a `Final_ASE_sites` sheet containing `SNP_id` and, for panel f, `DLC_site_max_MAF`. Panel b also exports a workbook with this schema. To substitute it for the study's final workbook, verify that the variant IDs and MAFs agree.

The PigGTEx gzip archive contains `chr`, `position`, `refAllele`, `altAllele`, `totalCount`, `FDR`, `SampleID`, `Breed` and `Tissue`. The script retains muscle records with FDR < 0.05 and total read count ≥ 10, excluding the `Cross` breed label. Matching uses chromosome, position and an unordered allele pair. Outputs include `Breed_ASE_overlap_counts.pdf`, a single-sheet Source Data workbook and the overlap-status TSV consumed by panel f.

```bash
python scripts/plot_fig4e_piggtex_breed_overlap.py \
  --dlc-ase data/ASE_final_FDR005_Mratio060_DLCmaxMAF005.xlsx \
  --piggtex-ase data/6020.ASE.site.sum.sort.pCADD.txt.gz \
  --outdir results/fig4e
```

## Fig. 4f — DLC and Duroc MAF heatmap for PigGTEx-unreported ASE sites

Run panel e first. The heatmap requires the original `plot_all_DLC_MAF005.py` helper, which provides `get_keys`, `scan_du`, `draw` and `log`, and defines the two category labels in `LABELS`. Obtain that helper from the study's original heatmap workflow and place it under `scripts/`. It is an outstanding dependency of this module. The input schemas for the DU site list and frequency table follow the helper's `scan_du()` implementation.

The wrapper uses the DLC maximum site MAF and the DU maximum MAF returned by the helper. Unavailable DU MAFs are displayed in gray. Outputs include the heatmap from the helper, site-level TSVs and `ASE_new_SourceData.xlsx`.

```bash
python scripts/plot_fig4f_novel_ase_maf_heatmap.py \
  --function-file scripts/plot_all_DLC_MAF005.py \
  --ase-status results/fig4e/05_DLC_sites_overlap_status.tsv \
  --du-list data/DU_site_list.tsv \
  --du-frequency data/DU_frequency.tsv \
  --outdir results/fig4f
```
