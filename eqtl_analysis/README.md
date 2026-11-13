# DLC muscle cis-eQTL analysis

Upstream analysis scripts for the Duroc × Luchuan F1 study. Install the R packages MatrixEQTL and TreeQTL. Run the commands from this directory and substitute the actual study input paths.

## 1. Map associations with MatrixEQTL

Provide tab-separated genotype dosage, rank-normalized gene expression and covariate matrices. The first column holds row IDs. The remaining columns must list the same sample IDs in the same order. Provide SNP and gene position tables with headers in MatrixEQTL format (SNP ID, chromosome, position and gene ID, chromosome, start, end, respectively).

```bash
Rscript scripts/01_matrixeqtl.R data/genotypes.tsv data/snp_map.tsv \
  data/rank_normalized_expression.tsv data/gene_map.tsv data/covariates.tsv results/mapping
```

The script uses the original linear model, 1 Mb cis distance, and nominal output thresholds of 0.01 for cis and 0.00001 for trans. It exports the original association tables, versions with signed partial correlation `r` and a Q–Q plot. Genotype filtering, expression rank normalization and covariate preparation are upstream.

## 2. Apply TreeQTL

Use the same coordinate maps and the cis association table from step 1 without the added `r` column.

```bash
Rscript scripts/02_treeqtl.R data/snp_map.tsv data/gene_map.tsv \
  results/mapping/out_cis_SNP_gene_modelLINEAR_Covariates_1e-2 results/treeqtl
```

This preserves the original BH calls and writes `eAssoc_cis_by_snp.txt` and `eAssoc_cis_by_gene.txt`.

## 3. Annotate cis-eQTLs and filter by effect size

The original annotation script reads a prepared cis-eQTL TSV. Prepare it by joining the gene-level TreeQTL associations to gene metadata, SNP coordinates and the step 1 partial correlations by SNP–gene pair. Use the TreeQTL BBFDR < 0.05 study set. The original strand-aware TSS distance is POS − Start for forward-strand genes and End − POS for reverse-strand genes. These preparation operations appeared as commented blocks in the supplied script. This module starts step 3 from the resulting table.

```bash
Rscript scripts/03_annotate_cis_eqtl.R data/prepared_cis_eqtl_with_r.tsv \
  data/pcadd.tsv data/chromatin_states.tsv results/annotation
```

Keep the original column order for the prepared TSV and pCADD table because the supplied script uses positional column selections. The prepared TSV includes `SNP_id`, `CHROM`, `POS`, `BBFDR`, `r` and `pos_to_TSS`. The pCADD table uses `Chrom` and `Position`, with the original annotation fields in columns 3–5. The chromatin-state table uses `SNP_id` and `type`. The script retains the original position-based pCADD join and SNP-ID-based chromatin-state join. It writes `cis_eqtl_annotated.tsv` for |r| ≥ 0.3 and a separate `cis_eqtl_annotated_2kb.tsv` for the subset within 2 kb of the TSS.

The Fig. 3 plotting module uses additional processed inputs, including the lead-eQTL scatterplot table and `eQTL_MAF005_kept.tsv`. Their generation is outside the three supplied scripts.
