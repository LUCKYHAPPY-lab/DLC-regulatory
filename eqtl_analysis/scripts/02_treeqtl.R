#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: Rscript 02_treeqtl.R SNP_MAP GENE_MAP MATRIXEQTL_CIS OUTDIR"
if (identical(args, "--help")) {
  cat(usage, "\n")
  quit(save = "no")
}
if (length(args) != 4L) stop(usage, call. = FALSE)

outdir <- args[4]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

library(MatrixEQTL)
library(TreeQTL)


snps_location_file_name <- args[1]
gene_location_file_name <- args[2]

snp_map = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
gene_map = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

dist <- 1e6

# Run TreeQTL separately for local vs. distal analysis ------------------------
# Identify eSNPs which affect nearby genes ------------------------------------
# Get number of genes nearby to each SNP
n_tests_per_SNP <- get_n_tests_per_SNP(snp_map, gene_map, nearby = TRUE, dist = dist)
# Get list of eSNPs
eSNPs <- get_eSNPs(n_tests_per_SNP, args[3],method = "BH")
# Generate txt output file with full set of eAssociations
get_eAssociations(eSNPs, n_tests_per_SNP, args[3], file.path(outdir, "eAssoc_cis_by_snp.txt"),by_snp = TRUE)

# Identify eGenes regulated by nearby SNPs -------------------------------------
# Get number of SNPs nearby to each gene
n_tests_per_gene <- get_n_tests_per_gene(snp_map, gene_map, nearby = TRUE, dist = dist)
# Get list of eGenes
eGenes <- get_eGenes(n_tests_per_gene, args[3],method = "BH")
# Generate txt output file with full set of eAssociations
get_eAssociations(eGenes, n_tests_per_gene, args[3], file.path(outdir, "eAssoc_cis_by_gene.txt"),by_snp = FALSE)
