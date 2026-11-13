#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: Rscript 01_matrixeqtl.R GENOTYPES SNP_MAP EXPRESSION GENE_MAP COVARIATES OUTDIR"
if (identical(args, "--help")) {
  cat(usage, "\n")
  quit(save = "no")
}
if (length(args) != 6L) stop(usage, call. = FALSE)

outdir <- args[6]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# MatrixEQTL expects the same sample order in all three matrices.
sample_header <- function(path) {
  strsplit(readLines(path, n = 1L), "\t", fixed = TRUE)[[1]][-1L]
}
samples <- sample_header(args[1])
if (!length(samples) || anyDuplicated(samples) ||
    !identical(samples, sample_header(args[3])) ||
    !identical(samples, sample_header(args[5]))) {
  stop("Genotype, expression and covariate headers must list the same unique samples in the same order.")
}

library(MatrixEQTL)


# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name <- args[1]
snps_location_file_name <- args[2]

# Gene expression file name
expression_file_name <- args[3]
gene_location_file_name <- args[4]

# Covariates file name
# Set to character() for no covariates
#covariates_file_name = character();
covariates_file_name <- args[5]



# Output file name
output_file_name_cis = file.path(outdir, "out_cis_SNP_gene_modelLINEAR_Covariates_1e-2");
output_file_name_tra = file.path(outdir, "out_tra_SNP_gene_modelLINEAR_Covariates_1e-5");

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-2;
pvOutputThreshold_tra = 1e-5;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
snps = snps, 
gene = gene, 
cvrt = cvrt,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel, 
errorCovariance = errorCovariance, 
verbose = TRUE, 
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos, 
genepos = genepos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

#unlink(output_file_name_tra);
#unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values
pdf(file.path(outdir, "qqplot_linear_1e-2_1e-5.pdf"))
plot(me)
dev.off()




dfFull = me$param$dfFull;
##cis
tstat = me$cis$eqtls$statistic;
r = tstat / sqrt( dfFull + tstat^2 );
cis_result<-me$cis$eqtls
cis_result$r<-r
##trans
tstat = me$trans$eqtls$statistic;
r = tstat / sqrt( dfFull + tstat^2 );
trans_result<-me$trans$eqtls
trans_result$r<-r

write.table(cis_result,file.path(outdir, "out_cis_SNP_gene_modelLINEAR_Covariates_1e-2_add_R"),quote=FALSE,sep="\t",row.names=F)
write.table(trans_result,file.path(outdir, "out_tra_SNP_gene_modelLINEAR_Covariates_1e-5_add_R"),quote=FALSE,sep="\t",row.names=F)
