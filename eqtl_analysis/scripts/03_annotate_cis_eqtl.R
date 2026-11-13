#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
usage <- "Usage: Rscript 03_annotate_cis_eqtl.R PREPARED_CIS_TSV PCADD CHROMATIN_STATES OUTDIR"
if (identical(args, "--help")) {
  cat(usage, "\n")
  quit(save = "no")
}
if (length(args) != 4L) stop(usage, call. = FALSE)

outdir <- args[4]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Preserve the original coordinate-based pCADD and SNP-ID-based chromatin joins.
eQTL<-read.table(args[1],sep="\t",header=T,stringsAsFactors=F)

SNP_pCADD<-read.table(args[2],sep="\t",header=T,stringsAsFactors=F)
SNP_pCADD$eQTL_id2<-paste(SNP_pCADD$Chrom,SNP_pCADD$Position,sep=":")

SNP_encode<-read.table(args[3],sep="\t",header=T,stringsAsFactors=F)


eQTL1<-merge(eQTL,SNP_encode[,c("SNP_id","type")],by="SNP_id",all.x=TRUE)
eQTL1$eQTL_id2<-paste(eQTL1$CHROM,eQTL1$POS,sep=":")
eQTL2<-merge(eQTL1,SNP_pCADD[,3:6],by="eQTL_id2",all.x=TRUE)
eQTL2<-eQTL2[order(eQTL2$BBFDR,decreasing=F),]
eQTL2<-eQTL2[,c(2,4:22)]
eQTL2<-eQTL2[abs(eQTL2$r)>=0.3,]

write.table(eQTL2,file.path(outdir, "cis_eqtl_annotated.tsv"),sep="\t",quote=F,row.names=F)



write.table(eQTL2[abs(eQTL2$pos_to_TSS)<=2000,],file.path(outdir, "cis_eqtl_annotated_2kb.tsv"),sep="\t",quote=F,row.names=F)
