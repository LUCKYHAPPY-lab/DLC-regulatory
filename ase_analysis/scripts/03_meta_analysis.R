#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript 03_meta_analysis.R HET_MASKED_TSV SAMPLE_MAP_CSV ANNOTATION_TSV OUTDIR")
}
data <- read.delim(args[1], stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA")
sample_map <- read.csv(args[2], colClasses = "character", check.names = FALSE)
if (!"Ear_ID" %in% names(sample_map)) stop("Sample map must contain Ear_ID.")
ear_id <- sample_map$Ear_ID
info_columns <- paste0(ear_id, "_info")
if (!length(ear_id) || anyNA(ear_id) || anyDuplicated(ear_id) ||
    !all(c("SNP_id", info_columns) %in% names(data))) {
  stop("Check the sample map and the matching <Ear_ID>_info columns.")
}
data <- data[, c("SNP_id", info_columns), drop = FALSE]
annotation_file <- args[3]
out_dir <- args[4]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library("poolr")

binom<-function(hapone){
hapone<-as.character(hapone)
ref<-as.numeric(hapone[3])
all<-as.numeric(hapone[3])+as.numeric(hapone[4])
p<-binom.test(ref,all, p = 0.5)$p.value
return(p)
}


effect_ASE<-function(hapone){
hapone<-as.character(hapone)
SNP_id<-hapone[1]
df_info <- hapone[-1L]

df1<-data.frame(sam=ear_id,info=df_info)
df1$info<-as.character(df1$info)

df1<-df1[!is.na(df1$info),]
if (nrow(df1)>0){
df1$ref<-unlist(strsplit(df1$info,split="_"))[seq(1,nrow(df1)*3,3)]
df1$ref<-as.numeric(df1$ref)
df1$alt<-unlist(strsplit(df1$info,split="_"))[seq(2,nrow(df1)*3,3)]
df1$alt<-as.numeric(df1$alt)
df1$fdr<-unlist(strsplit(df1$info,split="_"))[seq(3,nrow(df1)*3,3)]
df1$fdr<-as.numeric(df1$fdr)
df1$ratio<-round(df1$ref/(df1$ref+df1$alt),3)
df1$p<-apply(df1,1,binom)

if(nrow(df1)>=19){
binomtest_p =poolr::binomtest(df1$p)$p
}else{
binomtest_p= "NA"
}

df2<-df1[df1$fdr<=0.05,]
num_all<-nrow(df1)
num_fdr<-nrow(df2)
trend_ref_num_all<-sum(df1$ratio>0.5)
trend_alt_num_all<-sum(df1$ratio<0.5)
median_FC_all<-median(df1$ratio)

median_Mratio<-median(c(df1$ratio[df1$ratio>=0.5],1-(df1$ratio[df1$ratio<0.5])))

trend_ref_num_fdr<-sum(df2$ratio>0.5)
trend_alt_num_fdr<-sum(df2$ratio<0.5)
median_FC_fdr<-median(df2$ratio)

aa<-c(SNP_id,num_all,num_fdr,trend_ref_num_all,trend_alt_num_all,median_FC_all,median_Mratio,trend_ref_num_fdr,trend_alt_num_fdr,median_FC_fdr,binomtest_p)

}else{
aa<-c(SNP_id,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
}

return(aa)
}

myresult<-apply(data, 1, effect_ASE)
myresult1<-data.frame(t(myresult))
colnames(myresult1)<-c("SNP_id","num_all","num_fdr","trend_ref_num_all","trend_alt_num_all","median_FC_all","median_Mratio","trend_ref_num_fdr","trend_alt_num_fdr","median_FC_fdr","binomtest_p")

write.table(myresult1,file.path(out_dir, "ase_meta_analysis_unannotated.tsv"), col.names=T,row.names=F,sep="\t",quote=F)


data1=read.table(annotation_file,sep="\t",header=T,stringsAsFactors = F)
# Retain the same annotation columns and apply BH after the left join.
annotation_columns <- c(1, 13, 19, 4, 5, 7, 22)
if (ncol(data1) < max(annotation_columns) ||
    !"Uploaded_variation" %in% names(data1)[annotation_columns]) {
  stop("The annotation table does not match the original column layout.")
}
kk <- merge(myresult1, data1[, annotation_columns], by.x = "SNP_id",
            by.y = "Uploaded_variation", all.x = TRUE)


kk$num_all<-as.numeric(as.character(kk$num_all))
kk$num_fdr<-as.numeric(as.character(kk$num_fdr))
kk$ratio<-round(kk$num_fdr/kk$num_all,3)
kk2<-kk

kk2$binomtest_p<-as.numeric(as.character(kk2$binomtest_p))
kk2$FDR_binomtest<-p.adjust(kk2$binomtest_p,method="fdr")

write.table(kk2,file.path(out_dir, "ase_meta_analysis_annotated.tsv"), col.names=T,row.names=F,sep="\t",quote=F)
