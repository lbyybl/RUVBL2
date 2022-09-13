#--- this script is used Deseq2 to analysis the differential of Ruvbl's RNA-seq
#----deseq2
library(DESeq2)
library(dplyr)
library(gplots)
library(clusterProfiler)
library(org.Mm.eg.db)
rm(list=ls())
#---- read feature count data and calculate the TPM value
setwd("") # set work dir

countdata<-read.table("gene_matrix.txt",sep="\t",header = T,row.names = 1)
head(countdata)

countdata_tpm<-read.table("gene_matrix_tpm.txt",sep="\t",header = T,row.names = 1)
head(countdata_tpm)

#--------------------
sub_fun <- function(string){
  a <- sub('.genes.results','',string)
  # b <- sub('_rRNA.ex.bam','',a)
  c <- gsub('\\.','_',a)
  return(c)
  print(c)
}
sub_fun2 <- function(string){
  a <- sub('R2_','',string)
  # b <- sub('_rRNA.ex.bam','',a)
  c <- gsub('_1','',a)
  return(c)
  print(c)
}
for (i in 1:length(names(countdata_tpm))){
  names(countdata_tpm)[i] <- sub_fun(names(countdata_tpm)[i])
}

for (i in 1:length(names(countdata))){
  names(countdata)[i] <- sub_fun(names(countdata)[i])
}
# pass <- rowSums(countdata_tpm > 1)
# pass <- pass[pass > (ncol(countdata_tpm) / 2)]
# countdata_tpm <- as.data.frame(subset(countdata_tpm, rownames(countdata_tpm) %in% names(pass)))
# countdata <- as.data.frame(subset(countdata, rownames(countdata) %in% names(pass)))

comp <- t(combn(seq(1,15,2),2))

setwd('') # set dir
for (i in 1:dim(comp)[1]){
  col_one <- comp[i,]
  col <- c(col_one[1],col_one[1]+1,col_one[2],col_one[2]+1)
  # head(countdata[,col])
  names <- sub_fun2(names(countdata)[col_one])
  
  hour <- names[2]
  hour_control <- names[1]
  fold <- log2(1)
  countdata_tpm2 <- countdata_tpm[,col]
  mycounts <- countdata[,col]
  # names(countdata)
  
  merge_data <- countdata_tpm2
  for (i in 1:ncol(merge_data)){
    merge_data[,i] <- log(merge_data[,i]+1)
  }
  #------------------------------------------------------------
  
  
  head(mycounts)
  # colname <- colnames(mycounts)
  mycounts <- apply(mycounts,2,round)
  # colnames(mycounts) <- colname
  
  condition <- factor(c(rep(hour_control,2),rep(paste0("deg",hour),2)), levels = c(hour_control,paste0("deg",hour)))
  condition
  colData <- data.frame(row.names=colnames(mycounts), condition)
  
  colData
  # mycounts <- as.data.frame(mycounts)
  # rownames(mycounts) <- mycounts$gene
  # mycounts <- mycounts[,2:5]
  head(mycounts)
  #----build dds object and DEseq；
  dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
  dds <- DESeq(dds)
  dds
  
  res= results(dds)
  res = res[order(res$pvalue),]

  write.csv(mycounts,file=paste0('RAW_count',hour_control,'vs',hour,'.csv'))
  write.csv(res,file=paste0("All_results",hour_control,'vs',hour,".csv"))
  diff_gene_deseq2 <- subset(res, padj < 0.05 & abs(log2FoldChange) > fold)
  up_gene_deseq2 <- subset(res, padj < 0.05 & log2FoldChange > fold)
  down_gene_deseq2 <- subset(res, padj < 0.05 & log2FoldChange < -fold)
  data.table::fwrite(as.data.frame(down_gene_deseq2),paste0('down',hour_control,'vs',hour,'.tsv'),quote = F,sep = '\t',col.names = F,row.names = T)
  data.table::fwrite(as.data.frame(up_gene_deseq2),paste0('up',hour_control,'vs',hour,'.tsv'),quote = F,sep = '\t',col.names = F,row.names = T)
 }



