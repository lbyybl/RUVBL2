#--- this script is used Deseq2 to analysis the differential of Ruvbl's RNA-seq
#----deseq2
library(DESeq2)
library(dplyr)
library(gplots)
library(clusterProfiler)
library(org.Mm.eg.db)
rm(list=ls())
#---- read feature count data and calculate the TPM value
setwd("/WORK/lbyybl/WH/rvb/RNAseq/sample20200330/result/RSEM")#设置工作目录

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

setwd('/WORK/lbyybl/WH/rvb/RNAseq/sample20200330/result/RSEM/diff_result/')
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
  #----构建dds对象，开始DEseq流程；
  dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
  # keep <- rowSums(counts(dds)) >= 10
  # dds <- dds[keep,]
  dds <- DESeq(dds)
  dds
  #---总体结果查看
  
  res= results(dds)
  res = res[order(res$pvalue),]
  head(res)
  summary(res) 
  #所有结果先进行输出
  # setwd('diff_result')
  write.csv(mycounts,file=paste0('RAW_count',hour_control,'vs',hour,'.csv'))
  write.csv(res,file=paste0("All_results",hour_control,'vs',hour,".csv"))
  table(res$padj<0.05)
  summary(res)
  diff_gene_deseq2 <- subset(res, padj < 0.05 & abs(log2FoldChange) > fold)
  up_gene_deseq2 <- subset(res, padj < 0.05 & log2FoldChange > fold)
  down_gene_deseq2 <- subset(res, padj < 0.05 & log2FoldChange < -fold)
  data.table::fwrite(as.data.frame(down_gene_deseq2),paste0('down',hour_control,'vs',hour,'.tsv'),quote = F,sep = '\t',col.names = F,row.names = T)
  data.table::fwrite(as.data.frame(up_gene_deseq2),paste0('up',hour_control,'vs',hour,'.tsv'),quote = F,sep = '\t',col.names = F,row.names = T)
  
  
  diff_name <- c(rownames(up_gene_deseq2),rownames(down_gene_deseq2))
  
  
  dim(diff_gene_deseq2)
  dim(up_gene_deseq2)
  dim(down_gene_deseq2)
  i <- which(rownames(mycounts) %in% diff_name)
  mycounts[i,]
  mycol <- colorpanel(1000,"blue","white","red")
  
  library(pheatmap)
  breaksList = seq(-3, 3, by = 1)
  #pdf('AD38noDOXAD38withDOX_heatmap.pdf',width = 6,height = 12)
  pheatmap(as.matrix(merge_data[i,]),color = mycol,clustering_method = 'ward',
           scale = 'row',show_rownames = F ,filename = paste0('heatmap',hour_control,'vs',hour,'.pdf'))
  
}

#--- enrichmen
# setwd('/DATA2/work/lbyybl/coorlaborate/YB/YB_RNA_seq/total/output/graph/38nodaoxvs38jiadox')
keytypes(org.Hs.eg.db)
down <- read.table('down.tsv',col.names = c('gene','baseMean','log2FoldChange',
                                            'lfcSE','stat','pvalue','padj'))
conv_id <- read.table('/DATA/work/lbyybl/genomes/hg19/anno/ensemble_genesymbol.txt',
                      col.names = c('gene','symbol'))
down_conv <- merge(down,conv_id,by='gene')
fwrite(down_conv,'down2.tsv',sep = '\t')
ego2 <- enrichGO(gene         = down_conv$symbol,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1)

pdf('down_BP_enrichment.pdf',width = 10,height = 5)
dotplot(ego2, showCategory=30)
dev.off()

up <- read.table('up.tsv',col.names = c('gene','baseMean','log2FoldChange',
                                        'lfcSE','stat','pvalue','padj'))
up_conv <- merge(up,conv_id,by='gene')
fwrite(up_conv,'up2.tsv',sep = '\t')
ego2 <- enrichGO(gene         = up_conv$symbol,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1)

pdf('up_BP_enrichment.pdf',width = 10,height = 5)
dotplot(ego2, showCategory=30)
dev.off()

