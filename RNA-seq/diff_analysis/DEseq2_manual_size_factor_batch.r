#--- this script is used Deseq2 to analysis the differential of Ruvbl's RNA-seq
#----deseq2
library(DESeq2)
library(dplyr)
library(gplots)
library(clusterProfiler)
library(org.Mm.eg.db)
rm(list=ls())
#---- read feature count data and calculate the TPM value
setwd("/WORK/lbyybl/ypj/cross_talk/EUseq/sample20211107")#设置工作目录
dir()
# countdata<-read.table("genes_matrix.txt",sep="\t",header = F,skip=1,row.names = 1)
countdata <- read.table('Euseq.txt',stringsAsFactors = F, row.names = 1, header = T, skip = 1)
countdata <- as.data.frame(countdata)
name.c <- lapply(names(countdata), function(x){
  a <- gsub('_rmdup_unique.ex.bam','',x)
  a <- gsub('rm_rRNA.','',a)
  return(a)
}) %>% unlist()
names(countdata) <- name.c
# colnames(countdata) <- c('Control_rep1','Control_rep2','Experiment_rep1','Experiment_rep2')
head(countdata)

#------------------------------
metadata <- countdata[,1:5]#提取基因信息count数据前的几列
head(metadata)
countdata <- countdata[,6:ncol(countdata)]#提取counts数，counts数据主题部分

name<-"TPM"#设置输出文件前缀名
#-----TPM Calculation------

kb <- metadata$Length / 1000

rpk <- countdata / kb

calcu_tpm <- function(data){
  total_data <- data %>%
    summarise_all(sum)
  for (i in 1:ncol(data)){
    data[,i] <- (data[,i])/total_data[,i]*1e6
  }
  # data <- log(data)
  return(data)
}
countdata_tpm <- calcu_tpm(rpk)

colSums(countdata_tpm)
fwrite(countdata_tpm,'gene_tmp_count.txt',sep = '\t',row.names = T,col.names = T)
setwd('diff_result')
#--------------------
countdata_bk <- countdata
countdata_tpm_bk <- countdata_tpm
sizeft_bk <- c(47183,55667,72497,65564,104945,110945,176441,184673)
sub_fun2 <- function(string){
  a <- sub('R2_','',string)
  # b <- sub('_rRNA.ex.bam','',a)
  c <- gsub('rep1','',a)
  return(c)
  print(c)
}


comp <- t(combn(seq(1,8,2),2))[c(1,6),]
for (i in 1:dim(comp)[1]){
  col_one <- comp[i,]
  col <- c(col_one[1],col_one[1]+1,col_one[2],col_one[2]+1)
  names <- sub_fun2(names(countdata_bk)[col_one])
  hour <- names[2]
  hour_control <- names[1]
  fold <- log2(2)
  countdata <- countdata_bk[,col]
  countdata_tpm <- countdata_tpm_bk[,col]
  # head(countdata_tpm)
  # pass <- rowSums(countdata_tpm > 1)
  # pass <- pass[pass > (ncol(countdata_tpm) / 2)]
  # countdata_tpm <- as.data.frame(subset(countdata_tpm, rownames(countdata_tpm) %in% names(pass)))
  # countdata <- as.data.frame(subset(countdata, rownames(countdata) %in% rownames(countdata_tpm)))
  countdata_tpm[rownames(countdata_tpm)=='Mipepos',]
  metadata[rownames(metadata)=='Mipepos',]
  #---------------------------------------------
  countdata_tpm2 <- countdata_tpm_bk[,col]
  mycounts <- countdata_bk[,col]
  # names(countdata)
  
  merge_data <- countdata_tpm2
  for (i in 1:ncol(merge_data)){
    merge_data[,i] <- log(merge_data[,i]+1)
  }
  
  #---------------------------------------------
  # mycounts <- countdata
  # head(mycounts)
  condition <- factor(c(rep(hour_control,2),rep(paste0("deg",hour),2)), levels = c(hour_control,paste0("deg",hour)))
  condition
  colData <- data.frame(row.names=colnames(mycounts), condition)
  
  colData
  
  #----构建dds对象，开始DEseq流程；
  dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
  
  dds <- estimateSizeFactors(dds)
  sizeft <- sizeft_bk[col]
  sizeft <- (1/sizeft)*sizeft[1]
  sizeFactors(dds) <- 1/sizeft
  mycounts2 <- mycounts
  for (i in 1:ncol(mycounts)){
    mycounts2[i] <- mycounts[i]*sizeft[i]
  }
  
  mycounts[rownames(mycounts)=='Myc',]
  mycounts2[rownames(mycounts2)=='Myc',]
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  
  dds
  #---总体结果查看
  vsd <- vst(dds, blind=FALSE)
  head(assay(vsd),3)
  mycounts3 <- counts(dds, normalized=T)
  mycounts3[rownames(mycounts3)=='Myc',]
  res= results(dds)
  res = res[order(res$pvalue),]
  head(res)
  summary(res) 
  #所有结果先进行输出
  #setwd('diff_result')
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
  d_i <- which(rownames(mycounts) %in% diff_name)
  # mycounts[i,]
  mycol <- colorpanel(1000,"blue","white","red")
  
  library(pheatmap)
  breaksList = seq(-3, 3, by = 1)
  #pdf('AD38noDOXAD38withDOX_heatmap.pdf',width = 6,height = 12)
  pheatmap(as.matrix(merge_data[d_i,]),color = mycol,clustering_method = 'ward',
           scale = 'row',show_rownames = F ,filename = paste0('heatmap',hour_control,'vs',hour,'.pdf'))
  
  dev.new()
}




