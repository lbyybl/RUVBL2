#--- this script is used Deseq2 to analysis the differential of Ruvbl's RNA-seq
#----deseq2
library(DESeq2)
library(dplyr)
library(gplots)
library(clusterProfiler)
library(org.Mm.eg.db)
rm(list=ls())
#---- read feature count data and calculate the TPM value
setwd("/WORK/lbyybl/zhn/RPB7/enrich_site/Proseq_count")#设置工作目录
dir()
# countdata<-read.table("genes_matrix.txt",sep="\t",header = F,skip=1,row.names = 1)
countdata <- read.table('ProSeqcounts2.txt',stringsAsFactors = F, row.names = 1, header = T, skip = 1)
countdata <- countdata[,c(1:5,30:33)]
colnames(countdata) <- c('chr','st','en','sd','length','RPB7_un_rep1','RPB7_un_rep2','RPB7_IAA_rep1','RPB7_IAA_rep2')
countdata <- as.data.frame(countdata)
# colnames(countdata) <- c('Control_rep1','Control_rep2','Experiment_rep1','Experiment_rep2')
head(countdata)

#------------------------------
metadata <- countdata[,1:5]#提取基因信息count数据前的几列
head(metadata)
countdata <- countdata[,6:ncol(countdata)]#提取counts数，counts数据主题部分

name<-"TPM"#设置输出文件前缀名
#-----TPM Calculation------

kb <- metadata$length / 1000

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

#--------------------
# countdata <- countdata[,c(5,6,1,2)]
# countdata_tpm <- countdata_tpm[,c(5,6,1,2)]
head(countdata_tpm)
pass <- rowSums(countdata_tpm > 1)
pass <- pass[pass > (ncol(countdata_tpm) / 2)]
countdata_tpm <- as.data.frame(subset(countdata_tpm, rownames(countdata_tpm) %in% names(pass)))
countdata <- as.data.frame(subset(countdata, rownames(countdata) %in% rownames(countdata_tpm)))
countdata_tpm[rownames(countdata_tpm)=='Mipepos',]
metadata[rownames(metadata)=='Mipepos',]
#---------------------------------------------
# #--- select expression gene
# coverage <- countdata[,c(5,6,7,10,11)]
# head(coverage)
# bg <- colSums(countdata[,c(6,7,10,11)])/sum(countdata$length)
# coverage$exp <- 'no'
# head(coverage)
# for (i in 1:nrow(coverage)){
#   for (j in 2:5){
#     exp <- bg[j-1]*coverage$length[i]
#     pvalue <- poisson.test(coverage[i,j],round(exp),conf.level = 0.95,alternative = 'greater')$p.value
#     if (pvalue <= 0.05){
#       coverage[i,6] <- 'yes'
#     }
#   }
# }
# 
# sum(coverage$exp=='yes')
# coverage <- coverage[coverage$exp=='yes',]
# countdata <- countdata[,c(10,11,6,7)]
# # countdata_tpm <- countdata_tpm[,1:4]
# head(countdata)
# countdata <- as.data.frame(subset(countdata, rownames(countdata) %in% rownames(coverage)))
#--------------------
# 
# # pass <- rowSums(countdata_tpm > 1)
# # pass <- pass[pass > (ncol(countdata_tpm) / 2)]
# # countdata_tpm <- as.data.frame(subset(countdata_tpm, rownames(countdata_tpm) %in% names(pass)))
# # countdata <- as.data.frame(subset(countdata, rownames(countdata) %in% names(pass)))
# 
# setwd('/DATA3/lbyybl/download/RNAseq/upf3a/RSEM/diff_result/')
# 
hour <- 'RPB7'
hour_control <- 'Untreated'
fold <- log2(1.5)
countdata_tpm2 <- countdata_tpm
mycounts <- countdata
# # names(countdata)
# 
merge_data <- countdata_tpm2
for (i in 1:ncol(merge_data)){
  merge_data[,i] <- log(merge_data[,i]+1)
}
# #------------------------------------------------------------


head(mycounts)
# colname <- colnames(mycounts)
# mycounts <- apply(mycounts,2,round)
# colnames(mycounts) <- colname

condition <- factor(c(rep(hour_control,2),rep(paste0("deg",hour),2)), levels = c(hour_control,paste0("deg",hour)))
condition
colData <- data.frame(row.names=colnames(mycounts), condition)

colData

head(mycounts)
#----构建dds对象，开始DEseq流程；
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)

dds <- estimateSizeFactors(dds)
sizeft <- c(125746,137318,133340,119436)
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
i <- which(rownames(mycounts) %in% diff_name)
mycounts[i,]
mycol <- colorpanel(1000,"blue","white","red")

library(pheatmap)
breaksList = seq(-3, 3, by = 1)
#pdf('AD38noDOXAD38withDOX_heatmap.pdf',width = 6,height = 12)
pheatmap(as.matrix(merge_data[i,]),color = mycol,clustering_method = 'ward',
         scale = 'row',show_rownames = F ,filename = paste0('heatmap',hour_control,'vs',hour,'.pdf'))

dev.new()

# #--- enrichmen
# # setwd('/DATA2/work/lbyybl/coorlaborate/YB/YB_RNA_seq/total/output/graph/38nodaoxvs38jiadox')
# keytypes(org.Hs.eg.db)
# down <- read.table('down.tsv',col.names = c('gene','baseMean','log2FoldChange',
#                                             'lfcSE','stat','pvalue','padj'))
# conv_id <- read.table('/DATA/work/lbyybl/genomes/hg19/anno/ensemble_genesymbol.txt',
#                       col.names = c('gene','symbol'))
# down_conv <- merge(down,conv_id,by='gene')
# fwrite(down_conv,'down2.tsv',sep = '\t')
# ego2 <- enrichGO(gene         = down_conv$symbol,
#                  OrgDb         = org.Hs.eg.db,
#                  keyType       = 'SYMBOL',
#                  ont           = "BP",
#                  pAdjustMethod = "BH",
#                  pvalueCutoff  = 1,
#                  qvalueCutoff  = 1)
# 
# pdf('down_BP_enrichment.pdf',width = 10,height = 5)
# dotplot(ego2, showCategory=30)
# dev.off()
# 
# up <- read.table('up.tsv',col.names = c('gene','baseMean','log2FoldChange',
#                                         'lfcSE','stat','pvalue','padj'))
# up_conv <- merge(up,conv_id,by='gene')
# fwrite(up_conv,'up2.tsv',sep = '\t')
# ego2 <- enrichGO(gene         = up_conv$symbol,
#                  OrgDb         = org.Hs.eg.db,
#                  keyType       = 'SYMBOL',
#                  ont           = "BP",
#                  pAdjustMethod = "BH",
#                  pvalueCutoff  = 1,
#                  qvalueCutoff  = 1)
# 
# pdf('up_BP_enrichment.pdf',width = 10,height = 5)
# dotplot(ego2, showCategory=30)
# dev.off()

