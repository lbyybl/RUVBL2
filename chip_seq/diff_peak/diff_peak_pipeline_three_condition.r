#!/usr/bin/env R
# Thu Nov 14 11:35:07 2019
# Boyuan_Li


command=matrix(c(
	"dir", "d", "1", "character", 
	"sample_csv_file", "f", "1", "character" 
	),byrow=T,ncol=4)


args=getopt::getopt(command)


if (is.null(args$dir) || is.null(args$sample_csv_file)) {
	cat(paste(getopt::getopt(command, usage = T), "\n"))
	q()
}

#--QC for chip-seq
require(DiffBind)
library(BiocParallel)
library(DiffBind)
library(data.table)
library(dplyr)
register(SerialParam())
data_blacklist <- read.table('/DATA/work/lbyybl/wh/ruvb2/chip-seq/sample190312/local_result/R2-IAA-Pol2-C1_FKDL190724841-1a/bigwig/R2-IAA-Pol2-C1_blacklist')
data_blacklist <- GRanges(
  seqnames = Rle(data_blacklist$V1),
  ranges = IRanges(data_blacklist$V2,data_blacklist$V3)
)
#--- diffbind
#---一旦读入了peaksets，合并函数就找到所有重叠的peaks，并导出一致性的peaksets???
diffbind <- function(file,name){
  dbObj <- dba(sampleSheet=file)
  #---计算每个peaks/regions的count信息。先对一致性的peaks数据集进行标准化，然后根据他们的峰值（point of greatest read overlap???
  #---再次中心化并修剪一些peaks，最终得到更加标准的peak间隔。使用函数dba.count()，参数bUseSummarizeOverlaps可以得到更加标准的计算流程???
  dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE,filter = 10)
  #---差异分析
  # Establishing a contrast 
  dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION,minMembers = 2)
  dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
  comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
  comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
  
  # EdgeR
  # out <- as.data.frame(comp1.edgeR)
  # write.table(out, file=paste0("diff_peak/",name,"_edgeR.txt"), sep="\t", quote=F, col.names = NA)
  # DESeq2
  out <- as.data.frame(comp1.deseq)
  write.table(out, file=paste0("diff_peak/",name,"_deseq2.txt"), sep="\t", quote=F, col.names = NA)
  
  # Create bed files for each keeping only significant peaks (p < 0.05)
  # EdgeR
  # out <- as.data.frame(comp1.edgeR)
  # edge.bed <- out[ which(out$FDR < 0.05), 
                   # c("seqnames", "start", "end", "strand", "Fold")]
  # write.table(edge.bed, file=paste0("diff_peak/",name,"_edgeR_sig.bed"), sep="\t", quote=F, row.names=F, col.names=F)
  
  # DESeq2
  out <- as.data.frame(comp1.deseq)
  deseq.bed <- out[ which(out$FDR < 0.05), 
                    c("seqnames", "start", "end", "strand", "Fold")]
  write.table(deseq.bed, file=paste0("diff_peak/",name,"_deseq2_sig.bed"), sep="\t", quote=F, row.names=F, col.names=F)
  
  fold <- out$Fold[out$Fold>-200 & out$Fold < 200]
  max_fold <- round(as.numeric(quantile(abs(fold),1)))
  conc <- out$Conc[out$Conc>-0 & out$Conc < 10000000]
  max_conc <- round(as.numeric(quantile(conc,1)))
  min_conc <- floor(min(as.numeric(conc)))
  #setwd('/DATA/work/lbyybl/wh/ruvb2/chip-seq/Ruvbl1/sample20190509/graph/compare/diff_peak')
  pdf(paste0(name,'MA_plot.pdf'),width = 5,height = 5)
  dba.plotMA(dbObj,bXY=F,xrange = c(min_conc,max_conc),yrange = c(-max_fold,max_fold),bSmooth = T,fold = log2(2))
  dev.off()
  pdf(paste0(name,'scatter_plot.pdf'),width = 5,height = 5.5)
  dba.plotMA(dbObj,bXY=T,xrange = c(min_conc,max_conc),yrange = c(min_conc,max_conc),bSmooth = T,fold = log2(2))
  dev.off()
  # pdf(paste0(name,'volcano_plot.pdf'),width = 7,height = 5.5)
  # dba.plotVolcano(dbObj,bLabels = F)
  # dev.off()
  # pdf(paste0(name,'heatmap.pdf'))
  # dba.plotHeatmap(dbObj,contrast=1,correlations=FALSE)
  # dev.off()
  return(dbObj)
}
setwd(args$dir)
dir.create('diff_peak')
sample <- fread(args$sample_csv_file)
dox05h <- sample[c(3,4,1,2),]
fwrite(dox05h,'dox05hsample.csv',sep = ',')
R05h1 <- sample[c(5,6,3,4),]
fwrite(R05h1,'R05h1sample.csv',sep = ',')
Rdox1h <- sample[c(5:6,1:2),]
fwrite(Rdox1h,'Rdox1hsample.csv',sep = ',')

dox05h<-diffbind('dox05hsample.csv','H2O2_0vs2')
R05h1<-diffbind('R05h1sample.csv','H2O2_2vs10')
Rdox1h<-diffbind('Rdox1hsample.csv','H2O2_0vs10')


