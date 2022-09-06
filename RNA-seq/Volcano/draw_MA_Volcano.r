#--- this script is used to draw the MAplot of output by diffbind
rm(list = ls())
library(ggrastr,lib.loc = "/home/boyuanli/tools/anaconda2/envs/R40/lib/R/library")
source('/home/boyuanli/bashscript/bin/RNA-seq/MAplot/draw_smoothMA_fun.r')
setwd('/WORK/lbyybl/ypj/cross_talk/EUseq/sample20211107')
files <- system('ls diff_result/All_resultsR3*',intern = T)
name <- lapply(files, function(x){
  a <- basename(x)
  a <- gsub('All_results','',a)
  b <- gsub('.csv','',a)
  return(b)
}) %>% unlist()
data <- lapply(files, function(x){
  # browser()
  f <- read.csv(x,stringsAsFactors = F,header = T,row.names = 1)
  f$padj[is.na(f$padj)] <- 1
  f$log2FoldChange[is.na(f$log2FoldChange)] <- 0
  return(f)
})
names(data) <- name

head(data[[1]])

lapply(name, function(x){
  tmp <- data[[x]]
  pdf(paste0(x,'.smooth.pdf'),width = 5,height = 5)
  draw_smooth_MA(data = tmp,xrange = c(0,13),yrange = c(-5,5),Fold.cutoff = log2(1.5))
  dev.off()
})


fd <- log2(1.5)
xrange <- c(0,13)
xmin <- min(xrange)
xmax <- max(xrange)
yrange <- c(-5,5)
ymin <- min(yrange)
ymax <- max(yrange)
draw_MA <- function(data,name){
  ggplot(data,aes(log2(baseMean),log2FoldChange)) + rasterise(geom_point(
    color = ifelse(data$log2FoldChange > fd & data$padj <= 0.05, '#ca3e47',ifelse(
      data$log2FoldChange < -fd & data$padj <= 0.05, '#005792','grey'
    ))
  ),dpi=300) +scale_y_continuous(limits = yrange)+
    scale_x_continuous(limits = xrange)+
    annotate("text",x=xmax*0.9,y=ymax*0.9,label=paste0('n = ',sum(data$log2FoldChange > fd & data$padj <= 0.05))) + 
    annotate("text",x=xmax*0.9,y=-ymax*0.9,label=paste0('n = ',sum(data$log2FoldChange < -fd & data$padj <= 0.05))) + 
    geom_hline(yintercept = 0,linetype='dashed')+
    theme_bw()+
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank()
    ) + xlab('log2 mean reads') + ylab(paste0('log2(',name,')'))
  ggsave(paste0(name,'_MA.pdf'),width = 5,height = 5)
}
lapply(name, function(x){
  tmp <- data[[x]]
  draw_MA(tmp,x)
})

#---
xrange <- c(-5,5)
xmin <- min(xrange)
xmax <- max(xrange)
yrange <- c(0.001,20)
ymin <- min(yrange)
ymax <- max(yrange)
draw_valcano <- function(data,name){
  ggplot(data,aes(log2FoldChange,-log10(padj))) + rasterise(geom_point(
    color = ifelse(data$log2FoldChange > fd & data$padj <= 0.05, '#ca3e47',ifelse(
      data$log2FoldChange < -fd & data$padj <= 0.05, '#005792','grey'
    ))
  ),dpi=300) +scale_y_continuous(limits = yrange)+
    scale_x_continuous(limits = xrange)+
    geom_vline(xintercept = c(-fd,fd),linetype='dashed',
               color='darkslateblue',size=1)+
    geom_hline(yintercept = -log10(0.05),linetype='dashed',
               color='darkslateblue',size=1)+
    annotate("text",x=xmax*0.9,y=ymax*0.9,label=paste0('n = ',sum(data$log2FoldChange > fd & data$padj <= 0.05))) + 
    annotate("text",x=-xmax*0.9,y=ymax*0.9,label=paste0('n = ',sum(data$log2FoldChange < -fd & data$padj <= 0.05))) + 
    geom_hline(yintercept = 0,linetype='dashed')+
    theme_bw()+
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank()
    ) + xlab('log2 mean reads') + ylab(paste0('log2(',name,')'))
  ggsave(paste0(name,'_Valcano.pdf'),width = 5,height = 5)
}

lapply(name, function(x){
  tmp <- data[[x]]
  draw_valcano(tmp,x)
})
