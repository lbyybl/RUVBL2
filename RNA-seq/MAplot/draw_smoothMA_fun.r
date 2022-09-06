#data <- gene_test

draw_smooth_MA <- function(data,signi.type='padj',signi.cutoff=0.05,Fold.cutoff=log2(1.5),
                           pch=20,cex=.15,xlab='log concentration',ylab='log2 foldchange',xrange=10,yrange=3){
  if (signi.type=='padj'){
    up <- data[data$log2FoldChange >= Fold.cutoff & data$padj <= signi.cutoff, ]
    down <- data[data$log2FoldChange <= -Fold.cutoff & data$padj <= signi.cutoff, ]
    unchange <- data #%>% dplyr::filter(!(gene %in% c(up$gene,down$gene)))
  } else if(signi.type=='pvalue'){
    up <- data[data$log2FoldChange >= Fold.cutoff & data$pvalue <= signi.cutoff, ]
    down <- data[data$log2FoldChange <= -Fold.cutoff & data$pvalue <= signi.cutoff, ]
    unchange <- data # %>% dplyr::filter(!(gene %in% c(up$gene,down$gene)))
  } else{
    stop("the signi.type should be 'padj' or 'pvalue' !!!",call.=FALSE)
  }
  
  if (yrange!=3){
    max_fold <- max(yrange)
  }else{
    fold <- data$log2FoldChange[data$log2FoldChange>-200 & data$log2FoldChange < 200]
    max_fold <- round(as.numeric(quantile(abs(fold),1)))
  }
  
  if (xrange!=10){
    max <- max(xrange)
    min <- min(xrange)
  }else{
    basemean <- log2(data$baseMean[data$baseMean >-0 & data$baseMean < 10000000])
    max <- round(as.numeric(quantile(basemean,1)))
    min <- floor(min(as.numeric(basemean)))
  }
  
  
  smoothScatter(log2(unchange$baseMean),unchange$log2FoldChange,pch=20,cex=0.15,col='grey',
                xlim=c(min,max),
                xlab=xlab,
                ylim=c(-max_fold,max_fold),
                ylab=ylab)  
  points(log2(up$baseMean),up$log2FoldChange,pch=20,cex=cex,col='red')
  points(log2(down$baseMean),down$log2FoldChange,pch=20,cex=cex,col='blue')
  abline(h=0,col='dodgerblue')
  # browser()
  text(max*0.95,max_fold*0.95,paste0('n = ',nrow(up)))
  text(max*0.95,-max_fold*0.95,paste0('n = ',nrow(down)))
}
#data$log2FoldChange <- -data$log2FoldChange
#pdf('MA_smooth_plot.pdf',width = 5,height = 5)
#draw_smooth_MA(data,cex = 0.5,xrange = c(1,15.5),yrange = 4)
#dev.off()
