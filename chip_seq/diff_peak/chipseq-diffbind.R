suppressMessages(library(dplyr))
suppressMessages(library(DiffBind))
suppressMessages(library(optparse))

option_list <- list(
  make_option(c("-c", "--csv"), type = "character", default=NULL,
              help=".csv file read by DiffBind"),
  make_option(c("-o", "--outputdir"), type="character", default=".",
              help="output directory, default is current directory"),
  make_option(c("-p", "--prefix"), type="character", default="result",
              help="output prefix, default is 'result'"),
  make_option(c("-m", "--mergeoverlap"), type="numeric", default=1,
              help="mergeOverlap for dba config"),
  make_option(c("-g", "--isgene"), type="logical", default=FALSE, action="store_true",
              help="assign gene name or not, default is FALSE"),
  make_option(c("-q", "--fdrthreshold"), type="numeric", default=0.05,
              help="FDR threshold, default is 0.05"),
  make_option(c("-f", "--log2fc"), type="numeric", default=0,
              help="log2 fold change, default is 0, and this number will be log2 transform"),
  make_option(c("-u", "--usepvalue"), type="logical", default=FALSE, action="store_true",
              help="use p value instead of q value, threshold remains to be set by '-q'"),
  make_option(c("-l", "--xyplot"), type="logical", default=FALSE, action="store_true",
              help="plot MAplot in x-y mode, default is FALSE"),
  make_option(c("-x", "--xrange"), type="character", default=FALSE,
              help="xrange for MAplot, separate by ','"),
  make_option(c("-y", "--yrange"), type="character", default=FALSE,
              help="yrange for MAplot, separate by ','")
)

opt <- parse_args(OptionParser(option_list=option_list))
#browser()
mergeoverlap <- opt$mergeoverlap
isgene <- opt$isgene
fdrtr <-  opt$fdrthreshold
fctr <-  log2(opt$log2fc)
usepvalue <- opt$usepvalue
xyplot <-  opt$xyplot
xrange <- opt$xrange
yrange <- opt$yrange

# R version 4.0.3; DiffBind version 3.2.2
# 此脚本可能有参数需要调整
# 注意: 在3.2.2版本中wt和ko位置需要wt在前, 同时需要让bed文件的score>0

### analyze
message("Start analyzing by DiffBind ...")
# 步骤: dbs -> count -> contrast -> analyze -> report
# mergeOverlap: 只有overlap多少bp的才会被merge起来
dbObj <- dba(sampleSheet = opt$csv, config=data.frame(mergeOverlap=mergeoverlap, doGreylist=FALSE))
# dbObj <- dba(sampleSheet = opt$csv)
# 该步骤中summit设置两边的flank, 在3.2.2版本中需要关闭summits
# binding 数据中默认的是DBA_SCORE_TMM_MINUS_FULL, 此值会输出, 但最终的FC仍然会用raw重新建模
# 要求有overlap的地方都数一次, 此参数会被传入 summarizeOverlaps(mode=?) 中 [注意:源代码中此步出错, 并没有传入参数, counts.R LINE 610]
dbObj$config$intersectMode <- 'Union'
dbObj$config$mode <- 'Union'
dbObj$config$inter.feature <- FALSE
dbObj$config$singleEnd <- FALSE
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, summits=FALSE, mapQCth = 0, bParallel = TRUE,
                   score = 'DBA_SCORE_NORMALIZED')
#--- normalize
dbObj <- dba.normalize(dbObj)
# 此步中如果不加minMembers=2一个group只有2个rep会报错
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION, minMembers = 2)
# DBA_ALL_METHODS会算DESeq2以及EdgeR, 这里就只算一下DESeq2
# 如果有多个组进行比较, 不需要每个写一遍, 下面的步骤中写清楚contrast是几号
# 如果blacklist和greylist不是FALSE会报错
# dbObj <- dba.normalize()
dbObj <- dba.analyze(dbObj, method=DBA_DESEQ2, bBlacklist = FALSE, bGreylist = FALSE)

message("Assign contrasts ...")
for (i in 1:length(dbObj$contrasts)){
  contrastname<-paste0(dbObj$contrasts[[i]]$name1, "_vs_", dbObj$contrasts[[i]]$name2)
  outputname<-paste0(opt$outputdir, "/", opt$prefix, "-", contrastname)
  
  # th: 输出FDR <= th的结果, th=1表示全部输出
  dbObj.DB <- dba.report(dbObj, method=DBA_DESEQ2, contrast = i, th=1)
  
  ### write data
  deseq2res <- dbObj.DB %>% as.data.frame()
  deseq2res$reg <- "unchanged"
  deseq2res$reg[deseq2res$FDR<=fdrtr & abs(deseq2res$Fold)>=fctr] <- "DEG"
  
  #--- get the change region
  if (opt$yrange == FALSE){
	  out <- deseq2res
	  fold <- out$Fold[out$Fold>-200 & out$Fold < 200]
	  max_fold <- round(as.numeric(quantile(abs(fold),1)))
	  conc <- out$Conc[out$Conc>-0 & out$Conc < 10000000]
	  max_conc <- round(as.numeric(quantile(conc,1)))
	  min_conc <- floor(min(as.numeric(conc)))
	  yrange <-c(-max_fold,max_fold)
	}
  
  # 连接表达数据
  contrastgroup <- as.vector(dbObj$contrasts[[i]]$group1) | as.vector(dbObj$contrasts[[i]]$group2)
  contrastscore <- dbObj$binding[, -1:-3][rownames(deseq2res), contrastgroup]
  
  # 连接基因名
  if (isgene){
    configdf <- read.csv(opt$csv, header = T, stringsAsFactors = F)
    peakfile <- read.delim(configdf$Peaks[1], header = F, stringsAsFactors=F)
    colnames(peakfile) <- c("seqnames", "start", "end", "name", "score", "strand")
    
    output <- cbind(deseq2res, contrastscore) %>% .[,c(-4,-5)] %>% left_join(., peakfile, by = c("seqnames", "start", "end"), all.x = T) %>% .[,-16]
    write.table(output, paste0(outputname, "-diffbind.txt"), quote = F, col.names = T, row.names = F, sep = "\t")
  }
  else{
    output <- cbind(deseq2res, contrastscore)
    write.table(output, paste0(outputname, "-diffbind.txt"), quote = F, col.names = T, row.names = F, sep = "\t")
  }
  
  # source('/home/lbyybl//bashscript/bin/ruvb2/chip-seq/diff_peak/draw_smoothMA_fun.r')
  # output$Conc <- output$Conc + 1
  # output$FDR[is.na(output$FDR)] <- 1 
  # output$Fold[is.na(output$Fold)] <- 0 
  # output$Fold <- log2((2^output[,5]-1)/(2^output[,6]-1))
  # output$Fold <- output$Conc_INO80_IAA-output$Conc_INO80_untreated
  # output$fold <- log2(rowSums(output[,13:14])/rowSums(output[,11:12]))
  # output$fold2 <- log2((rowSums(output[,13:14])+1)/(rowSums(output[,11:12])+1))
  # draw_smooth_MA(output,xrange = c(1, 14))
  
  ### plot
  if (! usepvalue){
    # MAplot: FDR <= th; Fold >= fold
    pdf(paste0(outputname, "-diffbind-maplot.pdf"), width = 5, height = 5)
    if (xrange != FALSE & yrange != FALSE){
      if (grepl(',',paste(xrange,''))) {xrange <- strsplit(xrange, split = ",")[[1]] %>% as.numeric()}
      if (grepl(',',paste(yrange,''))) {yrange <- strsplit(yrange, split = ",")[[1]] %>% as.numeric()}
      dba.plotMA(dbObj, method = DBA_DESEQ2, contrast = i,bNormalized = T, 
                 bXY = F, bSmooth = T,  th = fdrtr, fold = fctr, bLoess = F, 
                 xrange = xrange, yrange = yrange)
    }else if (xrange != FALSE & yrange == FALSE){
      if (grepl(',',paste(xrange,''))) {xrange <- strsplit(xrange, split = ",")[[1]] %>% as.numeric()}
	  yrange <- c(-max_fold,max_fold)
      dba.plotMA(dbObj, method = DBA_DESEQ2, contrast = i, bXY = F, 
                 bSmooth = T, th = fdrtr, fold = fctr, bLoess = F, 
                 xrange = xrange, yrange = yrange)
    }else if (xrange == FALSE & yrange != FALSE){
      if (grepl(',',paste(yrange,''))) {yrange <- strsplit(yrange, split = ",")[[1]] %>% as.numeric()}
      dba.plotMA(dbObj, method = DBA_DESEQ2, contrast = i, bXY = F, 
                 bSmooth = T, th = fdrtr, fold = fctr, bLoess = F, 
                 yrange = yrange)
    }else{
	  yrange <- c(-max_fold,max_fold)
      dba.plotMA(dbObj, method = DBA_DESEQ2, contrast = i, bXY = F, 
                 bSmooth = T, th = fdrtr, fold = fctr, bLoess = F, 
                 yrange = yrange)
    }
    dev.off()
    
    if (xyplot){
      pdf(paste0(outputname, "-diffbind-maplot-xy.pdf"), width = 5, height = 5)
      dba.plotMA(dbObj, method = DBA_DESEQ2, contrast = i, bXY = T, bSmooth = T, th = fdrtr, fold = fctr, bLoess = F)
      dev.off()
    }
  }
  else{    
    # MAplot: p <= th; Fold >= fold
    pdf(paste0(outputname, "-diffbind-maplot.pdf"), width = 5, height = 5)
    if (xrange != FALSE & yrange != FALSE){
      if (grepl(',',paste(xrange,''))) {xrange <- strsplit(xrange, split = ",")[[1]] %>% as.numeric()}
      if (grepl(',',paste(yrange,''))) {yrange <- strsplit(yrange, split = ",")[[1]] %>% as.numeric()}
      dba.plotMA(dbObj, method = DBA_DESEQ2, contrast = i, bXY = F, 
                 bSmooth = T, bUsePval = T, th = fdrtr, fold = fctr, 
                 bLoess = F, xrange = xrange, yrange = yrange)
    }else if (xrange != FALSE & yrange == FALSE){
      if (grepl(',',paste(xrange,''))) {xrange <- strsplit(xrange, split = ",")[[1]] %>% as.numeric()}
	  yrange <- c(-max_fold,max_fold)
      dba.plotMA(dbObj, method = DBA_DESEQ2, contrast = i, bXY = F, 
                 bSmooth = T, bUsePval = T, th = fdrtr, fold = fctr, 
                 bLoess = F, xrange = xrange, yrange = yrange)
    }else if (xrange == FALSE & yrange != FALSE){
      if (grepl(',',paste(yrange,''))) {yrange <- strsplit(yrange, split = ",")[[1]] %>% as.numeric()}
      dba.plotMA(dbObj, method = DBA_DESEQ2, contrast = i, bXY = F, 
                 bSmooth = T, bUsePval = T, th = fdrtr, fold = fctr, 
                 bLoess = F, yrange = yrange)
    }else{
	  yrange <- c(-max_fold,max_fold)
      dba.plotMA(dbObj, method = DBA_DESEQ2, contrast = i, bXY = F, 
                 bSmooth = T, bUsePval = T, th = fdrtr, fold = fctr, 
                 bLoess = F, yrange = yrange)
    }
    dev.off()
    
    if (xyplot){
      pdf(paste0(outputname, "-diffbind-maplot-xy.pdf"), width = 5, height = 5)
      dba.plotMA(dbObj, method = DBA_DESEQ2, contrast = i, bXY = T, bSmooth = T, bUsePval = T, th = fdrtr, fold = fctr, bLoess = F)
      dev.off()
    }
  }
  
}



