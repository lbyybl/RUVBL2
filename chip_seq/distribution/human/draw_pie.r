#--- this script is used to draw the distribution of peaks and diffrential peak sites
setwd('/WORK/lbyybl/qff/ATACseq20220322/graph/distribution/')
rm(list=ls())
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(ChIPseeker)
# txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#--- mRNA gene promoter
gene <- read.table('/WORK/lbyybl/qff/chipseq_H3K56bhb/graph/distribution/genecode_hg19_gene.bed',stringsAsFactors = F,
                   col.names = c('chr','st','en','ensmbl','symbol','strand'))
# protein_coding <- read.table('/DATA/work/lbyybl/genomes/mm10/gene_code_protein_coding.gene.txt',stringsAsFactors = F,
#                              col.names = c('symbol','ensmbl'))
# protein_coding <- gene[gene$ensmbl %in% protein_coding$ensmbl,]

#-- promoter
gene_promoter <- gene %>% dplyr::mutate(st=ifelse(strand=="+",st-2000,
                                                         en-2000))
gene_promoter$en <- gene_promoter$st + 4000

#-- gene body
gene_GeneBody <- gene %>% dplyr::filter(en-st>2000) %>% dplyr::mutate(st=ifelse(strand=="+",st+2000,
                                                                                st),
                                                                      en=ifelse(strand=="+",en,
                                                                                en-2000))

#--- Super Enhancer
SE <- fread('/WORK/lbyybl/qff/chipseq_H3K56bhb/graph/distribution/ele/SE.bed',stringsAsFactors = F,
            col.names = c('chr','st','en'))

#--- Typical Enhancer
TE <- fread('/WORK/lbyybl/qff/chipseq_H3K56bhb/graph/distribution/ele/TE.bed',stringsAsFactors = F,
            col.names = c('chr','st','en'))

enhancer <- rbind(SE,TE)
#--- insulator
IN <- fread('/WORK/lbyybl/qff/chipseq20210131/peak/broadPeak/enhancer//venn/distribution/insulator.bed',stringsAsFactors = F,
            col.names = c('chr','st','en','name','score'))

#--- convert dataframe to Granges
name1 <- c('gene_promoter','gene_GeneBody','enhancer','IN')
for (i in name1[1:length(name1)]){
  data <- get0(i)
  assign(paste0(i,'_G'),
         makeGRangesFromDataFrame(data,ignore.strand = T,seqnames.field = 'chr',start.field = 'st',end.field = 'en'))
}

#---function overlapping peak with each element
peakovelement <- function(data){
  # browser()
  data <- data[,1:3]
  colnames(data) <- c('chr','st','en')
  
  data_G <- makeGRangesFromDataFrame(data,ignore.strand = T,seqnames.field = 'chr',start.field = 'st',end.field = 'en')
  
  data[,name1[1:length(name1)]] <- 0
  
  for (i in name1[1:length(name1)]){
    ele_G <- get0(paste0(i,'_G'))
    ov <- findOverlaps(data_G,ele_G) %>% as.data.frame()
    data[ov$queryHits,i] <- 1
  }
  data$ele1 <- 'Intergetic'
  ele1_name <- c('enhancer','gene_promoter','IN','gene_GeneBody')
  mu1 <- function(line){
    line <- line*c(length(line):1)
    return(line)
  }
  data <- as.data.frame(data)
  # match(ele1_name,names(data))
  data[,ele1_name] <- t(apply(data[,ele1_name],1,mu1))
  # data[,ele2_name] <- t(apply(data[,ele2_name],1,mu1))
  for (i in 1:nrow(data)){
    loc1 <- which.max(data[i,ele1_name])
    if (data[i,ele1_name][loc1]>0){
      data$ele1[i] <- ele1_name[loc1]
    }
  }
  return(data)
}


#--- annotation peak
# H3K9bhb <- '/DATA2/work/lbyybl/mouse_liver/H3K9bhb/H3K9bhb_ST_GSM1704768_peak.bed'

sts_data <- function(data){
  line <- nrow(data)
  sts <- data %>%
    dplyr::group_by(ele1) %>%
    summarise(n=n())
  sts$ratio <- sts$n/line
  return(as.data.frame(sts))
}
draw_pie <- function(data_ratio,name){
  order <- c('enhancer','gene_promoter','gene_GeneBody','IN','Intergetic')
  mc <- match(order,data_ratio$ele1)
  data_ratio <- data_ratio[mc,]
  
  pie(data_ratio$ratio,labels = paste0(round(data_ratio$ratio*100,2),'%'),radius=0.7,
      cex=0.8,col=c('#1DB8D3','#F47A45','#C0C1C3','#F4B432','#4171B3'),init.angle=90,clockwise=T,
      border="white")
  lenghd <- c('enhancer','promoter','gene body','insulator','others')
  legend("bottomright",legend=lenghd,cex=0.6,bty="n",
         fill=c('#1DB8D3','#F47A45','#C0C1C3','#F4B432','#4171B3'))
  title(paste0(name,'\nn=',sum(data_ratio$n)))
}

peaks <- system('ls /WORK/lbyybl/qff/ATACseq20220322/peak/*_idr.narrowPeak',intern = T)
name <- lapply(peaks,function(x){
  a <- basename(x)
  a <- gsub('_idr.narrowPeak','',a)
  return(a)
}) %>% unlist()
names(peaks) <- name
lapply(name,function(x){
  tmp <- fread(peaks[[x]],stringsAsFactors = F,
                select = c(1:3),col.names = c('chr','st','en'))
  tmp <- peakovelement(tmp)
  tmp.ratio <- sts_data(tmp)
  pdf(paste0(x,'_pie.pdf'))
  draw_pie(tmp.ratio,x)
  dev.off()
})

#--- for diff peak result

difffiles <- system('ls /WORK/lbyybl/qff/ATACseq20220322/diff_peak/*WT-diffbind.txt',intern = T)
fold <- 0
name <- lapply(difffiles, function(x){
  a <- basename(x)
  a <- gsub('-diffbind.txt','',a)
}) %>% unlist()

difffiles <- lapply(name, function(x){
  data <- fread(difffiles[[x]],stringsAsFactors = F)
  if (names(data)[1]=='seqnames'){
    data <- data[,1:ncol(data)]
  }else{
    data <- data[,2:ncol(data)]
  }
  data$Fold[is.na(data$Fold)] <- 0
  data$FDR[is.na(data$FDR)] <- 1
  return(data)
})
names(difffiles) <- name
lapply(name,function(x){
  # browser()
  tmp <- difffiles[[x]]
  tmp.up <- tmp[tmp$Fold > fold & tmp$FDR < 0.05,]
  tmp.down <- tmp[tmp$Fold < -fold & tmp$FDR < 0.05,]
  tmp.up <- peakovelement(tmp.up)
  tmp.up.ratio <- sts_data(tmp.up)
  pdf(paste0(x,'_pie_up.pdf'))
  draw_pie(tmp.up.ratio,x)
  dev.off()
  tmp.down <- peakovelement(tmp.down)
  tmp.down.ratio <- sts_data(tmp.down)
  pdf(paste0(x,'_pie_down.pdf'))
  draw_pie(tmp.down.ratio,x)
  dev.off()
})


