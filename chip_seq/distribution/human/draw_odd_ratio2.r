#--- this script is used to draw the distribution of peaks and diffrential peak sites
setwd('/WORK/lbyybl/qff/chipseq_H3K56bhb/graph/distribution')
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

#--- insulator
IN <- fread('/WORK/lbyybl/qff/chipseq20210131/peak/broadPeak/enhancer//venn/distribution/insulator.bed',stringsAsFactors = F,
            col.names = c('chr','st','en','name','score'))

#--- convert dataframe to Granges
name1 <- ls()
for (i in name1[2:length(name1)]){
  data <- get0(i)
  assign(paste0(i,'_G'),
         makeGRangesFromDataFrame(data,ignore.strand = T,seqnames.field = 'chr',start.field = 'st',end.field = 'en'))
}

#---function overlapping peak with each element
peakovelement <- function(file){
  # browser()
  data <- fread(file,stringsAsFactors = F,
                select = c(1:3),col.names = c('chr','st','en'))
  
  data_G <- makeGRangesFromDataFrame(data,ignore.strand = T,seqnames.field = 'chr',start.field = 'st',end.field = 'en')
  
  data[,name1[2:length(name1)]] <- 0
  
  for (i in name1[2:length(name1)]){
    ele_G <- get0(paste0(i,'_G'))
    ov <- findOverlaps(data_G,ele_G) %>% as.data.frame()
    data[ov$queryHits,i] <- 1
  }
  data$ele1 <- 'Intergetic'
  ele1_name <- c('SE','TE','gene_promoter','IN','gene_GeneBody')
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
H3K56bhb <- '/WORK/lbyybl/qff/chipseq_H3K56bhb/peak/H3K56bhb_myc_idr_nocut.broadPeak'

H3K56bhb_ov <- peakovelement(H3K56bhb)

sts_data <- function(data){
  line <- nrow(data)
  sts <- data %>%
    dplyr::group_by(ele1) %>%
    summarise(n=n())
  sts$ratio <- sts$n/line
  return(as.data.frame(sts))
}

data_ratio <- sts_data(H3K56bhb_ov)

#--- calculate the ratio of each elements
for (i in name1[2:length(name1)]){
  data <- get0(i)
  data <- unique(data)
  length <- nrow(data)
  assign(paste0('length_',i),length)
}

len_name <- ls(pattern = 'length_')
sum_lenthg <- do.call(sum,lapply(len_name, get0))

list_len_name <- lapply(len_name, get0)
list_len_name <- lapply(list_len_name, function(x) x/sum_lenthg)

ele_ratio <- data.frame('ele1'=len_name,
                        'ratio'=unlist(list_len_name))
odds_ratio <- data_ratio[c(1,2,3,5,6),]
odds_ratio$ratio <- odds_ratio$ratio/ele_ratio$ratio

draw_bar <- function(data,name){
  ymax <- ceiling(max(data$ratio)/2)*2
  ggplot(data,aes(factor(ele1,levels = c('gene_promoter','TE','SE','gene_GeneBody','IN')),ratio))+
    geom_col(width = 0.4,color='black',fill = '#29ACE4') +
    scale_y_continuous(limits = c(0,ymax),
                       breaks = seq(0,ymax,by=1),
                       labels = as.character(seq(0,ymax,by=1)),
                       expand = c(0,0))+
    scale_x_discrete(#limits = c('Prom','Tyen','Supe','Geby','insu'),
      labels = c('promoter','typical enhancer','super enhancer','gene body','insulator'),
      expand = c(0.2,0.2)) +
    theme_bw()+theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = 'black'),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45,vjust = 1,hjust=1),
      legend.position = 'bottom',
      legend.direction = 'vertical',
      # plot.margin = unit(c(6, 6, 6, 6), "lines"),
      legend.title = element_blank()
    ) + ylab('Relative distribution\n\t(odds ratio)') +
    xlab('')
}

draw_bar(odds_ratio,'H3K56bhb')

ggsave('H3K56bhb_ST.pdf',width = 4, height = 3.5)

