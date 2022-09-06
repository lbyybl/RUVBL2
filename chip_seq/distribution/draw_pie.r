#--- this script is used to draw the distribution of peaks and diffrential peak sites
setwd('/WORK/lbyybl/ypj/cross_talk/ChIPseq/graph/distribution')
rm(list=ls())
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(plotrix)
# txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#--- mRNA gene promoter
gene <- read.table('/DATA/work/lbyybl/genomes/mm10/genecode_ensemble_gene.bed',stringsAsFactors = F,
                   col.names = c('chr','st','en','ensmbl','symbol','strand'))
# protein_coding <- read.table('/DATA/work/lbyybl/genomes/mm10/gene_code_protein_coding.gene.txt',stringsAsFactors = F,
#                              col.names = c('symbol','ensmbl'))
# protein_coding <- gene[gene$ensmbl %in% protein_coding$ensmbl,]

#-- promoter
gene_promoter <- gene %>% dplyr::mutate(st=ifelse(strand=="+",st-2000,
                                                         en-2000))
gene_promoter$en <- gene_promoter$st + 4000

#-- gene body
gene_GeneBody <- gene %>% dplyr::filter(en-st>4000) %>% dplyr::mutate(st=st+2000,en=en-2000)

#--- TES 
gene_TES <- gene %>% dplyr::mutate(st=ifelse(strand=="+",en-2000,st-2000))
gene_TES$en <- gene_TES$st + 4000

#--- tRNA
tRNA <- fread('/DATA/work/lbyybl/genomes/mm10/tRNA_mm10.bed',stringsAsFactors = F,
              col.names = c('chr','st','en','id','score','strand'))
#--- lincRNA
lincRNA <- fread('/DATA/work/lbyybl/genomes/mm10/anno_from_JH/genecode_ensemble_lnkRNA.bed',stringsAsFactors = F,
                 col.names = c('chr','st','en','ensmbl','symbol','strand'))
#--- LINE
LINE <- fread('/DATA/work/lbyybl/genomes/mm10/anno_from_JH/LINE.bed',stringsAsFactors = F,
                   col.names = c('chr','st','en','id','type','strand'))
#--- SINE 
SINE <- fread('/DATA/work/lbyybl/genomes/mm10/anno_from_JH/SINE.bed',stringsAsFactors = F,
              col.names = c('chr','st','en','id','type','strand'))
#--- LTR
LTR <- fread('/DATA/work/lbyybl/genomes/mm10/anno_from_JH/LTR.bed',stringsAsFactors = F,
             col.names = c('chr','st','en','id','type','strand'))

#--- Super Enhancer
SE <- fread('/DATA/work/lbyybl/genomes/mm10/DNA_elements_made_by_Boyuan/enhancer/SE.bed',stringsAsFactors = F,
            col.names = c('chr','st','en','strand','element'))

#--- Typical Enhancer
TE <- fread('/DATA/work/lbyybl/genomes/mm10/DNA_elements_made_by_Boyuan/enhancer/TE.bed',stringsAsFactors = F,
            col.names = c('chr','st','en','strand','element'))
#--- convert dataframe to Granges
name1 <- ls()
for (i in name1[2:length(name1)]){
  data <- get0(i)
  assign(paste0(i,'_G'),
         makeGRangesFromDataFrame(data,ignore.strand = T,seqnames.field = 'chr',start.field = 'st',end.field = 'en'))
}

#---function overlapping peak with each element
peakovelement <- function(file){
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
  data$ele2 <- 'Unknown'
  ele1_name <- c('tRNA','SE','TE','gene_promoter','lincRNA','gene_TES','gene_GeneBody')
  ele2_name <- c('LINE','SINE','LTR')
  mu1 <- function(line){
    line <- line*c(length(line):1)
    return(line)
  }
  data <- as.data.frame(data)
  # match(ele1_name,names(data))
  data[,ele1_name] <- t(apply(data[,ele1_name],1,mu1))
  data[,ele2_name] <- t(apply(data[,ele2_name],1,mu1))
  for (i in 1:nrow(data)){
    loc1 <- which.max(data[i,ele1_name])
    loc2 <- which.max(data[i,ele2_name])
    if (data[i,ele1_name][loc1]>0){
      data$ele1[i] <- ele1_name[loc1]
    }
    if (data[i,ele2_name][loc2]>0 & data$ele1[i]=='Intergetic'){
      data$ele2[i] <- ele2_name[loc2]
    }else if (data$ele1[i]!='Intergetic') {
      data$ele2[i] <- 'Others'
    }
  }
  return(data)
}

diffpeakovelement <- function(data){
  # browser()
  data_G <- makeGRangesFromDataFrame(data,ignore.strand = T,seqnames.field = 'seqnames',start.field = 'start',end.field = 'end')
  
  data[,name1[2:length(name1)]] <- 0
  
  for (i in name1[2:length(name1)]){
    ele_G <- get0(paste0(i,'_G'))
    ov <- findOverlaps(data_G,ele_G) %>% as.data.frame()
    data[ov$queryHits,i] <- 1
  }
  data$ele1 <- 'Intergetic'
  data$ele2 <- 'Unknown'
  ele1_name <- c('tRNA','SE','TE','gene_promoter','lincRNA','gene_TES','gene_GeneBody')
  ele2_name <- c('LINE','SINE','LTR')
  mu1 <- function(line){
    line <- line*c(length(line):1)
    return(line)
  }
  data <- as.data.frame(data)
  # match(ele1_name,names(data))
  data[,ele1_name] <- t(apply(data[,ele1_name],1,mu1))
  data[,ele2_name] <- t(apply(data[,ele2_name],1,mu1))
  for (i in 1:nrow(data)){
    loc1 <- which.max(data[i,ele1_name])
    loc2 <- which.max(data[i,ele2_name])
    if (data[i,ele1_name][loc1]>0){
      data$ele1[i] <- ele1_name[loc1]
    }
    if (data[i,ele2_name][loc2]>0 & data$ele1[i]=='Intergetic'){
      data$ele2[i] <- ele2_name[loc2]
    }else if (data$ele1[i]!='Intergetic') {
      data$ele2[i] <- 'Others'
    }
  }
  return(data)
}
#--- draw_Pie plot for data
draw_multi_pie <- function(data,name){
  # browser()
  data_sts_ele1 <- data %>% dplyr::group_by(ele1) %>% dplyr::summarise(n=n(),per=n()/nrow(data)) %>% as.data.frame()
  names(data_sts_ele1)[1] <- 'ele'
  data_sts_ele1$type <- 'ele1'
  data2 <- data %>% dplyr::filter(!(ele2=='Others'))
  data_sts_ele2 <- data2 %>% dplyr::group_by(ele2) %>% dplyr::summarise(n=n(),per=n()/nrow(data2)) %>% as.data.frame()
  names(data_sts_ele2)[1] <- 'ele'
  if (nrow(data_sts_ele2)>0){
    data_sts_ele2$type <- 'ele2'
  }else{
    data_sts_ele2$type <- NULL
  }

  
  sts_df <- rbind(data_sts_ele1,data_sts_ele2)
  sts_df2 <- data.frame('ele'=c('gene_GeneBody','gene_promoter','gene_TES','Intergetic','lincRNA','SE','TE','tRNA','LINE','LTR','SINE','Unknown'),
                        'n'=1e-100,'per'=1e-100,'type'=c(rep('ele1',8),rep('ele2',4)))
  for (i in 1:nrow(sts_df2)){
    if (sts_df2$ele[i] %in% sts_df$ele){
      loc <- which(sts_df$ele==sts_df2$ele[i])
      sts_df2$n[i] <- sts_df$n[loc]
      sts_df2$per[i] <- sts_df$per[loc]
    }
  }

  sts_df <- sts_df2

  iniR=0.2
  colors=list(NO='white',gene_GeneBody='#EB6437',gene_promoter='#3182bd',gene_TES='#DEB21A',Intergetic='#999999',
              lincRNA='#E68A6F',SE='#F2E29C',TE='#99C0D4',tRNA='#B97198',LINE='#d95f0e',LTR='LightGreen',SINE='#11A5B0',Unknown='#F8DCCE')
  name = gsub("_"," ",name)
  #0 circle: blank
  pie(1, radius=iniR, init.angle=90, col=c('white'), border = NA, labels='',main = paste0(name,'\nn = ',nrow(data)))
  
  sf <- sts_df$per[4]
  if (sf > 1e-50){
    #3 circle: show genic:introns and Intergetic:not_near_genes | upstream
    floating.pie(0,0,c(sum(sts_df$per[1:3]), sts_df$per[9]*sf, sts_df$per[10]*sf, 
                       sts_df$per[11]*sf, sts_df$per[12]*sf,
                       sum(sts_df$per[5:8])),
                 radius=3.3*iniR, startpos=pi/2, col=as.character(colors[c('NO','LINE','LTR','SINE','Unknown','NO')]),border=NA)
  }

  
  floating.pie(0,0,c(1,100),
               radius=2.3*iniR, startpos=pi/2, 
               col=as.character(colors[c('NO','NO')]),border=NA)
  
  #4 circle: show genic:exons and Intergetic:downstream
  floating.pie(0,0,c(sts_df$per[1], sts_df$per[2], sts_df$per[3], 
                     sts_df$per[4],sts_df$per[5], sts_df$per[6], sts_df$per[7],sts_df$per[8]),
               radius=2.1*iniR, startpos=pi/2, 
               col=as.character(colors[c('gene_GeneBody','gene_promoter','gene_TES','Intergetic',
                                         'lincRNA','SE','TE','tRNA')]),border=NA)
  
  
  legend(0, 5*iniR, gsub("_"," ",paste0(sts_df$ele,' ',round(sts_df$per*100,2),'%')), col=as.character(colors[-1]), pch=19,bty='n', ncol=2)
}

#--- annotation peak
Polr1c_H2O2_0.2 <- '/WORK/lbyybl/ypj/cross_talk/ChIPseq/Plor1c_GFP/peak_withinput/Polr1c_02_GFP_overlap_peak_idr.narrowPeak'
Polr1c_H2O2_0 <- '/WORK/lbyybl/ypj/cross_talk/ChIPseq/Plor1c_GFP/peak_withinput/Polr1c_0_GFP_overlap_peak_idr.narrowPeak'
Polr1c_H2O2_2 <- '/WORK/lbyybl/ypj/cross_talk/ChIPseq/Plor1c_GFP/peak_withinput/Polr1c_10_GFP_overlap_peak_idr.narrowPeak'

Polr1c_degron_RPB3 <- '/WORK/lbyybl/ypj/cross_talk/ChIPseq/Polr1c_degron_Rpb3/peak_withinput/Polr1c_degron_Rpb3-Polr1c_1h_Rpb3_ChIP_overlap_peak_idr.narrowPeak'
Polr1c_untreated_RPB3 <- '/WORK/lbyybl/ypj/cross_talk/ChIPseq/Polr1c_degron_Rpb3/peak_withinput/Polr1c_degron_Rpb3-Polr1c_untreated_Rpb3_ChIP_overlap_peak_idr.narrowPeak'

RPB3_degron_Polr1c <- '/WORK/lbyybl/ypj/cross_talk/ChIPseq/Rpb3_degron_polr1c/peak_withinput/0113-Rpb3-degron-Polr1c-ChIP-Rpb3-IAA-1h_polr1c-ChIP_overlap_peak_idr.narrowPeak'
RPB3_untreated_Polr1c <- '/WORK/lbyybl/ypj/cross_talk/ChIPseq/Rpb3_degron_polr1c/peak_withinput/0113-Rpb3-degron-Polr1c-ChIP-Rpb3-untreated_polr1c-ChIP_overlap_peak_idr.narrowPeak'

RPB3_H2O2_0.2 <- '/WORK/lbyybl/ypj/cross_talk/ChIPseq/Rpb3_GFP1/peak_withinput/Rpb3-0-2mM-GFP_overlap_peak_idr.narrowPeak'
RPB3_H2O2_0 <- '/WORK/lbyybl/ypj/cross_talk/ChIPseq/Rpb3_GFP1/peak_withinput/Rpb3-0mM-GFP_overlap_peak_idr.narrowPeak'
RPB3_H2O2_2 <- '/WORK/lbyybl/ypj/cross_talk/ChIPseq/Rpb3_GFP1/peak_withinput/Rpb3-10mM-GFP_overlap_peak_idr.narrowPeak'

peaks <- c('Polr1c_H2O2_0.2','Polr1c_H2O2_0','Polr1c_H2O2_2',
           'Polr1c_degron_RPB3','Polr1c_untreated_RPB3',
           'RPB3_degron_Polr1c','RPB3_untreated_Polr1c',
           'RPB3_H2O2_0.2','RPB3_H2O2_0','RPB3_H2O2_2')
for (i in peaks){
  data <- peakovelement(get0(i))
  pdf(paste0(i,'.pdf'),width = 15,height = 10)
  draw_multi_pie(data,paste0(i,'mM'))
  dev.off()
}

#--- annotation diffpeak
Polr1c_0vs2 <- "/WORK/lbyybl/ypj/cross_talk/ChIPseq/Plor1c_GFP/diff_peak/H2O2_0vs10_deseq2.txt"  
Polr1c_0vs0.2 <- "/WORK/lbyybl/ypj/cross_talk/ChIPseq/Plor1c_GFP/diff_peak/H2O2_0vs2_deseq2.txt"  
Polr1c_0.2vs2 <- "/WORK/lbyybl/ypj/cross_talk/ChIPseq/Plor1c_GFP/diff_peak/H2O2_2vs10_deseq2.txt"   
Polr1c_degron_RPB3 <- "/WORK/lbyybl/ypj/cross_talk/ChIPseq/Polr1c_degron_Rpb3/diff_peak/Ruvbl_deseq2.txt"
RPB3_degron_Polr1c <- "/WORK/lbyybl/ypj/cross_talk/ChIPseq/Rpb3_degron_polr1c/diff_peak/Ruvbl_deseq2.txt"
RPB3_0vs2 <- "/WORK/lbyybl/ypj/cross_talk/ChIPseq/Rpb3_GFP1/diff_peak/H2O2_0vs10_deseq2.txt"
RPB3_0vs0.2 <- "/WORK/lbyybl/ypj/cross_talk/ChIPseq/Rpb3_GFP1/diff_peak/H2O2_0vs2_deseq2.txt" 
RPB3_0.2vs2 <- "/WORK/lbyybl/ypj/cross_talk/ChIPseq/Rpb3_GFP1/diff_peak/H2O2_2vs10_deseq2.txt" 

diff_peaks <- c('Polr1c_0vs2','Polr1c_0vs0.2','Polr1c_0.2vs2','Polr1c_degron_RPB3','RPB3_degron_Polr1c','RPB3_0vs2','RPB3_0vs0.2','RPB3_0.2vs2')
for (i in diff_peaks){
  # if (i=='Polr1c_0vs2') browser()
  data <- fread(get0(i),stringsAsFactors = F, header = T)
  data_up <- unique(data[data$Fold>1 & data$FDR < 0.05,c(2:12)]) %>% as.data.frame()
  data_down <- unique(data[data$Fold < -1 & data$FDR < 0.05,c(2:12)]) %>% as.data.frame()
  if (nrow(data_up)>0){
    data_up <- diffpeakovelement(data_up)
    fwrite(data_up,paste0('diff_peak/',i,'_up.txt'),sep = '\t')
    pdf(paste0('diff_peak/',i,'_up.pdf'),width = 15,height = 10)
    draw_multi_pie(data_up,paste0(i,'_UP'))
    dev.off()
  }
  if (nrow(data_down)>0){
    data_down <- diffpeakovelement(data_down)
    fwrite(data_down,paste0('diff_peak/',i,'_down.txt'),sep = '\t')
    pdf(paste0('diff_peak/',i,'_down.pdf'),width = 15,height = 10)
    draw_multi_pie(data_down,paste0(i,'_DOWN'))
    dev.off()
  }
}

