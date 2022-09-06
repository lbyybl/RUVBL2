#---- this script is used to draw the H3K56bhb overlap map with BRD4/MED1/H3K27ac peaks
setwd('/WORK/lbyybl/qff/chipseq20210131/peak/broadPeak/enhancer/venn')
library(ggsignif)
rm(list=ls())
library(ChIPpeakAnno)
H3K27ac <- list('/WORK/lbyybl/qff/chipseq20210131/peak/Nabhb_H3K27ac_idr.broadPeak',
                '/WORK/lbyybl/qff/chipseq20210131/peak/NC_H3K27ac_idr.broadPeak')
BRD4 <- list('/WORK/lbyybl/qff/chipseq20210419/peak/Nabhb_BRD4_idr_nocut.broadPeak',
             '/WORK/lbyybl/qff/chipseq20210419/peak/NC_BRD4_idr_nocut.broadPeak')
MED1 <- list('/WORK/lbyybl/qff/chipseq20210329/peak/Nabhb_MED1_idr_nocut.broadPeak',
             '/WORK/lbyybl/qff/chipseq20210329/peak/NC_MED1_idr_nocut.broadPeak')

H3K56bhb_H3K27ac <- list('/WORK/lbyybl/qff/chipseq_H3K56bhb/peak/FQ_13_106_H3K56bhb_H3K27ac_idr_nocut.broadPeak',
                         '/WORK/lbyybl/qff/chipseq_H3K56bhb/peak/FQ_13_135_WT_H3K27ac_idr_nocut.broadPeak')
H3K56bhb_BRD4 <- list('/WORK/lbyybl/qff/chipseq_H3K56bhb/peak/FQ_13_106_H3K56bhb_BRD4_idr_nocut.broadPeak',
                      '/WORK/lbyybl/qff/chipseq_H3K56bhb/peak/FQ_13_130_WT_BRD4_idr_nocut.broadPeak')
H3K56bhb_MED1 <- list('/WORK/lbyybl/qff/chipseq_H3K56bhb/peak/FQ_13_106_H3K56bhb_MED1_1_idr_nocut.broadPeak',
                      '/WORK/lbyybl/qff/chipseq_H3K56bhb/peak/FQ_13_130_WT_MED1_idr_nocut.broadPeak')

sample_names2 <- ls()
sample_names1 <- c('H3K27ac','BRD4','MED1')
#--- give name for list elements
give_name_forele <- function(list){
  for (i in 1:length(list)){
    name <- basename(list[[i]])
    pref <- gsub('_idr.*.broadPeak','',name)
    pref <- gsub('_H3K27ac_2_rmdup_uniqe_sample_peaks.broadPeak','',pref)
    pref <- gsub('FQ_.*_H3K56bhb','H3K56bhb',pref)
    pref <- gsub('FQ_.*_WT','WT',pref)
    pref <- gsub('H3K56bhb.*','H3K56bhb',pref)
    pref <- gsub('WT.*','WT',pref)
    names(list)[i] <- pref
  }
  return(list)
}
for (i in sample_names2){
  data <- get(i)
  data <- give_name_forele(data)
  assign(i,data)
}

#--- read_data 
read_data <- function(file){
  data <- fread(file,stringsAsFactors = F,header = F)
  data <- as.data.frame(data)
  data <- data[,1:3]
  names(data) <- c('chr','st','en')
  data <- makeGRangesFromDataFrame(data,start.field = 'st',end.field = 'en')
  return(data)
}

for (i in sample_names2){
  file <- get(i)
  data <- lapply(file, read_data)
  assign(i,data)
}
#--- overlap
for (i in sample_names1){
  file <- get(i)
  data1 <-file[[1]]
  data2 <- file[[2]]
  assign(paste0(i,'_ov'), findOverlapsOfPeaks(data1,data2))
  # pdf(paste0(i,'overlap.pdf'))
  # makeVennDiagram(ov,fill=c("#F4C6C3", "#E1E1E2"),
  #                 col=c("#BDBCBD", "#BDBCBD"))
  # dev.off()
}

#--- overlap with H3K56bhb
length(H3K27ac_ov$peaklist$data1)

length(unique(findOverlaps(H3K27ac$Nabhb_H3K27ac,H3K27ac$NC_H3K27ac)))
data <- H3K27ac_ov
# data2 <- H3K56bhb
get_overlap <- function(data,data2=H3K56bhb){
  nabhb_ovp <- length(unique(findOverlaps(data2[[1]],data$peaklist$data1,select = 'first')))
  nabhb_oth <- length(data2[[1]])-nabhb_ovp
  
  wt_ovp <- length(unique(findOverlaps(data2[[2]],data$peaklist$data1,select = 'first')))
  wt_oth <- length(data2[[2]])-wt_ovp
  
  df <- data.frame('sample'=c(rep(names(data2)[1],2),rep(names(data2)[2],2)),
                   'num' = c(nabhb_ovp,nabhb_oth,wt_ovp,wt_oth),
                   'overratio'=c(nabhb_ovp/length(data2[[1]]),nabhb_oth/length(data2[[1]]),wt_ovp/length(data2[[2]]),wt_oth/length(data2[[2]]))*100,
                   'class' = c(rep(c('overlap','nooverlap'),2)))
}

K27ac_ovK56 <- get_overlap(H3K27ac_ov,H3K56bhb_H3K27ac)
BRD4_ovK56 <- get_overlap(BRD4_ov,H3K56bhb_BRD4)
MED1_ovK56 <- get_overlap(MED1_ov,H3K56bhb_MED1)

draw_bar <- function(data,name,ymax){
  pvalue <- fisher.test(matrix(data$num,nrow = 2))$p.value
  ggplot(data,aes(sample,num,fill=class, label = num))+
    geom_col(position = position_dodge2(width = 0.4,reverse = T),width = 0.4,color='black') +
    geom_text(position = position_dodge2(0.4,reverse = T),vjust=-0.2)+
    geom_signif(comparisons = list(c("H3K56bhb", "WT")), 
                map_signif_level=function(p)sprintf('p = %.2g',pvalue)) + 
    scale_y_continuous(limits = c(0,ymax),
                       breaks = seq(0,ymax,by=(ymax/5)),
                       labels = as.character(seq(0,ymax,by=(ymax/5))),
                       expand = c(0,0))+
    scale_x_discrete(expand = c(0.2,0.2)) +
    scale_fill_manual(values = c('#C4C5C7','#29ACE4'),
                      labels=c(paste0('non-overlap with Nabhb-',name,' unique peak'),
                               paste0('overlap with Nabhb-',name,' unique peak')))+
    # scale_fill_discrete() +
    theme_bw()+theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = 'black'),
      axis.text = element_text(size = 12),
      legend.position = 'bottom',
      legend.direction = 'vertical',
      # plot.margin = unit(c(6, 6, 6, 6), "lines"),
      legend.title = element_blank()
    ) + ylab('overlap percentage (%)') +
    xlab('') #+ coord_cartesian(clip = "off")
}

draw_bar(K27ac_ovK56,'H3K27ac',85000)
ggsave('H3K27ac_overlap_H3K56bhb_num.pdf',width = 4, height = 5)
draw_bar(BRD4_ovK56,'BRD4',15000)
ggsave('BRD4_overlap_H3K56bhb_num.pdf',width = 4, height = 5)
draw_bar(MED1_ovK56,'MED1',12500)
ggsave('MED1_overlap_H3K56bhb_num.pdf',width = 4, height = 5)

pvalue <- fisher.test(matrix(BRD4_ovK56$num,nrow = 2))$p.value
pvalue2 <- fisher.test(matrix(MED1_ovK56$num,nrow = 2))$p.value

# pvalue <- fisher.test(matrix(K27ac_ovK56$num,nrow = 2))$p.value
# ggplot(K27ac_ovK56,aes(sample,overratio,fill=class, label = num))+
#   geom_col(position = position_dodge2(width = 0.4,reverse = T),width = 0.4,color='black') +
#   geom_text(position = position_dodge2(0.4,reverse = T),vjust=-0.2)+
#   geom_signif(comparisons = list(c("H3K56bhb_PolII", "WT_PolII")), 
#               map_signif_level=function(p)sprintf('p = %.2g',pvalue)) + 
#   scale_y_continuous(limits = c(0,100),
#                      breaks = seq(0,100,by=10),
#                      labels = as.character(seq(0,100,by=10)),
#                      expand = c(0,0))+
#   scale_x_discrete(labels=c('H3K56bhb','WT'),expand = c(0.2,0.2)) +
#   scale_fill_manual(values = c('#29ACE4','#C4C5C7'),
#                     labels=c('non-overlap with Nabhb-H3K27ac unique peak',
#                              'overlap with Nabhb-H3K27ac unique peak'))+
#   # scale_fill_discrete() +
#   theme_bw()+theme(
#     panel.grid = element_blank(),
#     panel.background = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = 'black'),
#     axis.text = element_text(size = 12),
#     legend.position = 'bottom',
#     legend.direction = 'vertical',
#     legend.title = element_blank()
#   ) + ylab('overlap percentage (%)') +
#   xlab('') 



# ggplot(K27ac_ovK56,aes(sample,num,fill=class))+geom_col(position = 'dodge2',width = 0.4) +
#   # scale_fill_discrete(limits = c('overlap','nonoverlap'))+
#   scale_y_continuous(limits = c(0,ceiling(max(K27ac_ovK56$num)/1000)*1000),
#                      breaks = seq(0,ceiling(max(K27ac_ovK56$num)/1000)*1000,length.out=5),
#                      labels = as.character(seq(0,ceiling(max(K27ac_ovK56$num)/1000)*1000,length.out=5)),
#                      expand = c(0,0))+
#   theme_bw()+theme(
#     panel.grid = element_blank(),
#     panel.background = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = 'black')
#   )
