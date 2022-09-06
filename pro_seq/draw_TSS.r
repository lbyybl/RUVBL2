#--- this script is used to draw the meta of pro-seq TSS
setwd('/WORK/lbyybl/ypj/cross_talk/EUseq/graph/heatmap')
rm(list = ls())
name <- 'TSS'
sample_num <- 5
sample_name <- c('r1c_1h','r1c_E','r1c_U','r1c_W','RPB1_Untreated')
step <- 600
data <- fread('genecode_ensemble_gene_peak.gz',stringsAsFactors = F,skip = 1)
bk <- c('ENSMUSG00000092612.2','ENSMUSG00000101756.1','ENSMUSG00000064347.1','ENSMUSG00000064348.1','ENSMUSG00000064336.1','ENSMUSG00000064337.1','ENSMUSG00000064349.1')
data <- data %>%
  dplyr::filter(!(V4 %in% bk)) %>%
  dplyr::filter(V1 %in% paste0('chr',c(1:21,'X','Y')))

# data <- data %>% as.data.frame()
data_plus <- data[data$V6=='+',]
data_minus <- data[data$V6=='-',]
# head(data[1:10,1:10])
data_plus[,7:ncol(data_plus)] <- abs(data_plus[,7:ncol(data_plus)])
for (i in seq(1,2*sample_num,by = 2)){
  # browser()
  data_minus[,c((7+step*(i-1)):(6+step*i),(7+step*i):(6+step*(i+1)))] <- abs(data_minus[,c((7+step*i):(6+step*(i+1)),(7+step*(i-1)):(6+step*i))])
  # data_minus[c((7+(i-1)*step):(6+i*step),(7+(i)*step):(6+(i+1)*step))] <- data_minus[c((7+(i)*step):(6+(i+1)*step),(7+(i-1)*step):(6+i*step))]
  # (7+step*1):(6+step*2),(7+step*0):(6+step*1)
}

data <- rbind(data_plus,data_minus)
data[is.na(data)] <- 0
data$sum <- rowSums(data[,7:ncol(data)])
data <- data[order(-data$sum),]
data <- data[round(0.05*nrow(data)):round(0.95*nrow(data)),1:(ncol(data)-1)]
# data$sum <- rowSums(data[,1.5*step:(1.59*step)])
# data <- data[order(-data$sum),]
# data <- data[round(0.1*nrow(data)):round(0.8*nrow(data)),1:(ncol(data)-1)]

data_sts <- data[,7:ncol(data)] %>% dplyr::summarise_all(mean) %>% as.data.frame()
data_sts <- as.data.frame(t(data_sts))
colnames(data_sts) <- 'value'
data_sts$type <- rep(paste0(rep(sample_name,each=2),rep(c('_F','_R'),sample_num)),each=step)
data_sts$sample <- rep(sample_name,each=2*step)
data_sts$coor <- rep(seq(-3000,3000,length.out = step),2*sample_num)
data_sts$value <- data_sts$value * rep(rep(c(1,-1),sample_num),each=step)
# data_sts$value <- -data_sts$value

ggplot(data_sts,aes(coor,value,fill=type,color=sample))+geom_line()+
  scale_x_continuous(limits=c(-3000,3000),breaks = c(-3000,0,3000),
                     labels = c('-3000','TSS','3000'))+xlab('')+
  geom_hline(yintercept = 0)+
  ylab('Read density(normalized)')+
  #scale_y_continuous(limits = c(0.4,1),expand=c(0,0))+
  theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(), panel.background=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2))
ggsave(paste0(name,'_divergent_rep.pdf'),width=4,height = 3)
