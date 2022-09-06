#--- this script is used to draw the MAplot of output by diffbind
rm(list = ls())
library(ggrastr,lib.loc = "/home/boyuanli/tools/anaconda2/envs/R40/lib/R/library")
source('/home/boyuanli/bashscript/bin/ruvb2/chip-seq/diff_peak/draw_smoothMA_fun.r')
setwd('/WORK/lbyybl/ypj/cross_talk/ChIPseq/Rpb3_degron_Rpb1Rpc1')
data <- read.table('diff_peak/Ruvbl_deseq2.txt',row.names = 1,stringsAsFactors = F,
                   header = T)
head(data)
data$FDR[is.na(data$FDR)] <- 0
pdf('Rpb3_degron_Rpc1.pdf',width = 5,height = 5)
draw_smooth_MA(data = data,xrange = c(2,4),yrange = c(-1.5,1.5),Fold.cutoff = log2(2))
dev.off()

fd <- log2(2)
xrange <- c(2,4)
xmin <- min(xrange)
xmax <- max(xrange)
yrange <- c(-1.5,1.5)
ymin <- min(yrange)
ymax <- max(yrange)
ggplot(data,aes(log2(Conc),Fold)) + rasterise(geom_point(
  color = ifelse(data$Fold > fd & data$FDR <= 0.05, '#ca3e47',ifelse(
    data$Fold < -fd & data$FDR <= 0.05, '#005792','grey'
  ))
),dpi=300) +scale_y_continuous(limits = yrange)+
  scale_x_continuous(limits = xrange)+
  annotate("text",x=xmax*0.9,y=ymax*0.9,label=paste0('n = ',sum(data$Fold > fd & data$FDR <= 0.05))) + 
  annotate("text",x=xmax*0.9,y=-ymax*0.9,label=paste0('n = ',sum(data$Fold < -fd & data$FDR <= 0.05))) + 
  geom_hline(yintercept = 0,linetype='dashed')+
  theme_bw()+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank()
  ) + xlab('log2 mean reads') + ylab('log2(IAA 1hr/untreated)')
ggsave('Polr1c_degron_RPB1_MA.pdf',width = 5,height = 5)
