#---- Pro-seq profile 
rm(list=ls())
# setwd('/WORK/lbyybl/WH/rvb/pro-seq/Ruvbl2/sample20190921/graph/TSS_profile/')
setwd('/WORK/lbyybl/WH/rvb/RNAseq/sample20200330/result/gene_count/GO_result/kmeans5/noisoform/Proseq')
trans_df <- function(data,name){
  data <- abs(data)
  colnames(data) <- paste0('c',0:10999)
  data <- data %>%
    summarise_all(mean)
  data <- t(data)
  colnames(data) <- "value"
  data <- as.data.frame(data)
  data$sample <- name
  data$coor <- 1:11000
  return(data)
}
prefix <- 'unchange'
tss_profile <- fread(paste0(prefix,'_noiso_gene.gz'),skip = 1,stringsAsFactors = F)
tss_profile[,7:ncol(tss_profile)] <- abs(tss_profile[,7:ncol(tss_profile)])
TSS_plus <- tss_profile %>%
  dplyr::filter(V6=="+")
colnames(TSS_plus) <- paste0('V',1:ncol(TSS_plus))
TSS_minus <- tss_profile %>%
  dplyr::filter(V6=="-")
names(tss_profile)
step <- 11000
TSS_minus <- TSS_minus[,c((1+step*1):(6+step*2),(7+step*0):(6+step*1),
                          (7+step*3):(6+step*4),(7+step*2):(6+step*3),
                          (7+step*5):(6+step*6),(7+step*4):(6+step*5))]
colnames(TSS_minus) <- paste0('V',1:ncol(TSS_minus))
tss_profile <- rbind(TSS_plus,TSS_minus)


dox1 <- tss_profile[,(7+step*0):(6+step*1)]
dox2 <- tss_profile[,(7+step*1):(6+step*2)]
h051 <- tss_profile[,(7+step*2):(6+step*3)]
h052 <- tss_profile[,(7+step*3):(6+step*4)]
h11 <- tss_profile[,(7+step*4):(6+step*5)]
h12 <- tss_profile[,(7+step*5):(6+step*6)]

dox_mean<- trans_df(dox1,'dox')
h05_mean <- trans_df(h051,'h05')
h1_mean <- trans_df(h11,'h1')

doxr_mean<- trans_df(dox2,'dox')
h05r_mean <- trans_df(h052,'h05')
h1r_mean <- trans_df(h12,'h1')

plot(1:11000,tss_profile[5,7:11006],type='l')
plot(1:11000,tss_profile[5,11007:22006],type='l')
# plot(dox$coor,dox$value,type='l',xlim=c(-3000,3000),ylim=c(0,100))

merge <- rbind(dox_mean,h05_mean)
merge <- rbind(merge,h1_mean)

merger <- rbind(doxr_mean,h05r_mean)
merger <- rbind(merger,h1r_mean)

xrange <- c(0,11000)
yrange <- c(-25,65)

ggplot(merge,aes(coor,value,color=sample))+geom_line()+
  geom_line(data=merger,aes(coor,-value,color=sample))+
  xlab('')+ 
  ylab('Read density(normalized)')+ 
  scale_y_continuous(limits = yrange)+
  scale_x_continuous(breaks = c(0,3000,8000,10999), 
               labels = c('-3kb','TSS','TES','3kb'))+
  theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(), panel.background=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2)) #+
  # geom_vline(xintercept =c(90,330),color='blue')
ggsave(paste0(prefix,'_gene.pdf'),width=4,height = 3)
# tss_min <- c(-305,-495)
# tss_max <- c(120,330)
# # ctcf_max <-c(130,330,520,700,880,-130,-330,-520,-700,-880)
# # ctcf_min <- c(-230,-420,-610,-790,230,420,610,790)
# #plot(1:400,dox_mean[1,],type='l')
# #---quantify
# max_local <- function(data,loc){
#   data <- data %>%
#     filter(coor>(loc-50) & coor <(loc+50))
#   return(max(data$value))
# }
# min_local <- function(data,loc){
#   data <- data %>%
#     filter(coor>(loc-50) & coor <(loc+50))
#   return(min(data$value))
# }
# find_real_coor <- function(data,coor_df,type){
#   if (type=='max'){
#     for (i in 1:nrow(coor_df)){
#       va <- max_local(data,coor_df$cand[i])
#       coor_df$value[i]<- va
#       loc<-which(data$value==va)
#       coor_df$real_coor[i]<-data$coor[loc]
#     }
#   }else if (type=='min'){
#     for (i in 1:nrow(coor_df)){
#       va <- min_local(data,coor_df$cand[i])
#       coor_df$value[i]<- va
#       loc<-which(data$value==va)
#       coor_df$real_coor[i]<-data$coor[loc]
#     }
#   }
#   return(coor_df)
# }
# find_value <- function(data,coor){
#   data2 <- data.frame('coor'=coor,
#                       'value'=0)
#   for (i in 1:nrow(data2)){
#     loc <- which(data$coor==data2$coor[i])
#     data2$value[i] <- data$value[loc]
#   }
#   data2<-data2%>%
#     arrange(data2$coor)
#   return(data2)
# }
# #----find coor to quantify
# all <- dox_mean
# all$value <- dox_mean$value+h05_mean$value+h1_mean$value
# cand_max <- data.frame('cand'=c(120,330),
#                        'value'=0,'real_coor'=0)
# real_max <- find_real_coor(dox_mean,cand_max,'max')
# 
# cand_min <- data.frame('cand'=c(-305,-495),
#                        'value'=0,'real_coor'=0)
# real_min <- find_real_coor(dox_mean,cand_min,'min')
# 
# dox_max <- find_value(dox_mean,real_max$real_coor)
# dox_min <- find_value(dox_mean,real_min$real_coor)
# h05_max <- find_value(h05_mean,real_max$real_coor)
# h05_min <- find_value(h05_mean,real_min$real_coor)
# h1_max <- find_value(h1_mean,real_max$real_coor)
# h1_min <- find_value(h1_mean,real_min$real_coor)
# find_change <- function(dox_max,h05_max,h1_max){
#   names(dox_max) <- c('coor','dox')
#   names(h05_max) <- c('coor','h05')
#   names(h1_max) <- c('coor','h1')
#   merge <- merge(dox_max,h05_max)
#   merge_max <- merge(merge,h1_max)
#   merge_max$dox_05 <- merge_max$h05/merge_max$dox
#   merge_max$dox_1 <- merge_max$h1/merge_max$dox
#   merge_max$h05_1 <- merge_max$h1/merge_max$h05
#   return(merge_max)
# }
# 
# merge_max <- find_change(dox_max,h05_max,h1_max)
# merge_min <- find_change(dox_min,h05_min,h1_min)
# fwrite(merge_max,'tss_unbinding_max.csv')
# fwrite(merge_min,'tss_unbinding_min.csv')
# 
