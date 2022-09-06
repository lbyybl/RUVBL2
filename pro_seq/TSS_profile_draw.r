#---- Pro-seq profile 
rm(list=ls())
setwd('/WORK/lbyybl/WH/rvb/pro-seq/Ruvbl2/sample20190921/graph/TSS_profile/rep/')
trans_df <- function(data,name){
  data <- abs(data)
  colnames(data) <- paste0('c',0:5999)
  data <- data %>%
    summarise_all(mean)
  data <- t(data)
  colnames(data) <- "value"
  data <- as.data.frame(data)
  data$sample <- name
  data$coor <- ((1:6000)-3000)
  return(data)
}

tss_profile <- fread('Ruvbl_unbinding_TSS_peak.gz',skip = 1,stringsAsFactors = F)
name <- 'unbind'
TSS_plus <- tss_profile %>%
  dplyr::filter(V6=="+")
colnames(TSS_plus) <- paste0('V',1:ncol(TSS_plus))
TSS_minus <- tss_profile %>%
  dplyr::filter(V6=="-")
names(tss_profile)
step <- 6000
TSS_minus <- TSS_minus[,c((1+step*1):(6+step*2),(7+step*0):(6+step*1),
                          (7+step*3):(6+step*4),(7+step*2):(6+step*3),
                          (7+step*5):(6+step*6),(7+step*4):(6+step*5),
                          (7+step*7):(6+step*8),(7+step*6):(6+step*7),
                          (7+step*9):(6+step*10),(7+step*8):(6+step*9),
                          (7+step*11):(6+step*12),(7+step*10):(6+step*11))]
colnames(TSS_minus) <- paste0('V',1:ncol(TSS_minus))
tss_profile <- rbind(TSS_plus,TSS_minus)

dox11 <- tss_profile[,(7+step*0):(6+step*1)]
dox12 <- tss_profile[,(7+step*1):(6+step*2)]
dox21 <- tss_profile[,(7+step*2):(6+step*3)]
dox22 <- tss_profile[,(7+step*3):(6+step*4)]
h0511 <- tss_profile[,(7+step*4):(6+step*5)]
h0512 <- tss_profile[,(7+step*5):(6+step*6)]
h0521 <- tss_profile[,(7+step*6):(6+step*7)]
h0522 <- tss_profile[,(7+step*7):(6+step*8)]
h111 <- tss_profile[,(7+step*8):(6+step*9)]
h112 <- tss_profile[,(7+step*9):(6+step*10)]
h121 <- tss_profile[,(7+step*10):(6+step*11)]
h122 <- tss_profile[,(7+step*11):(6+step*12)]
dox1_mean<- trans_df(dox11,'dox1')
dox2_mean <- trans_df(dox21,'dox2')
h051_mean <- trans_df(h0511,'h051')
h052_mean <- trans_df(h0521,'h052')
h11_mean <- trans_df(h111,'h11')
h12_mean <- trans_df(h121,'h12')

dox1r_mean<- trans_df(dox12,'dox1')
dox2r_mean <- trans_df(dox22,'dox2')
h051r_mean <- trans_df(h0512,'h051')
h052r_mean <- trans_df(h0522,'h052')
h11r_mean <- trans_df(h112,'h11')
h12r_mean <- trans_df(h122,'h12')
plot(-dox1_mean$coor,-dox1_mean$value,type='l')
# plot(dox$coor,dox$value,type='l',xlim=c(-3000,3000),ylim=c(0,100))

merge1 <- rbind(dox1_mean,h051_mean)
merge1 <- rbind(merge1,h11_mean)

merge2 <- rbind(dox2_mean,h052_mean)
merge2 <- rbind(merge2,h12_mean)

merge1r <- rbind(dox1r_mean,h051r_mean)
merge1r <- rbind(merge1r,h11r_mean)

merge2r <- rbind(dox2r_mean,h052r_mean)
merge2r <- rbind(merge2r,h12r_mean)

ggplot(merge1,aes(coor,value,color=sample))+geom_line()+
  geom_line(data=merge1r,aes(coor,-value,color=sample))+
  scale_x_continuous(limits=c(-3000,3000))+xlab('')+
  ylab('Read density(normalized)')+ 
  #scale_y_continuous(limits = c(0.4,1),expand=c(0,0))+
  theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(), panel.background=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2))
ggsave(paste0(name,'_divergent_rep1.pdf'),width=4,height = 3)
ggplot(merge2,aes(coor,value,color=sample))+geom_line()+
  geom_line(data=merge2r,aes(coor,-value,color=sample))+
  scale_x_continuous(limits=c(-3000,3000))+xlab('')+
  ylab('Read density(normalized)')+ 
  #scale_y_continuous(limits = c(0.4,1),expand=c(0,0))+
  theme(plot.title=element_text(hjust=0.5),
        panel.grid=element_blank(), panel.background=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2))
ggsave(paste0(name,'_divergent_rep2.pdf'),width=4,height = 3)


tss_min <- c(-305,-495)
tss_max <- c(120,330)
# ctcf_max <-c(130,330,520,700,880,-130,-330,-520,-700,-880)
# ctcf_min <- c(-230,-420,-610,-790,230,420,610,790)
#plot(1:400,dox_mean[1,],type='l')
#---quantify
max_local <- function(data,loc){
  data <- data %>%
    filter(coor>(loc-50) & coor <(loc+50))
  return(max(data$value))
}
min_local <- function(data,loc){
  data <- data %>%
    filter(coor>(loc-50) & coor <(loc+50))
  return(min(data$value))
}
find_real_coor <- function(data,coor_df,type){
  if (type=='max'){
    for (i in 1:nrow(coor_df)){
      va <- max_local(data,coor_df$cand[i])
      coor_df$value[i]<- va
      loc<-which(data$value==va)
      coor_df$real_coor[i]<-data$coor[loc]
    }
  }else if (type=='min'){
    for (i in 1:nrow(coor_df)){
      va <- min_local(data,coor_df$cand[i])
      coor_df$value[i]<- va
      loc<-which(data$value==va)
      coor_df$real_coor[i]<-data$coor[loc]
    }
  }
  return(coor_df)
}
find_value <- function(data,coor){
  data2 <- data.frame('coor'=coor,
                      'value'=0)
  for (i in 1:nrow(data2)){
    loc <- which(data$coor==data2$coor[i])
    data2$value[i] <- data$value[loc]
  }
  data2<-data2%>%
    arrange(data2$coor)
  return(data2)
}
#----find coor to quantify
all <- dox_mean
all$value <- dox_mean$value+h05_mean$value+h1_mean$value
cand_max <- data.frame('cand'=c(120,330),
                       'value'=0,'real_coor'=0)
real_max <- find_real_coor(dox_mean,cand_max,'max')

cand_min <- data.frame('cand'=c(-305,-495),
                       'value'=0,'real_coor'=0)
real_min <- find_real_coor(dox_mean,cand_min,'min')

dox_max <- find_value(dox_mean,real_max$real_coor)
dox_min <- find_value(dox_mean,real_min$real_coor)
h05_max <- find_value(h05_mean,real_max$real_coor)
h05_min <- find_value(h05_mean,real_min$real_coor)
h1_max <- find_value(h1_mean,real_max$real_coor)
h1_min <- find_value(h1_mean,real_min$real_coor)
find_change <- function(dox_max,h05_max,h1_max){
  names(dox_max) <- c('coor','dox')
  names(h05_max) <- c('coor','h05')
  names(h1_max) <- c('coor','h1')
  merge <- merge(dox_max,h05_max)
  merge_max <- merge(merge,h1_max)
  merge_max$dox_05 <- merge_max$h05/merge_max$dox
  merge_max$dox_1 <- merge_max$h1/merge_max$dox
  merge_max$h05_1 <- merge_max$h1/merge_max$h05
  return(merge_max)
}

merge_max <- find_change(dox_max,h05_max,h1_max)
merge_min <- find_change(dox_min,h05_min,h1_min)
fwrite(merge_max,'tss_unbinding_max.csv')
fwrite(merge_min,'tss_unbinding_min.csv')

# c <- rbind(dox_max,dox_min)
# ggplot(c,aes(coor,value))+geom_line()+geom_line(data = dox_mean,aes(coor,value))
