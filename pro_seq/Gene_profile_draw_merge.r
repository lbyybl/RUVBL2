#---- Pro-seq profile 
rm(list=ls())
setwd('') # set work dir
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
        panel.border = element_rect(colour = "black", fill=NA, size=0.2)) 
ggsave(paste0(prefix,'_gene.pdf'),width=4,height = 3)

