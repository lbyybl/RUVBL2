#--- this script is used to draw the venn plot of RPB7 RPB1 RPB7 big compolnent and RPB7 small compolnent
rm(list=ls())
setwd('/WORK/lbyybl/zhn/RPB7/ChIP-Seq/sample20210814/graph/venn')
RPB1 <- fread('/WORK/lbyybl/zhn/RPB7/ChIP-Seq/sample20210814/graph/venn/peak/RPB1_peaks.narrowPeak',stringsAsFactors = F,
              select = c(1:3),col.names = c('chr','st','en')) %>% unique()
RPB1 <- makeGRangesFromDataFrame(RPB1,seqnames.field = 'chr',start.field = 'st',end.field = 'en')

RPB7 <- fread('/WORK/lbyybl/zhn/RPB7/ChIP-Seq/sample20210814/graph/venn/peak/RPB7_peaks.narrowPeak',stringsAsFactors = F,
              select = c(1:3),col.names = c('chr','st','en')) %>% unique()
RPB7 <- makeGRangesFromDataFrame(RPB7,seqnames.field = 'chr',start.field = 'st',end.field = 'en')

RPB7_big <- fread('/WORK/lbyybl/zhn/RPB7/ChIP-Seq/sample20210814/peak/R7H_S10_12_GFPIP_peaks.narrowPeak',stringsAsFactors = F,
                  select = c(1:3),col.names = c('chr','st','en')) %>% unique()
RPB7_big <- makeGRangesFromDataFrame(RPB7_big,seqnames.field = 'chr',start.field = 'st',end.field = 'en')

RPB7_small <- fread('/WORK/lbyybl/zhn/RPB7/ChIP-Seq/sample20210814/peak/R7H_S15_17_GFPIP_peaks.narrowPeak',stringsAsFactors = F,
                    select = c(1:3),col.names = c('chr','st','en')) %>% unique()
RPB7_small <- makeGRangesFromDataFrame(RPB7_small,seqnames.field = 'chr',start.field = 'st',end.field = 'en')

ov <- findOverlapsOfPeaks(RPB1,RPB7,RPB7_big,
                          RPB7_small,connectedPeaks = 'merge',minoverlap = 1)
pdf('overlapall.pdf')
makeVennDiagram(ov,minoverlap = 1,connectedPeaks = 'merge',
                fill=c("#CC79A7", "#56B4E9", "#F0E442",'purple'),
                col=c("#D55E00", "#0072B2", "#E69F00",'red'),
                cat.col=c("#D55E00", "#0072B2", "#E69F00",'black'))
dev.off()
ov2 <- findOverlapsOfPeaks(RPB7,RPB7_big,
                           RPB7_small,connectedPeaks = 'merge',minoverlap = 1)
pdf('overlapRPB7.pdf')
makeVennDiagram(ov2,minoverlap = 1,connectedPeaks = 'merge',
                fill=c("#CC79A7", "#56B4E9", 'purple'),
                col=c("#D55E00", "#0072B2", 'red'),
                cat.col=c("#D55E00", "#0072B2", 'black'))
dev.off()

peak_list <- ov2$peaklist
names(peak_list) <- lapply(names(peak_list),function(x){gsub('///','_ov_',x)})

peak_list <- lapply(peak_list, as.data.frame)
setwd('overlap_txt/')
for (i in names(peak_list)){
  data <- peak_list[[i]]
  fwrite(data,paste0(i,'.txt'),sep = '\t')
}
