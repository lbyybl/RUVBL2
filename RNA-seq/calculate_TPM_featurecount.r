#--- this script is used to calculate the TPM value from feature count result
# file <- '/WORK/lbyybl/ypj/cross_talk/ChIPseq/graph/RPB1_polycomb/data/Proseq/NTD_Proseq_count.txt'
calculatete_TPM_featurecount <- function(file){
  countdata<-read.table(file,skip = 1,sep="\t",header = T,row.names = 1)
  
  metadata <- countdata[,1:5]#提取基因信息count数据前的几列
  
  countdata <- countdata[,6:ncol(countdata)]#提取counts数，counts数据主题部分
  
  name<-"TPM"#设置输出文件前缀名
  #-----TPM Calculation------
  
  kb <- metadata$Length / 1000
  
  rpk <- countdata / kb
  
  calcu_tpm <- function(data){
    total_data <- data %>%
      summarise_all(sum)
    for (i in 1:ncol(data)){
      data[,i] <- (data[,i])/total_data[,i]*1e6
    }
    # data <- log(data)
    return(data)
  }
  countdata_tpm <- calcu_tpm(rpk)
  
  
  colSums(countdata_tpm)
  return(countdata_tpm)
}


