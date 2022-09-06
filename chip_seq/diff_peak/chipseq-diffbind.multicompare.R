#!/bin/bash
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
option_list <- list(
  make_option(c("-c", "--csv"), type = "character", default=NULL,
              help=".csv file read by DiffBind"),
  make_option(c("-p", "--parameter"), type="character", default="-f 1.5 -o diff_peak -p Peak",
              help="the parameter you want to use")
)

opt <- parse_args(OptionParser(option_list=option_list))


sampl_csv <- read.csv('/WORK/lbyybl/qff/ATACseq20220322/sample.csv',stringsAsFactors = F)
condition_num <- length(unique(sampl_csv$Condition))
csv_condiation <- unique(sampl_csv$Condition)
compare <- combn(condition_num,2) %>% as.data.frame()
compare <- apply(compare,2,function(x){
  csv_condiation[as.vector(x)]
}) %>% t() %>% as.data.frame()
colnames(compare) <- c('cond1','cond2')

apply(compare,1,function(x){
  # browser()
  tmp <- sampl_csv %>% dplyr::filter(Condition %in% x)
  tmp_name <- paste(x,collapse = '_vs_')
  write.csv(tmp,file = paste0(tmp_name,'sample.csv'),row.names = F,quote = F)
})

files <- apply(compare,1,function(x){
  tmp_name <- paste(x,collapse = '_vs_')
  return(paste0(tmp_name,'sample.csv'))
}) 

cmd1 <- paste0('ls ',paste(files,collapse = ' '))
cmd2 <- ' | parallel -j3 "Rscript /home/boyuanli/bashscript/bin/chip_seq/diff_peak/chipseq-diffbind.R -c {} "'
par <- opt$parameter
cmd <- paste0("system('",cmd1,cmd2,par,"')")
eval(parse(text = cmd))

