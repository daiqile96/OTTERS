#!/usr/bin/env Rscript

Sys.setlocale("LC_ALL", "C")
.libPaths("/mnt/YangFSS/data/qdai/.local/R/4.0.4/lib64/R/library")
options(stringsAsFactors=F)
library(dplyr)
options(scipen=999)
### load in arguments
args=(commandArgs(TRUE))
print(args)
if(length(args)==0) {
  stop("Error: No arguments supplied!")
} else {
  
  chr = as.numeric(args[[1]])
  medianN_dir = as.character(args[[2]])
  Lassosum_dir = as.character(args[[3]])
  out_dir = as.character(args[[4]])
  
}

gene_list_tab = read.table(medianN_dir, header = F)
colnames(gene_list_tab) = c("TargetID", "N")
gene_list = unlist(gene_list_tab$TargetID)
N = unlist(gene_list_tab$N)

Lassosum_out = file.path(Lassosum_dir, paste0("chr", chr, "_lassosum_", gene_list, 
                                              ".txt"))

format_lassosum_out = function(i){
  
    if (file.exists(Lassosum_out[i])){
      
      # extract SDPR results
      temp_out = read.table(Lassosum_out[i], header = T)
      colnames(temp_out) = c("CHR", "POS", "A1", "A2", "TargetID", "beta")
      
      # merge bim file with SDPR results 
      temp = temp_out %>% 
        mutate(TargetID = gene_list[i]) %>% 
        mutate(POS = as.character(POS)) %>% 
        select(CHR, POS, A1, A2, TargetID, beta)
      
      write.table(temp, out_dir, append = T, 
                  row.names=F, col.names = F, quote = F, sep = "\t")
      
    } 

}

sapply(1:length(gene_list), function(i) format_lassosum_out(i))
