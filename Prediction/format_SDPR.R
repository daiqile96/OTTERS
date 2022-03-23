#!/usr/bin/env Rscript

Sys.setlocale("LC_ALL", "C")
options(stringsAsFactors=F)

### Load in arguments
library(dplyr)
args=(commandArgs(TRUE))
print(args)
if(length(args)==0) {
  stop("Error: No arguments supplied!")
} else {
  
  chr = as.character(args[[1]])
  medianN_dir = as.character(args[[2]])
  bim_dir = as.character(args[[3]])
  out_dir = as.character(args[[4]])
}

gene_list = unlist(read.table(medianN_dir, header = F)[, "V1"])                      

genes_bim_dir = file.path(bim_dir, gene_list)
genes_out_dir = file.path(genes_bim_dir, paste0(gene_list, "_out.txt"))
genes_bim_dir = file.path(genes_bim_dir, paste0(gene_list, ".bim"))      

out_all = NULL

for (i in 1:length(gene_list)){
  
  if (file.exists(genes_out_dir[i])){
    
    # extract SDPR results
    temp_out = read.table(genes_out_dir[i], header = T)
    colnames(temp_out) = c("SNP", "A1", "beta")
    
    # read bim file to get POS
    temp_bim = read.table(genes_bim_dir[i], header = F)
    colnames(temp_bim) = c("CHR", "SNP", "bp", "POS", "A1", "A2")
    
    # merge bim file with SDPR results 
    temp = merge(temp_bim, temp_out, by = c("SNP", "A1")) %>% 
      mutate(TargetID = gene_list[i]) %>% 
      select(CHR, POS, A1, A2, TargetID, beta)
    
    out_all = rbind(out_all, temp)
    
    print(paste("FINISH THE", i, "GENE:", gene_list[i]))
  
  } 
  
}

write.table(out_all, file = out_dir,
           row.names= F, col.names = F, quote = F, sep = "\t")


