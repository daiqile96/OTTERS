#!/usr/bin/env Rscript

# This script will convert summary statistics with Z score to 
# summary statisics that only containing clumped SNPs in three versions:
# with P-value and Beta that will be used in P+T
# with Zscore that will be used in SDPR

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
  clump = as.character(args[[2]])
  clumping_dir = as.character(args[[3]])
  medianN_dir = as.character(args[[4]])
  Zscore_dir = as.character(args[[5]])
  p_cutoff = as.numeric(args[[6]])
  out_dir = as.character(args[[7]])
}

gene_list_tab = read.table(medianN_dir, header = F)
colnames(gene_list_tab) = c("TargetID", "N")
gene_list = unlist(gene_list_tab$TargetID)
N = unlist(gene_list_tab$N)

clumping_dirs = file.path(clumping_dir, paste0(gene_list, "_clump", clump, ".valid.snp"))
sst_dirs = paste0(Zscore_dir, paste0("chr", chr, "_", gene_list, "_sst_Zscore.txt"))

format_PT_out = function(i){
  
  if (file.exists(clumping_dirs[i]) & file.exists(sst_dirs[i])){
    
    clumping_tab = as.vector(unlist(read.table(clumping_dirs[i])))

    sst_tab = read.table(sst_dirs[i], sep = "\t", header = T)
    sst_clumped = sst_tab[sst_tab$SNP %in% clumping_tab, ]
    
    sst_PT = sst_clumped %>% 
      mutate(P = pchisq(Z^2, df=1, lower.tail=FALSE),
             ES = Z/sqrt(N[i]),
             TargetID = gene_list[i],
             CHROM = chr,
             POS = sapply(1:length(SNP), function(i) strsplit(SNP[i], "_")[[1]][2])) %>% 
      mutate(POS = as.character(POS)) %>% 
      filter(P < p_cutoff) %>% 
      select(CHROM, POS, A1, A2, TargetID, ES)
    
    if (nrow(sst_PT) > 0){
      
      write.table(sst_PT, out_dir, append = T, 
                  row.names=F, col.names = F, quote = F, sep = "\t")
      
    }

  }
  
}

sapply(1:length(gene_list), function(i) format_PT_out(i))
