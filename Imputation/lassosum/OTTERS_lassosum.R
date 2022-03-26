#!/usr/bin/env Rscript

###################################################################
# Import packages needed
library(foreach)
library(bigstatsr)


###############################################################
# parse input arguments
Sys.setlocale("LC_ALL", "C")
options(stringsAsFactors = F)
# args <- (commandArgs(TRUE))
# print(args)
# if (length(args) == 0) {
#   stop("Error: No arguments supplied!")
# } else {
#   chr <- as.numeric(args[[1]])
#   medianN_path <- as.character(args[[2]])
#   bim_dir <- as.character(args[[5]])
#   sst_dir <- as.character(args[[6]])
#   out_dir <- as.character(args[[7]])
#   LD_blocks <- as.character(args[[8]])
#   n_thread <- as.character(args[[8]])
# }

## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      The R Script
 
      Arguments:
      --arg1=someValue   - numeric, blah blah
      --arg2=someValue   - character, blah blah
      --arg3=someValue   - logical, blah blah
      --help              - print this text \n")
  
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

## Check mandatory parameters
if (is.null(argsL$chr)) {
  cat('* Please specify the chromosome --chr\n')
  q(save="no")
} else if(is.null(argsL$medianN_path)) {
  cat('* Please specify the path to the median sample size file using --medianN_path\n')
  q(save="no")
} else if (is.null(argsL$bim_dir)) {
  cat('* Please specify the directory to the reference PLINK file --bim_dir\n')
  q(save="no")
} else if (is.null(argsL$sst_dir)) {
  cat('* Please specify the path to the summary statistics with standardized beta using --sst_dir\n')
  q(save="no")
} else if (is.null(argsL$LDblocks)) {
  cat('* Please specify the name of LDblocks --LDblocks\n')
  q(save="no")
} 

## Check optional parameters and assign default values
if (is.null(argsL$n_thread)){
  argsL$n_thread <- 1
}

print(argsL)

###############################################################
# time calculation
start_time <- Sys.time()

# Create the output file
chr = argsL$chr
out_path <- file.path(argsL$out_dir, paste0("chr", chr, "_lassosum.txt"))
print(paste0('Creating output file: ',out_path))
out_cols = c('CHROM','POS','REF','ALT','TargetID', 'ES')
write.table(data.frame(t(out_cols)),
            out_path,
            quote = F,
            row.names= F,
            col.names= F,
            sep = "\t")

# function to train lassosum for one gene
train_lassosum <- function(num) {
  
  gene_name = gene_names[num]
  
  library(data.table)
  library(lassosum)
  
  r2 = NULL
  
  # Specify the PLINK file of the reference panel
  bfile <- file.path(argsL$bim_dir, gene_name, gene_name)
  # Specify the summary statistics of standardized beta in single variant test
  sst_path <- file.path(argsL$sst_dir, paste0("chr", chr,
                                              "_", gene_name,
                                              "_clumped_sst_Beta.txt"))
  
  # check if we have reference data for the target gene
  if (file.exists(paste0(bfile, ".bim"))){
    
    # the standardized beta in single variant test = correlation
    ss <- fread(sst_path)
    cor <- ss$Beta
    
    # check if we have summary statistics for the target gene
    if (nrow(ss) > 0){
      
      ss$SNPPos <- sapply(1:length(ss$SNP), function(i) strsplit(ss$SNP[i], "_")[[1]][2])
      ss$Chrom <- chr
      
      # train lassosum
      out <- lassosum.pipeline(cor=cor,
                               chr=as.numeric(ss$Chrom),
                               pos=as.numeric(ss$SNPPos),
                               A1=ss$A1,
                               A2=ss$A2, # A2 is not required but advised
                               s = c(0.2, 0.5, 0.9, 1),
                               lambda = exp(seq(log(0.0001), log(0.1), length.out = 20)),
                               ref.bfile = bfile, # The reference panel dataset
                               test.bfile = bfile, # We don't have test data here
                               LDblocks = argsL$LDblocks,
                               exclude.ambiguous = T,
                               destandardize = F,
                               trace = 0)
      
      # perform pseudovalidation
      v <- pseudovalidate(out)
      lassosum_out <- subset(out, s=v$best.s, lambda=v$best.lambda)
      
      # save estimated beta
      sumstats = lassosum_out$sumstats[, c("chr", "pos", "A1", "A2")]
      beta = unlist(lassosum_out$beta)
      results = data.frame(sumstats, TargetID = gene_name, beta = beta)
      
      
      
      write.table(results,
                  out_path,
                  quote = F,
                  row.names= F,
                  col.names= F,
                  sep = "\t",
                  append = T)
      
      print(paste("Done train gene", i, ":", gene_name))
      
    } else {
      
      print("No summary statistics")
    }
    
  } else {
    
    print(paste("No", gene_name))
    
  }
  
}


###############################################################
### Read in gene annotation
medianN <- read.table(argsL$medianN_path, header = F)
colnames(medianN) <- c("TargetID", "N")
gene_names = unlist(medianN$TargetID)

###############################################################
# thread process
print(paste0('Starting lassosum for ', length(gene_names), ' target genes'))
cl <- parallel::makeCluster(as.numeric(argsL$n_thread))
doParallel::registerDoParallel(cl)

foreach(i = 1:10, .combine = 'c') %dopar% train_lassosum(num = i)
print(paste('Done Training lassosum on Chromosome', chr))

# ###############################################################
# # time calculation
# end_time <- Sys.time()
# print(paste("Computation time: ", end_time - start_time))

# # Install lassosum
# install.packages(c("RcppArmadillo", "data.table", "Matrix"), dependencies=TRUE)
# library(devtools)
# install_github("tshmak/lassosum")
# # Data of lassosum:
# # cd /Library/Frameworks/R.framework/Versions/4.0/Resources/library/lassosum/data
