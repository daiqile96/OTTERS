#!/usr/bin/env Rscript

###################################################################
# Import packages needed
library(foreach)
library(bigstatsr)
library(data.table)
library(lassosum)

###############################################################
# parse input arguments
Sys.setlocale("LC_ALL", "C")
options(stringsAsFactors = F)

## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
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
} else if(is.null(argsL$gene_name)) {
  cat('* Please specify the path to the median sample size file using --medianN_path\n')
  q(save="no")
} else if(is.null(argsL$medianN)) {
  cat('* Please specify the path to the median sample size file using --medianN_path\n')
  q(save="no")
} else if (is.null(argsL$bim_file)) {
  cat('* Please specify the directory to the reference PLINK file --bim_dir\n')
  q(save="no")
} else if (is.null(argsL$sst_file)) {
  cat('* Please specify the path to the summary statistics with standardized beta using --sst_dir\n')
  q(save="no")
} else if (is.null(argsL$LDblocks)) {
  cat('* Please specify the name of LDblocks --LDblocks\n')
  q(save="no")
} else if (is.null(argsL$out_path)) {
  cat('* Please specify the output path\n')
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
gene_name=argsL$gene_name
chr = argsL$chr

# Specify the PLINK file of the reference panel
bfile <- argsL$bim_file
# Read the summary statistics of standardized beta in single variant test
# the standardized beta in single variant test = correlation
ss <- fread(argsL$sst_file)
cor <- ss$Beta
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
                         exclude.ambiguous = F,
                         destandardize = F,
                         trace = 0)

# perform pseudovalidation
v <- pseudovalidate(out)
lassosum_out <- subset(out, s=v$best.s, lambda=v$best.lambda)

# save estimated beta
sumstats = lassosum_out$sumstats[, c("chr", "pos", "A1", "A2")]
beta = unlist(lassosum_out$beta)
results = data.frame(sumstats, TargetID = gene_name, ES = beta) 
results = results[, c("chr", "pos", "A1", "A2", "TargetID", "ES")]

write.table(results,
            argsL$out_path,
            quote = F,
            row.names= F,
            col.names= F,
            sep = "\t",
            append = T)

