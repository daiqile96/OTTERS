###################################################################
# Import packages needed
library(ACAT)

###################################################################
# load in arguments
args=(commandArgs(TRUE))

if(length(args)==0) {
  stop("Error: No arguments supplied!")
} else {
  twas_dir = args[[1]]
  methods = strsplit(args[[2]], spilt = ',')
}

###################################################################

get_pvalues = function(method){
  
  twas_file = file.path(twas_dir, paste0(method, '.txt'))
  twas = read.table(twas.file, head = T)
  
  out = twas[, c('CHROM', 'GeneStart', 'GeneEnd', 
                 'TargetID',  
                 'FUSION_PVAL')]
  
  return(twas)
  
}


