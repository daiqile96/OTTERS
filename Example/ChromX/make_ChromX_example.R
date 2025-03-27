
# List all files in the Example folder
files <- list.files("../", full.names = TRUE)
files <- files[!file.info(files)$isdir]

# Copy each file with the new name
for (file in files) {
  new_name <- paste0("ChromX_", basename(file))
  file.copy(
    from = file,
    to = file.path(new_name),
    overwrite = TRUE  # Set to TRUE to overwrite existing files
  )
}

anno = read.table('ChromX_exp_anno.txt',
                  header = TRUE,
                  sep = "\t")
anno$CHROM = 'X'
write.table(anno,
            file = "ChromX_exp_anno.txt",
            col.names = TRUE,  
            row.names = FALSE,  
            sep = "\t",
            quote = FALSE)

bim_file = read.table('ChromX_Exp_geno.bim')
bim_file$V1 = 'X'
write.table(bim_file,
            file = "ChromX_Exp_geno.bim",
            col.names = FALSE,  
            row.names = FALSE,  
            sep = "\t",
            quote = FALSE)


eQTLSumStas <- read.table("ChromX_Exp_eQTLSumStats.txt", 
                          header = TRUE,
                          sep = "\t")
eQTLSumStas$CHROM = 'X'
write.table(eQTLSumStas,
            file = "ChromX_Exp_eQTLSumStats.txt",
            col.names = TRUE,  
            row.names = FALSE,  
            sep = "\t",
            quote = FALSE)

GWASSumStas <- read.table("ChromX_Exp_GWASSumStats.txt", 
                          header = TRUE,
                          sep = "\t")
GWASSumStas$CHROM = 'X'
write.table(GWASSumStas,
            file = "ChromX_Exp_GWASSumStats.txt",
            col.names = TRUE,  
            row.names = FALSE,  
            sep = "\t",
            quote = FALSE)


