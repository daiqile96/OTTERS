# lassosum: A GReX imputation model infers posterior eQTL effect sizes using LASSO regression

## Usage:

Rscript OTTERS_lassosum.R --medianN_path=PATH_TO_medianN --LDblocks=LD_BLOCKS --bim_dir=PATH_TO_BIM --sst_dir=PATH_TO_SUM_STATS out_dir=OUTPUT_DIR  --chr=CHR [--n_thread=THREAD]

  - PATH_TO_medianN: The full path of the file showing the median sample size of all eQTLs for the target gene in the eQTL summary statistics. The genes in this file **should be on the same chromosome**. The file is assumed to be in this format:

    | TargetID 	|  MedianN 	|
    |----------	|-----------|
    | ENSG0001 	|    3000 	|
    | ENSG0002 	|    4000 	|

  - LDblocks: Either (1) one of "EUR.hg19", "AFR.hg19", "ASN.hg19", "EUR.hg38", "AFR.hg38", "ASN.hg38", to use blocks defined by Berisa and Pickrell (2015) based on the 1000 Genome data, or (2) a vector to define LD blocks, or (3) a data.frame of regions in bed format. This is the same parameter as the LDblocks paramter in the [lassosum package](https://github.com/tshmak/lassosum).
  
  - CHR: An integer for the chromosome of the genes in the MedianN file. 

  - PATH_TO_BIM: The directory of binary files for reference LD panels. The reference panel is assumed to be in [PLINK 1 format](https://www.cog-genomics.org/plink/1.9/input#bed). Under the BIM_DIR, the path of the binary file for each target gene is assumed to be ***TargetID/TargetID.bim/bed/fam***

  - SST_DIR: The directory of summary statistics of standardized beta for the target gene. The summary statitistics is assumed to be in this format:
   
    | SNP       	| A1 	| A2 	| Beta    |
    |-----------	|----	|----	|--------	|
    | rs737657  	| A  	| G  	| -2.044 	|
    | rs7086391 	| T  	| C  	| -2.257 	|
    | rs1983865 	| T  	| C  	| 3.652  	|
    
    with the name ```chr${chr}_${TargetID}_clumped_sst_Beta.txt```. See [here]() for guidance to generate standardized beta for each gene. 
  
  - OUTPUT_DIR: Output directory
  
  - THREAD (optional): Number of simultaneous processes to use for parallel computation. Default is 1.

Rscript /home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/2_Final_RDA/Lassosum/OTTERS_lassosum.R \
--medianN_path=${medianN_path} \
--LDblocks=EUR.hg38 \
--bim_dir=${bim_dir} \
--sst_dir=${sst_dir} \
--out_dir=${out_dir} \
--chr=4 \
--n_thread=2 \



 