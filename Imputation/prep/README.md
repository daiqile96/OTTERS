# Prepare Inputs to run P+T, lassosum, SDPR, and PRS-CS as the GReX imputation model.

## Usage:

```
python prep.py \
 --OTTERS_dir=PATH_TO_OTTERS \
 --anno_dir=PATH_TO_ANNO \
 --geno_dir=PATH_TO_GENO \
 --sst_file=PATH_TO_SST \
 --out_dir=OUTPUT_DIR 
 --chrom=CHROM 
 [--r2=R2 --window=WINDOW --thread=THREAD --help]
```

 - PATH_TO_OTTERS: The directory of OTTERS source code

 - PATH_TO_ANNO: Full path and the file name of the gene annotation file. The annotation file is assumed to be in this format:

    | CHROM | GeneStart | GeneEnd |     TargetID    | 
    |:-----:|:---------:|:-------:|:---------------:|
    |   1   |    100    |   200   |     ENSG0000    | 

   See [Example](Example/exp_anno.txt) of annotation file. 

 - PATH_TO_GENO:  The directory and preflix of PLINK binary files(.bim/.bed/.fam) for genotype data from LD reference panel.

 - PATH_TO_SST: Full path and the file name of the summary statistics. 

    | CHROM | POS | A1 | A2 |   Z    |  TargetID  |   N  |
    |:-----:|:---:|:--:|:--:|:------:|:----------:|:----:|
    |   1   | 100 |  C |  T |   3    |  ENSG0000  |  0.2 |

 - OUTPUT_DIR: The directory to save output files. These output files will be used as input files for GReX Imputation models. 

 - CHROM: An integer for the chromosome of tested genes. 

 - R2 (optional): The R square threshold used to perform LD-clumping. (default: 2, which means will not perform LD-clumping])

 - WINDOW (optional): Window size (in base pairs) around gene region from which to include SNPs (default: 1000000 [+- 1MB region around gene region])
  
 - THREAD (optional): Number of simultaneous processes to use for parallel computation. Default is 1.


## Example:

```bash
cd ${OTTERS_DIR}/Example

exp_anno=exp_anno.txt
geno_dir=Exp_geno
sst_dir=Exp_SumStats.txt
input_to_imputation=Inputs
chr=4

# Step 1: Prepare inputs with LD-clumping R^2 = 0.05
python3 ${OTTERS_DIR}/Imputation/prep/prep.py \
--OTTERS_dir=${OTTERS_DIR} \
--anno_dir=${exp_anno} \
--geno_dir=${geno_dir} \
--sst_file=${sst_dir}.gz \
--out_dir=${input_to_imputation} \
--chrom=${chr} \
--r2=0.99 \
--thread=2
```

## Outputs of the Example: 

   - Example/Inputs
     - ENSG000001
       - ENSG000001(.bim/.bed/.fam): PLINK binary genotype data from LD reference panel for gene ENSG000001.
       - ENSG000001.clumped: LD-clumping results for this gene.
       - Summary Statistics:
         - ENSG000001_Zscore.txt: eQTL summary statistics (Zscore for each eQTL in the single variant regression). The `SNP` column is the ID for each eQTL, generated as CHROM_POS_A2_A1. 

            |    SNP   | A1 | A2 | Z |
            |:--------:|:--:|:--:|:-:|
            | 1_100_T_C| C  | T  | 1 |
           
         - ENSG000001_beta.txt: eQTL summary statistics (standardized regression coefficient for each eQTL in the single variant). 

            |    SNP   | A1 | A2 | Beta |
            |:--------:|:--:|:--:|:-:|
            | 1_100_T_C| C  | T  | 0.5 |

     - ENSG000002
       - ENSG000002(.bim/.bed/.fam)
       - ENSG000002.clumped
       - Summary Statistics:
         - ENSG000002_Zscore.txt: eQTL summary statistics (Zscore for each eQTL in the single variant regression). The `SNP` column is the ID for each eQTL, generated as CHROM_POS_A2_A1. 

            |    SNP   | A1 | A2 | Z |
            |:--------:|:--:|:--:|:-:|
            | 1_100_T_C| C  | T  | 1 |
           
         - ENSG000002_beta.txt: eQTL summary statistics (standardized regression coefficient for each eQTL in the single variant). 

            |    SNP   | A1 | A2 | Beta |
            |:--------:|:--:|:--:|:-:|
            | 1_100_T_C| C  | T  | 0.5 |

     - medianN.txt: the median sample size of eQTLs for target genes. 

          |   TargetID   |  N  | 
          |:------------:|:---:|
          |  ENSG000001  | 1000| 
          |  ENSG000002  | 2000| 
