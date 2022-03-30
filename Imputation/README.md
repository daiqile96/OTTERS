
# An integrated tool to run P+T, lassosum, SDPR, and PRS-CS as the GReX imputation model.

## Usage:
python OTTERS_imputation.py --OTTERS_dir=PATH_TO_OTTERS --SDPR_dir=PATH_TO_SDPR --anno_dir=PATH_TO_ANNO --b_dir=PATH_TO_BINARY --sst_dir=SUM_STATS_FILE --out_dir=OUTPUT_DIR --chrom=CHROM 
       [--r2=R2 --window=WINDOW --thread=THREAD --help]
                
 - OTTERS_DIR: The directory of OTTERS source code

 - SDPR_DIR: The directory of SDPR source code

 - PATH_TO_ANNO: Full path and the file name of the gene annotation file. The annotation file is assumed to be in this format:

    | CHROM | GeneStart | GeneEnd |     TargetID    | GeneName | 
    |:-----:|:---------:|:-------:|:---------------:|:--------:|
    |   1   |    100    |   200   |     ENSG0000    |     X    |

   ***We will run imputation models for all the genes in this annotation file. ***

 - PATH_TO_BINARY:  The directory of PLINK binary files for genotype data.

 - SUM_STATS_FILE: Full path and the file name of the summary statistics. 

    | CHROM | POS | REF | ALT | Zscore |  TargetID   |   N  |
    |:-----:|:---:|:---:|:---:|:------:|:-----------:|:----:|
    |   1   | 100 |  C  |  T  |   3    |   ENSG0000  |  0.2 |

 - OUTPUT_DIR: Output directory

 - CHROM: An integer for the chromosome of tested genes.  

 - WINDOW (optional): Window size (in base pairs) around gene region from which to include SNPs (default: 1000000 [+- 1MB region around gene region])
  
 - THREAD (optional): Number of simultaneous processes to use for parallel computation. Default is 1.


## Example:

```bash
cd ${OTTERS_DIR}/Example

python3 ${OTTERS_dir}/imputation.py \
--OTTERS_dir=${OTTERS_DIR} \
--SDPR_dir=${SDPR_DIR} \
--anno_dir=${PATH_TO_ANNO} \
--sst_dir=${SUM_STATS_FILE} \
--b_dir=${PATH_TO_BINARY} \
--out_dir=${OUTPUT_DIR} \
--chrom=${chr} 
```

## Outputs:

The estimated eQTL weights for SDPR, P+T, lassosum will be saved under OUTPUT_DIR with the format:

  - | CHROM | POS | REF | ALT |     TargetID    |  ES  |
    |:-----:|:---:|:---:|:---:|:---------------:|:----:|
    |   1   | 100 |  C  |  T  |     ENSG0000    |  0.2 |

and name as
 
  - lassosum.txt
  - SDPR.txt
  - P0.05_R0.99.txt for P+T with p-value thresholding 0.05 and LD-clumping $R^2$ 0.99.