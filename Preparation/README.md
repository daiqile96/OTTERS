# Perform LD clumping and Prepare clumped summary statistics

## Usage:
python prepare_input.py --OTTERS_dir=PATH_TO_OTTERS --anno_dir=PATH_TO_ANNO --b_dir=PATH_TO_BINARY --sst_dir=SUM_STATS_FILE --chrom=CHROM 
       [--r2=R2 --window=WINDOW --thread=THREAD --help]
                
 - PATH_TO_OTTERS: The directory of OTTERS source code

 - PATH_TO_ANNO: Full path and the file name of the gene annotation file. The annotation file is assumed to be in this format:

    | CHROM | GeneStart | GeneEnd |     TargetID    | GeneName | 
    |:-----:|:---------:|:-------:|:---------------:|:--------:|
    |   1   |    100    |   200   |     ENSG0000    |     X    |

 - PATH_TO_BINARY:  The directory of PLINK binary files for genotype data.

 - SUM_STATS_FILE: Full path and the file name of the summary statistics. 

    | CHROM | POS | REF | ALT | Zscore |  TargetID   |   N  |
    |:-----:|:---:|:---:|:---:|:------:|:-----------:|:----:|
    |   1   | 100 |  C  |  T  |   3    |   ENSG0000  |  0.2 |

 - OUTPUT_DIR: Output directory and output filename prefix.

 - CHROM: An integer for the chromosome of tested genes.  

 - WINDOW (optional): Window size (in base pairs) around gene region from which to include SNPs (default: 1000000 [+- 1MB region around gene region])
  
 - THREAD (optional): Number of simultaneous processes to use for parallel computation. Default is 1.

 - SEED (optional): Non-negative integer which seeds the random number generator.


## Example:

```bash
cd ${OTTERS_dir}/Example

python3 ${OTTERS_dir}/prepare_input.py \
--anno_dir=${anno_dir} \
--b_dir=${bim_dir} \
--sst_dir=${sst_dir} \
--chrom=${chr} 
```

## Output:

For a gene with ID ENSG0001 on chrom k: 

The median sample size of all eQTLs for the target gene in the eQTL summary statistics will be saved under ```${bim_dir}/CHRk/```

  - MedianN file:
   
    | TargetID 	|  MedianN 	|
    |----------	|-----------|
    | ENSG0001 	|    3000 	|
    | ENSG0002 	|    4000 	|

Other results for gene ENSG0001 will be saved under ```${bim_dir}/CHRk/ENSG0001/```.
  
 - Genotype data from the reference panel in PLINK binary file format for the target gene.
    - ENSG0001.bed
    - ENSG0001.bim
    - ENSG0001.fam

 - Clumping results for the target gene 
   - ENSG0001.clumped
    
 - Summary statistics of Zscore that will be used in SDPR:
   
    | SNP       	| A1 	| A2 	| Z      	|
    |-----------	|----	|----	|--------	|
    | rs737657  	| A  	| G  	| -2.044 	|
    | rs7086391 	| T  	| C  	| -2.257 	|
    | rs1983865 	| T  	| C  	| 3.652  	|

- Summary statistics of standardized beta that will be used in lassosum and P+T:
    
    | SNP       	| A1 	| A2 	| Beta      |
    |-----------	|----	|----	|-----------|
    | rs737657  	| A  	| G  	| -2.044 	|
    | rs7086391 	| T  	| C  	| -2.257 	|
    | rs1983865 	| T  	| C  	| 3.652  	|

