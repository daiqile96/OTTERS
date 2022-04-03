# A tool to generate LD reference matrix

## Usage:

python ${OTTERS_dir}/LD/OTTERS_LD.sh --OTTERS_dir=PATH_TO_OTTERS --genome_block=PATH_TO_BLK --b_path=PATH_TO_BINARY --out_dir=OUTPUT_DIR --chrom=CHROM [--thread=THREAD]

 - OTTERS_DIR: The directory of OTTERS source code

 - PATH_TO_BINARY: The full path and prefix of PLINK binary files for genotype data from LD reference panel.

 - PATH_TO_BLK: The full path and name of the genome block annotation file used for generating the reference LD covariance files.

    |     CHROM    |   Start  | End | 
    |:------------:|:--------:|:---:|
    |    ENSG0000  |   1000   | 2000|
   
   The genome block annotation file is a tab-delimited text file with 3 columns CHROM Start End, denoting the chromosome number, block start position, and block ending position.

 - OUTPUT_DIR: Output directory

 - CHROM: An integer for the chromosome of tested genes.  
  
 - THREAD (optional): Number of simultaneous processes to use for parallel computation. Default is 1.


