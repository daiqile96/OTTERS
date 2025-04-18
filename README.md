# OTTERS: **O**mnibus **T**ranscriptome **T**est using **E**xpression **R**eference **S**ummary data

![OTTERS Framework](F1.png)

A powerful TWAS framework leveraging summary-level reference data. 

We provide an integrated tool that can 
 -  perform LD-clumping 
 -  prepare the required inputs (LD reference and eQTL summary statistics) in the required formats for P+T, lassosum, SDPR, and PRS-CS
 -  perform TWAS Stage I: run P+T, lassosum, SDPR, and PRS-CS to train eQTL weights.
 -  perform TWAS Stage II: gene-based association test using GWAS summary statistics.

In this tool, we 
 - integrate [TABIX](http://www.htslib.org/doc/tabix.html) and [PLINK 1.9](https://www.cog-genomics.org/plink) tools to extract input data per target gene more efficiently
 - integrate [PRS-CS software](https://github.com/getian107/PRScs), [lassosum R package](https://github.com/tshmak/lassosum), and [SDPR software](https://github.com/eldronzhou/SDPR) to train eQTL weights
 - enable parallel computation to train GReX imputation models / perform gene-based association test simultaneously for multiple genes. 


## Latest Updates

### **January 14, 2025: Added `--lassosum_ld_blocks` argument**
  * Purpose: Select LD blocks for lassosum analysis  
  * Available options (from [Berisa & Pickrell (2015)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4731402/)):  
    - `EUR.hg19`, `ASN.hg19`, `AFR.hg19`  
    - `EUR.hg38`, `ASN.hg38`, `AFR.hg38`  
  * Default: `EUR.hg38`  

### **March 29, 2025: Added support for sex chromosomes (`--chrom=X` and `--chrom=Y`)**  
* Note: In OTTERS outputs:  
  - Chromosome X → labeled as 23  
  - Chromosome Y → labeled as 24  
* Example resources:  
  - [ChromX Analysis Script](https://raw.githubusercontent.com/daiqile96/OTTERS/main/Example/ChromX/ChromX_example.sh)  
    - [data to run this example analysis](https://github.com/daiqile96/OTTERS/tree/development/Example/ChromX)  

### **March 31, 2025: Added VCF/VCF.gz support for LD reference panels**
  * New `--geno_type` options: `vcf` 
  * Default remains: `plink` (binary PLINK format)
  * Note: 
    * Require bgzip and tabixed vcf files, e.g. [Exp_geno.vcf.gz](https://github.com/daiqile96/OTTERS/tree/development/Example/VCF/Exp_geno.vcf.gz) and [Exp_geno.vcf.gz.tbi](https://github.com/daiqile96/OTTERS/tree/development/Example/VCF/Exp_geno.vcf.gz.tbi)
    * For each gene, OTTERS uses the Tabix tool to extract the relevant region from the VCF, and then calls PLINK to convert the extracted vcf into binary PLINK format. See the functions [extract_vcf_region_with_tabix](https://github.com/daiqile96/OTTERS/blob/main/ottersutils.py#L233) and [call_PLINK_extract](https://github.com/daiqile96/OTTERS/blob/main/ottersutils.py#L267) for details.
  * Example resources:  
    - [Analysis with VCF.gz genotype data](https://raw.githubusercontent.com/daiqile96/OTTERS/main/Example/VCF/vcf_example.sh)  
      - `geno_dir=Exp_geno.vcf.gz` 
      - `--geno_type=vcf`
      - [data to run this example analysis](https://github.com/daiqile96/OTTERS/tree/development/Example/VCF)


## Getting Started

Download and install following required tools, modules, and packages:

  - Download OTTERS using:
    
    ```bash
    git clone https://github.com/daiqile96/OTTERS.git 
    ```

  - [BGZIP](http://www.htslib.org/doc/bgzip.html) and [TABIX](http://www.htslib.org/doc/tabix.html)

     *Here is my code to download and install BGZIP and TABIX:*

     ```bash
     wget https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2
     tar jxvf tabix-0.2.6.tar.bz2
     cd tabix-0.2.6
     make

     # copy the binary bgzip and tabix to my path: ~/projects/bin
     cp bgzip tabix ~/projects/bin

     # test if BGZIP and TABIX are successfully installed
     export PATH=~/projects/bin:$PATH 
     bgzip
     tabix
     ```

  - [PLINK 1.9](https://www.cog-genomics.org/plink/)
    
    *Here is my code to download and install PLINK 1.9:*

    ```bash
    # Download the latest binary of PLINK 1.9
    wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220402.zip

    # Unzip the archive
    unzip plink_linux_x86_64_20220402.zip -d ~/projects/bin

    # Remove archive
    rm plink2_linux_x86_64_latest.zip

    # Test if PLINK 1.9 is successfully installed
    export PATH=~/projects/bin:$PATH 
    plink
    ```
  - To permanently store the path, please add `export PATH=~/projects/bin:$PATH` to your ~/.bash_profile (or ~/.bashrc).

    Then reload file .bash_profile from the command line:
    ```bash
    source ~/.bash_profile
    ```

  - Python modules/libraries:

    - [pandas 1.4.4](https://pandas.pydata.org)
    - [scipy 1.7.3](https://scipy.org)
    - [numpy 1.21.5](https://numpy.org)
    - [pysam 0.19.1](https://pysam.readthedocs.io/en/latest/api.html) 

    If you're new to python, [here](Example/Exp_otters_env.sh) is an example to set up the Python Environment to use OTTERS.

  - To apply SDPR and lassosum as imputation models, please install:
    - [SDPR](https://github.com/eldronzhou/SDPR) to perform SDPR
      
      *Here is my code to download SDPR*
      
      ```bash
      cd ~/projects/bin
      git clone https://github.com/eldronzhou/SDPR.git

      # make sure that dynamic libraries gsl/lib/libgsl.so and MKL/lib/libmkl_rt.so are not changed
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/projects/bin/SDPR/MKL/lib
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/projects/bin/SDPR/gsl/lib
      ```
      **Please make sure that dynamic libraries are not changed every time when use SDPR.**
    
    - R packages:
      - [lassosum](https://github.com/tshmak/lassosum) to perform lassosum (lassosum requires several R packages, please follow the lassosum installation guidance to install all of them).
      - [fdrtool](https://cran.r-project.org/web/packages/fdrtool/index.html) to perform pseudo-validation implemented in lassosum. 
  
  
## Example: 

  ```bash
  # activate the environment
  conda activate otters

  # set number of threads to be used
  N_THREADS=1

  # set up my OTTERS directory and SDPR directory
  OTTERS_DIR=/home/qdai8/projects/bin/OTTERS
  SDPR_DIR=/home/qdai8/projects/bin/SDPR

  # make sure the dynamic libraries of SDPR are not changed (For SDPR)
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/MKL/lib
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/gsl/lib

  # prevent automatically using  all available cores on a compute node (For SDPR and PRS-CS)
  export MKL_NUM_THREADS=$N_THREADS
  export NUMEXPR_NUM_THREADS=$N_THREADS
  export OMP_NUM_THREADS=$N_THREADS

  # Load R to perform lassosum
  module load R

  # Start to run OTTERS
  cd ${OTTERS_DIR}/Example

  # Input for OTTERS STAGE I 
  # Annotation File 
  exp_anno=exp_anno.txt
  # Genotype data from LD reference panel
  geno_dir=Exp_geno
  # eQTL summary statistics 
  sst_file=Exp_eQTLSumStats.txt
  # Input for OTTERS STAGE II
  # GWAS summary statistics 
  gwas_sst_file=Exp_GWASSumStats.txt

  # Set chromosome number (The example is for chromosome 4)
  chr=4
  # Set LD-clumping threshold in STAGE I
  clump_r2=0.99
  # Set output directory for STAGE I
  out_dir=Results
  
  # STAGE I
  # train eQTL weights using P+T, lassosum, SDPR and PRS-CS. 
  # It may take several minutes to complete.
  python3 ${OTTERS_DIR}/training.py \
  --OTTERS_dir=${OTTERS_DIR} \
  --SDPR_dir=${SDPR_DIR} \
  --anno_dir=${exp_anno} \
  --geno_dir=${geno_dir} \
  --sst_file=${sst_file} \
  --out_dir=${out_dir} \
  --chrom=${chr} \
  --r2=${clump_r2} \
  --models=PT,lassosum,SDPR,PRScs \
  --lassosum_ld_blocks=EUR.hg38 \
  --thread=$N_THREADS
  
  # Set output directory for STAGE II
  twas_dir=TWAS

  # STAGE II
  # gene-based association test using eQTL-weight trained from P+T, lassosum, SDPR and PRS-CS.
  python3 ${OTTERS_DIR}/testing.py \
  --OTTERS_dir=${OTTERS_DIR} \
  --weight_dir=${OTTERS_DIR}/Example/Results \
  --models=P0.001,P0.05,lassosum,SDPR,PRScs \
  --anno_dir=${exp_anno} \
  --geno_dir=${geno_dir} \
  --out_dir=${twas_dir} \
  --gwas_file=${gwas_sst_file} \
  --chrom=${chr} \
  --thread=$N_THREADS

  # get imputed genetically regulated gene expression
  impute_dir=GReX
  # samples to perform imputation
  samples=Exp_samples.txt
  # imputation
  python3 ${OTTERS_DIR}/imputing.py \
  --OTTERS_dir=${OTTERS_DIR} \
  --weight_dir=${OTTERS_DIR}/Example/Results \
  --models=P0.001,P0.05,lassosum,SDPR,PRScs \
  --anno_dir=${exp_anno} \
  --geno_dir=${geno_dir} \
  --out_dir=${impute_dir} \
  --chrom=${chr} \
  --samples=${samples} \
  --thread=$N_THREADS
  ```

## Required Inputs 

 - Annotation File should contain following columns with the same name in same order (text file): 

    | CHROM | GeneStart | GeneEnd |     TargetID    | 
    |:-----:|:---------:|:-------:|:---------------:|
    |   1   |    100    |   200   |     ENSG0000    |

    Example: Example/exp_anno.txt

 - Genotype data from LD reference panel in PLINK binary files :

    Example: Example/Exp_geno.bed; Example/Exp_geno.bim; Example/Exp_geno.fam

 - eQTL summary statistics should contain following columns with the same name in same order (text file):

    *Please sort by chromosome and then by SNP position in ascending order*

    | CHROM | POS | A1 | A2 | Zscore |  TargetID   |   N  |
    |:-----:|:---:|:--:|:--:|:------:|:-----------:|:----:|
    |   1   | 100 |  C |  T |   3    |   ENSG0000  |  200 |

    To convert a two-sided pvalue and beta (coefficients) to Zscore, in python, we can do:
    ```python
    import numpy as np
    from scipy.stats import norm
    Z = np.sign(beta) * abs(norm.ppf(pvalue/2.0))
    ```
    Example: Example/Exp_eQTLSumStats.txt

 - GWAS summary statistics should contain following columns with the same name in same order (text file):
  
    *Not required to be sorted*

    | CHROM | POS | A1 | A2 | Zscore |  
    |:-----:|:---:|:--:|:--:|:------:|
    |   1   | 100 |  C |  T |   2    |

    Example: Example/Exp_GWASSumStats.txt

## Notes:

- OTTERS is based on [PLINK 1.9](https://www.cog-genomics.org/plink/), [PRS-CS software](https://github.com/getian107/PRScs), [lassosum R package](https://github.com/tshmak/lassosum), and [SDPR software](https://github.com/eldronzhou/SDPR). We really appreciate the authors of these softwares. 
