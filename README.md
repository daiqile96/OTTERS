# OTTERS: **O**mnibus **T**ranscriptome **T**est using **E**xpression **R**eference **S**ummary data

![OTTERS Framework](Manuscript/F1.png)

A powerful TWAS framework leveraging summary-level reference data. 

## Stage I

We provide an [integrated tool](Imputation/README.md) that can 
 -  perform LD-clumping 
 -  prepare the required inputs (LD reference and eQTL summary statistics) in the required formats for P+T, lassosum, SDPR, and PRS-CS
 -  perform P+T, lassosum, SDPR, and PRS-CS.

In this tool, we integrate [TABIX](http://www.htslib.org/doc/tabix.html) and [PLINK 1.9](https://www.cog-genomics.org/plink) tools to extract input data per target gene more efficiently, and enable parallel computation to train imputation models simultaneously for multiple genes. 

The example of using OTTERS to train eQTL weights using eQTL summary statistics and LD reference data can be found [here](Imputation/README.md).

## Stage II

In Stage II, we impute GReX and perform gene-based association analysis in the test GWAS dataset. In this stage, we applied the [TIGAR tool](https://github.com/yanglab-emory/TIGAR) to perform GReX imputation and gene-based association analysis. 
  
 - [GRex Imputation](https://github.com/daiqile96/OTTERS/tree/main/Testing#grex-imputation)
  
 - [Gene-based association analysis](https://github.com/daiqile96/OTTERS/tree/main/Testing#gene-based-association-test)

 - [Omnibus Test](https://github.com/daiqile96/OTTERS/tree/main/Testing#omnibus-test)

Notes:

- OTTERS is based on [PLINK 1.9](https://www.cog-genomics.org/plink/), [PRS-CS software](https://github.com/getian107/PRScs), [lassosum R package](https://github.com/tshmak/lassosum), and [SDPR software](https://github.com/eldronzhou/SDPR). We really appreciate the authors of these softwares. 

- When run SDPR and PRS-CS, please use the following commands to prevent automatically using  all available cores on a compute node:

  ```bash
  export MKL_NUM_THREADS=$N_THREADS
  export NUMEXPR_NUM_THREADS=$N_THREADS
  export OMP_NUM_THREADS=$N_THREADS
  ```

- The scripts for generating all the Tables and Figures of the Manuscript can be found [here](https://htmlpreview.github.io/?https://github.com/daiqile96/OTTERS/blob/main/Manuscript/FiguresAndTables.html). 
