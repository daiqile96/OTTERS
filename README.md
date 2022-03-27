# OTTERS: Omnibus Transcriptome Test using Expression Reference Summary data

A powerful TWAS framework leveraging summary-level reference data

## Stage 0

Please refer to [here](Preparation/README.md) to perform LD-clumping and generate inputs that required in the Stage I imputation models for target genes. 

## Stage I

Stage I estimates eQTL weights from eQTL summary data and reference LD panel using four imputation models (P+T, lassosum, SDPR, PRS-CS)

Please refer to the following pages for manual of using each imputation model:

- [P+T](Imputation/P+T/README.md)

- [lassosum](Imputation/lassosum/README.md)

- [SDPR](Imputation/SDPR/README.md)

- [PRS-CS](Imputation/PRScs/README.md)


## Stage II

In Stage II, we impute the respective GReX using each method in Stage I and perform the respective gene-based association analysis in the test GWAS dataset. In this stage, we applied the [TIGAR tool](https://github.com/yanglab-emory/TIGAR) to perform GReX imputation and gene-based association analysis. 

 - [Format output eQTL weights](https://github.com/daiqile96/OTTERS/tree/main/Testing#format)
  
 - [GRex Imputation](https://github.com/daiqile96/OTTERS/tree/main/Testing#grex-imputation)
  
 - [Gene-based association analysis](https://github.com/daiqile96/OTTERS/tree/main/Testing#gene-based-association-test)

 - [Omnibus Test](https://github.com/daiqile96/OTTERS/tree/main/Testing#omnibus-test)

Notes:

- The PRS-CS implemented in OTTERS was based on the [PRS-CS software](https://github.com/getian107/PRScs) but prioritized for TWAS. In OTTERS-PRS-CS, we use [tabix](http://www.htslib.org/doc/tabix.html) tool to efficiently extract the summary statistics and the LD reference for each target gene and allow parallel computation for multiple genes. For P+T, SDPR, and lassosum, currently we provide the manual of using [PLINK 1.9](https://www.cog-genomics.org/plink), [lassosum](https://github.com/tshmak/lassosum), and [SDPR](https://github.com/eldronzhou/SDPR) to train eQTL weights. The python script prioritizing these softwares for TWAS will be available soon. We appreciate the authors of these softwares. 

- When run ***SDPR*** and ***PRS-CS***, please use the following commands to prevent automatically using  all available cores on a compute node:

  ```bash
  export MKL_NUM_THREADS=$N_THREADS
  export NUMEXPR_NUM_THREADS=$N_THREADS
  export OMP_NUM_THREADS=$N_THREADS
  ```

- The scripts for generating all the Tables and Figures of the Manuscript can be found [here](https://htmlpreview.github.io/?https://github.com/daiqile96/OTTERS/blob/main/Manuscript/FiguresAndTables.html). 
