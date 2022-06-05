# OTTERS: **O**mnibus **T**ranscriptome **T**est using **E**xpression **R**eference **S**ummary data

![OTTERS Framework](Manuscript/F1.png)

A powerful TWAS framework leveraging summary-level reference data. 

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

  - Python modules/libraries:
    - pandas
    - scipy
    - numpy
  
      *Here is my code to set up the Python Environment to use OTTERS:*
      ```bash
      # install miniconda
      wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
      bash Miniconda3-latest-Linux-x86_64.sh
      source ~/.bashrc
      # remove sh file 
      rm Miniconda3-latest-Linux-x86_64.sh

      # create the environment 
      conda create --name otters python=3.9 pandas numpy scipy 
      # deactivate the conda environment
      conda deactivate
      ```

      *To use this environment to run OTTERS*:

      ```bash
        # activate the environment
      conda activate otters

      # RUN OTTERS HERE

      # deactivate the environment
      conda deactivate
      ```

  - If you want to apply SDPR and lassosum as imputation models, please install:
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

    - [lassosum](https://github.com/tshmak/lassosum) to perform lassosum. 
    - [fdrtool](https://github.com/tshmak/lassosum) to perform pseudo-validation implemented in lassosum. 
  
  

## Stage I

We provide an integrated tool that can 
 -  perform LD-clumping 
 -  prepare the required inputs (LD reference and eQTL summary statistics) in the required formats for P+T, lassosum, SDPR, and PRS-CS
 -  perform P+T, lassosum, SDPR, and PRS-CS.

In this tool, we 
 - integrate [TABIX](http://www.htslib.org/doc/tabix.html) and [PLINK 1.9](https://www.cog-genomics.org/plink) tools to extract input data per target gene more efficiently
 - enable parallel computation to train imputation models simultaneously for multiple genes. 

More details about using OTTERS to perform Stage I TWAS can be found [here](Imputation/README.md).

### *The example of using OTTERS to train eQTL weights using eQTL summary statistics and LD reference data:*

  ```bash
  # activate the environment
  conda activate otters

  # set number of threads to be used
  N_THREADS=1
 
  # set up my OTTERS directory and SDPR directory
  OTTERS_DIR=/home/qdai8/projects/bin/OTTERS
  SDPR_DIR=/home/qdai8/projects/bin/SDPR
  # make sure the dynamic libraries are not changed
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/MKL/lib
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/gsl/lib
  # prevent automatically using  all available cores on a compute node
  export MKL_NUM_THREADS=$N_THREADS
  export NUMEXPR_NUM_THREADS=$N_THREADS
  export OMP_NUM_THREADS=$N_THREADS

  # Start to run OTTERS
  cd ${OTTERS_DIR}/Example
  
  # set up input files
  exp_anno=exp_anno.txt
  geno_dir=Exp_geno
  sst_file=Exp_SumStats.txt
  input_to_imputation=Inputs
  chr=4

  # Step 1: Prepare inputs with LD-clumping R^2 = 0.99
  python3 ${OTTERS_DIR}/Imputation/prep/prep.py \
  --OTTERS_dir=${OTTERS_DIR} \
  --anno_dir=${exp_anno} \
  --geno_dir=${geno_dir} \
  --sst_file=${sst_file}.gz \
  --out_dir=${input_to_imputation} \
  --chrom=${chr} \
  --r2=0.99 \
  --thread=$N_THREADS

  # Step 2: Run Imputation models
  output_dir=Outputs

  # P+T
  python3 ${OTTERS_DIR}/Imputation/PT/OTTERS_PT.py \
  --OTTERS_dir=${OTTERS_DIR} \
  --in_dir=${input_to_imputation} \
  --out_dir=${output_dir} \
  --chrom=${chr} \
  --pt=0.1,0.05 \
  --thread=$N_THREADS

  # lassosum
  python3 ${OTTERS_DIR}/Imputation/lassosum/OTTERS_lassosum.py \
  --OTTERS_dir=${OTTERS_DIR} \
  --in_dir=${input_to_imputation} \
  --out_dir=${output_dir} \
  --chrom=${chr} \
  --thread=$N_THREADS

  # SDPR
  SDPR_dir=${SDPR_DIR}
  python3 ${OTTERS_DIR}/Imputation/SDPR/OTTERS_SDPR.py \
  --OTTERS_dir=${OTTERS_DIR} \
  --SDPR_dir=${SDPR_DIR} \
  --in_dir=${input_to_imputation} \
  --out_dir=${output_dir} \
  --chrom=${chr} \
  --r2=0.1 \
  --thread=$N_THREADS

  # PRS-CS
  # This may take several minutes to run
  python3 ${OTTERS_DIR}/Imputation/PRScs/OTTERS_PRScs.py \
  --OTTERS_dir=${OTTERS_DIR} \
  --in_dir=${input_to_imputation} \
  --out_dir=${output_dir} \
  --chrom=${chr} \
  --n_iter=500 \
  --n_burnin=200 \
  --phi=1e-4 \
  --thread=$N_THREADS
  ```

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
