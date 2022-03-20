# OTTERS
A powerful TWAS framework leveraging summary-level reference data

## Example Data

The example LD reference can be downloaded [here](https://www.dropbox.com/sh/7ubnuzamh45pwgs/AACKCL2CsTXIkbynLozVAzXna?dl=0).
Other example input files provided under `./Example/`.

## Stage I: estimates eQTL weights from eQTL summary data and reference LD panel using four imputation models (P+T, lassosum, SDPR, PRS-CS)

### ***PRS-CS***

A GReX imputation model infers posterior eQTL effect sizes under continuous shrinkage (CS) priors
using eQTL summary statistics and an external LD reference panel.

Usage:
python OTTERS_PRScs.py --anno_dir=PATH_TO_ANNO --ld_dir=PATH_TO_LD --clumped_dir=PATH_TO_CLUMPED --sst_file=SUM_STATS_FILE --out_dir=OUTPUT_DIR --chrom=CHROM
                [--window=WINDOW_SIZE --a=PARAM_A --b=PARAM_B --phi=PARAM_PHI --n_iter=MCMC_ITERATIONS --n_burnin=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR --thread=THREAD --seed=SEED]
                
 - PATH_TO_OTTERS: The directory of OTTERS source code

 - PATH_TO_ANNO: Full path and the file name of the gene annotation file. 

 - PATH_TO_LD:  Full path and the file name of the LD reference file.

 - PATH_TO_CLUMPED: Full path and the file name of the clumped eQTL file. 
                    If not be provided, all the eQTLs in the eQTL summary statistics for the target gene will be used, 
                    i.e, no clumping will be performed.

 - SUM_STATS_FILE: Full path and the file name of the summary statistics. 

 - OUTPUT_DIR: Output directory and output filename prefix of the posterior effect size estimates.

 - CHROM: An integer for the chromosome of tested genes.  

 - WINDOW (optional): Window size (in base pairs) around gene region from which to include SNPs (default: 1000000 [+- 1MB region around gene region])

 - PARAM_A (optional): Parameter a in the gamma-gamma prior in the continous shrinkage prior proposed by PRS-CS. Default is 1. 

 - PARAM_B (optional): Parameter b in the gamma-gamma prior in the continous shrinkage prior proposed by PRS-CS. Default is 0.5.

 - PARAM_PHI (optional): Global shrinkage parameter phi. If phi is not specified, it will be learnt from the data using a fully Bayesian approach.
                         This usually works well for polygenic traits with large GWAS sample sizes (hundreds of thousands of subjects).
                         For GWAS with limited sample sizes (including most of the current disease GWAS), fixing phi to 1e-4 or 1e-2,
                         or doing a small-scale grid search (e.g., phi=1e-6, 1e-4, 1e-2, 1) to find the optimal phi value often improves perdictive performance.

 - MCMC_ITERATIONS (optional): Total number of MCMC iterations. Default is 1,000.

 - MCMC_BURNIN (optional): Number of burnin iterations. Default is 500.

 - MCMC_THINNING_FACTOR (optional): Thinning of the Markov chain. Default is 5.
 
 - THREAD (optional): Number of simultaneous processes to use for parallel computation. Default is 1.

 - SEED (optional): Non-negative integer which seeds the random number generator.


Example:

```bash
OTTERS_dir=/home/qdai/YangFSSdata2/qdai/tool/OTTERS
cd ${OTTERS_dir}/Example

python3 ${OTTERS_dir}/PRScs/OTTERS_PRScs.py \
--OTTERS_dir=${OTTERS_dir} \
--anno_dir=anno_exp.txt \
--ld_dir=${path_to_the_downloaded_example_ld}/ld_exp.txt.gz \
--clumped_dir=exp_clumped.txt.gz \
--sst_dir=sst_exp.txt.gz \
--out_dir=${OTTERS_dir}/Example/ \
--chrom=4 \
--phi=1e-4 \
--window=1000000 \
--thread=1 \
--seed=20210522  \
```
