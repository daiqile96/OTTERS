## SDPR: A GReX imputation model infers posterior eQTL effect sizes under DPR regression model. 

## Usage:

```bash
python OTTERS_SDPR.py \
--OTTERS_dir=PATH_TO_OTTERS \
--in_dir=INPUR_DIR \
--out_dir=OUTPUT_DIR \
--chrom=CHROM \
[--out_file=OUTPUT_FILE --M=PARAM_M --opt_llk=LIKLIHOOD_EQ --iter=MCMC_ITERATIONS --burn=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR 
 --r2=R2_FOR_BLOCK --a=PARAM_A --c=PARAM_C --a0k=PARAM_A0K --b0k=PARAM_B0K --thread=THREAD]
```

 - PATH_TO_OTTERS: The directory of OTTERS source code

 - INPUT_DIR:  Full path and the file name of the LD reference file.

 - OUTPUT_DIR: Output directory of the posterior effect size estimates.

 - CHROM: An integer for the chromosome of tested genes.

 - OUTPUT_FILE (optional): Output filename of the posterior effect size estimates. Default is SDPR.txt.

 - PARAM_M (optional): Max number of variance components. M must be greater than 4. Default is 1000.

 - LIKLIHOOD_EQ (optional): Which likelihood to evaluate. 1 for equation 6 (slightly shrink the correlation of SNPs) and 2 for equation 5 (SNPs genotyped on different arrays in a separate cohort). Please refer to manuscript or manual (3.2.3-2.3.4) for more details. Default is 1.

 - MCMC_ITERATIONS (optional): number of iterations for MCMC. Default is 1000.

 - MCMC_BURNIN  (optional): number of burn-in for MCMC. Default is 200.

 - MCMC_THINNING_FACTOR (optional): Thinning for MCMC. Default is 1 (no thin).

 - R2_FOR_BLOCK (optional): r2 cut-off for partition of independent blocks. Default is 0.1.

 - PARAM_A (optional): factor to shrink the reference LD matrix. Default is 0.1. Please refer to the manual for more information.

 - PARAM_C (optional): factor to correct for the deflation. Default is 1. Please refer to the manual for more information.

 - PARAM_A0K (optional): hyperparameter for inverse gamma distribution. Default is 0.5.

 - PARAM_B0K(optional): hyperparameter for inverse gamma distribution. Default is 0.5.

 - THREAD (optional): number of threads to use. Default is 1.

## Example:

```bash
input_to_imputation=Example/Inputs
output_dir=Example/Outputs
chr=4

python3 ${OTTERS_dir}/Imputation/SDPR/OTTERS_SDPR.py \
--OTTERS_dir=${OTTERS_dir} \
--SDPR_dir=${SDPR_dir} \
--in_dir=${input_to_imputation} \
--out_dir=${output_dir} \
--chrom=${chr} \
--r2=0.1 \
--thread=2
```

## Output:

Example/Outputs/SDPR.txt: the estimated eQTL weights

  - | CHROM | POS | A1 | A2 |     TargetID    |  ES  |
    |:-----:|:---:|:--:|:--:|:---------------:|:----:|
    |   1   | 100 |  C |  T |     ENSG0000    |  0.2 |