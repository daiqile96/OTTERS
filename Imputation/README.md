# An integrated tool to run P+T, lassosum, SDPR, and PRS-CS as the GReX imputation model.

## Input Files
                
 - Annotation File : 

    | CHROM | GeneStart | GeneEnd |     TargetID    | GeneName | 
    |:-----:|:---------:|:-------:|:---------------:|:--------:|
    |   1   |    100    |   200   |     ENSG0000    |     X    |

   Example: Example/exp_anno.txt

 - PLINK binary files(.bim/.bed/.fam) for genotype data from LD reference panel.

   Example: Example/Exp_geno.bed; Example/Exp_geno.bim; Example/Exp_geno.fam

 - eQTL summary statistics (bgzipped and tabixed)

    | CHROM | POS | A1 | A2 | Zscore |  TargetID   |   N  |
    |:-----:|:---:|:--:|:--:|:------:|:-----------:|:----:|
    |   1   | 100 |  C |  T |   3    |   ENSG0000  |  200 |

   Example: Example/Exp_SumStats.txt.gz; Example/Exp_SumStats.txt.gz.tbi


## Example:

```bash
cd ${OTTERS_DIR}/Example

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
--thread=2

# Step 2: Run Imputation models
output_dir=Outputs

# P+T
python3 ${OTTERS_DIR}/Imputation/PT/OTTERS_PT.py \
--OTTERS_dir=${OTTERS_DIR} \
--in_dir=${input_to_imputation} \
--out_dir=${output_dir} \
--chrom=${chr} \
--pt=0.1,0.05 \
--thread=2

# lassosum
python3 ${OTTERS_DIR}/Imputation/lassosum/OTTERS_lassosum.py \
--OTTERS_dir=${OTTERS_DIR} \
--in_dir=${input_to_imputation} \
--out_dir=${output_dir} \
--chrom=${chr} \
--thread=2

# SDPR
SDPR_dir=${SDPR_DIR}
python3 ${OTTERS_DIR}/Imputation/SDPR/OTTERS_SDPR.py \
--OTTERS_dir=${OTTERS_DIR} \
--SDPR_dir=${SDPR_DIR} \
--in_dir=${input_to_imputation} \
--out_dir=${output_dir} \
--chrom=${chr} \
--r2=0.1 \
--thread=2

# PRS-CS
# This may take several minutes to run
python3 ${OTTERS_dir}/Imputation/PRScs/OTTERS_PRScs.py \
--OTTERS_dir=${OTTERS_dir} \
--in_dir=${input_to_imputation} \
--out_dir=${output_dir} \
--chrom=${chr} \
--n_iter=500 \
--n_burnin=200 \
--phi=1e-4 \
--thread=1

```

## Outputs of the Example: 

The estimated eQTL weights for SDPR, P+T, lassosum will be saved under OUTPUT_DIR with the format:

  - | CHROM | POS | A1 | A2 |     TargetID    |  ES  |
    |:-----:|:---:|:--:|:--:|:---------------:|:----:|
    |   1   | 100 |  C |  T |     ENSG0000    |  0.2 |
 
  - Example/Outputs/P0.05.txt & Exp/Outputs/P0.1.txt for P+T with p-value thresholding 0.1 and 0.05.
  - Example/Outputs/lassosum.txt
  - Example/Outputs/SDPR.txt
  - Example/Outputs/PRScs.txt

## Detailed usage


  - [Step 1: Prepare inputs and perform LD-clumping](prep/README.md)

  - Step 2: Run Imputation models

    - [P+T](PT/README.md) 

    - [lassosum](lassosum/README.md) 

    - [SDPR](SDPR/README.md) 

    - [PRS-CS](PRScs/README.md) 

