## Usage:

```bash
python OTTERS_PT.py \
--OTTERS_dir=PATH_TO_OTTERS \
--in_dir=INPUR_DIR \
--out_dir=OUTPUT_DIR \
--chrom=CHROM
[--thread --pt=PT]

 - PATH_TO_OTTERS: The directory of OTTERS source code

 - INPUT_DIR:  Full path and the file name of the LD reference file.

 - OUTPUT_DIR: Output directory of the effect size estimates.

 - CHROM: An integer for the chromosome of tested genes.

 - PT(optional): The p-value threshold used in P+T, separated by comma, e.g., --pt=0.1,0.05,0.001. Default is 0.05.

 - THREAD (optional): Number of simultaneous processes to use for parallel computation. Default is 1.
```

## Example:

```bash
input_to_imputation=Exp/Inputs
output_dir=Exp/Outputs
chr=4

python3 ${OTTERS_dir}/Imputation/PT/OTTERS_PT.py \
--OTTERS_dir=${OTTERS_dir} \
--in_dir=${input_to_imputation} \
--out_dir=${output_dir} \
--chrom=${chr} \
--pt=0.1,0.05 \
--thread=2
```

## Output:

Exp/Outputs/P0.05.txt: the estimated eQTL weights

  - | CHROM | POS | A1 | A2 |     TargetID    |  ES  |
    |:-----:|:---:|:--:|:--:|:---------------:|:----:|
    |   1   | 100 |  C |  T |     ENSG0000    |  0.2 |
 