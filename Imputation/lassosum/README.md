# lassosum: A GReX imputation model infers posterior eQTL effect sizes using lasso regression.

## Usage:

```bash
Usage:
python OTTERS_PRScs.py --OTTERS_dir=PATH_TO_OTTERS --in_dir=INPUR_DIR --out_dir=OUTPUT_DIR --chrom=CHROM
[--out_file=OUTPUT_FILE --thread=THREAD]

 - PATH_TO_OTTERS: The directory of OTTERS source code

 - INPUT_DIR:  Full path and the file name of the LD reference file.

 - OUTPUT_DIR: Output directory of the posterior effect size estimates.

 - CHROM: An integer for the chromosome of tested genes.

 - OUTPUT_FILE (optional): Output filename of the posterior effect size estimates. Default is lassosum.txt

 - THREAD (optional): Number of simultaneous processes to use for parallel computation. Default is 1.
```

## Example:

```bash
input_to_imputation=Example/Inputs
output_dir=Example/Outputs
chr=4

python3 ${OTTERS_dir}/Imputation/lassosum/OTTERS_lassosum.py \
--OTTERS_dir=${OTTERS_dir} \
--in_dir=${input_to_imputation} \
--out_dir=${output_dir} \
--chrom=${chr} \
--thread=2
```

## Output:

Example/Outputs/lassosum.txt: the estimated eQTL weights

  - | CHROM | POS | A1 | A2 |     TargetID    |  ES  |
    |:-----:|:---:|:--:|:--:|:---------------:|:----:|
    |   1   | 100 |  C |  T |     ENSG0000    |  0.2 |
 

