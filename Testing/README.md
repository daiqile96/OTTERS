Stage II: GReX imputation and gene-based association analysis 

- [Format](#format)
  - [PRS-CS](#prs-cs)
  - [lassosum](#lassosum)
  - [SDPR](#sdpr)
- [GReX imputation](#grex-imputation)
- [Gene-based association test](#gene-based-association-test)


# Format 

We first outputs from PRS-CS, SDPR, P+T and lassosum to TIGAR eQTL weights format.

We used the [TIGAR](https://github.com/yanglab-emory/TIGAR) tool to impute GReX. The TIGAR tool requires the eQTL weights to be in the format as:

    | CHROM | POS | REF | ALT |     TargetID    |  ES  |
    |:---:|:--:|:--:|:--:|:--------:|:--:|
    |   1   | 100 |  C  |  T  |     ENSG0000    |  0.2 |

## PRS-CS

The output eQTL weights of PRS-CS implemented in OTTERS are already in the TIGAR format. 

## lassosum

The output eQTL weights of lassosum implemented in OTTERS are already in the TIGAR format.

## SDPR

 - Input:

   - MedianN

   - BIM_DIR(SDPR_OUTPUT_DIR)

   - OUT_DIR

 - Example 

    ```bash
    # Format SDPR results.
    FORMAT=${OUT_DIR}/CHR${chr}_SDPR_formatted.txt
    Rscript ${OTTERS_DIR}/Testing/format_SDPR.R ${CHR} ${MedianN} ${BIM} ${FORMAT}

    # SORT, BGZIP, and TABIX the formatted results
    SORT=${OUT_DIR}/CHR${chr}_SDPR_sort.txt

    echo -e "CHROM\tPOS\tALT\tREF\tTargetID\tES" > ${SORT}
    echo Sorting.
    sort -n -k1 -k2 ${FORMAT} >> ${SORT}
    echo Done sorting.

    echo Bgzipping.
    #bgzip and tabix
    bgzip -f ${SORT}

    echo tabixing.
    tabix -f -b2 -e2 -S1 ${SORT}.gz
    echo Done.
    ```

# GReX imputation 

With the formatted, sorted, tabixed outputs from Step 1 (SORT.gz), we can directly use the TIGAR tool to predict the GReX of test samples. The SORT.gz will be used as the ```-- weight``` option in the TIGAR tool, i.e.:

```bash
${TIGAR_dir}/TIGAR_GReX_Pred.sh \
  --other options\
  --weight ${SORT}.gz
```
[Here](https://github.com/yanglab-emory/TIGAR#2-predict-grex) provides the detailed usage of using TIGAR to predict GReX. 

# Gene-based association test

Similar to the GReX imputation, with SORT.gz, we can directly use the TIGAR tool to perform TWAS with individual-level/summary-level GWAS data. The SORT.gz will be used as the ```-- weight``` option in the TIGAR tool, i.e.:

```bash
${TIGAR_dir}/TIGAR_TWAS.sh \
  --other options\
  --weight ${SORT}.gz
```

[Here](https://github.com/yanglab-emory/TIGAR#3-twas) provides the detailed usage of using TIGAR to perform TWAS.


 


    