Stage II: GReX imputation and gene-based association analysis 

- [Sort, Bgzip and Tabix the estimated eQTL weights](#sort-bgzip-and-tabix-the-estimated-eqtl-weights)
- [GReX imputation](#grex-imputation)
- [Gene-based association test](#gene-based-association-test)
- [Omnibus Test](#omnibus-test)


# Sort, Bgzip and Tabix the estimated eQTL weights

We take the estimated eQTL weights from SDPR as an example:

 - Example 

    ```bash
    # Format SDPR results.
    RAW=${OUT_DIR}/SDPR.txt

    # SORT, BGZIP, and TABIX the formatted results
    SORT=${OUT_DIR}/SDPR_sort.txt

    echo -e "CHROM\tPOS\tALT\tREF\tTargetID\tES" > ${SORT}
    echo Sorting.
    sort -n -k1 -k2 ${RAW} >> ${SORT}
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


 # Omnibus Test

After obtaining the TWAS p-values for the eQTL weights trained by each imputation model, we apply [ACAT-O](https://github.com/yaowuliu/ACAT) to combine the p-values based on different imputation models.  

    