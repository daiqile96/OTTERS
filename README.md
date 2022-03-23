# OTTERS: A powerful TWAS framework leveraging summary-level reference data

## Stage I: estimates eQTL weights from eQTL summary data and reference LD panel using four imputation models (P+T, lassosum, SDPR, PRS-CS)

Please refer to the following pages for manual of using each imputation model:

- [P+T](OTTERS/Imputation/P+T/README.md)

- [lassosum](Imputation/lassosum/README.md)

- [SDPR](Imputation/SDPR/README.md)

- [PRS-CS](PRScs/README.md)


Notes: 

- The PRS-CS implemented in OTTERS was based on the [PRS-CS software](https://github.com/getian107/PRScs) but prioritized for TWAS. In OTTERS-PRS-CS, we use [tabix](http://www.htslib.org/doc/tabix.html) tool to efficiently extract the summary statistics and the LD reference for each target gene and allow parallel computation for multiple genes. For P+T, SDPR, and lassosum, currently we provide the manual of using [PLINK 1.9](https://www.cog-genomics.org/plink), [lassosum](https://github.com/tshmak/lassosum), and [SDPR](https://github.com/eldronzhou/SDPR) to train eQTL weights. The python script prioritizing these softwares for TWAS will be available soon.

- When run ***SDPR*** and ***PRS-CS***, please use the following commands to prevent automatically using  all available cores on a compute node:

```bash
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS
```

## Stage 2: 

