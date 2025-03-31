# activate the environment
conda activate otters

# set number of threads to be used
N_THREADS=1

# set up my OTTERS directory and SDPR directory
OTTERS_DIR=/home/qdai8/projects/Temp/OTTERS
SDPR_DIR=/home/qdai8/projects/bin/SDPR

# make sure the dynamic libraries of SDPR are not changed (For SDPR)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/MKL/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/gsl/lib

# prevent automatically using  all available cores on a compute node (For SDPR and PRS-CS)
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

# Load R to perform lassosum
module load R

# Start to run OTTERS
cd ${OTTERS_DIR}/Example

# Input for OTTERS STAGE I 
# Annotation File 
exp_anno=exp_anno.txt
# Genotype data from LD reference panel
geno_dir=Exp_geno
# eQTL summary statistics 
sst_file=Exp_eQTLSumStats.txt
# Input for OTTERS STAGE II
# GWAS summary statistics 
gwas_sst_file=Exp_GWASSumStats.txt

# Set chromosome number (The example is for chromosome X)
chr=4
# Set LD-clumping threshold in STAGE I
clump_r2=0.99
# Set output directory for STAGE I
out_dir=Results

# STAGE I
# train eQTL weights using P+T, lassosum, SDPR and PRS-CS. 
# It may take several minutes to complete.
python3 ${OTTERS_DIR}/training.py \
--OTTERS_dir=${OTTERS_DIR} \
--SDPR_dir=${SDPR_DIR} \
--anno_dir=${exp_anno} \
--geno_dir=${geno_dir} \
--sst_file=${sst_file} \
--out_dir=${out_dir} \
--chrom=${chr} \
--r2=${clump_r2} \
--models=PT,lassosum,SDPR,PRScs \
--lassosum_ld_blocks=EUR.hg38 \
--thread=$N_THREADS \
--geno_type=vcf

# Set output directory for STAGE II
twas_dir=TWAS

# STAGE II
# gene-based association test using eQTL-weight trained from P+T, lassosum, SDPR and PRS-CS.
python3 ${OTTERS_DIR}/testing.py \
--OTTERS_dir=${OTTERS_DIR} \
--weight_dir=${OTTERS_DIR}/Example/Results \
--models=P0.001,P0.05,lassosum,SDPR,PRScs \
--anno_dir=${exp_anno} \
--geno_dir=${geno_dir} \
--out_dir=${twas_dir} \
--gwas_file=${gwas_sst_file} \
--chrom=${chr} \
--thread=$N_THREADS

# get imputed genetically regulated gene expression
impute_dir=GReX
# samples to perform imputation
samples=ChromX_Exp_samples.txt
# imputation
python3 ${OTTERS_DIR}/imputing.py \
--OTTERS_dir=${OTTERS_DIR} \
--weight_dir=${OTTERS_DIR}/Example/ChromX/Results \
--models=P0.001,P0.05,lassosum,SDPR,PRScs \
--anno_dir=${exp_anno} \
--geno_dir=${geno_dir} \
--out_dir=${impute_dir} \
--chrom=${chr} \
--samples=${samples} \
--thread=$N_THREADS