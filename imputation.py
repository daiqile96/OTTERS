#!/usr/bin/env python

"""

Usage:
python imputation.py --OTTERS_dir=PATH_TO_OTTERS --anno_dir=PATH_TO_ANNO --geno_dir=PATH_TO_GENO --sst_file=PATH_TO_SST --out_dir=OUTPUT_DIR --chrom=CHROM 
       [--r2=R2 --window=WINDOW --thread=THREAD 

        --prscs_a=prscs_PARAM_A --prscs_b=prscs_PARAM_B --prscs_phi=prscs_PARAM_PHI 
        --prscs_n_iter=prscs_MCMC_ITERATIONS --prscs_n_burnin=prscs_MCMC_BURNIN 
        --prscs_thin=prscs_MCMC_THINNING_FACTOR

        --SDPR_dir=PATH_TO_SDPR --SDPR_M=SDPR_PARAM_M
        --SDPR_opt_llk=SDPR_LIKLIHOOD_EQ --SDPR_iter=SDPR_MCMC_ITERATIONS 
        --SDPR_burn=SDPR_MCMC_BURNIN --SDPR_thin=SDPR_MCMC_THINNING_FACTOR 
        --SDPR_r2=SDPR_R2_FOR_BLOCK
        --SDPR_a=SDPR_PARAM_A --SDPR_c=SDPR_PARAM_C
        --SDPR_a0k=SDPR_PARAM_A0K --SDPR_b0k=SDPR_PARAM_B0K 

        --help]

 - PATH_TO_OTTERS: The directory of OTTERS source code

 - PATH_TO_ANNO: Full path and the file name of the gene annotation file. The annotation file is assumed to be in this format:

    | CHROM | GeneStart | GeneEnd |     TargetID    | 
    |:-----:|:---------:|:-------:|:---------------:|
    |   1   |    100    |   200   |     ENSG0000    | 

 - PATH_TO_GENO:  The directory and preflix of PLINK binary files for genotype data.

 - PATH_TO_SST: Full path and the file name of the summary statistics. 

    | CHROM | POS | A1 | A2 | Zscore |  TargetID  |   N  |
    |:-----:|:---:|:--:|:--:|:------:|:----------:|:----:|
    |   1   | 100 |  C |  T |   3    |  ENSG0000  |  0.2 |

 - OUTPUT_DIR: Output directory

 - CHROM: An integer for the chromosome of tested genes. ***We will parpare inputs for all the genes in the annotation file that are on this chromosome. ***

 - R2: The R squre threshold used to perform LD-clumping.

 - WINDOW (optional): Window size (in base pairs) around gene region from which to include SNPs (default: 1000000 [+- 1MB region around gene region])

 - THREAD (optional): Number of simultaneous processes to use for parallel computation. Default is 1.

 
 - prscs_PARAM_A (optional): Parameter a in the gamma-gamma prior in the continous shrinkage prior proposed by PRS-CS. Default is 1. 

 - prscs_PARAM_B (optional): Parameter b in the gamma-gamma prior in the continous shrinkage prior proposed by PRS-CS. Default is 0.5.

 - prscs_PARAM_PHI (optional): Global shrinkage parameter phi. If phi is not specified, it will be learnt from the data using a fully Bayesian approach.
                         This usually works well for polygenic traits with large GWAS sample sizes (hundreds of thousands of subjects).
                         For GWAS with limited sample sizes (including most of the current disease GWAS), fixing phi to 1e-4 or 1e-2,
                         or doing a small-scale grid search (e.g., phi=1e-6, 1e-4, 1e-2, 1) to find the optimal phi value often improves perdictive performance.

 - prscs_MCMC_ITERATIONS (optional): Total number of MCMC iterations. Default is 1,000.

 - prscs_MCMC_BURNIN (optional): Number of burnin iterations. Default is 500.

 - prscs_MCMC_THINNING_FACTOR (optional): Thinning of the Markov chain. Default is 5.

 - PATH_TO_SDPR (optional): The directory of SDPR binary tool.

 - SDPR_PARAM_M (optional): Max number of variance components. M must be greater than 4. Default is 1000.

 - SDPR_LIKLIHOOD_EQ (optional): Which likelihood to evaluate. 1 for equation 6 (slightly shrink the correlation of SNPs) and 2 for equation 5 (SNPs genotyped on different arrays in a separate cohort). Please refer to manuscript or manual (3.2.3-2.3.4) for more details. Default is 1.

 - SDPR_MCMC_ITERATIONS (optional): number of iterations for MCMC. Default is 1000.

 - SDPR_MCMC_BURNIN  (optional): number of burn-in for MCMC. Default is 200.

 - SDPR_MCMC_THINNING_FACTOR (optional): Thinning for MCMC. Default is 1 (no thin).

 - SDPR_R2_FOR_BLOCK (optional): r2 cut-off for partition of independent blocks. Default is 0.1.

 - SDPR_PARAM_A (optional): factor to shrink the reference LD matrix. Default is 0.1. Please refer to the manual for more information.

 - SDPR_PARAM_C (optional): factor to correct for the deflation. Default is 1. Please refer to the manual for more information.

 - SDPR_PARAM_A0K (optional): hyperparameter for inverse gamma distribution. Default is 0.5.

 - SDPR_PARAM_B0K(optional): hyperparameter for inverse gamma distribution. Default is 0.5.

"""

import multiprocessing
import numpy as np
import subprocess
import shutil
import sys
import os

from time import time
from scipy.stats.distributions import chi2

import ottersutils as ots
import prepare as prep

############################################################
# time calculation
start_time = time()

############################################################

# Parse input arguments
param_dict = ots.parse_param()

# Read in gene annotation file
GeneAnno, TargetID, n_targets = ots.read_anno(anno_dir=param_dict['anno_dir'],
                                              chrom=param_dict['chrom'])
# Create directory for output files
out_dir = param_dict['out_dir']
ots.check_path(out_dir)

# Create output files and set tool paths
out_cols = ['CHROM', 'POS', 'A1', 'A2', 'TargetID', 'ES']
for model in param_dict['models']:
    if model == 'PT':
        for p in param_dict['pt']:
            ots.create_file_title(out_cols, out_dir, 'P' + str(p) + '.txt')
    elif model in ['SDPR', 'PRScs', 'lassosum']:
        ots.create_file_title(out_cols, out_dir, model + '.txt')
        if model == 'SDPR':
            SDPR_path = os.path.join(param_dict['SDPR_dir'], 'SDPR')
        elif model == 'PRScs':
            sys.path.append(os.path.join(param_dict['OTTERS_dir'], 'PRSmodels'))
            import mcmc_gtb
        elif model == 'lassosum':
            lassosum_path = os.path.join(param_dict['OTTERS_dir'], 'PRSmodels/lassosum.R')
    else:
        print('Please specify models in PT,PRScs,SDPR,lassosum')


############################################################
@ots.error_handler
def thread_process(num):

    target = TargetID[num]
    print('num=' + str(num) + ' TargetID=' + target)
    target_anno = GeneAnno.iloc[[num]]

    target_dir, target_sst, median_N = prep.prepare(target=target,
                                                    target_anno=target_anno,
                                                    chrom=param_dict['chrom'],
                                                    window=param_dict['window'],
                                                    geno_dir=param_dict['geno_dir'],
                                                    out_dir=out_dir,
                                                    sst_dir=param_dict['sst_file'],
                                                    clump_r2=param_dict['r2'])

    ################# Start Imputation Models #############################

    print('...Train eQTL weights...')

    # P+T #
    if 'PT' in param_dict['models']:
        # print('*Start P+T...*')
        target_sst['P'] = chi2.sf(np.power(target_sst['Z'], 2), 1)

        for p in param_dict['pt']:

            PT_out = target_sst[target_sst.P < float(p)]

            ots.save_results_PT(p=p, out_df=PT_out, out_dir=out_dir)

    # PRS-CS #
    if 'PRScs' in param_dict['models']:

        # generate LD reference for PRS-CS
        V, blk_size = ots.PRScs_LD_cmd(bim_dir=target,
                                       sst_df=target_sst,
                                       out_file=target,
                                       work_dir=target_dir)

        # start to run PRS-CS MCMC
        if blk_size > 0:
            try:
                target_sst['PRScs'] = mcmc_gtb.mcmc(a=param_dict['prscs_a'],
                                                    b=param_dict['prscs_b'],
                                                    phi=param_dict['prscs_phi'],
                                                    sst_dict=target_sst[['snpID', 'Beta']],
                                                    n=median_N,
                                                    ld_blk=V,
                                                    blk_size=blk_size,
                                                    n_iter=param_dict['prscs_n_iter'],
                                                    n_burnin=param_dict['prscs_n_burnin'],
                                                    thin=param_dict['prscs_thin'],
                                                    seed=param_dict['seed'])
            except subprocess.CalledProcessError:
                print('PRS-CS failed.')
                return None

        # save PRS-CS results
        ots.save_results(model='PRScs', out_df=target_sst, out_dir=out_dir)

    # ############################# SDPR ###############################
    if 'SDPR' in param_dict['models']:

        # generate LD reference for SDPR
        ots.SDPR_LD_args(SDPR_path=SDPR_path,
                         chrom=param_dict['chrom'],
                         r2=param_dict['SDPR_r2'],
                         bim_dir=target,
                         work_dir=target_dir)

        # start to run SDPR MCMC
        try:
            SDPR_mcmc_args = ots.SDPR_args(SDPR_path=SDPR_path,
                                           sst_file=target + '_Zscore.txt',
                                           N=median_N,
                                           chrom=param_dict['chrom'],
                                           M=param_dict['SDPR_M'],
                                           opt_llk=param_dict['SDPR_opt_llk'],
                                           iter=param_dict['SDPR_iter'],
                                           burn=param_dict['SDPR_burn'],
                                           thin=param_dict['SDPR_thin'],
                                           a=param_dict['SDPR_a'],
                                           c=param_dict['SDPR_c'],
                                           a0k=param_dict['SDPR_a0k'],
                                           b0k=param_dict['SDPR_b0k'])

            proc = subprocess.check_call(SDPR_mcmc_args,
                                         stdout=subprocess.PIPE,
                                         cwd=target_dir,
                                         shell=True)
        except subprocess.CalledProcessError:
            print('SDPR failed')
            return None

        # save SDPR results
        ots.format_save_results(work_dir=target_dir,
                                out_dir=out_dir,
                                model='SDPR',
                                sst_df=target_sst)

    ############################# Lassosum ###############################
    if 'lassosum' in param_dict['models']:
        # print("*Start lassosum...*")
        try:
            lassosum_arg = ots.lassosum_cmd(param_dict['chrom'], target, target + '_beta.txt', 'lassosum.txt',
                                            lassosum_path, int(median_N))
            proc = subprocess.check_call(lassosum_arg,
                                         stdout=subprocess.PIPE,
                                         cwd=target_dir,
                                         shell=True)
        except subprocess.CalledProcessError:
                print('lassosum failed for TargetID: ' + target)
                return None

        # save lassosum results
        ots.format_save_results(work_dir=target_dir,
                                out_dir=out_dir,
                                model='lassosum',
                                sst_df=target_sst)

    ############################ Clean temporary #########################
    print('Remove temporary files. \n')
    shutil.rmtree(target_dir)

############################################################
if __name__ == '__main__':
    print('Start train eQTL weights for ' + str(n_targets) + ' target genes.\n')
    pool = multiprocessing.Pool(param_dict['thread'])
    pool.imap(thread_process, [num for num in range(n_targets)])
    pool.close()
    pool.join()
    print('Done.')

############################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = ots.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)
