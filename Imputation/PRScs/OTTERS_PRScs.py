#!/usr/bin/env python

"""
A TWAS method that infers posterior eQTL effect sizes under continuous shrinkage (CS) priors
using eQTL summary statistics and an external LD reference panel.

Usage:
python OTTERS_PRScs.py --OTTERS_dir=PATH_TO_OTTERS --N_path=PATH_TO_N --ld_dir=PATH_TO_LD --sst_file=SUM_STATS_FILE --out_dir=OUTPUT_DIR --chrom=CHROM
                [--window=WINDOW_SIZE --a=PARAM_A --b=PARAM_B --phi=PARAM_PHI --n_iter=MCMC_ITERATIONS --n_burnin=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR --thread=THREAD --seed=SEED]

 - PATH_TO_OTTERS: The directory of OTTERS source code

 - PATH_TO_N: Full path and the file name of the gene sample size file.

 - PATH_TO_LD:  Full path and the file name of the LD reference file.

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

"""

import multiprocessing
import subprocess
import getopt
import sys

from io import StringIO
from time import time
from numpy import linalg

import os
import numpy as np
import pandas as pd


import mcmc_gtb

############################################################
# time calculation
start_time = time()

############################################################


def parse_param():
    long_opts_list = ['OTTERS_dir=', 'N_path=', 'ld_dir=', 'sst_dir=',
                      'out_dir=', 'chrom=', 'window=', 'thread=',
                      'a=', 'b=', 'phi=', 'n_iter=', 'n_burnin=', 'thin=',
                      'seed=', 'help']

    param_dict = {'OTTERS_dir': None, 'N_path': None, 'ld_dir': None,
                  'sst_dir': None, 'out_dir': None, 'chrom': None,
                  'window': 1000000, 'thread': 1, 'a': 1, 'b': 0.5,
                  'phi': None, 'n_iter': 1000, 'n_burnin': 500, 'thin': 5,
                  'seed': None}

    print('\n')

    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_opts_list)
        except:
            print('Option not recognized.')
            print('Use --help for usage information.\n')
            sys.exit(2)

        for opt, arg in opts:
            if opt == "-h" or opt == "--help":
                print(__doc__)
                sys.exit(0)
            elif opt == "--OTTERS_dir":
                param_dict['OTTERS_dir'] = arg
            elif opt == "--N_path":
                param_dict['N_path'] = arg
            elif opt == "--ld_dir":
                param_dict['ld_dir'] = arg
            elif opt == "--sst_dir":
                param_dict['sst_dir'] = arg
            elif opt == "--out_dir":
                param_dict['out_dir'] = arg
            elif opt == "--chrom":
                param_dict['chrom'] = str(arg)
            elif opt == "--window":
                param_dict['window'] = int(arg)
            elif opt == "--thread":
                param_dict['thread'] = int(arg)
            elif opt == "--a":
                param_dict['a'] = float(arg)
            elif opt == "--b":
                param_dict['b'] = float(arg)
            elif opt == "--phi":
                param_dict['phi'] = float(arg)
            elif opt == "--n_iter":
                param_dict['n_iter'] = int(arg)
            elif opt == "--n_burnin":
                param_dict['n_burnin'] = int(arg)
            elif opt == "--thin":
                param_dict['thin'] = int(arg)
            elif opt == "--thread":
                param_dict['thread'] = int(arg)
            elif opt == "--seed":
                param_dict['seed'] = int(arg)
    else:
        print(__doc__)
        sys.exit(0)

    if param_dict['OTTERS_dir'] is None:
        print('* Please specify the directory to OTTERS --OTTERS_dir\n')
        sys.exit(2)
    elif param_dict['N_path'] is None:
        print('* Please specify the full path to the sample size file using --anno_dir\n')
        sys.exit(2)
    elif param_dict['ld_dir'] is None:
        print('* Please specify the directory to the LD reference file using --ld_dir\n')
        sys.exit(2)
    elif param_dict['sst_dir'] is None:
        print('* Please specify the eQTL summary statistics file using --sst_dir\n')
        sys.exit(2)
    elif param_dict['chrom'] is None:
        print('* Please specify the chromosome --chrom\n')
        sys.exit(2)
    elif param_dict['out_dir'] is None:
        print('* Please specify the output directory using --out_dir\n')
        sys.exit(2)

    for key in param_dict:
        print('--%s=%s' % (key, param_dict[key]))

    print('\n')
    return param_dict


############################################################

param_dict = parse_param()
sys.path.append(param_dict['OTTERS_dir'])
import OTTERSutils as ots

print('Reading sample size file.\n')
N_chunks = pd.read_csv(
    param_dict['N_path'],
    sep='\t',
    header=0,
    chunksize=10000,
    iterator=True,
    dtype={'TargetID': object, 'N': np.int32})

GeneN = pd.concat([x[x['CHROM'] == param_dict['chrom']] for x in N_chunks]).reset_index(drop=True)

if GeneN.empty:
    raise SystemExit('There are no valid sample size data for chromosome ' + param_dict['chrom'] + '\n')

TargetID = GeneN.TargetID
n_targets = TargetID.size
N = GeneN.N


############################################################

@ots.error_handler
def thread_process(num):
    print('Reading eQTL summary statistics data.')
    target = TargetID[num]
    target_n = N[num]
    print('num=' + str(num) + '\nTargetID=' + target)

    # extract summary statistics and reference LD matrix for the gene
    print('Extracting eQTL summary statistics data.')

    sst_path = os.path.join(param_dict['sst_dir'], target, target+'_beta.txt')

    sst_chunks = pd.read_csv(sst_path, sep='\t',
                             low_memory=False,
                             header=None,
                             names=['SNP', "A1", "A2", "Beta"],
                             iterator=True, chunksize=1000)

    target_sst = pd.concat([chunk for chunk in sst_chunks]).reset_index(drop=True)

    if len(target_sst) == 0:
        print('There is no summary statistics data for TargetID: ' + target + '\n')
        return None

    print('Extracting reference LD matrix.')

    # construct covariance matrix
    snp_search_ids = target_sst['SNP']
    ld_blk = ots.get_ld_data(param_dict['ld_dir'], snp_search_ids)

    if len(target_sst['SNP']) > len(ld_blk['SNP']):
        print('%d eQTLs does not have genotype data in LD reference panel.' % (len(target_sst['SNP'] not in (ld_blk['SNP']))))
        print('Please consider use OTTERS to prepare summary statistics')
        target_sst = target_sst[target_sst['SNP'].isin(ld_blk['SNP'])]

    blk_size, V_cov = ots.get_ld_matrix(ld_blk)

    if blk_size > 0:
        try:
            target_sst['ES'] = mcmc_gtb.mcmc(a=param_dict['a'], b=param_dict['b'],
                                             phi=param_dict['phi'],
                                             sst_dict=target_sst, n=target_n,
                                             ld_blk=V_cov,
                                             blk_size=blk_size,
                                             n_iter=param_dict['n_iter'],
                                             n_burnin=param_dict['n_burnin'],
                                             thin=param_dict['thin'],
                                             seed=param_dict['seed'])
        except subprocess.CalledProcessError as err:
            print('PRS-CS failed for TargetID: ' + target + '\n')
            return None

    if param_dict['phi'] is None:
        target_weights_path = param_dict['out_dir'] + 'eQTL_weights_a%d_b%.1f_phiauto_chr%d.txt' % (param_dict['a'], param_dict['b'], int(param_dict['chrom']))
    else:
        target_weights_path = param_dict['out_dir'] + 'eQTL_weights_a%d_b%.1f_phi%1.0e_chr%d.txt' % (param_dict['a'], param_dict['b'], param_dict['phi'], int(param_dict['chrom']))

    with open(target_weights_path, 'a') as ff:
        for pos, a1, a2, beta in zip(target_sst['SNPPos'], target_sst['A1'], target_sst['A2'], target_sst['ES']):
            ff.write('%d\t%d\t%s\t%s\t%s\t%.6e\n' % (int(param_dict['chrom']), int(pos), a1, a2, target, beta))

    print('Done ' + str(num) + ':' + target + '\n')

############################################################
if __name__ == '__main__':
    print('Starting summary for ' + str(n_targets) + ' target genes.\n')
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