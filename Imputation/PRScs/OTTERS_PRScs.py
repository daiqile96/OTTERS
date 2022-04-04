#!/usr/bin/env python

"""
A TWAS method that infers posterior eQTL effect sizes under continuous shrinkage (CS) priors
using eQTL summary statistics and an external LD reference panel.

Usage:
python OTTERS_PRScs.py --OTTERS_dir=PATH_TO_OTTERS --in_dir=INPUR_DIR --out_dir=OUTPUT_DIR --chrom=CHROM
                [--window=WINDOW_SIZE --a=PARAM_A --b=PARAM_B --phi=PARAM_PHI --n_iter=MCMC_ITERATIONS --n_burnin=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR --thread=THREAD --seed=SEED]

 - PATH_TO_OTTERS: The directory of OTTERS source code

 - INPUT_DIR:  Full path and the file name of the LD reference file.

 - OUTPUT_DIR: Output directory and output filename of the posterior effect size estimates.

 - CHROM: An integer for the chromosome of tested genes.

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
import os

import numpy as np
import pandas as pd

from time import time
from numpy import linalg

import mcmc_gtb

############################################################
# time calculation
start_time = time()

############################################################


def parse_param():
    long_opts_list = ['OTTERS_dir=', 'in_dir=',
                      'out_dir=', 'chrom=', 'window=', 'thread=',
                      'a=', 'b=', 'phi=', 'n_iter=', 'n_burnin=', 'thin=',
                      'seed=', 'help']

    param_dict = {'OTTERS_dir': None, 'in_dir': None,
                  'out_dir': None, 'chrom': None,
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
            elif opt == "--in_dir":
                param_dict['in_dir'] = arg
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
    elif param_dict['in_dir'] is None:
        print('* Please specify the full path to the sample size file using --anno_dir\n')
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

chrom = param_dict['chrom']
input_dir = param_dict['in_dir']

print('Reading sample size file.\n')
N_chunks = pd.read_csv(
    os.path.join(input_dir, 'medianN.txt'),
    sep='\t',
    header=0,
    chunksize=10000,
    iterator=True,
    dtype={'CHROM': object, 'TargetID': object, 'N': np.int32})

GeneN = pd.concat([x[x['CHROM'] == chrom] for x in N_chunks]).reset_index(drop=True)

if GeneN.empty:
    raise SystemExit('There are no valid sample size data for chromosome ' + chrom + '\n')

TargetID = GeneN.TargetID
n_targets = TargetID.size
N = GeneN.N

# create output file for PRS-CS
out_cols = ['CHROM', 'POS', 'A1', 'A2', 'TargetID', 'ES']
pd.DataFrame(columns=out_cols).to_csv(
    param_dict['out_dir'],
    sep='\t',
    index=None,
    header=True,
    mode='w')

# ############################################################

@ots.error_handler
def thread_process(num):
    print('Reading eQTL summary statistics data.')

    target = TargetID[num]
    print('num=' + str(num) + '\nTargetID=' + target)
    target_n = N[num]
    target_dir = os.path.join(input_dir, target)

    # extract summary statistics 
    print('Extracting eQTL summary statistics data.')
    sst_path = os.path.join(target_dir, target + '_beta.txt')
    sst_chunks = pd.read_csv(sst_path, sep='\t',
                             low_memory=False,
                             header=0,
                             names=['SNP', "A1", "A2", "BETA"],
                             iterator=True, chunksize=1000,
                             dtype={'SNP': object, 'A1': object, 'A2': object, 'BETA': np.float32})
    target_sst = pd.concat([chunk for chunk in sst_chunks]).reset_index(drop=True)

    if len(target_sst) == 0:
        print('There is no summary statistics data for TargetID: ' + target)
        print('Please use OTTERS to prepare input files \n')
        return None

    # Save the list of SNPs that are kept in the summary statistics 
    target_snplist = target_sst['SNP']
    target_snplist_path = os.path.join(target_dir, target + '.snplist')

    target_snplist.to_csv(target_snplist_path,
                          sep='\t',
                          index=None,
                          header=None,
                          mode='w')

    # Calculate LD for these SNPs.
    print('Generating reference LD matrix.')
    snplist = target + '.snplist'
    get_ld_cmd = ["plink --bfile " + target + " --keep-allele-order --extract " +
                  snplist + " --r square --out " + target + " --memory 2000 "]

    try:
        ld_proc = subprocess.check_call(get_ld_cmd,
                                    stdout=subprocess.PIPE,
                                    cwd=target_dir,
                                    shell=True,
                                    bufsize=1)
        print('LD calculation completed')
    except subprocess.CalledProcessError as err:
        print('calculation failed for TargetID: ' + target + '\n')
        return None

    target_ld = os.path.join(target_dir, target + '.ld')

    ld_chunks = pd.read_csv(target_ld, sep='\t',
                            low_memory=False,
                            header=None,
                            iterator=True,
                            chunksize=1000)

    ld = pd.concat([chunk for chunk in ld_chunks]).reset_index(drop=True)

    # Make sure the ld is positive definite matrix
    # PLINK rounds to 6 decimal points, which sometimes makes the correlation matrix become not positive definite
    V = ld.to_numpy()
    # PLINK return nan when all mutually-nonmissing set is homozygous major
    # We set nan = 0 here 
    V = np.nan_to_num(V)
    _, s, v = linalg.svd(V)
    h = np.dot(v.T, np.dot(np.diag(s), v))
    V = (V+h)/2
    blk_size = len(V)

    if blk_size > 0:
        try:
            target_sst['ES'] = mcmc_gtb.mcmc(a=param_dict['a'], b=param_dict['b'],
                                             phi=param_dict['phi'],
                                             sst_dict=target_sst, n=target_n,
                                             ld_blk=V,
                                             blk_size=blk_size,
                                             n_iter=param_dict['n_iter'],
                                             n_burnin=param_dict['n_burnin'],
                                             thin=param_dict['thin'],
                                             seed=param_dict['seed'])
        except subprocess.CalledProcessError as err:
            print('PRS-CS failed for TargetID: ' + target + '\n')
            return None

    target_sst['POS'] = [int(snp.split('_')[1]) for snp in target_sst['SNP']]
    with open(param_dict['out_dir'], 'a') as ff:
        for pos, a1, a2, beta in zip(target_sst['POS'], target_sst['A1'], target_sst['A2'], target_sst['ES']):
            ff.write('%d\t%d\t%s\t%s\t%s\t%.6e\n' % (int(chrom), int(pos), a1, a2, target, beta))

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
