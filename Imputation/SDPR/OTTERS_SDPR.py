#!/usr/bin/env python

"""

Usage:
python OTTERS_SDPR.py --OTTERS_dir=PATH_TO_OTTERS --in_dir=INPUR_DIR --out_dir=OUTPUT_DIR --chrom=CHROM
                [--out_file=OUTPUT_FILE --M=PARAM_M --opt_llk=LIKLIHOOD_EQ --iter=MCMC_ITERATIONS --burn=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR 
                 --r2=R2_FOR_BLOCK --a=PARAM_A --c=PARAM_C --a0k=PARAM_A0K --b0k=PARAM_B0K --thread=THREAD]

 - PATH_TO_OTTERS: The directory of OTTERS source code

 - INPUT_DIR:  Full path and the file name of the LD reference file.

 - OUTPUT_DIR: Output directory of the posterior effect size estimates.

 - OUTPUT_FILE (optional): Output filename of the posterior effect size estimates. Default is SDPR.txt.

 - CHROM: An integer for the chromosome of tested genes.

 - PARAM_M (optional): Max number of variance components. M must be greater than 4. Default is 1000.

 - LIKLIHOOD_EQ (optional): Which likelihood to evaluate. 1 for equation 6 (slightly shrink the correlation of SNPs) and 2 for equation 5 (SNPs genotyped on different arrays in a separate cohort). Please refer to manuscript or manual (3.2.3-2.3.4) for more details. Default is 1.

 - MCMC_ITERATIONS (optional): number of iterations for MCMC. Default is 1000.

 - MCMC_BURNIN  (optional): number of burn-in for MCMC. Default is 200.

 - MCMC_THINNING_FACTOR (optional): Thinning for MCMC. Default is 1 (no thin).

 - R2_FOR_BLOCK (optional): r2 cut-off for partition of independent blocks. Default is 0.1.

 - PARAM_A (optional): factor to shrink the reference LD matrix. Default is 0.1. Please refer to the manual for more information.

 - PARAM_C (optional): factor to correct for the deflation. Default is 1. Please refer to the manual for more information.

 - PARAM_A0K (optional): hyperparameter for inverse gamma distribution. Default is 0.5.

 - PARAM_B0K(optional): hyperparameter for inverse gamma distribution. Default is 0.5.

 - THREAD (optional): number of threads to use. Default is 1.

"""

import multiprocessing
import subprocess
import getopt
import sys
import os

import numpy as np
import pandas as pd

from time import time

############################################################
# time calculation
start_time = time()

############################################################


def parse_param():
    long_opts_list = ['OTTERS_dir=', 'SDPR_dir=', 'in_dir=',
                      'out_dir=', 'out_file=', 'chrom=', 'thread=', 'r2=',
                      'M=', 'opt_llk=', 'iter=', 'burn=', 'thin=',
                      'a=', 'c=', 'a0k=', 'b0k=']

    param_dict = {'OTTERS_dir': None, 'SDPR_dir': None,'in_dir': None,
                  'out_dir': None, 'out_file': 'SDPR.txt', 'chrom': None, 'thread': 1,
                  'r2': 0.1, 'M': 1000, 'opt_llk': 1, 'iter': 1000,
                  'burn': 200, 'thin': 1, 'a': 0.1, 'c': 1,
                  'a0k': 0.5, 'b0k': 0.5}

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
            elif opt == "--SDPR_dir":
                param_dict['SDPR_dir'] = arg
            elif opt == "--in_dir":
                param_dict['in_dir'] = arg
            elif opt == "--out_dir":
                param_dict['out_dir'] = arg
            elif opt == "--out_file":
                param_dict['out_file'] = arg
            elif opt == "--chrom":
                param_dict['chrom'] = str(arg)
            elif opt == "--thread":
                param_dict['thread'] = int(arg)
            elif opt == "--r2":
                param_dict['r2'] = float(arg)
            elif opt == "--M":
                param_dict['M'] = int(arg)
            elif opt == "--opt_llk":
                param_dict['opt_llk'] = int(arg)
            elif opt == "--iter":
                param_dict['iter'] = int(arg)
            elif opt == "--burn":
                param_dict['burn'] = int(arg)
            elif opt == "--thin":
                param_dict['thin'] = int(arg)
            elif opt == "--a":
                param_dict['a'] = float(arg)
            elif opt == "--c":
                param_dict['c'] = float(arg)
            elif opt == "--a0k":
                param_dict['a0k'] = float(arg)
            elif opt == "--b0k":
                param_dict['b0k'] = float(arg)
    else:
        print(__doc__)
        sys.exit(0)

    if param_dict['OTTERS_dir'] is None:
        print('* Please specify the directory to OTTERS --OTTERS_dir\n')
        sys.exit(2)
    elif param_dict['SDPR_dir'] is None:
        print('* Please specify the directory to SDPR --SDPR_dir\n')
        sys.exit(2)
    elif param_dict['in_dir'] is None:
        print('* Please specify the full path to input files --in_dir\n')
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
out_dir = os.path.join(param_dict['out_dir'], param_dict['out_file'])
ots.check_path(param_dict['out_dir'])

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

# create output file for SDPR
out_cols = ['CHROM', 'POS', 'A1', 'A2', 'TargetID', 'ES']
pd.DataFrame(columns=out_cols).to_csv(
    out_dir,
    sep='\t',
    index=None,
    header=True,
    mode='w')

# ############################################################

@ots.error_handler
def thread_process(num):
    target = TargetID[num]
    print('num=' + str(num) + '\nTargetID=' + target)
    target_n = N[num]
    target_dir = os.path.join(input_dir, target)

    # extract summary statistics and reference LD matrix for the gene
    print('Extracting eQTL summary statistics data...')
    sst_path = os.path.join(target_dir, target + '_Zscore.txt')
    sst_chunks = pd.read_csv(sst_path, sep='\t',
                             low_memory=False,
                             header=0,
                             names=['SNP', "A1", "A2", "Z"],
                             iterator=True, chunksize=1000,
                             dtype={'SNP': object, 'A1': object, 'A2': object, 'Z': np.float32})
    target_sst = pd.concat([chunk for chunk in sst_chunks]).reset_index(drop=True)

    if len(target_sst) == 0:
        print('There is no summary statistics data for TargetID: ' + target)
        print('Please use OTTERS to prepare summary statistics \n')
        return None

    print('Generating reference LD matrix...')
    target_ref_dir = os.path.join(target_dir, "ref")
    ots.check_path(target_ref_dir)

    SDPR_path = os.path.join(param_dict['SDPR_dir'], 'SDPR')

    try:
        SDPR_LD_args = [SDPR_path +
                        ' -make_ref -ref_prefix ' + target +
                        ' -chr ' + str(chrom) +
                        ' -r2 ' + str(param_dict['r2']) +
                        ' -ref_dir ref/']
        proc = subprocess.check_call(SDPR_LD_args,
                                     stdout=subprocess.PIPE,
                                     cwd=target_dir,
                                     shell=True)
        print('Done generating LD.')
    except subprocess.CalledProcessError as err:
        print('SDPR failed to generate LD for TargetID: ' + target + '\n')
        return None

    print("Start running SDPR...")

    try:
        SDPR_mcmc_args = [SDPR_path +
                          ' -mcmc -ref_dir ref/ ' +
                          ' -ss ' + target+'_Zscore.txt' +
                          ' -N ' + str(int(target_n)) +
                          ' -chr ' + str(chrom) +
                          ' -M ' + str(param_dict['M']) +
                          ' -opt_llk ' + str(param_dict['opt_llk']) +
                          ' -iter ' + str(param_dict['iter']) +
                          ' -burn ' + str(param_dict['burn']) +
                          ' -thin ' + str(param_dict['thin']) +
                          ' -a ' + str(param_dict['a']) +
                          ' -c ' + str(param_dict['c']) +
                          ' -a0k ' + str(param_dict['a0k']) +
                          ' -b0k ' + str(param_dict['b0k']) +
                          ' -out ' + target + '_SDPR.txt']
        proc = subprocess.check_call(SDPR_mcmc_args,
                                     stdout=subprocess.PIPE,
                                     cwd=target_dir,
                                     shell=True)
        print('Done training SDPR for TargetID: ' + target + '.')
    except subprocess.CalledProcessError:
        print('SDPR failed for TargetID: ' + target + '\n')
        return None

    # save the SDPR results.
    SDPR_out_dir = os.path.join(target_dir, target + '_SDPR.txt')
    SDPR_chunks = pd.read_csv(SDPR_out_dir, sep='\t',
                              low_memory=False,
                              header=0,
                              names=["SNP", "A1", "ES"],
                              iterator=True,
                              chunksize=1000)

    SDPR_target = pd.concat([chunk for chunk in SDPR_chunks]).reset_index(drop=True)

    # remove the orginal SDPR outputs and LD
    os.remove(SDPR_out_dir)
    os.remove(os.path.join(target_ref_dir, 'chr' + chrom + '.dat'))

    SDPR_out = SDPR_target.merge(target_sst, left_on=['SNP', 'A1'],
                                 right_on=['SNP', 'A1'], how="inner")

    SDPR_out['CHROM'] = chrom
    SDPR_out['TargetID'] = target
    SDPR_out['POS'] = [int(snp.split('_')[1]) for snp in target_sst['SNP']]
    SDPR_out[['CHROM', 'POS', 'A1', 'A2', 'TargetID', 'ES']].to_csv(
        out_dir,
        sep='\t',
        index=None,
        header=None,
        mode='a')

    print('Done saving SDPR results. \n')

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

