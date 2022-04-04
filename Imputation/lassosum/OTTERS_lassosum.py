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

 - THREAD (optional): Number of simultaneous processes to use for parallel computation. Default is 1.

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
    long_opts_list = ['OTTERS_dir=', 'in_dir=',
                      'out_dir=', 'chrom=', 'thread=', 'help']

    param_dict = {'OTTERS_dir': None, 'in_dir': None,
                  'out_dir': None, 'chrom': None,
                  'thread': 1}

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
            elif opt == "--thread":
                param_dict['thread'] = int(arg)
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
out_dir = param_dict['out_dir']
lassosum_path = os.path.join(param_dict['OTTERS_dir'], 'Imputation', 'lassosum', 'OTTERS_lassosum.R')

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
    out_dir,
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

    print("Start lassosum")

    try:
        lassosum_arg = ['Rscript ' + lassosum_path +
                        ' --gene_name=' + target +
                        ' --medianN=' + str(int(target_n)) +
                        ' --bim_file=' + target +
                        ' --sst_file=' + target + '_beta.txt' +
                        ' --out_path=' + out_dir +
                        ' --chr=' + str(chrom) +
                        ' --LDblocks=EUR.hg38']

        proc = subprocess.check_call(lassosum_arg,
                                     stdout=subprocess.PIPE,
                                     cwd=target_dir,
                                     shell=True,
                                     bufsize=1)
    except subprocess.CalledProcessError:
            print('SDPR failed for TargetID: ' + target)
            return None

    print('Done training lassosum for ' + str(num) + ':' + target)

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
