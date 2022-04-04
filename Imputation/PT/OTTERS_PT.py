#!/usr/bin/env python

"""
A TWAS method that infers posterior eQTL effect sizes under continuous shrinkage (CS) priors
using eQTL summary statistics and an external LD reference panel.

Usage:
python OTTERS_PRScs.py --OTTERS_dir=PATH_TO_OTTERS --in_dir=INPUR_DIR --out_dir=OUTPUT_DIR --chrom=CHROM
                [--pt=PT]

 - PATH_TO_OTTERS: The directory of OTTERS source code

 - INPUT_DIR:  Full path and the file name of the LD reference file.

 - OUTPUT_DIR: Output directory of the effect size estimates.

 - CHROM: An integer for the chromosome of tested genes.

 - PT(optional): The p-value threshold used in P+T, separated by comma, e.g., --pt=0.1,0.05,0.001. Default is 0.05.

 - THREAD (optional): Number of simultaneous processes to use for parallel computation. Default is 1.

"""

import multiprocessing
import getopt
import sys
import os

import numpy as np
import pandas as pd

from time import time
from scipy.stats.distributions import chi2

############################################################
# time calculation
start_time = time()

############################################################


def parse_param():
    long_opts_list = ['OTTERS_dir=', 'SDPR_dir=', 'in_dir=',
                      'out_dir=', 'chrom=', 'thread=', 'pt=']

    param_dict = {'OTTERS_dir': None, 'SDPR_dir': None,'in_dir': None,
                  'out_dir': None, 'chrom': None, 'thread': 1,
                  'pt': 0.05}

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
            elif opt == "--chrom":
                param_dict['chrom'] = str(arg)
            elif opt == "--thread":
                param_dict['thread'] = int(arg)
            elif opt == "--pt":
                param_dict['pt'] = arg.split(',')
    else:
        print(__doc__)
        sys.exit(0)

    if param_dict['OTTERS_dir'] is None:
        print('* Please specify the directory to OTTERS --OTTERS_dir\n')
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
out_dir = param_dict['out_dir']

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
for p in param_dict['pt']:
    PT_res = os.path.join(out_dir, 'P' + str(p) + '.txt')
    PT_cols = ['CHROM', 'POS', 'A1', 'A2', 'TargetID', 'ES']
    pd.DataFrame(columns=PT_cols).to_csv(
        PT_res,
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
        print('Please use OTTERS to prepare input files \n')
        return None

    target_sst['CHROM'] = chrom
    target_sst['TargetID'] = target
    target_sst['POS'] = [int(snp.split('_')[1]) for snp in target_sst['SNP']]
    target_sst['ES'] = target_sst['Z']/np.sqrt(target_n)
    target_sst['P'] = chi2.sf(np.power(target_sst['Z'], 2), 1)

    print("Start P+T")

    for p in param_dict['pt']:

        PT_res = os.path.join(param_dict['out_dir'], 'P' + str(p) + '.txt')

        PT = target_sst[target_sst.P < float(p)]

        PT[['CHROM', 'POS', 'A1', 'A2', 'TargetID', 'ES']].to_csv(
            PT_res,
            sep='\t',
            index=None,
            header=None,
            mode='a')

        print('Done training P+T for ' + str(num) + ':' + target + 'with p-value threshold: ' + p)

    print("Done P+T \n")

############################################################
if __name__ == '__main__':
    print('Starting summary for ' + str(n_targets) + ' target genes.\n')
    pool = multiprocessing.Pool(param_dict['thread'])
    pool.imap(thread_process, [num for num in range(3)])
    pool.close()
    pool.join()
    print('Done.')

############################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = ots.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)
