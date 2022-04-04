#!/usr/bin/env python

"""

Usage:
python OTTERS_imputation.py --OTTERS_dir=PATH_TO_OTTERS --anno_dir=PATH_TO_ANNO --geno_dir=PATH_TO_GENO --sst_dir=PATH_TO_SST --out_dir=OUTPUT_DIR --chrom=CHROM 
       [--r2=R2 --window=WINDOW --thread=THREAD --help]
                
 - OTTERS_DIR: The directory of OTTERS source code

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

"""

import multiprocessing
import subprocess
import getopt
import sys
import os

from io import StringIO
from time import time

import numpy as np
import pandas as pd
import shutil


############################################################
# time calculation
start_time = time()

############################################################


def parse_param():

    long_opts_list = ['OTTERS_dir=', 'anno_dir=', 'b_dir=',
                      'sst_dir=', 'out_dir=', 'chrom=', 'r2=',
                      'window=', 'thread=', 'help']

    param_dict = {'OTTERS_dir=': None, 'anno_dir=': None,
                  'b_dir': None, 'sst_dir': None, 'out_dir': None, 'chrom': None,
                  'r2': 2, 'window': 1000000, 'thread': 1}

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
            elif opt == "--anno_dir":
                param_dict['anno_dir'] = arg
            elif opt == "--b_dir":
                param_dict['b_dir'] = arg
            elif opt == "--sst_dir":
                param_dict['sst_dir'] = arg
            elif opt == "--out_dir":
                param_dict['out_dir'] = arg
            elif opt == "--chrom":
                param_dict['chrom'] = str(arg)
            elif opt == "--r2":
                param_dict['r2'] = float(arg)
            elif opt == "--window":
                param_dict['window'] = int(arg)
            elif opt == "--thread":
                param_dict['thread'] = int(arg)
    else:
        print(__doc__)
        sys.exit(0)

    if param_dict['OTTERS_dir'] is None:
        print('* Please specify the directory to OTTERS --OTTERS_dir\n')
        sys.exit(2)
    elif param_dict['anno_dir'] is None:
        print('* Please specify the directory to the gene annotation file using --anno_dir\n')
        sys.exit(2)
    elif param_dict['b_dir'] is None:
        print('* Please specify the directory to the binary file of LD reference panel --bim_dir\n')
        sys.exit(2)
    elif param_dict['sst_dir'] is None:
        print('* Please specify the eQTL summary statistics file using --sst_dir\n')
        sys.exit(2)
    elif param_dict['out_dir'] is None:
        print('* Please specify the output directory\n')
        sys.exit(2)
    elif param_dict['chrom'] is None:
        print('* Please specify the chromosome --chrom\n')
        sys.exit(2)

    for key in param_dict:
        print('--%s=%s' % (key, param_dict[key]))

    print('\n')
    return param_dict

############################################################

param_dict = parse_param()
sys.path.append(param_dict['OTTERS_dir'])
import OTTERSutils as ots

chr = param_dict['chrom']

print('Reading gene annotation data on chromosome ' + chr)
anno_chunks = pd.read_csv(
    param_dict['anno_dir'],
    sep='\t',
    header=0,
    chunksize=10000,
    iterator=True,
    dtype={'CHROM': object,
        'GeneEnd': np.int64,
        'GeneName': object,
        'GeneStart': np.int64,
        'TargetID': object})

GeneAnno = pd.concat([x[x['CHROM'] == chr] for x in anno_chunks]).reset_index(drop=True)

if GeneAnno.empty:
    raise SystemExit('There are no valid gene annotation data for chromosome ' + chr + '\n')

GeneAnno = ots.optimize_cols(GeneAnno)
TargetID = GeneAnno.TargetID
n_targets = TargetID.size

# create output file for median eQTL sample size for each gene
target_medianN = os.path.join(param_dict['out_dir'], 'medianN.txt')
medianN_cols = ['CHROM', 'TargetID', 'N']
pd.DataFrame(columns=medianN_cols).to_csv(
    target_medianN,
    sep='\t',
    index=None,
    header=True,
    mode='w')


############################################################

@ots.error_handler
def thread_process(num):

    target = TargetID[num]
    print('num=' + str(num) + '\nTargetID=' + target)
    target_exp = GeneAnno.iloc[[num]]

    start = str(max(int(target_exp.GeneStart) - param_dict['window'],0))
    end = str(int(target_exp.GeneEnd) + param_dict['window'])

    ################# PLINK Binary Files #####################

    print('Making PLINK binary files for the target gene')

    # set output path of the target gene
    target_b_dir = os.path.join(param_dict['out_dir'], target)
    ots.check_path(target_b_dir)

    # generate command to call PLINK to extract the binary file for the target gene
    extract_cmd = ots.call_PLINK_extract(param_dict['b_dir'], target_b_dir, target, chr, start, end)

    try:
        proc = subprocess.check_call(extract_cmd,
                                     stdout=subprocess.PIPE,
                                     shell=True,
                                     bufsize=1)
        print('Done making binary files')
    except subprocess.CalledProcessError:
        print('There is no genotype data for TargetID: ' + target)
        shutil.rmtree(target_b_dir)
        return None

    ################# eQTL summary statistics #####################

    print('Reading eQTL summary statistics data.')

    sst_proc_out = ots.call_tabix(param_dict['sst_dir'], chr, start, end)

    if not sst_proc_out:
        print('There is no summary statistics data for the range of ' + target + '\n')
        print('Remove binary files for TargetID: ' + target + '.\n')
        shutil.rmtree(target_b_dir)
        return None

    sst_chunks = pd.read_csv(StringIO(sst_proc_out.decode('utf-8')), sep='\t',
                             low_memory=False,
                             header=None,
                             names=["CHROM", "POS", "SNP", "A1", "A2", "Z", "P", "TargetID", "Gene", "NrSamples"],
                             iterator=True,
                             chunksize=1000,
                             dtype={'P': np.float64, 'NrSamples': np.int32})

    target_sst = pd.concat([chunk[chunk.TargetID == target] for chunk in sst_chunks]).reset_index(drop=True)

    if len(target_sst) == 0:
        print('There is no summary statistics data for TargetID: ' + target)
        print('Remove binary files for TargetID: ' + target + '.\n')
        shutil.rmtree(target_b_dir)
        return None

    print('Check overlap between eQTL summary stat and LD reference.')

    # read in snps in LD reference panel
    target_bim = os.path.join(target_b_dir, target + '.bim')
    ref_chunks = pd.read_csv(target_bim, sep='\t',
                             low_memory=False,
                             header=None,
                             names=["CHROM", "SNP", "bp", "POS", "A1", "A2"],
                             iterator=True,
                             chunksize=1000,
                             usecols=["CHROM", "POS", "A1", "A2"])

    target_ref = pd.concat([chunk for chunk in ref_chunks]).reset_index(drop=True)

    target_ref['snpID'] = ots.get_snpIDs(target_ref, flip=False)
    target_sst['snpID'] = ots.get_snpIDs(target_sst, flip=False)
    target_sst['snpIDflip'] = ots.get_snpIDs(target_sst, flip=True)

    snp_overlap = np.intersect1d(target_ref.snpID, target_sst[['snpID', 'snpIDflip']])

    if not snp_overlap.size:
        print('No overlapping test eQTLs between eQTL summary statistics and LD reference panel for TargetID: ' + target)
        print('Remove binary files for TargetID: ' + target + '.\n')
        shutil.rmtree('target_b_dir')
        return None
    else:
        print('%s overlapped eQTLs' % (len(snp_overlap)))

    # filter out non-matching snpID rows
    target_sst = target_sst[np.any(target_sst[['snpID', 'snpIDflip']].isin(snp_overlap), axis=1)].reset_index(drop=True)

    # if not in target_ref.snpID, assumed flipped; if flipped, flip Zscore sign
    target_sst['flip'] = np.where(target_sst.snpID.isin(target_ref.snpID.values), 1, -1)
    if not np.all(target_sst['flip'] == 1):
        idx = (target_sst['flip'] == -1)
        target_sst.loc[idx, ['snpID']] = target_sst.loc[idx, ['snpIDflip']].values
        target_sst.loc[idx, ['A1', 'A2']] = target_sst.loc[idx, ['A2', 'A1']].values
        target_sst['Z'] = target_sst['flip'] * target_sst['Z']

    ################# Calculate median sample size ####################

    print('Calculate median sample size of eQTLs')
    median_N = np.nanmedian(target_sst['NrSamples'])

    with open(target_medianN, 'a') as ff:
        ff.write('%s\t%s\t%d\n' % (chr, target, median_N))

    ################# LD clumping #############################
    print('Perform LD clumping')

    # generate summary statistics of p-value to perform LD-clumping
    target_p = os.path.join(target_b_dir, target+'.pvalue')
    target_sst[['snpID', 'A1', 'A2', 'P']].to_csv(
        target_p,
        sep='\t',
        index=None,
        header=True,
        mode='w')

    # set the path and prefix of the PLINK binary genotype files of the target gene
    target_b = os.path.join(target_b_dir, target)

    # use call_PLINK_clump to generate command to call PLINK to perform LD-clumping 
    # required input: 
    # path_to_binary_reference_genotype=target_b, clumping_r2_threshold=r2, path_to_summary_stat_pvalue=target_p
    # output: os.path.join(target_b_dir, target+'.clumped')
    clump_cmd = ots.call_PLINK_clump(target_b, param_dict['r2'], target_p)

    try:
        proc = subprocess.check_call(clump_cmd,
                                     stdout=subprocess.PIPE,
                                     shell=True,
                                     bufsize=1)
        print('Done Clumping')
    except subprocess.CalledProcessError:
        print('Clumping Failed for TargetID: ' + target + '\n')
        return None

    # read in remaining eQTLs after LD-clumping
    target_clumped = os.path.join(target_b_dir, target+'.clumped')
    with open(target_clumped) as file_in:
        lines = []
        clumped_snp = []
        for line in file_in:
            lines.append(line)
            res = line.split()
            if len(res) != 0:
                if res[0] == chr:
                    clumped_snp.append(res[2])

    # filter summary statistics by clumping results
    target_sst = target_sst[np.any(target_sst[['snpID']].isin(clumped_snp), axis=1)].reset_index(drop=True)

    ################# Prepare Input Summary Statistics for Imputation Models #############################

    print('Start prepare input summary statistics')

    # Prepare Zscore input for SDPR
    target_zscore = os.path.join(target_b_dir, target+'_Zscore.txt')
    target_sst[['snpID', 'A1', 'A2', 'Z']].to_csv(
        target_zscore,
        sep='\t',
        index=None,
        header=['SNP', "A1", "A2", "Z"],
        mode='w')
    print('Done generating Zscore')

    # Prepare the standardized beta input for lassosum and P+T
    target_beta = os.path.join(target_b_dir, target+'_beta.txt')
    target_sst['Beta'] = target_sst['Z']/np.sqrt(median_N)
    target_sst[['snpID', 'A1', 'A2', 'Beta']].to_csv(
        target_beta,
        sep='\t',
        index=None,
        header=['SNP', "A1", "A2", "Beta"],
        mode='w')
    print('Done generating standardized beta')
    print('Done prepare inputs for ' + str(num) + ':' + target + '\n')

   
############################################################
if __name__ == '__main__':
    print('Starting prepare inputs for ' + str(n_targets) + ' target genes.\n')
    pool = multiprocessing.Pool(param_dict['thread'])
    pool.imap(thread_process, [num for num in range(10)])
    pool.close()
    pool.join()
    print('Done.')

############################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = ots.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)
