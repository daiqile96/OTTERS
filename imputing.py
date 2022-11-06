#!/usr/bin/env python

"""

Usage:
python TWAS.py --OTTERS_dir=PATH_TO_OTTERS --anno_dir=PATH_TO_ANNO --geno_dir=PATH_TO_GENO --weight_file=PATH_TO_WEIGHT --gwas_file=PATH_TO_GWAS --out_dir=OUTPUT_DIR --chrom=CHROM
               [--window=WINDOW --thread=THREAD]

 - PATH_TO_OTTERS: The directory of OTTERS source code

 - PATH_TO_ANNO: Full path and the file name of the gene annotation file. The annotation file is assumed to be in this format:

    | CHROM | GeneStart | GeneEnd |     TargetID    | 
    |:-----:|:---------:|:-------:|:---------------:|
    |   1   |    100    |   200   |     ENSG0000    | 

 - PATH_TO_GENO:  The directory and preflix of PLINK binary files for genotype data.

 - PATH_TO_WEIGHT: Full path and the file name of estimated eQTL weights.

    | CHROM | POS | A1 | A2 | TargetID  |  ES   |
    |:-----:|:---:|:--:|:--:|:---------:|:-----:|
    |   1   | 100 |  C |  T | ENSG0000  |  0.2  |

 - PATH_TO_GWAS: Full path and the file name of GWAS summary statistics.

    | CHROM | POS | A1 | A2 | TargetID  |  Z   |
    |:-----:|:---:|:--:|:--:|:---------:|:-----:|
    |   1   | 100 |  C |  T | ENSG0000  |  2    |

 - OUTPUT_DIR: Output directory

 - CHROM: An integer for the chromosome of tested genes. ***We will parpare inputs for all the genes in the annotation file that are on this chromosome. ***

 - WINDOW (optional): Window size (in base pairs) around gene region from which to include SNPs (default: 1000000 [+- 1MB region around gene region])

 - THREAD (optional): Number of simultaneous processes to use for parallel computation. Default is 1.


"""
import multiprocessing
import subprocess
import getopt
import shutil
import pysam
import pandas as pd
import sys
import os


from time import time


############################################################
# time calculation
start_time = time()

############################################################


def parse_param():

    long_opts_list = ['OTTERS_dir=', 'anno_dir=', 'geno_dir=',
                      'weight_dir=', 'out_dir=', 'chrom=',
                      'window=', 'thread=', 'models=', 'samples=', 'help']

    param_dict = {'OTTERS_dir': None, 'anno_dir': None,
                  'geno_dir': None, 'weight_dir': None, 'out_dir': None, 'chrom': None,
                  'window': 1000000, 'thread': 1, 'models': None, 'samples': None}

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
            elif opt == "--geno_dir":
                param_dict['geno_dir'] = arg
            elif opt == "--weight_dir":
                param_dict['weight_dir'] = arg
            elif opt == "--models":
                param_dict['models'] = arg.split(',')
            elif opt == "--out_dir":
                param_dict['out_dir'] = arg
            elif opt == "--chrom":
                param_dict['chrom'] = str(arg)
            elif opt == "--window":
                param_dict['window'] = int(arg)
            elif opt == "--samples":
                param_dict['samples'] = arg
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
    elif param_dict['geno_dir'] is None:
        print('* Please specify the directory to the binary file of LD reference panel --geno_dir\n')
        sys.exit(2)
    elif param_dict['weight_dir'] is None:
        print('* Please specify the estimated eQTL weights file using --weight_dir\n')
        sys.exit(2)
    elif param_dict['out_dir'] is None:
        print('* Please specify the output directory\n')
        sys.exit(2)
    elif param_dict['chrom'] is None:
        print('* Please specify the chromosome --chrom\n')
        sys.exit(2)
    elif param_dict['models'] is None:
        print('* Please specify the imputation models --models\n')
        sys.exit(2)

    for key in param_dict:
        print('--%s=%s' % (key, param_dict[key]))

    print('\n')
    return param_dict

############################################################

import ottersutils as ots
param_dict = parse_param()

# set chrom and paths
out_dir = param_dict['out_dir']
ots.check_path(out_dir)

# Read in gene annotation file
GeneAnno, TargetID, n_targets = ots.read_anno(anno_dir=param_dict['anno_dir'],
                                              chrom=param_dict['chrom'])

############################################################
# Create directory for temporary files
tmp_dir = os.path.join(out_dir, 'tmp')
ots.check_path(tmp_dir)

# Sort and Tabix eQTL weights
models, tabix_files = ots.sort_tabix_weight(param_dict['weight_dir'], param_dict['models'], tmp_dir)

if not tabix_files:
    print('No weights file has been successfully sorted and tabixed. Please check the format of the weight file')
    sys.exit(2)


samples = pd.read_csv(param_dict['samples'],
                      sep='\t',
                      low_memory=False,
                      header=None,
                      names=['FID', 'IID'],
                      dtype={'FID': object, 'IID': object})

Result = samples.copy().T

# Create output files
for model in models:
    Result.to_csv(os.path.join(out_dir, model + '.txt'),
                  sep='\t',
                  index=True,
                  header=None,
                  mode='w')


############################################################
@ots.error_handler
def thread_process(num):

    target = TargetID[num]
    print('num=' + str(num) + '\nTargetID=' + target)
    target_info = GeneAnno.iloc[[num]]

    start, end = ots.extract_start_end_pos(target_info=target_info,
                                           window=param_dict['window'])

    ################# Read in genotype data #############

    # extract PLINK binary files for the target gene
    # set output path of the target gene
    target_dir = os.path.join(out_dir, target)
    ots.check_path(target_dir)

    # call PLINK to extract the binary file for the target gene
    extract_proc = ots.call_PLINK_extract(bim_path=param_dict['geno_dir'],
                                          out_path=target_dir,
                                          target=target,
                                          chrom=param_dict['chrom'],
                                          start_pos=start,
                                          end_pos=end)
    if not extract_proc:
        print('Remove temporary files. \n')
        shutil.rmtree(target_dir)
        return None

    # read in snps in reference panel
    target_ref = ots.read_format_ref_bim(ref_dir=target_dir,
                                         ref_file=target + '.bim')

    if target_ref is None:
        print('There is no reference bim file.')
        return None

    ################# GReX #####################

    for idx in range(len(models)):

        model = models[idx]
        w_file = tabix_files[idx]

        ################# Read in eQTL weights #####################

        target_w = ots.read_sst(sst_file=w_file,
                                sst='eQTL weights',
                                target=target,
                                chrom=param_dict['chrom'],
                                start_pos=start,
                                end_pos=end)

        if target_w is None:
            print('There is no estimated eQTL weights for ' + model)
            continue

        ########################### Find overlap snps #################

        # find overlap snps between reference panel, summary statistics, and estimated eQTL weights
        new_ref, target_w, snp_overlap = ots.match_snp_ID_double(df_ref=target_ref,
                                                                 df_1=target_w)

        if not snp_overlap.size:
            print('No overlapping test eQTLs between eQTL weights, GWAS summary statistics, and LD reference panel for' + model)
            continue

        target_w[['snpID', 'A1', 'ES']].to_csv(
            os.path.join(target_dir, model + '.txt'),
            sep='\t',
            index=None,
            header=True,
            mode='w')

        ################# Use PLINK to Calculate GReX #####################

        cmd = ["plink --bfile " + target +
               " --score " + model + ".txt" + " 1 2 3 header" +
               " --out GReX"]

        try:
            grex_proc = subprocess.check_call(cmd,
                                              stdout=subprocess.PIPE,
                                              cwd=target_dir,
                                              shell=True)

            # print('LD calculation completed.')
        except subprocess.CalledProcessError:
                print('GReX Imputation Failed For ' + model + '.\n')
                continue

        target_grex = pd.read_csv(os.path.join(target_dir, "GReX.profile"),
                                  low_memory=False,
                                  delim_whitespace=True,
                                  header=0,
                                  dtype={'FID': object, 'IID': object})

        tmp_merge = samples.merge(target_grex, how='left')

        tmp_merge[target] = tmp_merge[['SCORE']]

        tmp_merge.T.loc[[target]].to_csv(
            os.path.join(out_dir, model + '.txt'),
            index=True,
            sep='\t',
            header=None,
            mode='a')

    ############################ Clean temporary files ##################
    shutil.rmtree(target_dir)
    print('Done. \n')

############################################################
if __name__ == '__main__':
    print('Starting TWAS ' + str(n_targets) + ' target genes for ' + ','.join(models) + '.\n')
    pool = multiprocessing.Pool(param_dict['thread'])
    pool.imap(thread_process, [num for num in range(n_targets)])
    pool.close()
    pool.join()
    print('Remove temporary files.')
    shutil.rmtree(os.path.join(out_dir, 'tmp'))

############################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = ots.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)





