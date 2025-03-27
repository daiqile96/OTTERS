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
import sys
import os


from time import time


############################################################
# time calculation
start_time = time()

############################################################


def parse_param():

    long_opts_list = ['OTTERS_dir=', 'anno_dir=', 'geno_dir=',
                      'weight_dir=', 'gwas_file=', 'out_dir=', 'chrom=',
                      'window=', 'thread=', 'models=', 'help']

    param_dict = {'OTTERS_dir': None, 'anno_dir': None,
                  'geno_dir': None, 'weight_dir': None, 'gwas_file': None, 'out_dir': None, 'chrom': None,
                  'window': 1000000, 'thread': 1, 'models': None}

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
            elif opt == "--gwas_file":
                param_dict['gwas_file'] = arg
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
    elif param_dict['geno_dir'] is None:
        print('* Please specify the directory to the binary file of LD reference panel --geno_dir\n')
        sys.exit(2)
    elif param_dict['weight_dir'] is None:
        print('* Please specify the estimated eQTL weights file using --weight_dir\n')
    elif param_dict['gwas_file'] is None:
        print('* Please specify the GWAS summary statistics using --gwas_file\n')
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

# Create output files
for model in models:
    out_cols = ['CHROM', 'GeneStart', 'GeneEnd', 'TargetID', 'n_snps', 'FUSION_Z', 'FUSION_PVAL']
    ots.create_file_title(out_cols=out_cols,
                          out_dir=out_dir,
                          out_file=model + '.txt')

# Sort and Tabix GWAS summary statistics
tabix_gwas = ots.sort_tabix_sst(sst_dir=param_dict['gwas_file'],
                                out_dir=tmp_dir)
if not tabix_gwas:
    sys.exit(2)

############################################################
@ots.error_handler
def thread_process(num):

    target = TargetID[num]
    print('num=' + str(num) + '\nTargetID=' + target)
    target_info = GeneAnno.iloc[[num]]

    start, end = ots.extract_start_end_pos(target_info=target_info,
                                           window=param_dict['window'])

    ################# Read in GWAS summary statistics #####################

    target_g = ots.read_sst(sst_file=tabix_gwas,
                            sst='GWAS',
                            target=target,
                            chrom=param_dict['chrom'],
                            start_pos=start,
                            end_pos=end)

    if target_g is None:
        print('There is no GWAS summary statistics.')
        return None

    ################# Read in reference LD ##############

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

    ################# TWAS #####################

    Result = target_info.copy()
    # Result['n_snps'] = len(snp_overlap)

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
        sub_ref, sub_g, target_w, snp_overlap = ots.match_snp_ID_triple(df_ref=target_ref,
                                                                        df_1=target_g,
                                                                        df_2=target_w)

        if snp_overlap is None:
            print('No overlapping test eQTLs between eQTL weights, GWAS summary statistics, and LD reference panel for' + model)
            continue

        ################# Calculate LD Reference #####################

        # generate LD reference matrix
        V, blk_size = ots.plink_LD_cmd(bim_dir=target,
                                       sst_df=target_w,
                                       out_file=target,
                                       work_dir=target_dir)

        if blk_size is None:
            continue

        ################# Calculate FUSION TWAS statistics #####################

        Result['n_snps'] = len(snp_overlap)
        Result['Z'], Result['P'] = ots.get_fusion_zscore(V_cor=V,
                                                         w=target_w.ES.values,
                                                         Z_gwas=sub_g.Z.values)

        # write to file
        Result.to_csv(
            os.path.join(out_dir, model + '.txt'),
            sep='\t',
            index=None,
            header=None,
            mode='a')

        print('Finish ' + model + '.')

    ############################ Clean temporary files ##################
    # shutil.rmtree(target_dir)
    print('Done. \n')

############################################################
if __name__ == '__main__':
    print('Starting TWAS ' + str(n_targets) + ' target genes for ' + ','.join(models) + '.\n')
    pool = multiprocessing.Pool(param_dict['thread'])
    pool.imap(thread_process, [num for num in range(n_targets)])
    pool.close()
    pool.join()
    print('Remove temporary files.')
    # shutil.rmtree(os.path.join(out_dir, 'tmp'))

############################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = ots.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)





