import multiprocessing
import subprocess
import getopt
import shutil
import sys
import os

from io import StringIO
from time import time
from numpy import linalg
from functools import reduce
from scipy.stats.distributions import chi2

import numpy as np
import pandas as pd


############################################################
# time calculation
start_time = time()

############################################################


def parse_param():

    long_opts_list = ['OTTERS_dir=', 'anno_dir=', 'geno_dir=',
                      'weight_file=', 'gwas_file=', 'out_dir=', 'chrom=', 'r2=',
                      'window=', 'thread=', 'help']

    param_dict = {'OTTERS_dir': None, 'anno_dir': None,
                  'geno_dir': None, 'weight_file': None, 'gwas_file': None, 'out_dir': None, 'chrom': None,
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
            elif opt == "--geno_dir":
                param_dict['geno_dir'] = arg
            elif opt == "--weight_file":
                param_dict['weight_file'] = arg
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
    elif param_dict['weight_file'] is None:
        print('* Please specify the estimated eQTL weights file using --weight_file\n')
    elif param_dict['gwas_file'] is None:
        print('* Please specify the GWAS summary statistics using --gwas_file\n')
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

# set chrom and paths
chr = param_dict['chrom']
geno_dir = param_dict['geno_dir']
out_dir = param_dict['out_dir']
ots.check_path(out_dir)

print('Reading gene annotation data on chromosome ' + chr + '...')
anno_chunks = pd.read_csv(
    param_dict['anno_dir'],
    sep='\t',
    header=0,
    chunksize=10000,
    iterator=True,
    dtype={'CHROM': object,
        'GeneEnd': np.int64,
        'GeneStart': np.int64,
        'TargetID': object})

GeneAnno = pd.concat([x[x['CHROM'] == chr] for x in anno_chunks]).reset_index(drop=True)

if GeneAnno.empty:
    raise SystemExit('There are no valid gene annotation data for chromosome ' + chr + '\n')

GeneAnno = ots.optimize_cols(GeneAnno)
TargetID = GeneAnno.TargetID
n_targets = TargetID.size

out_twas_path = os.path.join(out_dir, 'TWAS.txt')
print('Creating file: ' + out_twas_path + '\n')
out_cols = ['CHROM', 'GeneStart', 'GeneEnd', 'TargetID', 'n_snps', 'FUSION_Z', 'FUSION_PVAL']
pd.DataFrame(columns=out_cols).to_csv(
    out_twas_path,
    sep='\t',
    index=None,
    header=True,
    mode='w')

############################################################


@ots.error_handler
def thread_process(num):

    target = TargetID[num]
    print('num=' + str(num) + '\nTargetID=' + target)
    target_info = GeneAnno.iloc[[num]]

    start = str(max(int(target_info.GeneStart) - param_dict['window'],0))
    end = str(int(target_info.GeneEnd) + param_dict['window'])

    ################# eQTL weights #####################
    print('*Extracting estimated eQTL effect sizes...*')

    w_proc_out = ots.call_tabix(param_dict['weight_file'], chr, start, end)

    if not w_proc_out:
        print('There is no estimated eQTLs for the range of ' + target)
        return None

    w_chunks = pd.read_csv(StringIO(w_proc_out.decode('utf-8')), sep='\t',
                           low_memory=False,
                           header=None,
                           names=["CHROM", "POS", "A1", "A2", "TargetID", "ES"],
                           iterator=True,
                           chunksize=1000)

    target_w = pd.concat([chunk[chunk.TargetID == target] for chunk in w_chunks]).reset_index(drop=True)

    if len(target_w) == 0:
        print('There is no estimated eQTLs for TargetID: ' + target)
        return None

    target_w['snpID'] = ots.get_snpIDs(target_w, flip=False)
    target_w['snpIDflip'] = ots.get_snpIDs(target_w, flip=True)

    ################# GWAS summary statistics #####################

    print('*Extracting GWAS summary statistics...*')

    sst_proc_out = ots.call_tabix(param_dict['gwas_file'], chr, start, end)

    if not sst_proc_out:
        print('There is no estimated eQTLs for the range of ' + target)
        return None

    sst_chunks = pd.read_csv(StringIO(sst_proc_out.decode('utf-8')), sep='\t',
                             low_memory=False,
                             header=None,
                             names=["CHROM", "POS", "A1", "A2", "TargetID", "Z"],
                             iterator=True,
                             chunksize=1000)

    target_sst = pd.concat([chunk[chunk.TargetID == target] for chunk in sst_chunks]).reset_index(drop=True)

    if len(target_sst) == 0:
        print('There is no GWAS summary statistics for TargetID: ' + target)
        return None

    target_sst['snpID'] = ots.get_snpIDs(target_sst, flip=False)
    target_sst['snpIDflip'] = ots.get_snpIDs(target_sst, flip=True)

    ################# Generate PLINK binary files for the target gene #####################

    print('*Generate PLINK binary files for the target gene...*')

    # set output path of the target gene
    target_dir = os.path.join(out_dir, target)
    ots.check_path(target_dir)

    # generate command to call PLINK to extract the binary file for the target gene
    extract_cmd = ots.call_PLINK_extract(geno_dir, target_dir, target, chr, start, end)

    try:
        proc = subprocess.check_call(extract_cmd,
                                     stdout=subprocess.PIPE,
                                     shell=True)
        print('Done making binary files.')
    except subprocess.CalledProcessError:
        print('There is no genotype data for TargetID: ' + target)
        print('Remove binary files for TargetID: ' + target + '.\n')
        shutil.rmtree(target_dir)
        return None

    # read in snps in LD reference panel
    target_bim = os.path.join(target_dir, target + '.bim')
    ref_chunks = pd.read_csv(target_bim, sep='\t',
                             low_memory=False,
                             header=None,
                             names=["CHROM", "SNP", "bp", "POS", "A1", "A2"],
                             iterator=True,
                             chunksize=1000,
                             usecols=["CHROM", "POS", 'bp', "A1", "A2"])
    target_ref = pd.concat([chunk for chunk in ref_chunks]).reset_index(drop=True)

    print('Formatting SNP ID in LD reference bim file...')
    target_ref['snpID'] = ots.get_snpIDs(target_ref, flip=False)
    target_ref[['CHROM', 'snpID', 'bp', 'POS', 'A1', 'A2']].to_csv(
        target_bim,
        sep='\t',
        index=None,
        header=None,
        mode='w')

    ########################### Find overlap snps #################

    snp_overlap = reduce(np.intersect1d, (target_w.snpID, target_sst[['snpID', 'snpIDflip']], target_w[['snpID', 'snpIDflip']]))

    if not snp_overlap.size:
        print('No overlapping test eQTLs between eQTL weights and LD reference panel for TargetID: ' + target)
        print('Remove binary files for TargetID: ' + target + '.\n')
        shutil.rmtree(target_dir)
        return None
    else:
        print('Find %s overlapped eQTLs' % (len(snp_overlap)))

    # filter out non-matching snpID rows
    target_sst = target_sst[np.any(target_sst[['snpID', 'snpIDflip']].isin(snp_overlap), axis=1)].reset_index(drop=True)
    target_w = target_w[np.any(target_w[['snpID', 'snpIDflip']].isin(snp_overlap), axis=1)].reset_index(drop=True)

    # if not in target_w.snpID, assumed flipped; if flipped, flip Zscore sign
    target_sst['flip'] = np.where(target_sst.snpID.isin(target_w.snpID.values), 1, -1)
    if not np.all(target_sst['flip'] == 1):
        idx = (target_sst['flip'] == -1)
        target_sst.loc[idx, ['snpID']] = target_sst.loc[idx, ['snpIDflip']].values
        target_sst.loc[idx, ['A1', 'A2']] = target_sst.loc[idx, ['A2', 'A1']].values
        target_sst['Z'] = target_sst['flip'] * target_sst['Z']

    target_w['flip'] = np.where(target_w.snpID.isin(target_ref.snpID.values), 1, -1)
    if not np.all(target_w['flip'] == 1):
        idx = (target_w['flip'] == -1)
        target_w.loc[idx, ['snpID']] = target_w.loc[idx, ['snpIDflip']].values
        target_w.loc[idx, ['A1', 'A2']] = target_w.loc[idx, ['A2', 'A1']].values
        target_w['ES'] = target_w['flip'] * target_w['ES']

    ################# Calculate LD Reference #####################

    # call PLINK to generate LD 
    target_snplist = target_w['snpID']
    target_snplist_path = os.path.join(target_dir, target + '.snplist')

    target_snplist.to_csv(target_snplist_path,
                          sep='\t',
                          index=None,
                          header=None,
                          mode='w')

    # Calculate LD for these SNPs.
    print('Generating reference LD matrix...')
    get_ld_cmd = ["plink --bfile " + target +
                  " --keep-allele-order" +
                  " --extract " + target + '.snplist' +
                  " --r square" +
                  " --out " + target +
                  " --memory 2000 "]

    try:
        ld_proc = subprocess.check_call(get_ld_cmd,
                                        stdout=subprocess.PIPE,
                                        cwd=target_dir,
                                        shell=True)
        print('LD calculation completed.')
    except subprocess.CalledProcessError as err:
        print('calculation failed for TargetID: ' + target + '\n')
        return None

    # read in LD covariance
    target_ld = os.path.join(target_dir, target + '.ld')

    ld_chunks = pd.read_csv(target_ld, sep='\t',
                            low_memory=False,
                            header=None,
                            iterator=True,
                            chunksize=1000)

    ld = pd.concat([chunk for chunk in ld_chunks]).reset_index(drop=True)
    # # remove the binary files and ld file 
    shutil.rmtree(target_dir)

    # format LD reference matrix calculated by PLINK
    V = ld.to_numpy()
    # PLINK return nan when all mutually-nonmissing set is homozygous major
    # We set nan = 0 here 
    V = np.nan_to_num(V)

    # Make sure the ld is positive definite matrix
    # PLINK rounds to 6 decimal points, which sometimes makes the correlation matrix become not positive definite
    _, s, v = linalg.svd(V)
    h = np.dot(v.T, np.dot(np.diag(s), v))
    V = (V+h)/2

    ################# Calculate TWAS statistics #####################

    def get_pval(z): return np.format_float_scientific(1-chi2.cdf(z**2, 1), precision=15, exp_digits=0)

    def get_z_denom(V, w):
        return np.sqrt(np.linalg.multi_dot([w, V, w]))

    def get_fusion_zscore(V_cor, w, Z_gwas, snp_sd=None):
        Z_twas = np.vdot(Z_gwas, w) / get_z_denom(V_cor, w)
        return Z_twas, get_pval(Z_twas)

    Result = target_info.copy()
    Result['n_snps'] = len(snp_overlap)
    Result['FUSION_Z'], Result['FUSION_PVAL'] = get_fusion_zscore(V, target_w.ES.values, target_sst.Z.values)

    # write to file
    Result.to_csv(
        out_twas_path,
        sep='\t',
        index=None,
        header=None,
        mode='a')

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