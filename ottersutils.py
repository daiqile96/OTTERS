#!/usr/bin/env python

#########################################################
import functools
import os
import subprocess
import sys
import traceback
from numpy import linalg
from functools import reduce
from io import StringIO
from scipy.stats.distributions import chi2
import pandas as pd
import numpy as np
import getopt
import pysam

#########################################################
## FUNCTIONS:

# error_handler
# format_elapsed_time
# check_path
# optimize_cols
# get_snpIDs

# call_tabix

# call_PLINK_clump
# call_PLINK_extract


def sort(in_file, out_file):

    sort_cmd = 'tail -n +2 ' + in_file + ' | sort -n -k1 -k2 >> ' + out_file

    try:
        proc = subprocess.check_call(sort_cmd,
                                     stdout=subprocess.PIPE,
                                     shell=True)
    except subprocess.CalledProcessError:
        print('Sorting failed for ' + in_file)
        return None

    return out_file


def tabix(in_file, format='bgzip', tmp_file=None):

    """ Tabix a txt or txt.gz file with option to keep the original file or not """

    # if the file is not compressed, firt compress and save to a temporary file
    if format == 'txt':
        s = 1
        try:
            pysam.tabix_compress(in_file, tmp_file, force=True)
        except:
            print('Compress failed for ' + in_file)
            return None
    else:
        s = 0
        tmp_file = in_file

    try:
        out_file = pysam.tabix_index(tmp_file, force=True,
                                     seq_col=0, start_col=1, end_col=1,
                                     line_skip=s)
    except:
        print('Tabix failed for ' + in_file)
        return None

    return out_file


def sort_tabix_sst(sst_dir, out_dir):

    print('Sort and Tabix summary statistics...')

    if not os.path.exists(sst_dir):
        print("Summary statistics doesn't exists.")
        return None

    # sort trained eQTL weights by chromosome and position of SNPs
    sort_sst = sort(sst_dir, os.path.join(out_dir, 'sort_sst.txt'))

    if not sort_sst:
        return None

    # tabix sorted file
    tabix_out = tabix(sort_sst)
    if not tabix_out:
        return None

    return tabix_out


def sort_tabix_weight(weight_dir, models, out_dir):

    print('Sort and Tabix trained eQTL weights...')
    # create a temp directory to save sorted and tabixed weights
    idx = 0
    fail_idx = 0
    success = []
    out = []

    for model in models:

        idx += 1
        w_file = os.path.join(weight_dir, model + '.txt')

        if not os.path.exists(w_file):
            fail_idx += 1
            print(w_file + 'does not exists')
            continue

        # sort trained eQTL weights by chromosome and position of SNPs
        sort_weights = sort(w_file, os.path.join(out_dir, 'sort_' + model + '.txt'))
        if not sort_weights:
            fail_idx += 1
            continue

        # tabix sorted file
        tabix_out = tabix(sort_weights)
        if not tabix_out:
            fail_idx += 1
            continue

        success.append(model)
        out.append(tabix_out)

    return success, out

#########################################################
# wrapper for thread_process functions; adds error catching/logging for failed targets
def error_handler(func):
    @functools.wraps(func)
    def wrapper(num, *args, **kwargs):
        try:
            return func(num, *args, **kwargs)

        except Exception as e:
            e_info = sys.exc_info()
            e_type = e_info[0].__name__

            # don't print traceback info for wrapper
            e_tracebk = ''.join(traceback.format_tb(e_info[2])[1:])

            print('Caught an exception for num={}:\n  {}: {}\nTraceback:\n{}'.format(num, e_type, e, e_tracebk))

        finally:
            sys.stdout.flush()

    return wrapper


# return human readable elapsed time string
def format_elapsed_time(time_secs):
    val = abs(int(time_secs))

    day = val // (3600*24)
    hour = val % (3600*24) // 3600
    mins = val % 3600 // 60
    secs = val % 60

    res = '%02d:%02d:%02d:%02d' % (day, hour, mins, secs)

    if int(time_secs) < 0:
        res = "-%s" % res

    return res


# check if a path exists. If not, create the path
def check_path(path):
    if not os.path.isdir(path):
        os.makedirs(path)


# return snpIDs
def get_snpIDs(df: pd.DataFrame, flip=False):
    chroms = df['CHROM'].astype('str').values
    pos = df['POS'].astype('str').values
    A1 = df['A1'].values
    A2 = df['A2'].values

    if flip:
        return ['_'.join(i) for i in zip(chroms, pos, A1, A2)]
    else:
        return ['_'.join(i) for i in zip(chroms, pos, A2, A1)]


# Decrease memory by downcasting 'CHROM' column to integer, integer and float columns to minimum size that will not lose info
def optimize_cols(df: pd.DataFrame):

    if 'CHROM' in df.columns:
        df['CHROM'] = df['CHROM'].astype(str)

    ints = df.select_dtypes(include=['int64']).columns.tolist()
    df[ints] = df[ints].apply(pd.to_numeric, downcast='integer')

    floats = df.select_dtypes(include=['float64']).columns.tolist()
    df[floats] = df[floats].apply(pd.to_numeric, downcast='float')

    return df


# determine indices of file cols to read in, dtype of each col
def get_cols_dtype(cols):

    dtype_dict = {
        'CHROM': object,
        'POS': np.int64,
        'A1': object,
        'A2': object,
        'SNP': object,
        'ES': np.float64,
        'Z': np.float64,
        'N': np.int64,
        'GeneEnd': np.int64,
        'GeneName': object,
        'GeneStart': np.int64,
        'snpID': object,
        'TargetID': object,
        'bp': object}

    out_dtype_dict = {x: dtype_dict[x] for x in cols}

    return out_dtype_dict

################ Functions related to PLINK ##################


def call_PLINK_extract(bim_path, out_path, target, chrom, start_pos, end_pos):

    # save the range of the gene
    range = os.path.join(out_path, 'range.txt')
    with open(range, 'w') as ff:
        ff.write('%s\t%s\t%s\t%s\n' % (chrom, start_pos, end_pos, target))

    # extract the genotype data for this range
    out_geno = os.path.join(out_path, target)

    cmd = ["plink --bfile "+bim_path+" --keep-allele-order --extract range " + range + " --make-bed --out " + out_geno]

    try:
        proc = subprocess.check_call(cmd,
                                     stdout=subprocess.PIPE,
                                     shell=True)
    except subprocess.CalledProcessError:
        print('There is no genotype reference data.')
        return None

    return True


def call_PLINK_clump(bim_path, r2, pvalue_path, work_dir, window=1000, p=1, snp_field="snpID", p_field="P"):

    cmd = ["plink --bfile " + bim_path+" --clump-p1 " + str(p) + " --clump-r2 " + str(r2) +
           " --clump-kb " + str(window)+" --clump " + pvalue_path +
           " --clump-snp-field " + snp_field + " --clump-field " + p_field +
           " --keep-allele-order --out "+bim_path]

    # perform LD-clumping
    try:
        proc = subprocess.check_call(cmd,
                                     stdout=subprocess.PIPE,
                                     cwd=work_dir,
                                     shell=True)
        print('Done LD Clumping.')
    except subprocess.CalledProcessError:
        print('LD Clumping Failed. \n')
        return None


def call_PLINK_predict(bim_path, score_path, out_path):
    cmd = ["plink --bfile " + bim_path +
           " --score " + score_path + " 1 2 5" +
           " --out " + out_path]
    return cmd


################ Function to read in and format reference bim file ##################

def read_format_ref_bim(ref_dir, ref_file):

    target_bim = os.path.join(ref_dir, ref_file)

    col_names = ["CHROM", "SNP", "bp", "POS", "A1", "A2"]
    dtypes = get_cols_dtype(col_names)

    ref_chunks = pd.read_csv(target_bim, sep='\t',
                             low_memory=False,
                             header=None,
                             names=col_names,
                             iterator=True,
                             chunksize=1000,
                             dtype=dtypes)

    target_ref = pd.concat([chunk for chunk in ref_chunks]).reset_index(drop=True)

    if len(target_ref) == 0:
        return None

    # format snp IDs in the reference bim file
    target_ref['snpID'] = get_snpIDs(target_ref, flip=False)

    target_ref[['CHROM', 'snpID', 'bp', 'POS', 'A1', 'A2']].to_csv(
        target_bim,
        sep='\t',
        index=None,
        header=None,
        mode='w')

    return target_ref


################ Functions to read in summary statistics ##################

def extract_start_end_pos(target_info, window):

    start = str(max(int(target_info.GeneStart) - window,0))
    end = str(int(target_info.GeneEnd) + window)

    return start, end


def read_sst_by_chunks(tabix_out, sst, target):

    if not tabix_out:
        return None

    if sst == 'eQTL weights':
        col_names = ["CHROM", "POS", "A1", "A2", "TargetID", "ES"]
    elif sst == 'GWAS':
        col_names = ["CHROM", "POS", "A1", "A2", "Z"]
    elif sst == 'eQTL sst':
        col_names = ["CHROM", "POS", "A1", "A2", "Z", "TargetID", "N"]

    dtypes = get_cols_dtype(col_names)

    chunks = pd.read_csv(StringIO(tabix_out.decode('utf-8')), sep='\t',
                         low_memory=False,
                         header=None,
                         names=col_names,
                         iterator=True,
                         chunksize=1000,
                         dtype=dtypes)

    if sst == 'GWAS':
        target_df = pd.concat([chunk for chunk in chunks]).drop_duplicates(["CHROM", "POS", "A1", "A2"]).reset_index(drop=True)
    else:
        target_df = pd.concat([chunk[chunk.TargetID == target] for chunk in chunks]).drop_duplicates(["CHROM", "POS", "A1", "A2"]).reset_index(drop=True)

    if target_df.empty:
        return None

    target_df = optimize_cols(target_df)

    return target_df


def read_sst(sst_file, sst, target, chrom, start_pos, end_pos):

    # call tabix to extract estimated eQTL effect sizes for target gene
    tabix_out = call_tabix(sst_file, chrom, start_pos, end_pos)

    # read in estimated eQTL effect sizes
    target_df = read_sst_by_chunks(tabix_out=tabix_out,
                                   sst=sst,
                                   target=target)

    return target_df


################ Function to read in annotation file ##################

def read_anno(anno_dir, chrom):

    col_names = ['CHROM', 'GeneEnd', 'GeneStart', 'TargetID']
    dtypes = get_cols_dtype(col_names)

    print('Read gene annotation data on chromosome ' + chrom + '...')
    anno_chunks = pd.read_csv(
        anno_dir,
        sep='\t',
        header=0,
        usecols=col_names,
        chunksize=10000,
        iterator=True,
        dtype=dtypes)

    GeneAnno = pd.concat([x[x['CHROM'] == chrom] for x in anno_chunks]).reset_index(drop=True)

    if GeneAnno.empty:
        print('There are no valid gene annotation data for chromosome ' + chrom + '\n')
        return None

    GeneAnno = optimize_cols(GeneAnno)
    TargetID = GeneAnno.TargetID
    n_targets = TargetID.size

    return GeneAnno, TargetID, n_targets


################ Functions related to TABIX ##################
# Call tabix, read in lines into byt array
def call_tabix(path, chrom, start, end):

    chrom = str(chrom)

    proc = subprocess.Popen(
        ["tabix "+path+" "+chrom+":"+start+"-"+end],
        shell=True,
        stdout=subprocess.PIPE)

    proc_out = bytearray()

    # process while subprocesses running
    while proc.poll() is None:
        line = proc.stdout.readline()
        if len(line) == 0:
            break
        proc_out += line

    # get any remaining lines
    for line in proc.stdout:
        proc_out += line

    return proc_out


################ Functions related to Imputation ##################
def lassosum_cmd(chrom, bim_dir, sst_dir, out_dir, lassosum_path, N, ld_blocks):

    print('lassosum using LD blocks ' + ld_blocks)

    cmd = ['Rscript ' + lassosum_path +
           ' --medianN=' + str(int(N)) +
           ' --bim_file=' + bim_dir +
           ' --sst_file=' + sst_dir +
           ' --out_path=' + out_dir +
           ' --chr=' + str(chrom) +
           ' --LDblocks=' + ld_blocks]

    return cmd


def SDPR_LD_args(SDPR_path, chrom, r2, bim_dir, work_dir):

    # create the directory under work_dir to save LD reference 
    check_path(os.path.join(work_dir, "ref"))

    cmd = [SDPR_path +
           ' -make_ref -ref_prefix ' + bim_dir +
           ' -chr ' + str(chrom) +
           ' -r2 ' + str(r2) +
           ' -ref_dir ref/']
    try:
        proc = subprocess.check_call(cmd,
                                     stdout=subprocess.PIPE,
                                     cwd=work_dir,
                                     shell=True)
    except subprocess.CalledProcessError:
        print('SDPR failed to generate LD. ')
        return None


def SDPR_args(SDPR_path, sst_file, N,
              chrom, M, opt_llk, iter, burn, thin,
              a, c, a0k, b0k):

    cmd = [SDPR_path +
           ' -mcmc -ref_dir ref/ ' +
           ' -ss ' + sst_file +
           ' -N ' + str(int(N)) +
           ' -chr ' + str(chrom) +
           ' -M ' + str(M) +
           ' -opt_llk ' + str(opt_llk) +
           ' -iter ' + str(iter) +
           ' -burn ' + str(burn) +
           ' -thin ' + str(thin) +
           ' -a ' + str(a) +
           ' -c ' + str(c) +
           ' -a0k ' + str(a0k) +
           ' -b0k ' + str(b0k) +
           ' -out SDPR.txt']

    return cmd


def save_results(model, out_df, out_dir):

    out_df[['CHROM', 'POS', 'A1', 'A2', 'TargetID', model]].to_csv(
        os.path.join(out_dir, model + '.txt'),
        sep='\t',
        index=None,
        header=None,
        mode='a')

    print('Finish ' + model + '.')


def save_results_PT(p, out_df, out_dir):

    out_df[['CHROM', 'POS', 'A1', 'A2', 'TargetID', 'Beta']].to_csv(
        os.path.join(out_dir, 'P' + str(p) + '.txt'),
        sep='\t',
        index=None,
        header=None,
        mode='a')

    print('Finish P+T with p-value < ' + str(p) + '.')


def format_save_results(work_dir, out_dir, model, sst_df):

    raw_out_dir = os.path.join(work_dir, model + '.txt')

    if model == 'SDPR':
        if not os.path.exists(raw_out_dir):
            print('SDPR failed. Please check the input arguments for SDPR')
            return None
        else:
            raw_names = ["snpID", "A1", model]
    elif model == 'lassosum':
        if not os.path.exists(raw_out_dir):
            print('lassosum failed. Please check the input arguments for lassosum')
            return None
        else:
            raw_names = ['CHROM', 'POS', 'A1', 'A2', model]

    raw_chunks = pd.read_csv(raw_out_dir, sep='\t',
                             low_memory=False,
                             header=0,
                             names=raw_names,
                             iterator=True,
                             chunksize=1000)

    raw_weights = pd.concat([chunk for chunk in raw_chunks]).reset_index(drop=True)

    if model == 'SDPR':
        out_df = raw_weights.merge(sst_df, left_on=['snpID', 'A1'],
                                   right_on=['snpID', 'A1'], how="inner")
    elif model == 'lassosum':
        out_df = raw_weights.merge(sst_df, on=['CHROM', 'POS', 'A1', 'A2'],
                                   how="inner")

    # remove the raw outputs
    os.remove(raw_out_dir)

    save_results(model, out_df, out_dir)


def pos_def_matrix(mat):
    """ convert the input matrix to the cloest positive definite matrix"""
    # Make sure the ld is positive definite matrix
    _, s, v = linalg.svd(mat)
    h = np.dot(v.T, np.dot(np.diag(s), v))
    mat_pos_def = (mat+h)/2

    return mat_pos_def


def plink_LD_cmd(bim_dir, sst_df, out_file, work_dir, convert=None):

    target_snplist = sst_df['snpID']
    target_snplist_path = os.path.join(work_dir, 'snplist.txt')

    target_snplist.to_csv(target_snplist_path,
                          sep='\t',
                          index=None,
                          header=None,
                          mode='w')

    cmd = ["plink --bfile " + bim_dir +
           " --keep-allele-order" +
           " --extract " + 'snplist.txt' +
           " --r square" +
           " --out " + out_file +
           " --memory 2000 "]

    try:
            ld_proc = subprocess.check_call(cmd,
                                            stdout=subprocess.PIPE,
                                            cwd=work_dir,
                                            shell=True)
            # print('LD calculation completed.')
    except subprocess.CalledProcessError:
            print('LD calculation failed. \n')
            return None, None

    ld_dir = os.path.join(work_dir, out_file + '.ld')

    ld_chunks = pd.read_csv(ld_dir, sep='\t',
                            low_memory=False,
                            header=None,
                            iterator=True,
                            chunksize=1000)

    ld = pd.concat([chunk for chunk in ld_chunks]).reset_index(drop=True)
    # The .ld file is large, delete it after using it.
    os.remove(ld_dir)

    # print('Start formatting PRS-CS LD.')
    V = ld.to_numpy()
    # PLINK return nan when all mutually-nonmissing set is homozygous major
    # We set nan = 0 here 
    V = np.nan_to_num(V)

    if convert:
        # PLINK rounds to 6 decimal points, which sometimes makes the correlation matrix become not positive definite
        # convert it to be positive definite
        V = pos_def_matrix(V)

    blk_size = len(V)

    return V, blk_size


################ Read In Functions ##################


def create_file_title(out_cols, out_dir, out_file):

    out_path = os.path.join(out_dir, out_file)
    pd.DataFrame(columns=out_cols).to_csv(
        out_path,
        sep='\t',
        index=None,
        header=True,
        mode='w')


def output_data(out_df, out_cols, out_path, write_mode, out_header=None):

    out_df[out_cols].to_csv(
        out_path,
        sep='\t',
        index=None,
        header=out_header,
        mode=write_mode)


def read_in_clumped(clumped_file, chrom):

    with open(clumped_file) as file_in:
        lines = []
        snps = []
        for line in file_in:
            lines.append(line)
            res = line.split()
            if len(res) != 0:
                if res[0] == chrom:
                    snps.append(res[2])
                elif chrom == 'X' and res[0] == '23':
                    snps.append(res[2])
                elif chrom == 'Y' and res[0] == '24':
                    snps.append(res[2])

    return snps


################ Match snpIDs of data frames ##################

def filter_df_rows(df, filter_by, filter_cols):

    df = df[np.any(df[filter_cols].isin(filter_by), axis=1)].reset_index(drop=True)

    return df


def handle_flip_snps(df_ref, df_1, statistics):
    """Handle flipped snp IDs"""

    # filter out non-matching snpID rows in df_ref and df_1
    df_1 = df_1[np.any(df_1[['snpID', 'snpIDflip']].isin(df_ref.snpID.values), axis=1)].reset_index(drop=True)

    # if snpID is not in df_ref.snpID, assumed flipped; if flipped, flip Zscore sign
    df_1['flip'] = np.where(df_1.snpID.isin(df_ref.snpID.values), 1, -1)
    if not np.all(df_1['flip'] == 1):
            idx = (df_1['flip'] == -1)
            df_1.loc[idx, ['snpID']] = df_1.loc[idx, ['snpIDflip']].values
            df_1.loc[idx, ['A1', 'A2']] = df_1.loc[idx, ['A2', 'A1']].values
            df_1[statistics] = df_1['flip'] * df_1[statistics]
    return df_1
def match_snp_ID_impute(df_ref, df_1):

    df_ref['snpID'] = get_snpIDs(df_ref, flip=False)
    df_1['snpID'] = get_snpIDs(df_1, flip=False)
    df_1['snpIDflip'] = get_snpIDs(df_1, flip=True)

    # find overlapped snpIDs
    overlap_list = np.intersect1d(df_ref['snpID'], df_1[['snpID', 'snpIDflip']])

    if overlap_list.size:
        df_ref = df_ref[df_ref.snpID.isin(overlap_list)]
        df_1 = handle_flip_snps(df_ref, df_1, 'ES')
    else:
        return None, None, None

    return df_ref, df_1, overlap_list

def match_snp_ID_double(df_ref, df_1):

    df_ref['snpID'] = get_snpIDs(df_ref, flip=False)
    df_1['snpID'] = get_snpIDs(df_1, flip=False)
    df_1['snpIDflip'] = get_snpIDs(df_1, flip=True)

    # find overlapped snpIDs
    overlap_list = np.intersect1d(df_ref['snpID'], df_1[['snpID', 'snpIDflip']])

    if overlap_list.size:
        df_ref = df_ref[df_ref.snpID.isin(overlap_list)]
        df_1 = handle_flip_snps(df_ref, df_1, 'Z')
    else:
        return None, None, None

    return df_ref, df_1, overlap_list


def match_snp_ID_triple(df_ref, df_1, df_2):
    """Match snps in df_1 and df_2 with snps in df_ref"""
    df_1['snpID'] = get_snpIDs(df_1, flip=False)
    df_1['snpIDflip'] = get_snpIDs(df_1, flip=True)
    df_2['snpID'] = get_snpIDs(df_2, flip=False)
    df_2['snpIDflip'] = get_snpIDs(df_2, flip=True)

    # find overlapped snpIDs
    overlap_list = reduce(np.intersect1d, (df_ref.snpID, df_1[['snpID', 'snpIDflip']], df_2[['snpID', 'snpIDflip']]))

    if overlap_list.size:
        df_ref = df_ref[df_ref.snpID.isin(overlap_list)]
        # filter out rows with non-matching snpIDs in df_1 and df_2
        df_1 = handle_flip_snps(df_ref, df_1, 'Z')
        df_2 = handle_flip_snps(df_ref, df_2, 'ES')
    else:
        return None, None, None, None

    return df_ref, df_1, df_2, overlap_list


################ TWAS with GWAS summary statistics ##################


def get_pval(z):
    return np.format_float_scientific(chi2.sf(np.power(z, 2), 1), precision=15, exp_digits=0)


def get_z_denom(V, w):
    return np.sqrt(np.linalg.multi_dot([w, V, w]))


def get_fusion_zscore(V_cor, w, Z_gwas, snp_sd=None):
    Z_twas = np.vdot(Z_gwas, w) / get_z_denom(V_cor, w)
    return Z_twas, get_pval(Z_twas)


