#!/usr/bin/env python

#########################################################
import functools
import os
import subprocess
import sys
import traceback
from itertools import groupby
from numpy import linalg
from io import StringIO

import pandas as pd
import numpy as np


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
        df['CHROM'] = df['CHROM'].astype(str).astype(int)

    ints = df.select_dtypes(include=['int64']).columns.tolist()
    df[ints] = df[ints].apply(pd.to_numeric, downcast='integer')

    floats = df.select_dtypes(include=['float64']).columns.tolist()
    df[floats] = df[floats].apply(pd.to_numeric, downcast='float')

    return df

################ Functions related to PLINK ##################


def call_PLINK_extract(bim_path, out_path, target, chr, start, end):

    # save the range of the gene
    range = os.path.join(out_path, 'range.txt')
    with open(range, 'w') as ff:
        ff.write('%s\t%s\t%s\t%s\n' % (chr, start, end, target))

    # extract the genotype data for this range
    out_geno = os.path.join(out_path, target)
    cmd = ["plink --bfile "+bim_path+" --extract range " + range + " --make-bed --out " + out_geno]
    return cmd


def call_PLINK_clump(bim_path, r2, pvalue_path, window=1000, p=1, snp_field="snpID", p_field="P"):

    cmd = ["plink --bfile " + bim_path+" --clump-p1 " + str(p) + " --clump-r2 " + str(r2) +
           " --clump-kb " + str(window)+" --clump " + pvalue_path +
           " --clump-snp-field " + snp_field + " --clump-field " + p_field +
           " --out "+bim_path]

    return cmd

################ Functions related to TABIX ##################

# Call tabix, read in lines into byt array
def call_tabix(path, chr, start, end):

    proc = subprocess.Popen(
        ["tabix "+path+" "+chr+":"+start+"-"+end],
        shell=True,
        stdout=subprocess.PIPE,
        bufsize=1)

    proc_out = bytearray()

    # process while subprocesses running
    while proc.poll() is None:
        line =  proc.stdout.readline()
        if len(line) == 0:
            break
        proc_out += line

    # get any remaining lines
    for line in proc.stdout:
        proc_out += line

    return proc_out

################ Functions related to LD ##################


# get header of ld file, get indices of columns to read in
def get_ld_cols(path):
    # get header
    file_cols = tuple(pd.read_csv(
        path,
        sep='\t',
        header=0,
        compression='gzip',
        low_memory=False,
        nrows=0).rename(columns={'#0': 'row'}))

    cols = ['row', 'SNP', 'COV']

    file_cols_ind = tuple([file_cols.index(x) for x in cols])

    return file_cols, file_cols_ind


# yields formatted tabix regions strings
def get_ld_regions_list(snp_ids):

    # 'chrm:' prefix for region string
    chrm = snp_ids[0].split('_')[0] + ':'

    # snp pos values as integers
    pos_vals = [int(snp.split('_')[1]) for snp in snp_ids]

    # get intervals of start,end positions; convert to tabix string
    for x, y in groupby(enumerate(pos_vals), lambda p: p[1]-p[0]):
        y = list(y)

        # chrm:start-end
        yield chrm + str(y[0][1]) + '-' + str(y[-1][1])


# yields formatted tabix regions strings
def get_regions_list(snp_ids):

    # split into groups by chromosome
    for chrm, grp in groupby(snp_ids, lambda x: x.split('_')[0]):
        # snp pos values as integers
        pos_vals = [int(snp.split('_')[1]) for snp in grp]

        # get intervals of start,end positions; convert to tabix string
        for x, y in groupby(enumerate(pos_vals), lambda p: p[1]-p[0]):
            y = list(y)

            # chrm:start-end
            yield chrm + ':' + str(y[0][1]) + '-' + str(y[-1][1])


# call tabix using regions string
def call_tabix_regions(path, regs_str, filter_line=lambda x: x):

    proc = subprocess.Popen(
        ['tabix '+path+' '+regs_str],
        shell=True,
        stdout=subprocess.PIPE,
        bufsize=1)
    proc_out = bytearray()

    # process while subprocesses running
    while proc.poll() is None:
        line = proc.stdout.readline()
        if len(line) == 0:
            break
        proc_out += filter_line(line)

    # leftover lines
    for line in proc.stdout:
        proc_out += filter_line(line)

    return proc_out


# get proc_out from function and parse data for regions
def get_ld_regions_data(regs_str, path, snp_ids, ld_cols, ld_cols_ind):

    proc_out = call_tabix_regions(path, regs_str)

    regs_data = pd.read_csv(
        StringIO(proc_out.decode('utf-8')),
        sep='\t',
        low_memory=False,
        header=None,
        names=ld_cols,
        usecols=ld_cols_ind, 
        dtype={
            'SNP': object,
            'row': np.int32,
            'COV': object}
        ).drop_duplicates(['SNP'], keep='first')

    regs_data = regs_data[regs_data.SNP.isin(snp_ids)]

    return regs_data


# read in covariance data for snps
def get_ld_data(path, snp_ids):

    # get columns names, indices for ld file
    ld_cols, ld_cols_ind = get_ld_cols(path)

    # format tabix regions from snp_ids; 'chrm:start-end'
    regs_lst = list(get_ld_regions_list(snp_ids))

    N = len(regs_lst)

    # arguments to pass
    regs_args = [path, snp_ids, ld_cols, ld_cols_ind]
    try:
        regs_str = ' '.join(regs_lst)
        cov_data = get_ld_regions_data(regs_str, *regs_args)

    except OSError:
        # argument may be too long for OS; if so try subset instead of getting all regions at once
        # print('Subseting regions to tabix.')
        n = 2500
        while n:
            try: 
                regs_str_lst = [' '.join(regs_lst[i:i+n]) for i in range(0, N, n)]
                cov_data = pd.concat([get_ld_regions_data(regs_str, *regs_args) for regs_str in regs_str_lst])
            except OSError:
                n -= 500
                pass
            else:
                n = 0

    return cov_data.set_index('row')


# get the ld matrix and size of the ld
def get_ld_matrix(MCOV):
    MCOV = MCOV.copy()

    MCOV['COV'] = MCOV['COV'].apply(lambda x: np.fromstring(x, dtype=np.float32, sep=','))

    inds = MCOV.index
    n_inds = inds.size
    V_upper = np.zeros((n_inds, n_inds))

    for i in range(n_inds):
        cov_i = MCOV.COV.loc[inds[i]]
        N = cov_i.size

        for j in range(i, n_inds):
            if inds[j] - inds[i] < N:
                V_upper[i, j] = cov_i[inds[j] - inds[i]]
            else:
                V_upper[i, j] = 0

    snp_Var = V_upper.diagonal()
    V = V_upper + V_upper.T - np.diag(snp_Var)
    V = np.nan_to_num(V)

    _, s, v = linalg.svd(V)
    h = np.dot(v.T, np.dot(np.diag(s), v))
    V = (V+h)/2
    size = len(V)

    return size, V
