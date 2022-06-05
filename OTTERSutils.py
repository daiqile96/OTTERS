#!/usr/bin/env python

#########################################################
import functools
import os
import subprocess
import sys
import traceback

import pandas as pd


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
    cmd = ["plink --bfile "+bim_path+" --keep-allele-order --extract range " + range + " --make-bed --out " + out_geno]
    return cmd


def call_PLINK_clump(bim_path, r2, pvalue_path, window, p=1, snp_field="snpID", p_field="P"):

    cmd = ["plink --bfile " + bim_path+" --clump-p1 " + str(p) + " --clump-r2 " + str(r2) +
           " --clump-kb " + str(window)+" --clump " + pvalue_path +
           " --clump-snp-field " + snp_field + " --clump-field " + p_field +
           " --keep-allele-order --out "+bim_path]

    return cmd


################ Functions related to TABIX ##################
# Call tabix, read in lines into byt array
def call_tabix(path, chr, start, end):

    proc = subprocess.Popen(
        ["tabix "+path+" "+chr+":"+start+"-"+end],
        shell=True,
        stdout=subprocess.PIPE)

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

