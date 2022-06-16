#!/usr/bin/env python

#########################################################
import functools
import os
import subprocess
import sys
import traceback
from numpy import linalg
import pandas as pd
import numpy as np
import getopt

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


################ Parse arguments ##################

def parse_param():

    long_opts_list = ['OTTERS_dir=', 'anno_dir=', 'geno_dir=',
                      'sst_file=', 'out_dir=', 'chrom=', 'r2=',
                      'window=', 'thread=', 'models=',
                      'pt=',
                      'prscs_a=', 'prscs_b=', 'prscs_phi=',
                      'prscs_n_iter=', 'prscs_n_burnin=', 'prscs_thin=',
                      'seed=',
                      'SDPR_dir=', 'SDPR_r2=', 'SDPR_M=', 'SDPR_opt_llk=',
                      'SDPR_iter=', 'SDPR_burn=', 'SDPR_thin=', 'SDPR_a=',
                      'SDPR_c=', 'SDPR_a0k=', 'SDPR_b0k=',
                      'help']

    param_dict = {'OTTERS_dir': None, 'anno_dir': None,
                  'geno_dir': None, 'sst_file': None, 'out_dir': None, 'chrom': None,
                  'r2': 2, 'window': 1000000, 'thread': 1, 'models': None,
                  'pt': [0.001, 0.05],
                  'prscs_a': 1, 'prscs_b': 0.5, 'prscs_phi': 0.0001,
                  'prscs_n_iter': 1000, 'prscs_n_burnin': 500, 'prscs_thin': 5,
                  'SDPR_dir': None, 'SDPR_r2': 0.1, 'SDPR_M': 1000, 'SDPR_opt_llk': 1,
                  'SDPR_iter': 1000, 'SDPR_burn': 200, 'SDPR_thin': 1, 'SDPR_a': 0.1,
                  'SDPR_c': 1, 'SDPR_a0k': 0.5, 'SDPR_b0k': 0.5,
                  'seed': None}

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
            elif opt == "--sst_file":
                param_dict['sst_file'] = arg
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
            elif opt == "--models":
                param_dict['models'] = arg.split(',')
            elif opt == "--pt":
                param_dict['pt'] = arg.split(',')
            elif opt == "--prscs_a":
                param_dict['prscs_a'] = float(arg)
            elif opt == "--prscs_b":
                param_dict['prscs_b'] = float(arg)
            elif opt == "--prscs_phi":
                param_dict['prscs_phi'] = float(arg)
            elif opt == "--prscs_n_iter":
                param_dict['prscs_n_iter'] = int(arg)
            elif opt == "--prscs_n_burnin":
                param_dict['prscs_n_burnin'] = int(arg)
            elif opt == "--prscs_thin":
                param_dict['prscs_thin'] = int(arg)
            elif opt == "--SDPR_dir":
                param_dict['SDPR_dir'] = arg
            elif opt == "--SDPR_r2":
                param_dict['SDPR_r2'] = float(arg)
            elif opt == "--SDPR_M":
                param_dict['SDPR_M'] = int(arg)
            elif opt == "--SDPR_opt_llk":
                param_dict['SDPR_opt_llk'] = int(arg)
            elif opt == "--SDPR_iter":
                param_dict['SDPR_iter'] = int(arg)
            elif opt == "--SDPR_burn":
                param_dict['SDPR_burn'] = int(arg)
            elif opt == "--SDPR_thin":
                param_dict['SDPR_thin'] = int(arg)
            elif opt == "--SDPR_a":
                param_dict['SDPR_a'] = float(arg)
            elif opt == "--SDPR_c":
                param_dict['SDPR_c'] = float(arg)
            elif opt == "--SDPR_a0k":
                param_dict['SDPR_a0k'] = float(arg)
            elif opt == "--SDPR_b0k":
                param_dict['SDPR_b0k'] = float(arg)
            elif opt == "--thread":
                param_dict['thread'] = int(arg)
            elif opt == "--seed":
                param_dict['seed'] = int(arg)
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
    elif param_dict['sst_file'] is None:
        print('* Please specify the eQTL summary statistics file using --sst_file\n')
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


################ Functions related to Imputation ##################
def lassosum_cmd(chrom, bim_dir, sst_dir, out_dir, lassosum_path, N):

    cmd = ['Rscript ' + lassosum_path +
           ' --medianN=' + str(int(N)) +
           ' --bim_file=' + bim_dir +
           ' --sst_file=' + sst_dir +
           ' --out_path=' + out_dir +
           ' --chr=' + str(chrom) +
           ' --LDblocks=EUR.hg38']

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


def PRScs_LD_cmd(bim_dir, sst_df, out_file, work_dir):

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
            return None

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

    # Make sure the ld is positive definite matrix
    # PLINK rounds to 6 decimal points, which sometimes makes the correlation matrix become not positive definite
    _, s, v = linalg.svd(V)
    h = np.dot(v.T, np.dot(np.diag(s), v))
    V = (V+h)/2
    blk_size = len(V)

    return V, blk_size



################ Read In Functions ##################

def read_anno(anno_dir, chrom, dtype={'CHROM': object, 'GeneEnd': np.int64, 'GeneStart': np.int64, 'TargetID': object}):

    print('Read gene annotation data on chromosome ' + chrom + '...')
    anno_chunks = pd.read_csv(
        anno_dir,
        sep='\t',
        header=0,
        chunksize=10000,
        iterator=True,
        dtype=dtype)

    GeneAnno = pd.concat([x[x['CHROM'] == chrom] for x in anno_chunks]).reset_index(drop=True)

    if GeneAnno.empty:
        print('There are no valid gene annotation data for chromosome ' + chrom + '\n')
        return None

    TargetID = GeneAnno.TargetID
    n_targets = TargetID.size

    return GeneAnno, TargetID, n_targets


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

    return snps


################ Match snpIDs of two data frames ##################

def filter_df_rows(df, filter_by, filter_cols):

    df = df[np.any(df[filter_cols].isin(filter_by), axis=1)].reset_index(drop=True)

    return df

def match_snp_ID(df_ref, df_sst):

    df_ref['snpID'] = get_snpIDs(df_ref, flip=False)
    df_sst['snpID'] = get_snpIDs(df_sst, flip=False)
    df_sst['snpIDflip'] = get_snpIDs(df_sst, flip=True)

    # find overlapped snpIDs
    overlap_list = np.intersect1d(df_ref['snpID'], df_sst[['snpID', 'snpIDflip']])

    if overlap_list.size:
        # filter out rows with non-matching snpIDs from the eQTL summary statistics
        df_sst = filter_df_rows(df=df_sst,
                                filter_by=overlap_list,
                                filter_cols=['snpID', 'snpIDflip'])

        # if not in df_ref.snpID, that means snpIDflip is in df_ref.snpID;
        # if snpIDflip is in df_ref.snpID, flip Zscore sign, flip A1 A2.
        df_sst['flip'] = np.where(df_sst.snpID.isin(df_ref.snpID.values), 1, -1)
        if not np.all(df_sst['flip'] == 1):
            idx = (df_sst['flip'] == -1)
            df_sst.loc[idx, ['snpID']] = df_sst.loc[idx, ['snpIDflip']].values
            df_sst.loc[idx, ['A1', 'A2']] = df_sst.loc[idx, ['A2', 'A1']].values
            df_sst['Z'] = df_sst['flip'] * df_sst['Z']

    return df_ref, df_sst, overlap_list

