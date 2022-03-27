#!/usr/bin/env python

#########################################################
import functools
import operator
import os
import subprocess
import sys
import traceback


from io import StringIO
from scipy.stats import norm

import pandas as pd
import numpy as np
#########################################################
## FUNCTIONS:

# error_handler

# calc_maf
# call_tabix
# call_tabix_header
# format_elapsed_time
# get_header

# exp_cols_dtype
# genofile_cols_dtype
# weight_cols_dtype
# zscore_cols_dtype

# get_snpIDs
# optimize_cols
# reformat_sample_vals
# reformat_vcf
# check_prep_vcf
# substr_in_strarray

#########################################################
# Testing functions

def check_path(path):
    if os.path.isdir(path):
        print(path + ' exists.')
    else:
        os.makedirs(path)    
        print(path + ' created.')

def call_PLINK_extract(bim_path, out_path, target, chr, start, end):

    # save the range of the gene 
    range = os.path.join(out_path, 'range.txt') 
    print(range)
    with open(range, 'w') as ff:
        ff.write('%s\t%s\t%s\t%s\n' % (chr, start, end, target))

    # extract the genotype data for this range
    out_geno = os.path.join(out_path, target) 
    cmd = ["plink --bfile "+bim_path+" --extract range "+range+" --make-bed --out "+ out_geno]
    return cmd

# extract_bim_cmd = call_PLINK(bim_path, "ZNF595", "4", "0", "188211")

def call_PLINK_clump(bim_path, r2, pvalue_path, window = 1000, p = 1, snp_field = "snpID", p_field = "P"):
    
    cmd = ["plink --bfile "+bim_path+" --clump-p1 "+str(p)+" --clump-r2 "+str(r2)+ 
           " --clump-kb "+str(window)+" --clump "+pvalue_path+
           " --clump-snp-field "+snp_field+" --clump-field "+p_field+
           " --out "+bim_path]

    return cmd


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


# returns absolute path
def get_abs_path(x): return os.path.abspath(os.path.expanduser(os.path.expandvars(x)))

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
    
# Call tabix to get header, read in lines into byt array
def call_tabix_header(path, out='tuple', rename={}):

    rename = {**{'#CHROM':'CHROM'}, **rename}

    proc = subprocess.Popen(
        ["tabix -H "+path],
        shell=True,
        stdout=subprocess.PIPE,
        bufsize=1)

    proc_out = bytearray()

    while proc.poll() is None:
        line =  proc.stdout.readline()
        if len(line) == 0:
            break
        if not line.startswith(b"##"):
            proc_out += line

    header = pd.read_csv(
        StringIO(proc_out.decode('utf-8')),
        sep='\t',
        error_bad_lines=False).rename(columns=rename)

    if out=='tuple':
        return tuple(header)

    elif out=='list':
        return list(header)

    return header

# RETURN HUMAN READABLE ELAPSED TIME STRING 
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


# get header from non-vcf file
def get_header(path, out='tuple', zipped=False, rename={}):

    compress_type = 'gzip' if zipped else None
    
    rename = {**{'#CHROM':'CHROM'}, **rename}

    header = pd.read_csv(
        path,
        sep='\t',
        header=0,
        compression=compress_type,
        low_memory=False,
        nrows=0).rename(columns=rename)

    if out=='tuple':
        return tuple(header)

    elif out=='list':
        return list(header)

    return header

# for testing on systems without tabix
def get_vcf_header(path, out='tuple'):
    
    proc = subprocess.Popen(
        ["zgrep -m1 -E 'CHROM' "+path],
        shell=True,
        stdout=subprocess.PIPE,
        bufsize=1)

    proc_out = bytearray()
    while proc.poll() is None:
        line =  proc.stdout.readline()
        if len(line) == 0:
            break
        proc_out += line

    header = pd.read_csv(
        StringIO(proc_out.decode('utf-8')),
        sep='\t',
        error_bad_lines=False).rename(columns={'#CHROM':'CHROM'})   

    if out=='tuple':
        return tuple(header)

    elif out=='list':
        return list(header)

    return header


# USED TO DETERMINE INDICES OF EXPRESSION FILE COLS TO READ IN, DTYPE OF EACH COL
def exp_cols_dtype(file_cols):
    cols = ['CHROM', 'GeneStart', 'GeneEnd', 'TargetID', 'GeneName']
    
    dtype_dict = {
        'CHROM': object,
        'GeneEnd': np.int64,
        'GeneName': object,
        'GeneStart': np.int64,
        'TargetID': object}

    file_cols_ind = tuple([file_cols.index(x) for x in cols])
    file_dtype = {file_cols.index(x):dtype_dict[x] for x in cols}

    return file_cols_ind, file_dtype

# USED TO DETERMINE INDICES OF GENOFILE COLS TO READ IN, DTYPE OF EACH COL
def genofile_cols_dtype(file_cols):
    cols = ['CHROM','POS','REF','ALT'] 
    dtype_dict = {
        'CHROM': object,
        'POS': np.int64,
        'A1': object,
        'A2': object,
        'Beta': np.int64,
        'P': np.int64}

    if type == 'vcf':
        cols.append('FORMAT')
        dtype_dict['FORMAT'] = object

    file_cols_ind = tuple([file_cols.index(x) for x in cols])
    file_dtype = {file_cols.index(x):dtype_dict[x] for x in cols}
    
    return file_cols_ind, file_dtype

# USED TO DETERMINE INDICES OF WEIGHT FILE COLS TO READ IN, DTYPE OF EACH COL
def weight_cols_dtype(file_cols, add_cols=[], drop_cols=[], get_id=True):
    cols = ['CHROM','POS','REF','ALT','TargetID','ES'] + add_cols
    dtype_dict = {
        'CHROM': object,
        'POS': np.int64,
        'REF': object,
        'ALT': object,
        'TargetID': object,
        'ES': np.float64,
        'MAF': np.float64,
        'snpID': object,
        'ID': object,
        'beta': np.float64,
        'b': np.float64}

    if get_id:
        if ('snpID' in file_cols):
            cols.append('snpID')

        elif ('ID' in file_cols):
            cols.append('ID')

    cols = [x for x in cols if (x not in drop_cols)]
    
    file_cols_ind = tuple([file_cols.index(x) for x in cols])
    file_dtype = {file_cols.index(x):dtype_dict[x] for x in cols}

    return file_cols_ind, file_dtype
# def weight_cols_dtype(file_cols, get_id=True, other_cols=[]):
#     cols = ['CHROM','POS','REF','ALT','TargetID','ES'] + other_cols
#     dtype_dict = {
#         'CHROM': object,
#         'POS': np.int64,
#         'REF': object,
#         'ALT': object,
#         'TargetID': object,
#         'ES': np.float64,
#         'MAF': np.float64,
#         'snpID': object,
#         'ID': object}

#     if get_maf:
#         cols.append('MAF')
#         dtype_dict['MAF'] = np.float64

#     if get_id:
#         if ('snpID' in file_cols):
#             cols.append('snpID')
#             dtype_dict['snpID'] = object

#         elif ('ID' in file_cols):
#             cols.append('ID')
#             dtype_dict['ID'] = object
    
#     file_cols_ind = tuple([file_cols.index(x) for x in cols])
#     file_dtype = {file_cols.index(x):dtype_dict[x] for x in cols}

#     return file_cols_ind, file_dtype

# USED TO DETERMINE INDICES OF ZSCORE FILE COLS TO READ IN, DTYPE OF EACH COL
def zscore_cols_dtype(file_cols):
    cols = ['CHROM','POS','REF','ALT','Zscore']
    dtype_dict = {
        'CHROM': object,
        'POS': np.int64,
        'REF': object,
        'ALT': object,
        'Zscore': np.float64}

    file_cols_ind = tuple([file_cols.index(x) for x in cols])
    file_dtype = {file_cols.index(x):dtype_dict[x] for x in cols}

    return file_cols_ind, file_dtype

# USED TO DETERMINE INDICES OF gwas FILE COLS TO READ IN, DTYPE OF EACH COL
def gwas_cols_dtype(file_cols):
    cols = ['CHROM','POS','REF','ALT','BETA','SE']
    dtype_dict = {
        'CHROM': object,
        'POS': np.int64,
        'REF': object,
        'ALT': object,
        'BETA': np.float64,
        'SE': np.float64}

    file_cols_ind = tuple([file_cols.index(x) for x in cols])
    file_dtype = {file_cols.index(x):dtype_dict[x] for x in cols}

    return file_cols_ind, file_dtype

# USED TO DETERMINE INDICES OF MCOV FILE COLS TO READ IN, DTYPE OF EACH COL
def MCOV_cols_dtype(file_cols, add_cols=[], drop_cols=[], get_id=True):
    cols = ['CHROM','POS','REF','ALT','COV'] + add_cols
    dtype_dict = {
        'CHROM': object,
        'POS': np.int64,
        'REF': object,
        'ALT': object,
        'COV': object,
        'snpID': object,
        'ID': object}

    if get_id:
        if ('snpID' in file_cols):
            cols.append('snpID')

        elif ('ID' in file_cols):
            cols.append('ID')

    cols = [x for x in cols if (x not in drop_cols)]
    
    file_cols_ind = tuple([file_cols.index(x) for x in cols])
    file_dtype = {file_cols.index(x):dtype_dict[x] for x in cols}
    return file_cols_ind, file_dtype

# RETURN SNP IDS
def get_snpIDs(df: pd.DataFrame, flip=False):
    chroms = df['chrom'].astype('str').values
    pos = df['SNPPos'].astype('str').values
    A1 = df['A1'].values
    A2 = df['A2'].values

    if flip:
        return ['_'.join(i) for i in zip(chroms,pos,A1,A2)]
    else:
        return ['_'.join(i) for i in zip(chroms,pos,A2,A1)]  

def get_mapped_dict(df: pd.DataFrame, flip = True, map = False):

    ATGC = ['A', 'T', 'G', 'C']
    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    snp = df['SNPPos'].astype('str').values
    a1 = df['A1'].values
    a2 = df['A2'].values
    
    all_snp = set(zip(snp, a1, a2)) 
    
    if flip:
        all_snp.update(set(zip(snp, a2, a1)))

    if map:
        for i in range(len(snp)):
            if a1[i] in ATGC and a2[i] in ATGC:
                all_snp.update({(snp[i], mapping[a1[i]], mapping[a2[i]])})
                all_snp.update({(snp[i], mapping[a2[i]], mapping[a1[i]])})

    return all_snp

# Decrease memory by downcasting 'CHROM' column to integer, integer and float columns to minimum size that will not lose info
def optimize_cols(df: pd.DataFrame):

    if 'CHROM' in df.columns:
        df['CHROM'] = df['CHROM'].astype(str).astype(int)

    ints = df.select_dtypes(include=['int64']).columns.tolist()
    df[ints] = df[ints].apply(pd.to_numeric, downcast='integer')

    floats = df.select_dtypes(include=['float64']).columns.tolist()
    df[floats] = df[floats].apply(pd.to_numeric, downcast='float')

    return df


# Reform vcf file
### input each sample genotype
### For GT Format:
###  code '0|0' or '0/0' as 0
###  code ('0|1' or '1|0')  or ('0/1' or '1/0') as 1
###  code '1|1' or '1/1' as 2
###  code '.|.' or './.' as nan(missing)

### For DS Format:
### code '.' as nan(missing)
def reformat_sample_vals(sample_col, Format):
    vals = sample_col.values
    if Format=='GT':
        vals[(vals=='0|0')|(vals=='0/0')] = 0
        vals[(vals=='1|0')|(vals=='1/0')|(vals=='0|1')|(vals=='0/1')] = 1
        vals[(vals=='1|1')|(vals=='1/1')] = 2
        vals[(vals=='.|.')|(vals=='./.')] = np.nan
        return vals.astype(np.float32)
    elif Format=='DS':
        vals[(vals==".")] = np.nan
        return vals.astype('float')


# reformats a vcf dataframe
def reformat_vcf(df: pd.DataFrame, Format, sampleID, uniqfrmts, singleformat=True):
    df = df.copy()
    if singleformat:
        val_ind = uniqfrmts[0].split(':').index(Format)
        df[sampleID]=df[sampleID].applymap(lambda x:x.split(":")[val_ind])
    else:
        # reformats sample values in row to include only specified format
        def vals_by_format(row):
            if row.needsplit:
                val_ind = row.FORMAT.split(':').index(Format)
                return row[sampleID].apply(lambda y: y.split(':')[val_ind])
            else:
                return row[sampleID]

        # specify which rows need to be reformatted
        df['needsplit'] = substr_in_strarray(':', uniqfrmts)        
        # apply to each row        
        df[sampleID] = df.apply(lambda x: vals_by_format(x), axis=1)
        df = df.drop(columns=['needsplit'])
    return df


def check_prep_vcf(df: pd.DataFrame, Format, sampleID):
    df = df.copy()

    # CHECK THAT ALL ROWS INCLUDE DATA IN THE ARGS.FORMAT FORMAT
    rowfrmts = np.unique(df.FORMAT.values)
    frmt_in_all = np.all(substr_in_strarray(Format,rowfrmts))

    if not frmt_in_all:
        raise Exception("Exception in check_prep_vcf(): Specified genotype format, format=" + Format + ", does not exist in all rows of the FORMAT column for this section of the input VCF file.")

    if rowfrmts.size > 1:
        #reformat multi
        df = reformat_vcf(df, Format, sampleID, rowfrmts, singleformat=False)
    # else assume rowfrmts.size == 1
    # if contains ':', needs to be reformatted
    elif (':' in rowfrmts[0]): 
        df = reformat_vcf(df, Format, sampleID, rowfrmts)
    # if doesnt contain ':' but isn't equivalent to Format then something's very wrong
    elif (rowfrmts[0] != Format):
        raise Exception('Exception in check_prep_vcf(): There is only one format in the FORMAT column for this section of the input VCF file, format_in_vcf='+str(rowfrmts[0]) +', which contains the specified genotype format, format=' + Format + ', but it cannot be parsed. ')

    df = df.drop(columns=['FORMAT'])

    return df


# returns a boolean array; whether substring is in a np object array
def substr_in_strarray(substr, strarray):
   return np.frompyfunc(lambda x: substr in x, 1,1)(strarray)


# CALCULATE MAF
def calc_maf(df: pd.DataFrame, sampleID, maf, filter=True, op=operator.gt):
    df = df.copy()

    # if all sample data for a row is NaN, drop the row
    drop_index = df.loc[df[sampleID].count(axis=1) == 0].index
    df = df.drop(index=drop_index)

    df['MAF'] = df[sampleID].apply(lambda x:sum(x)/(2*len(x.dropna())), axis=1)

    ### Dealing with NaN
    samp_maf_cols = np.append(sampleID,'MAF')
    df[samp_maf_cols] = df[samp_maf_cols].apply(lambda x: x.fillna(2*x.MAF),axis=1)

    # downcast floats
    df[samp_maf_cols] = df[samp_maf_cols].apply(pd.to_numeric, downcast='float')

    if filter:
        df = df[op(df.MAF, maf)].reset_index(drop=True)

    return df


# Calculate, filter p val of HWE
def calc_p_hwe(df: pd.DataFrame, sampleID, pval, filter=True, op=operator.gt):
    df = df.copy()

    df['p_HWE'] = df[sampleID].apply(lambda x:p_HWE(x.dropna()), axis=1)
    df['p_HWE'] = pd.to_numeric(df['p_HWE'], downcast='float')

    if filter:
        df = df[op(df.p_HWE, pval)].reset_index(drop=True)

    return df


### Prepare for HWE input
def p_HWE(sample_row):
    vals = sample_row.values
    if not vals.size:
        p_hwe = np.nan
    else:
        N_hets = vals[(vals>=0.5)&(vals<1.5)].size
        N_aa = vals[(vals>=0)&(vals<0.5)].size
        N_AA = vals[(vals>=1.5)&(vals<=2)].size

        p_hwe = calc_HWE(N_hets, N_AA, N_aa)

    return p_hwe


###### Calculating p-value for Hardy Weinberg Equilibrium exact test

### gij denote number of minor alleles for ith SNP and jth sample
### 0 <= gij< 0.5 denote as 0
### 0.5 <= gij < 1.5 denote as 1
### 1.5 <= gij <2 denote as 2

### Input value:
### 1.obs_hets: Observed heterozygosity = Number of 1 in each SNPs(i.e. 0.5 <= gij < 1.5)
### 2.obs_hom1: Observed AA homozygosity = Number of 0 in each SNPs(i.e. 0 <= gij< 0.5)
### 3.obs_hom2: Observed aa homozygosity = Number of 2 in each SNPs(i.e. 1.5 <= gij <= 2)

### Output: p-value for Hardy Weinberg Equilibrium exact test
def calc_HWE(obs_hets, obs_hom1, obs_hom2):
    if obs_hom1 < 0 or obs_hom2 < 0 or obs_hets < 0:
        raise Exception("FATAL ERROR - SNP-HWE: Current genotype configuration (%s  %s %s) includes negative count" % (obs_hets, obs_hom1, obs_hom2))

    obs_homc = obs_hom2 if obs_hom1 < obs_hom2 else obs_hom1
    obs_homr = obs_hom1 if obs_hom1 < obs_hom2 else obs_hom2

    rare_copies = 2*obs_homr + obs_hets
    genotypes   = obs_hets + obs_homc + obs_homr

    het_probs = [0.0]*(rare_copies + 1)

    #start at midpoint
    mid = int(rare_copies*(2*genotypes - rare_copies)/(2*genotypes))

    #check to ensure that midpoint and rare alleles have same parity
    if (rare_copies & 1)^(mid & 1):
        mid += 1

    curr_hets = mid
    curr_homr = (rare_copies - mid) / 2
    curr_homc = genotypes - curr_hets - curr_homr

    het_probs[mid] = 1.0
    sum_het_probs = float(het_probs[mid])

    for curr_hets in range(mid,1,-2):
        het_probs[curr_hets - 2] = het_probs[curr_hets]*curr_hets*(curr_hets - 1.0)/(4.0*(curr_homr + 1.0)*(curr_homc + 1.0))

        sum_het_probs += het_probs[curr_hets - 2];

        # 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
        curr_homr += 1
        curr_homc += 1

    curr_hets = mid
    curr_homr = (rare_copies - mid)/2
    curr_homc = genotypes - curr_hets - curr_homr

    for curr_hets in range(mid,rare_copies-1,2):
        het_probs[curr_hets + 2] = het_probs[curr_hets]*4.0*curr_homr*curr_homc/((curr_hets + 2.0)*(curr_hets + 1.0))

        sum_het_probs += het_probs[curr_hets + 2]

        #add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
        curr_homr -= 1
        curr_homc -= 1

    for i in range(0,rare_copies + 1):
        het_probs[i] /= sum_het_probs

    #alternate p-value calculation for p_hi/p_lo
    p_hi = float(het_probs[obs_hets])
    for i in range(obs_hets,rare_copies+1):
        p_hi += het_probs[i]

    p_lo = float(het_probs[obs_hets])
    for i in range(obs_hets-1,-1,-1):
        p_lo += het_probs[i]

    p_hi_lo = 2.0 * p_hi if p_hi < p_lo else 2.0 * p_lo

    p_hwe = 0.0
    #  p-value calculation for p_hwe
    for i in range(0,rare_copies + 1):
        if het_probs[i] > het_probs[obs_hets]:
            continue
        p_hwe += het_probs[i]

    p_hwe = 1.0 if p_hwe > 1.0 else p_hwe

    return p_hwe

# def get_cov_matrix(MCOV):

#     MCOV = MCOV.copy()
#     MCOV['COV'] =  MCOV['COV'].apply(lambda x:np.fromstring(x, dtype=np.float32, sep=','))
#     inds = MCOV.index
#     n_inds = inds.size
#     V_upper = np.zeros((n_inds, n_inds))
#     for i in range(n_inds):
#         cov_i = MCOV.COV.loc[inds[i]]
#         N = cov_i.size
#         for j in range(i,n_inds):
#             if inds[j] - inds[i] < N:
#                 V_upper[i,j] = cov_i[inds[j]-inds[i]]
#             else:
#                 V_upper[i,j] = 0
#     snp_Var = V_upper.diagonal()              
#     V = V_upper + V_upper.T - np.diag(snp_Var)

#     # Find the closest symmetric positive (semi)-definite matrix (see the proof of Theorem 2.1 in Higham's 1988 paper):
#     eignvals, eignvecs = np.linalg.eig(V)
#     print(np.min(eignvals[eignvals<0]))
#     eignvals[eignvals<0] = 0
#     pos_V = eignvecs.dot(np.diag(eignvals).dot(eignvecs.T))

#     return n_inds, pos_V

def get_psd_matrix(V):
    
    eignvals, eignvecs = np.linalg.eigh(V)
    print(np.min(eignvals[eignvals<0]))
    eignvals[eignvals<0] = 0
    pos_V = eignvecs.dot(np.diag(eignvals).dot(eignvecs.T))

    return pos_V

def get_cov_matrix(MCOV):
    MCOV = MCOV.copy()
    MCOV['COV'] =  MCOV['COV'].apply(lambda x:np.fromstring(x, dtype=np.float32, sep=','))
    inds = MCOV.index
    n_inds = inds.size
    V_upper = np.zeros((n_inds, n_inds))
    for i in range(n_inds):
        cov_i = MCOV.COV.loc[inds[i]]
        N = cov_i.size
        for j in range(i,n_inds):
            if inds[j] - inds[i] < N:
                V_upper[i,j] = cov_i[inds[j]-inds[i]]
            else:
                V_upper[i,j] = 0
    snp_Var = V_upper.diagonal()              
    V = V_upper + V_upper.T - np.diag(snp_Var)
    return n_inds, V

def get_corr_matrix(covariance):
    v = np.sqrt(np.diag(covariance))
    outer_v = np.outer(v, v)
    correlation = covariance / outer_v
    correlation[covariance == 0] = 0
    return correlation

def get_std_beta(beta, pvalue, sample_size, std_beta, snpid, idx, flip = False):
    beta_i = float(beta[idx])
    p_i = max(float(pvalue[idx]), 1e-323)
    sqrtn_i = np.sqrt(int(sample_size[idx]))
    if flip:
        beta_std_i = -1 * np.sign(beta_i)*abs(norm.ppf(p_i/2.0))/float(sqrtn_i) 
    else:
        beta_std_i =  np.sign(beta_i)*abs(norm.ppf(p_i/2.0))/float(sqrtn_i) 
    std_beta.update({snpid : beta_std_i})
    return(std_beta)

def get_std_beta_sameN(beta, pvalue, sample_size, std_beta, snpid, idx, flip = False):
    beta_i = float(beta[idx])
    p_i = max(float(pvalue[idx]), 1e-323)
    sqrtn_i = np.sqrt(sample_size)
    if flip:
        beta_std_i = -1 * np.sign(beta_i)*abs(norm.ppf(p_i/2.0))/float(sqrtn_i) 
    else:
        beta_std_i =  np.sign(beta_i)*abs(norm.ppf(p_i/2.0))/float(sqrtn_i) 
    std_beta.update({snpid : beta_std_i})
    return(std_beta)

def get_std_beta_sameN_fromZ(Z, sample_size, std_beta, snpid, idx, flip = False):
    Z_i = float(Z[idx])
    sqrtn_i = np.sqrt(sample_size)
    if flip:
        beta_std_i = -1 * (Z_i/float(sqrtn_i)) 
    else:
        beta_std_i =  (Z_i/float(sqrtn_i))
    std_beta.update({snpid : beta_std_i})
    return(std_beta)

def get_sst_dict(input_dict, pos, A1, A2, sst_beta, flip, snpid):
    input_dict['SNPPos'].append(pos)
    input_dict['A1'].append(A1)
    input_dict['A2'].append(A2)
    input_dict['BETA'].append(sst_beta[snpid])
    if flip:
        input_dict['FLP'].append(-1)
    else:
        input_dict['FLP'].append(1)
    return(input_dict)

def get_temp_dict(input_dict, pos, A1, A2, sst_beta, snpid):
    input_dict['SNPPos'].append(pos)
    input_dict['A1'].append(A1)
    input_dict['A2'].append(A2)
    input_dict['BETA'].append(sst_beta[snpid])
    input_dict['SNPID'].append(snpid)
    return(input_dict)