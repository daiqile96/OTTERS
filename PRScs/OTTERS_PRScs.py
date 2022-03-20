#!/usr/bin/env python

"""
A TWAS method that infers posterior eQTL effect sizes under continuous shrinkage (CS) priors
using eQTL summary statistics and an external LD reference panel.

Usage:
python OTTERS_PRScs.py --OTTERS_dir=PATH_TO_OTTERS --anno_dir=PATH_TO_ANNO --ld_dir=PATH_TO_LD --clumped_dir=PATH_TO_CLUMPED --sst_file=SUM_STATS_FILE --out_dir=OUTPUT_DIR --chrom=CHROM
                [--window=WINDOW_SIZE --a=PARAM_A --b=PARAM_B --phi=PARAM_PHI --n_iter=MCMC_ITERATIONS --n_burnin=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR --thread=THREAD --seed=SEED]
 
 - PATH_TO_OTTERS: The directory of OTTERS source code

 - PATH_TO_ANNO: Full path and the file name of the gene annotation file. 

 - PATH_TO_LD:  Full path and the file name of the LD reference file.

 - PATH_TO_CLUMPED: Full path and the file name of the clumped eQTL file. 
                    If not be provided, all the eQTLs in the eQTL summary statistics for the target gene will be used, 
                    i.e, no clumping will be performed.

 - SUM_STATS_FILE: Full path and the file name of the summary statistics. 

 - OUTPUT_DIR: Output directory and output filename prefix of the posterior effect size estimates.

 - CHROM: An integer for the chromosome of tested genes.  

 - WINDOW (optional): Window size (in base pairs) around gene region from which to include SNPs (default: 1000000 [+- 1MB region around gene region])

 - PARAM_A (optional): Parameter a in the gamma-gamma prior in the continous shrinkage prior proposed by PRS-CS. Default is 1. 

 - PARAM_B (optional): Parameter b in the gamma-gamma prior in the continous shrinkage prior proposed by PRS-CS. Default is 0.5.

 - PARAM_PHI (optional): Global shrinkage parameter phi. If phi is not specified, it will be learnt from the data using a fully Bayesian approach.
                         This usually works well for polygenic traits with large GWAS sample sizes (hundreds of thousands of subjects).
                         For GWAS with limited sample sizes (including most of the current disease GWAS), fixing phi to 1e-4 or 1e-2,
                         or doing a small-scale grid search (e.g., phi=1e-6, 1e-4, 1e-2, 1) to find the optimal phi value often improves perdictive performance.

 - MCMC_ITERATIONS (optional): Total number of MCMC iterations. Default is 1,000.

 - MCMC_BURNIN (optional): Number of burnin iterations. Default is 500.

 - MCMC_THINNING_FACTOR (optional): Thinning of the Markov chain. Default is 5.
 
 - THREAD (optional): Number of simultaneous processes to use for parallel computation. Default is 1.

 - SEED (optional): Non-negative integer which seeds the random number generator.

"""

import multiprocessing
import subprocess
import getopt
import sys

from io import StringIO
from time import time

import numpy as np
import pandas as pd


import mcmc_gtb

############################################################
# time calculation
start_time = time()

############################################################

def parse_param():
    long_opts_list = ['OTTERS_dir=', 'anno_dir=', 'ld_dir=', 'clumped_dir=', 'sst_dir=', 'out_dir=', 'chrom=', 
                      'window=', 'a=', 'b=', 'phi=', 'n_iter=', 'n_burnin=', 'thin=', 
                      'thread=', 'seed=', 'help']

    param_dict = {'OTTERS_dir=': None, 'anno_dir=': None,'ld_dir': None, 'clumped_dir': None, 'sst_dir': None, 'out_dir': None, 'chrom': None,  
                  'window': 1000000, 'a': 1, 'b': 0.5, 'phi': None, 'n_iter': 1000, 'n_burnin': 500, 'thin': 5, 
                  'thread': 1, 'seed': None}

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
            elif opt == "--OTTERS_dir": param_dict['OTTERS_dir'] = arg
            elif opt == "--anno_dir": param_dict['anno_dir'] = arg
            elif opt == "--ld_dir": param_dict['ld_dir'] = arg
            elif opt == "--clumped_dir": param_dict['clumped_dir'] = arg
            elif opt == "--sst_dir": param_dict['sst_dir'] = arg
            elif opt == "--out_dir": param_dict['out_dir'] = arg
            elif opt == "--chrom": param_dict['chrom'] = str(arg)
            elif opt == "--window": param_dict['window'] = int(arg)
            elif opt == "--a": param_dict['a'] = float(arg)
            elif opt == "--b": param_dict['b'] = float(arg)
            elif opt == "--phi": param_dict['phi'] = float(arg)
            elif opt == "--n_iter": param_dict['n_iter'] = int(arg)
            elif opt == "--n_burnin": param_dict['n_burnin'] = int(arg)
            elif opt == "--thin": param_dict['thin'] = int(arg)
            elif opt == "--thread": param_dict['thread'] = int(arg)
            elif opt == "--seed": param_dict['seed'] = int(arg)
    else:
        print(__doc__)
        sys.exit(0)
    
    if param_dict['OTTERS_dir'] == None:
        print('* Please specify the directory to OTTERS --OTTERS_dir\n')
        sys.exit(2)
    elif param_dict['anno_dir'] == None:
        print('* Please specify the directory to the gene annotation file using --anno_dir\n')
        sys.exit(2)
    elif param_dict['ld_dir'] == None:
        print('* Please specify the directory to the LD reference file using --ld_dir\n')
        sys.exit(2)
    elif param_dict['clumped_dir'] == None:
        print('* Please specify the directory to the clumped eQTL file using --clumped_dir\n')
        sys.exit(2)
    elif param_dict['sst_dir'] == None:
        print('* Please specify the eQTL summary statistics file using --sst_dir\n')
        sys.exit(2)
    elif param_dict['chrom'] == None:
        print('* Please specify the chromosome --chrom\n')
        sys.exit(2)
    elif param_dict['out_dir'] == None:
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

print('Reading header of gene annotation file.\n')
exp_cols = ots.get_header(param_dict['anno_dir'])
e_cols_ind, e_dtype = ots.exp_cols_dtype(exp_cols)

print('Reading gene annotation data.\n')
anno_chunks = pd.read_csv(
    param_dict['anno_dir'], 
    sep='\t',
    chunksize = 10000,
    iterator = True,    
    usecols=e_cols_ind,
    dtype=e_dtype)

GeneAnno = pd.concat([x[x['CHROM']==param_dict['chrom']] for x in anno_chunks]).reset_index(drop=True)

if GeneAnno.empty:
    raise SystemExit('There are no valid gene annotation data for chromosome ' + param_dict['chrom'] + '\n')

GeneAnno = ots.optimize_cols(GeneAnno) 

TargetID = GeneAnno.TargetID 
n_targets = TargetID.size 

############################################################
    
def call_PRScs(clumped_dict, sst_dict, ld_dict, targetid):

    median_ss = int(np.ceil(np.median(sst_dict['NrSamples'])))
    ATGC = ['A', 'T', 'G', 'C']
    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    clumped_snp = ots.get_mapped_dict(clumped_dict, flip = False, map = False)
    ld_snp =  ots.get_mapped_dict(ld_dict, flip = True, map = True)
    sst_snp = ots.get_mapped_dict(sst_dict, flip = True, map = True)
    comm_snp = clumped_snp & ld_snp & sst_snp

    print('%d common SNPs in the ld reference, eQTL sumstats, and validation set for TargetID: %s ...' % (len(comm_snp),targetid))
    
    # calculate standardized beta for summary statistics 
    sst_eff = {}
    sst_pos = sst_dict['SNPPos'].astype('str').values
    sst_a1 = sst_dict['A1'].values
    sst_a2 = sst_dict['A2'].values 
    sst_Z = sst_dict['Z'].values

    for i in range(len(sst_dict.SNPPos)):
        snp = sst_pos[i]; a1 = sst_a1[i]; a2 = sst_a2[i]
        if (snp, a1, a2) in comm_snp:
            snp_id = snp+":"+a1+":"+a2
            sst_eff = ots.get_std_beta_sameN_fromZ(Z = sst_Z, sample_size = median_ss, 
            std_beta = sst_eff, snpid = snp_id, idx = i, flip = False)
        elif (snp, a2, a1) in comm_snp:
            snp_id = snp+":"+a2+":"+a1 
            sst_eff = ots.get_std_beta_sameN_fromZ(sst_Z, median_ss, sst_eff, snp_id, i, flip = True)
        elif a1 in ATGC and a2 in ATGC and (snp, mapping[a1], mapping[a2]) in comm_snp:
            snp_id = snp+":"+mapping[a1]+":"+mapping[a2]
            sst_eff = ots.get_std_beta_sameN_fromZ(sst_Z, median_ss, sst_eff, snp_id, i)
        elif a1 in ATGC and a2 in ATGC and (snp, mapping[a2], mapping[a1]) in comm_snp:
            snp_id = snp+":"+mapping[a2]+":"+mapping[a1]
            sst_eff = ots.get_std_beta_sameN_fromZ(sst_Z, median_ss, sst_eff, snp_id, i, flip = True)

    # create temp dictionary for common snps
    temp_dict = {'BETA':[], 'SNPPos':[], 'A1':[], 'A2':[], 'SNPID':[]}
    ld_pos = ld_dict['SNPPos'].astype('str').values
    ld_a1 = ld_dict['A1'].astype('str').values
    ld_a2 = ld_dict['A2'].astype('str').values
    ld_idx = []
    for (ii, snp) in enumerate(ld_pos):
        snp = ld_pos[ii]; a1 = ld_a1[ii]; a2 = ld_a2[ii]
        if (snp, a1, a2) in comm_snp:
            snp_id = snp+":"+a1+":"+a2
            ots.get_temp_dict(temp_dict, snp, a1, a2, sst_eff, snp_id)
            ld_idx.append(ii)
        elif (snp, a2, a1) in comm_snp:
            snp_id = snp+":"+a2+":"+a1 
            ots.get_temp_dict(temp_dict, snp, a1, a2, sst_eff, snp_id)
            ld_idx.append(ii)
        elif a1 in ATGC and a2 in ATGC and (snp, mapping[a1], mapping[a2]) in comm_snp:
            snp_id = snp+":"+mapping[a1]+":"+mapping[a2]
            ots.get_temp_dict(temp_dict, snp, a1, a2, sst_eff, snp_id)
            ld_idx.append(ii)
        elif a1 in ATGC and a2 in ATGC and (snp, mapping[a2], mapping[a1]) in comm_snp:
            snp_id = snp+":"+mapping[a2]+":"+mapping[a1]
            ots.get_temp_dict(temp_dict, snp, a1, a2, sst_eff, snp_id)
            ld_idx.append(ii)

    # construct covariance matrix
    blk_size, V = ots.get_cov_matrix(ld_dict.iloc[ld_idx])
    temp_ld_blk = ots.get_corr_matrix(V)
    snps_kp_id = temp_dict['SNPID']

    # create final dictionary for common snps
    sst_dict = {'CHR':[], 'SNPPos':[], 'A1':[], 'A2':[], 'BETA':[], 'FLP':[]}
    for (ii, snp) in enumerate(ld_pos):
        a1 = ld_a1[ii]; a2 = ld_a2[ii]
        if ":".join((snp,a1,a2)) in snps_kp_id:
            snp_id = ":".join((snp,a1,a2))
            sst_dict['SNPPos'].append(str(snp))
            sst_dict['CHR'].append(str(param_dict['chrom']))
            sst_dict['A1'].append(a1)
            sst_dict['A2'].append(a2)
            sst_dict['BETA'].append(sst_eff[snp_id])
            sst_dict['FLP'].append(1)
        elif ":".join((snp,a2,a1)) in snps_kp_id:
            snp_id = ":".join((snp,a2,a1))
            sst_dict['SNPPos'].append(str(snp))
            sst_dict['CHR'].append(str(param_dict['chrom']))
            sst_dict['A1'].append(a2)
            sst_dict['A2'].append(a1)
            sst_dict['BETA'].append(sst_eff[snp_id])
            sst_dict['FLP'].append(-1)
        elif a1 in ATGC and a2 in ATGC and ":".join((snp, mapping[a1], mapping[a2])) in snps_kp_id:
            snp_id = ":".join((snp, mapping[a1], mapping[a2]))
            sst_dict['SNPPos'].append(str(snp))
            sst_dict['CHR'].append(str(param_dict['chrom']))
            sst_dict['A1'].append(mapping[a1])
            sst_dict['A2'].append(mapping[a2])
            sst_dict['BETA'].append(sst_eff[snp_id])
            sst_dict['FLP'].append(1)
        elif a1 in ATGC and a2 in ATGC and ":".join((snp, mapping[a2], mapping[a1])) in snps_kp_id:
            snp_id = ":".join((snp, mapping[a2], mapping[a1]))
            sst_dict['SNPPos'].append(str(snp))
            sst_dict['CHR'].append(str(param_dict['chrom']))
            sst_dict['A1'].append(mapping[a2])
            sst_dict['A2'].append(mapping[a1])
            sst_dict['BETA'].append(sst_eff[snp_id])
            sst_dict['FLP'].append(-1)

    # handle filpped SNPs 
    flip = sst_dict['FLP']
    ld_blk = temp_ld_blk*np.outer(flip,flip)
    blk_size = len(ld_blk)    

    # call PRS-CS to run MCMC 
    if blk_size != 0:
        try:
            sst_dict['ES'] = mcmc_gtb.mcmc(a = param_dict['a'], b = param_dict['b'], phi = param_dict['phi'], 
                sst_dict = sst_dict, n = median_ss, ld_blk = ld_blk, blk_size = blk_size, 
                n_iter = param_dict['n_iter'], n_burnin = param_dict['n_burnin'], 
                thin = param_dict['thin'], seed = param_dict['seed'])

        except subprocess.CalledProcessError as err:
            print('PRS-CS failed for TargetID: ' + targetid + '\n')
            return None

    return sst_dict

############################################################

@ots.error_handler
def thread_process(num):
    print('Reading eQTL summary statistics data.')
    target = TargetID[num]
    print('num=' + str(num) + '\nTargetID=' + target)
    target_exp = GeneAnno.iloc[[num]]

    start=str(max(int(target_exp.GeneStart) - param_dict['window'],0))
    end=str(int(target_exp.GeneEnd) + param_dict['window'])

    # extract summary statistics and reference LD matrix for the gene
    print('Extracting eQTL summary statistics data.')
    sst_proc_out = ots.call_tabix(param_dict['sst_dir'], param_dict['chrom'], start, end)

    if not sst_proc_out:
        print('There is no summary statistics data for range ' + start + 'to' + end + '\n')
        return None
    
    sst_chunks = pd.read_csv(StringIO(sst_proc_out.decode('utf-8')),sep='\t',
            low_memory=False,
            header=None,
            names=["chrom", "SNPPos", "SNP", "A1", "A2", "Z", "P", "TargetID", "Gene", "NrSamples"],
            iterator=True, chunksize=1000)
    target_sst = pd.concat([chunk[chunk.TargetID == target] for chunk in sst_chunks]).reset_index(drop=True)

    if len(target_sst) == 0:
        print('There is no summary statistics data for TargetID: ' + target + '\n')
        return None 

    print('Extracting reference LD matrix.')
    ld_proc_out = ots.call_tabix(param_dict['ld_dir'], param_dict['chrom'], start, end)
    ld_chunks = pd.read_csv(StringIO(ld_proc_out.decode('utf-8')),sep="\t|:",
        engine="python",
        header=None,
        names=["index", "chrom", "SNPPos", "A2", "A1", "V1", "V2", "COV"],
        iterator=True, chunksize=1000)
    target_ld = pd.concat([chunk for chunk in ld_chunks])

    if not ld_proc_out:
        print('There is no LD reference for range ' + start + 'to' + end + '\n')
        return None 

    print('Extracting clumped eQTLs')
    if param_dict['clumped_dir'] != None:
        clumped_proc_out = ots.call_tabix(param_dict['clumped_dir'], param_dict['chrom'], start, end)
        
        if not clumped_proc_out:
            print('There is no clumped eQTLs for TargetID: ' + target + '\n')
            return None 
        
        clumped_chunks = pd.read_csv(StringIO(clumped_proc_out.decode('utf-8')),sep='\t',
        low_memory=False,
        header=None,
        names=["chrom", "SNPPos", "A1", "A2", "TargetID"],
        iterator=True, chunksize=1000)
        target_clumped = pd.concat([chunk[chunk.TargetID == target] for chunk in clumped_chunks]).reset_index(drop=True)

        if len(target_clumped) == 0:
            print('There is no clumped eQTLs for TargetID: ' + target + '\n')
            return None 
    else:
        # if the path of clumped eQTLs is not provided, all the eQTLs in summary statistics will be used
        target_clumped = target_sst
    
    try:
        target_weights = call_PRScs(target_clumped, target_sst, target_ld, target)

    except subprocess.CalledProcessError as err:
        print('PRScs failed for TargetID: ' + target + '\n')
        return None
    
    target_weights['TargetID'] = target

    if param_dict['phi'] != None:
        target_weights_path = param_dict['out_dir'] + 'eQTL_weights_a%d_b%.1f_phi%1.0e_chr%d.txt' % (param_dict['a'], param_dict['b'], param_dict['phi'], int(param_dict['chrom']))
    else:
        target_weights_path = param_dict['out_dir'] + 'eQTL_weights_a%d_b%.1f_phiauto_chr%d.txt' % (param_dict['a'], param_dict['b'], int(param_dict['chrom']))
    
    with open(target_weights_path, 'a') as ff:
        for pos, a1, a2, beta in zip(target_weights['SNPPos'], target_weights['A1'], target_weights['A2'], target_weights['ES']):
            ff.write('%d\t%d\t%s\t%s\t%s\t%.6e\n' % (int(param_dict['chrom']), int(pos), a1, a2, target, beta))

    print('Done ' + str(num) + ':' + target + '\n')

############################################################
if __name__ == '__main__':
    print('Starting summary for ' + str(n_targets) + ' target genes.\n')
    pool = multiprocessing.Pool(param_dict['thread'])
    pool.imap(thread_process,[num for num in range(n_targets)])
    pool.close()
    pool.join()
    print('Done.')

############################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = ots.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)