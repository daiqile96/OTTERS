
from scipy.stats.distributions import chi2
from io import StringIO
import ottersutils as ots
import subprocess
import shutil
import pandas as pd
import numpy as np
import os


def prepare(target, target_anno, chrom, window,
            geno_dir, out_dir, sst_dir, clump_r2):

    ################# PLINK Binary Files #####################

    print('...Prepare Inputs...')
    start = str(max(int(target_anno.GeneStart) - window, 0))
    end = str(int(target_anno.GeneEnd) + window)

    # print('*Making PLINK binary files for the target gene...*')

    # set output path of the target gene
    target_dir = os.path.join(out_dir, target)
    ots.check_path(target_dir)

    # generate command to call PLINK to extract the binary file for the target gene
    extract_cmd = ots.call_PLINK_extract(geno_dir, target_dir, target, chrom, start, end)

    try:
        proc = subprocess.check_call(extract_cmd,
                                     stdout=subprocess.PIPE,
                                     shell=True)
        # print('Done making binary files.')
    except subprocess.CalledProcessError:
        print('There is no genotype data for TargetID: ' + target)
        print('Remove binary files for TargetID: ' + target + '.\n')
        shutil.rmtree(target_dir)
        return None

    ################# Read in eQTL summary statistics #####################

    # print('*Reading eQTL summary statistics data...*')

    sst_proc_out = ots.call_tabix(sst_dir, chrom, start, end)

    if not sst_proc_out:
        print('There is no summary statistics data for the range of ' + target)
        print('Remove binary files for TargetID: ' + target + '.\n')
        shutil.rmtree(target_dir)
        return None

    sst_chunks = pd.read_csv(StringIO(sst_proc_out.decode('utf-8')), sep='\t',
                             low_memory=False,
                             header=None,
                             names=["CHROM", "POS", "A1", "A2", "Z", "TargetID", "N"],
                             iterator=True,
                             chunksize=1000,
                             dtype={'P': np.float64, 'N': np.int32})

    target_sst = pd.concat([chunk[chunk.TargetID == target] for chunk in sst_chunks]).reset_index(drop=True)

    if len(target_sst) == 0:
        print('There is no summary statistics data for TargetID: ' + target)
        print('Remove binary files for TargetID: ' + target + '.\n')
        shutil.rmtree(target_dir)
        return None

    ################# Read in SNPs in LD reference #####################

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

    # here we format snpID in LD reference bim file...')
    target_ref['snpID'] = ots.get_snpIDs(target_ref, flip=False)
    ots.output_data(out_df=target_ref,
                    out_cols=['CHROM', 'snpID', 'bp', 'POS', 'A1', 'A2'],
                    out_path=target_bim,
                    write_mode='w')

    ########### Match SNPs in eQTL reference and LD reference #########
    # print('Check overlapped SNPs between eQTL reference and LD reference...')
    target_ref, target_sst, snp_overlap = ots.match_snp_ID(df_ref=target_ref,
                                                           df_sst=target_sst)

    if not snp_overlap.size:
        print('No overlapping test eQTLs')
        print('Remove binary files for TargetID: ' + target + '.\n')
        shutil.rmtree(target_dir)
        return None

    # else:
        # print('Find %s overlapped eQTLs' % (len(snp_overlap)))

    ################ Calculate median sample size ####################
    # print('*Calculate median sample size of eQTLs...*')
    median_N = np.nanmedian(target_sst['N'])

    ################# LD clumping #############################
    # print('*Perform LD clumping...*')

    # generate summary statistics of p-value to perform LD-clumping
    # chi.sf returns 0 for large Z scores
    # we divided it by 5 here to shrink Zscore to prevent p-values = 0
    target_sst['P'] = chi2.sf(np.power(target_sst['Z']/5, 2), 1)
    ots.output_data(out_df=target_sst,
                    out_cols=['snpID', 'A1', 'A2', 'P'],
                    out_path=os.path.join(target_dir, target + '.pvalue'),
                    write_mode='w',
                    out_header=True)

    # use PLINK to perform LD-clumping 
    ots.call_PLINK_clump(target, clump_r2, target + '.pvalue', target_dir, window/1000)

    # read in remaining eQTLs after LD-clumping
    target_clumped = os.path.join(target_dir, target + '.clumped')
    clumped_snp = ots.read_in_clumped(clumped_file=target_clumped,
                                      chrom=chrom)

    # filter summary statistics by clumping results
    target_sst = ots.filter_df_rows(df=target_sst,
                                    filter_by=clumped_snp,
                                    filter_cols=['snpID'])

    ################# Prepare Input Summary Statistics for Imputation Models #############################

    # print('*Start prepare input summary statistics...*')

    # Prepare Zscore input for SDPR
    # print('Done generating Zscore.')
    ots.output_data(out_df=target_sst,
                    out_cols=['snpID', 'A1', 'A2', 'Z'],
                    out_path=os.path.join(target_dir, target+'_Zscore.txt'),
                    write_mode='w',
                    out_header=['SNP', "A1", "A2", "Z"])

    # Prepare the standardized beta input for lassosum and P+T
    target_sst['Beta'] = target_sst['Z']/np.sqrt(median_N)
    ots.output_data(out_df=target_sst,
                    out_cols=['snpID', 'A1', 'A2', 'Beta'],
                    out_path=os.path.join(target_dir, target+'_beta.txt'),
                    write_mode='w',
                    out_header=['SNP', "A1", "A2", "Beta"])

    # print('Done generating standardized beta.')
    print('Done prepare inputs.')

    return target_dir, target_sst, median_N

