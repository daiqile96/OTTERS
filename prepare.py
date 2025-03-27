from scipy.stats.distributions import chi2
import ottersutils as ots
import shutil
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
    extract_proc = ots.call_PLINK_extract(bim_path=geno_dir,
                                          out_path=target_dir,
                                          target=target,
                                          chrom=chrom,
                                          start_pos=start,
                                          end_pos=end)

    if not extract_proc:
        print('Remove temporary files. \n')
        shutil.rmtree(target_dir)
        return None, None, None

    ################# Read in eQTL summary statistics #####################
    target_sst = ots.read_sst(sst_file=sst_dir,
                              sst='eQTL sst',
                              target=target,
                              chrom=chrom,
                              start_pos=start,
                              end_pos=end)

    if target_sst is None:
        print('There is no estimated eQTL summary statistics')
        print('Remove temporary files. \n')
        shutil.rmtree(target_dir)
        return None, None, None

    ################# Read in SNPs in LD reference #####################
    # read in snps in LD reference panel
    target_ref = ots.read_format_ref_bim(ref_dir=target_dir,
                                         ref_file=target + '.bim')

    if target_ref is None:
        print('There is no reference bim file.')
        return None

    ########### Match SNPs in eQTL reference and LD reference #########
    # print('Check overlapped SNPs between eQTL reference and LD reference...')
    target_ref, target_sst, snp_overlap = ots.match_snp_ID_double(df_ref=target_ref,
                                                                  df_1=target_sst)

    if not snp_overlap.size:
        print('No overlapping test eQTLs')
        print('Remove temporary files for TargetID: ' + target + '.\n')
        shutil.rmtree(target_dir)
        return None, None, None

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
    print(clumped_snp)

    # filter summary statistics by clumping results
    target_sst = ots.filter_df_rows(df=target_sst,
                                    filter_by=clumped_snp,
                                    filter_cols=['snpID'])
    # prevent duplicated snpIDs
    target_sst = target_sst.drop_duplicates(subset=['snpID']).reset_index(drop=True)

    ####### Prepare Inputs for Imputation Models ########

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

