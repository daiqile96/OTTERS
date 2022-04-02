import os
import sys
import getopt

import pandas as pd
import numpy as np
import subprocess
import multiprocessing

from time import time

############################################################
# time calculation
start_time = time()

############################################################

print("Load parameters")


def parse_param():

    long_opts_list = ['OTTERS_dir=', 'genome_block=', 'b_path=',
                      'out_dir=', 'chrom=', 'thread=', 'help']

    param_dict = {'OTTERS_dir': None, 'genome_block': None, 'b_path': None,
                  'out_dir': None, 'chrom': None, 'thread': 1}

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
            elif opt == "--genome_block":
                param_dict['genome_block'] = arg
            elif opt == "--b_path":
                param_dict['b_path'] = arg
            elif opt == "--out_dir":
                param_dict['out_dir'] = arg
            elif opt == "--chrom":
                param_dict['chrom'] = str(arg)
            elif opt == "--thread":
                param_dict['thread'] = int(arg)
    else:
        print(__doc__)
        sys.exit(0)

    if param_dict['OTTERS_dir'] is None:
        print('* Please specify the directory to OTTERS --OTTERS_dir\n')
        sys.exit(2)
    elif param_dict['genome_block'] is None:
        print('* Please specify the full path to block information --genome_block\n')
        sys.exit(2)
    elif param_dict['b_path'] is None:
        print('* Please specify the directory and prefix to the binary file of LD reference panel --b_path\n')
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

param_dict = parse_param()
sys.path.append(param_dict['OTTERS_dir'])
import OTTERSutils as ots


############################################################
# define functions to convert ld matrix


def cov_fmt(x):
    return ('%.4f' % x).rstrip('0').rstrip('.')

# trim array by positionin matrix (length should be rownumber:total for each row);
# format each element in each row, join all together separated by comma


def cov_str(cov_lst):
    return [','.join([cov_fmt(x) for x in row]) for row in [cov_lst[i][i:len(cov_lst)] for i in range(len(cov_lst))]]


def call_PLINK_snplist(b_path, out_dir, blk, chr, start, end):

    # save the range of the gene
    range = os.path.join(out_dir, 'range' + str(blk+1) + '.txt')
    print(range)

    with open(range, 'w') as ff:
        ff.write('%s\t%s\t%s\t%s\n' % (chr, start, end, blk))

    out_snplist = os.path.join(out_dir, 'ldblk' + str(blk+1))

    # extract the genotype data for this range
    cmd = ["plink --bfile " + b_path + " --extract range " + range + " --write-snplist --out " + out_snplist]
    return cmd


# Read in block information
print('Reading block annotation file.')

# read in block file
Blocks = pd.read_csv(
    param_dict['genome_block'],
    sep='\t',
    usecols=['CHROM', 'Start', 'End'],
    dtype={'CHROM': object, 'Start': object, 'End': object})
Blocks = Blocks[Blocks['CHROM'] == param_dict['chrom']].reset_index(drop=True)
Blocks = ots.optimize_cols(Blocks)
n_blocks = len(Blocks)


# set chrome, b_dir, and out_dir
chr = param_dict['chrom']
b_path = param_dict['b_path']
out_dir = param_dict['out_dir']

# write columns out to file
out_cols = ['#0', 'SNP', 'CHROM', 'POS', 'COV']
out_path = os.path.join(out_dir, 'chr' + chr + '_ld.txt')
pd.DataFrame(columns=out_cols).to_csv(
    out_path,
    sep='\t',
    index=None,
    header=True,
    mode='w')


###############################################################


def thread_process(num):
    Block = Blocks.loc[num]
    blk = num + 1
    print('BLOCK=' + str(blk))

    ## Use PLINK to generate snplist and calculate LD
    print('Start to generate snplist and calculate LD')
    out = os.path.join(out_dir, 'ldblk' + str(blk))

    # generate snplist
    out_range = os.path.join(out_dir, 'range' + str(blk) + '.txt')

    with open(out_range, 'w') as ff:
        ff.write('%s\t%s\t%s\t%s\n' % (chr, Block.Start, Block.End, blk))

    get_snplist_cmd = ["plink --bfile " + b_path + " --extract range " +
                       out_range + " --write-snplist --out " + out]

    try:
        snplist_proc = subprocess.check_call(get_snplist_cmd,
                                    stdout=subprocess.PIPE,
                                    shell = True,
                                    bufsize=1)
        print('Done generate snplist for block: ' + str(num+1))
    except subprocess.CalledProcessError:
        print('Snplist generation failed for block: ' + str(num+1))
        return None

    # calculate LD
    get_ld_cmd = ["plink --bfile " + param_dict['b_path'] + " --keep-allele-order --extract range " +
                  out_range + " --r square --out " + out + " --memory 2000 "]

    try:
        ld_proc = subprocess.check_call(get_ld_cmd,
                                    stdout=subprocess.PIPE,
                                    shell = True,
                                    bufsize=1)
        print('Block LD calculation completed for block: ' + str(num+1) + '\n')
    except subprocess.CalledProcessError:
        print('Block LD calculation failed for block: ' + str(num+1) + '\n')
        return None

    ## Format the calculated LD
    print('Start to save and format the calculated LD')
    snplist = out + '.snplist'
    ld = out + '.ld'

    # read in the generated snplist
    with open(snplist) as ff:
        snplist = [line.strip() for line in ff]

    if len(snplist) > 0:

        # with open(ld) as ff:
        #     ld = [[float(val) for val in (line.strip()).split()] for line in ff]

        ld_chunks = pd.read_csv(ld, sep='\t',
                                low_memory=False,
                                header=None,
                                iterator=True,
                                chunksize=1000)

        ld_blk = pd.concat([chunk for chunk in ld_chunks]).reset_index(drop=True)

        ld_blk = ld_blk.to_numpy()

        if len(snplist) == len(ld_blk):

            print('blk %d size %s' % (blk, np.shape(ld_blk)))

            pos = [str(i.split('_')[1]) for i in snplist]

            out_df = pd.DataFrame({'snpID': snplist, 'chrom': str(chr),
                                   'pos': pos, 'cov': cov_str(ld_blk)})

            out_ld = out + '.txt'

            out_df.to_csv(
                out_ld,
                sep='\t',
                index=None,
                header=None,
                mode='w')

            print('Block LD calculation completed for block.\n')

        else:
            print('The lengths of snplist and LD are different, please check')
            return None

    else:
        print('There is no SNP in block ' + str(n_blocks) + '.\n')
        return None


############################################################

if __name__ == '__main__':
    print('Starting LD for ' + str(n_blocks) + ' target blocks.\n')
    pool = multiprocessing.Pool(param_dict['thread'])
    pool.imap(thread_process, [num for num in range(n_blocks)])
    pool.close()
    pool.join()
    print('Done.')

############################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = ots.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)

