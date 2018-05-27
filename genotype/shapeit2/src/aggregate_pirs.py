#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
import argparse


def get_header_size(pir_file):
    """
    Get size of header (lines)
    """
    with open (pir_file, 'r') as f:
        line = f.readline()
    return int(line.split(' ')[1])+1


def get_header(pir_file):
    """
    Get header from PIR file
    """
    n = get_header_size(pir_file)
    with open(pir_file) as f:
        header = ''.join([next(f) for _ in range(n)])
    return header, n


def check_headers(pir_files):
    """
    Check whether headers have same size
    """
    s = [get_header_size(f) for f in pir_files]
    if len(np.unique(s))!=1:
        raise ValueError('Header sizes do not match: {}'.format(s))


def concatenate_pir_files(pir_files, output_pir):
    """
    Concatenate PIR files to 'output_pir'. All files must have the same header
    """
    check_headers(pir_files)
    header, header_lines = get_header(pir_files[0])

    with open(output_pir, 'w') as out:
        out.write(header)  # write header once

        for pir_file in pir_files:
            with open(pir_file) as f:
                for _ in range(header_lines):  # skip header
                    next(f)
                for line in f:
                    out.write(line)


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Combine output from extractPIRs run on individual samples. Equivalent to running extractPIRs on multiple samples.')
    parser.add_argument('pir_files_tsv',
        help='TSV containing paths to PIR files, in the format: [[bam1_chr1, ..., bam1_chrN], ..., [bamM_chr1, ..., bamM_chrN]].\
        \nPIR file names must be in the format <sample_id>.<chr>.pir')
    parser.add_argument('chr_list', help='File listing chromosomes to process.')
    parser.add_argument('prefix', help='Prefix for output files: <prefix>.<chr>.pir')
    parser.add_argument('-o', '--output_dir', help='Output directory')
    args = parser.parse_args()

    with open(args.chr_list) as f:
        chr_order = f.read().strip().split('\n')

    print('Sorting PIR files by chromosome.', flush=True)
    pir_files_df = pd.read_csv(args.pir_files_tsv, header=None, sep='\t')
    pir_files = pir_files_df.values.tolist()

    # sort by chromosome order
    sorted_pir_files = []
    for p in pir_files:
        pir_dict = {os.path.split(i)[1].split('.')[1]:i for i in p}
        sorted_pir_files.append([pir_dict[i] for i in chr_order])
    pir_files_df = pd.DataFrame(sorted_pir_files, columns=chr_order)

    print('Starting PIR aggregation for {} samples.'.format(pir_files_df.shape[0]), flush=True)
    for c in pir_files_df:
        print('  * processing chromosome {}'.format(c))
        chr_files = pir_files_df[c]

        chr_name = np.unique([os.path.split(i)[1].split('.')[1] for i in chr_files])
        if len(chr_name)!=1 or chr_name[0]!=c:
            raise ValueError('Chromosome names do not match for {}: {}'.format(c, chr_name))

        concatenate_pir_files(chr_files, os.path.join(args.output_dir, args.prefix+'.'+c+'.pir'))
    print('Done.', flush=True)
