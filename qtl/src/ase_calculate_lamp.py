#!/usr/bin/env python3
# Author: Francois Aguet
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import os


def lamp_from_readcount_list(readcount_files, coverage_cutoff=8, other_ratio_cutoff=0.05):
    """
    Determine ratio of foreign alleles across all samples from an individual
    """
    # concatenate samples
    merged_readcount_df = []
    for i,rfile in enumerate(readcount_files):
        readcount_df = pd.read_csv(rfile, sep='\t', usecols=['variantID', 'totalCount', 'otherBases'], index_col=0)
        readcount_df = readcount_df[readcount_df['totalCount']>=coverage_cutoff]
        merged_readcount_df.append(readcount_df)
    merged_readcount_df = pd.concat(merged_readcount_df, axis=0)
    
    merged_readcount_df['allCount'] = merged_readcount_df['totalCount'] + merged_readcount_df['otherBases']
    merged_readcount_df['otherRatio'] = merged_readcount_df['otherBases'] / merged_readcount_df['allCount']
    idx = merged_readcount_df['otherRatio'] < other_ratio_cutoff
    sumother = merged_readcount_df.loc[idx, 'otherBases'].sum()
    sumtotal = merged_readcount_df.loc[idx, 'allCount'].sum()
    
    return sumother / sumtotal / 2.0


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Calculate lamp')
    parser.add_argument('readcount_files', help='Output from ASE pipeline')
    parser.add_argument('prefix', help='Prefix for output files')
    parser.add_argument('-o', '--output_dir', default='.')
    args = parser.parse_args()
    
    with open(args.readcount_files) as f:
        file_list = f.read().strip().split('\n')
    sample_ids = [os.path.split(i)[1].split('.')[0] for i in file_list]
    individual_ids = ['-'.join(i.split('-')[:2]) for i in sample_ids]

    readcounts_df = pd.DataFrame(np.array([individual_ids, file_list]).T, index=sample_ids, columns=['individual_id', 'ase_read_count_file'])
    dfg = readcounts_df.groupby('individual_id')

    lamp_s = pd.Series(index=np.unique(individual_ids))
    for k,i in enumerate(lamp_s.index):
        print('Processing individual {0:d}/{1:d}'.format(k+1, lamp_s.shape[0]), end='\r')
        g = dfg.get_group(i)['ase_read_count_file']
        lamp_s.loc[i] = lamp_from_readcount_list(g.values)
    lamp_s.name = 'lamp'
    lamp_s = pd.DataFrame(lamp_s)
    lamp_s.index.name = 'individual_id'
    lamp_s.to_csv(os.path.join(args.output_dir, args.prefix+'.lamp_values.txt'), sep='\t', float_format='%.6g')

    # plot distribution
    fig = plt.figure()
    ax = fig.add_subplot(111)
    v = ax.hist(lamp_s['lamp'].values, 20, edgecolor='k', lw=0.5)
    ax.plot([lamp_s['lamp'].mean()]*2, [0, np.max(v[0])], 'k--')
    ax.plot([lamp_s['lamp'].median()]*2, [0, np.max(v[0])], '--', color=[0.5,0.5,0.5])
    ax.set_xlabel('Foreign allele read frequency ($\epsilon$)', fontsize=16)
    ax.text(lamp_s['lamp'].mean(), np.max(v[0]), r' $\bar\epsilon$ = {0:.6e}'.format(lamp_s['lamp'].mean()), fontsize=12, va='top')
    ax.text(lamp_s['lamp'].median(), 0.9*np.max(v[0]), r' median($\epsilon$) = {0:.6e}'.format(lamp_s['lamp'].median()), fontsize=12, va='top', color=[0.5,0.5,0.5])
    ax.set_ylabel('Frequency', fontsize=16)
    plt.savefig(os.path.join(args.output_dir, args.prefix+'.lamp_distribution.pdf'))
