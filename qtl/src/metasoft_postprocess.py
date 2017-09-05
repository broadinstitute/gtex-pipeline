#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import os
import gzip
import feather


def load_pair_data(path):
    if path.endswith('.txt.gz'):
        return pd.read_csv(path, sep='\t', usecols=['pair_id', 'pval_nominal'], index_col=0, dtype={'pair_id':str, 'pval_nominal':np.float64})
    elif path.endswith('.ft'):
        df = feather.read_dataframe(path, columns=['pair_id', 'pval_nominal'])
        df.set_index('pair_id', inplace=True)
        return df
    else:
        raise ValueError('Input format not recognized.')


parser = argparse.ArgumentParser(description='METASOFT post-processing.')
parser.add_argument('metasoft_output_chunks', help="List of metasoft outputs.")
parser.add_argument('metasoft_pairs', help="List of metasoft input pair files (containing nominal p-values).")
parser.add_argument('prefix', help='Prefix for output file: <prefix>.metasoft_input.[chunk000.]txt.gz')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

with open(args.metasoft_output_chunks) as f:
    chunk_paths = f.read().strip().split('\n')

with open(args.metasoft_pairs) as f:
    pair_paths = f.read().strip().split('\n')

# parse sample IDs from pairs
sample_ids = np.array([os.path.split(i)[1].split('.')[0] for i in pair_paths])
assert len(sample_ids)==len(np.unique(sample_ids))
# sort by sample ID
i = np.argsort(sample_ids)
sample_ids = sample_ids[i]
pair_paths = np.array(pair_paths)[i]

# prepare header
header = ['RSID', '#STUDY', 'PVALUE_FE', 'BETA_FE', 'STD_FE', 'PVALUE_RE', 'BETA_RE', 'STD_RE',
    'PVALUE_RE2', 'STAT1_RE2', 'STAT2_RE2', 'PVALUE_BE', 'I_SQUARE', 'Q', 'PVALUE_Q', 'TAU_SQUARE']
header += ['pval_'+i for i in sample_ids] + ['mval_'+i for i in sample_ids]

# concatenate chunks
print('Loading chunks')
output_df = []
for i,c in enumerate(chunk_paths):
    print('  * loading chunk {}/{}'.format(i+1, len(chunk_paths)), flush=True)
    output_df.append(
        pd.read_csv(c, sep='\t', header=None, skiprows=1, names=header, usecols=header, index_col=0)
    )
print('Concatenating chunks')
output_df = pd.concat(output_df, axis=0)

# sort chunks by input order
pair_df = load_pair_data(pair_paths[0])
if not np.all(output_df.index==pair_df.index):
    print('Sorting output')
    output_df = output_df.loc[pair_df.index]

print('Substituting p-values')
output_df['pval_'+sample_ids[0]] = pair_df['pval_nominal']
for i,p in zip(sample_ids[1:], pair_paths[1:]):
    pair_df = load_pair_data(p)
    output_df['pval_'+i] = pair_df['pval_nominal']

print('Writing output')
with gzip.open(os.path.join(args.output_dir, args.prefix+'.metasoft.txt.gz'), 'wt', compresslevel=6) as f:
    output_df.to_csv(f, sep='\t', float_format='%.6g', na_rep='NA')
