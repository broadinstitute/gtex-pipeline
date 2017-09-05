#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import subprocess
import os
import gzip
import feather


def load_pair_data(path):
    if path.endswith('.txt.gz'):
        return pd.read_csv(path, sep='\t', usecols=['pair_id', 'slope', 'slope_se'], index_col=0, dtype={'pair_id':str, 'slope':np.float32, 'slope_se':np.float32})
    elif path.endswith('.ft'):
        df = feather.read_dataframe(path, columns=['pair_id', 'slope', 'slope_se'])
        df.set_index('pair_id', inplace=True)
        return df
    else:
        raise ValueError('Input format not recognized.')


parser = argparse.ArgumentParser(description='Prepare METASOFT input.')
parser.add_argument('variant_gene_pair_files', help="List of variant-gene pair association result. Header must specify 'slope' and 'slope_se' columns.")
parser.add_argument('prefix', help='Prefix for output file: <prefix>.metasoft_input.[chunk000.]txt.gz')
parser.add_argument('--chunks', default=None, type=int, help='')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
parser.add_argument('--write_full', action='store_true', help='Write full input table')
args = parser.parse_args()

with open(args.variant_gene_pair_files) as f:
    paths = f.read().strip().split('\n')

sample_ids = np.array([os.path.split(i)[1].split('.')[0] for i in paths])
assert len(sample_ids)==len(np.unique(sample_ids))
# sort by sample ID
i = np.argsort(sample_ids)
sample_ids = sample_ids[i]
paths = np.array(paths)[i]

print('Reading input files')
df = load_pair_data(paths[0])

# input format: pair_id, tissue1_slope, tissue1_slope_se, tissue2_slope, tissue2_slope_s2, ...
metasoft_df = pd.DataFrame(0, index=df.index, columns=[j for i in sample_ids for j in [i+'_slope', i+'_slope_se']], dtype=np.float32)
metasoft_df[sample_ids[0]+'_slope'] = df['slope']
metasoft_df[sample_ids[0]+'_slope_se'] = df['slope_se']
for k,(i,p) in enumerate(zip(sample_ids[1:], paths[1:])):
    print('  * processing {}/{}'.format(k+2, len(paths)), flush=True)
    df = load_pair_data(p)
    metasoft_df[i+'_slope'] = df['slope']
    metasoft_df[i+'_slope_se'] = df['slope_se']
print()

print('Writing Metasoft input')
# split into chunks for parallelization
if args.chunks is not None:
    chunk_size = int(np.ceil(metasoft_df.shape[0] / args.chunks))
    for i in np.arange(args.chunks):
        print('  * writing chunk {}/{}'.format(i+1, args.chunks), flush=True)
        with gzip.open(os.path.join(args.output_dir, args.prefix+'.metasoft_input.chunk{:03d}.txt.gz'.format(i)), 'wt', compresslevel=1) as f:
            metasoft_df.iloc[i*chunk_size:(i+1)*chunk_size].to_csv(f, sep='\t', float_format='%.6g', na_rep='NA')
    print()

if args.write_full:
    print('Writing full table')
    with gzip.open(os.path.join(args.output_dir, args.prefix+'.metasoft_input.txt.gz'), 'wt', compresslevel=1) as f:
        metasoft_df.to_csv(f, sep='\t', float_format='%.6g', na_rep='NA')
