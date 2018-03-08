#!/usr/bin/env python3
# Author: Francois Aguet
import pandas as pd
import numpy as np
import argparse
import gzip
import os

parser = argparse.ArgumentParser(description='Combine GCT files')
parser.add_argument('input_files', nargs='+', help='List of GCT files, or file containing paths to GCTs.')
parser.add_argument('output_prefix', help='Prefix for output file: <output_prefix>.gct.gz')
parser.add_argument('--parquet', action='store_true', help='Write to parquet format instead of GCT')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

if len(args.input_files)==1 and '.gct' not in args.input_files[0]:
    with open(args.input_files[0]) as f:
        paths = f.read().strip().split('\n')
else:
    paths = args.input_files

sample_ids = np.array([os.path.split(i)[1].split('.')[0] for i in paths])
assert len(sample_ids)==len(np.unique(sample_ids))
# sort by sample_id
i = np.argsort(sample_ids)
sample_ids = sample_ids[i]
paths = np.array(paths)[i]

df = pd.read_csv(paths[0], sep='\t', skiprows=3, header=None, index_col=0, names=['Name','Description', sample_ids[0]])
if df[sample_ids[0]].dtype==np.float64:
    dtype = np.float32
elif df[sample_ids[0]].dtype==np.int64:
    dtype = np.int32
else:
    dtype = df[sample_ids[0]].dtype.type

gct_df = pd.DataFrame(0, index=df.index, columns=['Description']+list(sample_ids), dtype=dtype)
gct_df['Description'] = df['Description']
gct_df[sample_ids[0]] = df[sample_ids[0]].astype(dtype)
for k,(i,p) in enumerate(zip(sample_ids[1:], paths[1:])):
    print('\rProcessing {}/{}'.format(k+2, len(paths)), end='', flush=True)
    df = pd.read_csv(p, sep='\t', skiprows=3, header=None, usecols=[0,2], index_col=0, names=['Name',i], dtype={'Name':str, i:dtype})
    gct_df[i] = df[i]
print()

if args.parquet:
    gct_df.to_parquet(os.path.join(args.output_dir, args.output_prefix+'.parquet'))
else:
    with gzip.open(os.path.join(args.output_dir, args.output_prefix+'.gct.gz'), 'wt', compresslevel=6) as f:
        f.write('#1.2\n')
        f.write('{0:d}\t{1:d}\n'.format(gct_df.shape[0], gct_df.shape[1]-1))
        gct_df.to_csv(f, sep='\t', float_format='%.6g')
