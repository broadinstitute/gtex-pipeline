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
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

if len(args.input_files)==1 and '.gct' not in args.input_files[0]:
    with open(args.input_files[0]) as f:
        paths = f.read().strip().split('\n')
else:
    paths = args.input_files

sample_ids = np.array([os.path.split(i)[1].split('.')[0] for i in paths])
i = np.argsort(sample_ids)
sample_ids = sample_ids[i]
paths = np.array(paths)[i]

# sort by sample_id
gct_df = [pd.read_csv(paths[0], sep='\t', skiprows=3, header=None, index_col=0, names=['Name','Description', sample_ids[0]])]
gct_df += [pd.read_csv(i, sep='\t', skiprows=3, header=None, usecols=[0,2], index_col=0, names=['Name','Description',j]) for i,j in zip(paths[1:], sample_ids[1:])]
gct_df = pd.concat(gct_df, axis=1)

with gzip.open(os.path.join(args.output_dir, args.output_prefix+'.gct.gz'), 'wt') as f:
    f.write('#1.2\n')
    f.write('{0:d}\t{1:d}\n'.format(gct_df.shape[0], gct_df.shape[1]-1))
    gct_df.to_csv(f, sep='\t', float_format='%.6g')
