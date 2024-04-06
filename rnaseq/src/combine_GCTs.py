#!/usr/bin/env python3
# Author: Francois Aguet
import pandas as pd
import numpy as np
import argparse
import gzip
import os
import qtl.io

parser = argparse.ArgumentParser(description='Combine GCT files')
parser.add_argument('input_files', nargs='+', help='List of GCT files, or file containing paths to GCTs.')
parser.add_argument('prefix', help='Prefix for output file: <prefix>.gct.gz')
parser.add_argument('--parquet', action='store_true', help='Write to parquet format instead of GCT')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

if len(args.input_files) == 1 and '.gct' not in args.input_files[0]:
    with open(args.input_files[0]) as f:
        paths = f.read().strip().split('\n')
else:
    paths = args.input_files

path_dict = {os.path.split(i)[1].split('.')[0]:i for i in paths}
sample_ids = sorted(path_dict.keys())
assert len(sample_ids) == len(paths)

# detect format
if all([i.endswith('.parquet') for i in paths]):
    sample_id = sample_ids[0]
    df = pd.read_parquet(path_dict[sample_id])
    gct_df = pd.DataFrame(0, index=df.index, columns=['Description']+list(sample_ids), dtype=df[sample_id].dtype)
    gct_df['Description'] = df['Description']
    gct_df[sample_id] = df[sample_id]
    for k, sample_id in enumerate(sample_ids[1:], 2):
        print(f"\rProcessing {k}/{len(sample_ids)}", end='' if k < len(sample_ids) else None, flush=True)
        df = pd.read_parquet(path_dict[sample_id], columns=[sample_id])
        gct_df[sample_id] = df[sample_id]
else:  # .gct or .gct.gz
    sample_id = sample_ids[0]
    df = pd.read_csv(path_dict[sample_id], sep='\t', skiprows=3, header=None, index_col=0, names=['Name', 'Description', sample_id])
    if df[sample_id].dtype == np.float64:
        dtype = np.float32
    elif df[sample_id].dtype == np.int64:
        dtype = np.int32
    else:
        dtype = df[sample_id].dtype.type
    gct_df = pd.DataFrame(0, index=df.index, columns=['Description']+list(sample_ids), dtype=dtype)
    gct_df['Description'] = df['Description']
    gct_df[sample_id] = df[sample_id].astype(dtype)
    for k, sample_id in enumerate(sample_ids[1:], 2):
        print(f"\rProcessing {k}/{len(sample_ids)}", end='' if k < len(sample_ids) else None, flush=True)
        df = pd.read_csv(path_dict[sample_id], sep='\t', skiprows=3, header=None, usecols=[0,2],
                         index_col=0, names=['Name', sample_id], dtype={'Name':str, sample_id:dtype})
        gct_df[sample_id] = df[sample_id]

if args.parquet:
    gct_df.to_parquet(os.path.join(args.output_dir, f"{args.prefix}.gct.parquet"))
else:
    qtl.io.write_gct(gct_df, os.path.join(args.output_dir, f"{args.prefix}.gct.gz"), float_format='%.6g', compresslevel=6)
