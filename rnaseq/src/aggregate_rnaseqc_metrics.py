#!/usr/bin/env python3
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='Aggregate RNA-SeQC metrics from multiple samples')
parser.add_argument('input_files_tsv', help='TSV file with paths to metrics files')
parser.add_argument('prefix', help='Prefix for output file name')
parser.add_argument('--annotation_headers', default='', help='Comma-separate list')
parser.add_argument('--annotation_tsvs', nargs='+', default=[])
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

annotation_headers = args.annotation_headers.split(',')
if len(annotation_headers)==1 and annotation_headers[0]=='':
    annotation_headers = []
assert len(args.annotation_tsvs)==len(annotation_headers)

path_s = pd.read_csv(args.input_files_tsv, sep='\t', index_col=0, header=None, names=['sample_id','metrics_path'], squeeze=True)
if path_s.isnull().all():  # ID not provided
    path_s = pd.Series(path_s.index, index=[os.path.split(i)[1].split('.metrics.tsv')[0] for i in path_s.index])

# check format
df = pd.read_csv(path_s.iloc[0], sep='\t', header=None)
if df.shape[0]==2:  # RNA-SeQC v1.1.9
    dfs = [pd.read_csv(i, sep='\t') for i in path_s]
elif df.shape[1]==2:  # RNA-SeQC v2
    dfs = [pd.read_csv(i, sep='\t', header=None, index_col=0).T for i in path_s]
else:
    raise ValueError('Unrecognized input format (shape {}).'.format(df.shape))
metrics_df = pd.concat(dfs, axis=0)
metrics_df.index = metrics_df['Sample']

# add optional annotations as additional columns
for h,tsv in zip(annotation_headers, args.annotation_tsvs):
    annotation_df = pd.read_csv(tsv, sep='\t', index_col=0, header=None, names=['sample_id',h])
    metrics_df[h] = annotation_df.loc[metrics_df.index, h]

metrics_df.to_csv(os.path.join(args.output_dir, args.prefix+'.metrics.tsv'), sep='\t', index=False, float_format='%.8g')
