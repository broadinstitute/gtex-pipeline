#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import gzip

parser = argparse.ArgumentParser(description='Combine significant pairs from multiple eQTL mapping runs.')
parser.add_argument('signifpair_list_file', help="File listing of 'signifpairs' outputs from eQTL pipeline.")
parser.add_argument('prefix', help='Prefix for output file: <prefix>.combined_signifpairs.txt.gz')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

with open(args.signifpair_list_file) as f:
    file_paths = f.read().strip().split('\n')

print('Loading significant pairs.')
dfs = []
for f in file_paths:
    dfs.append(pd.read_csv(f, sep='\t', usecols=['variant_id', 'gene_id']))
dfs = pd.concat(dfs, axis=0)
dfs = dfs.drop_duplicates()

print('Sorting significant pairs.')
dfs['chr'] = dfs['variant_id'].apply(lambda x: x.split('_',1)[0])
dfs['pos'] = dfs['variant_id'].apply(lambda x: int(x.split('_',2)[1]))
dfs = dfs.sort_values(['chr', 'pos', 'gene_id'])

print('Writing output.')
with gzip.open(os.path.join(args.output_dir, args.prefix+'.combined_signifpairs.txt.gz'), 'wt', compresslevel=6) as f:
    dfs.to_csv(f, sep='\t', index=False)
