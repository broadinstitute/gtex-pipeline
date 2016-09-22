# Author: Francois Aguet
import pandas as pd
import argparse
import gzip
import os

parser = argparse.ArgumentParser(description='Combine GCT files')
parser.add_argument('input_files', nargs='+', help='GCT files')
parser.add_argument('col_id', choices=['expected_count', 'TPM', 'FPKM', 'IsoPct'], help='Column header')
parser.add_argument('output_prefix', help='Prefix for output file: ${prefix}.gct.gz')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

sample_ids = [os.path.split(i)[1].split('.')[0] for i in args.input_files]
ds = pd.read_csv(args.input_files[0], sep='\t',index_col=0, dtype=str)
ds = ds[[ds.columns[0], args.col_id]]
ds.rename(columns={args.col_id:sample_ids[0]}, inplace=True)
gct_df = [ds]
for i in range(1,len(args.input_files)):
    ds = pd.read_csv(args.input_files[i], sep='\t', index_col=0, dtype=str)[args.col_id]
    ds.name = sample_ids[i]
    gct_df.append(ds)
gct_df = pd.concat(gct_df, axis=1)

with gzip.open(os.path.join(args.output_dir, args.output_prefix+'.txt.gz'), 'wt') as f:
    gct_df.to_csv(f, sep='\t')
