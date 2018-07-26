#!/usr/bin/env python3
# Author: Francois Aguet
import pandas as pd
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description='Combine covariates into a single matrix')
parser.add_argument('expression_covariates', help='')
parser.add_argument('prefix', help='')
parser.add_argument('--genotype_pcs', default=None, help='Genotype PCs')
parser.add_argument('--add_covariates', default=[], nargs='+', help='Additional covariates')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

print('Combining covariates ... ', end='', flush=True)
expression_df = pd.read_csv(args.expression_covariates, sep='\t', index_col=0, dtype=str)
if args.genotype_pcs is not None:
    genotype_df = pd.read_csv(args.genotype_pcs, sep='\t', index_col=0, dtype=str)
    combined_df = pd.concat([genotype_df[expression_df.columns], expression_df], axis=0)
else:
    combined_df = expression_df
for c in args.add_covariates:
    additional_df = pd.read_csv(c, sep='\t', index_col=0, dtype=str)
    combined_df = pd.concat([combined_df, additional_df[expression_df.columns]], axis=0)

# identify and drop colinear covariates
C = combined_df.astype(np.float64).T
Q,R = np.linalg.qr(C-np.mean(C, axis=0))
colinear_ix = np.abs(np.diag(R)) < np.finfo(np.float64).eps * C.shape[1]
if np.any(colinear_ix):
    print('Colinear covariates detected:')
    for i in C.columns[colinear_ix]:
        print("  * dropped '{}'".format(i))
    combined_df = combined_df.loc[~colinear_ix]

combined_df.to_csv(os.path.join(args.output_dir, args.prefix+'.combined_covariates.txt'), sep='\t')#, float_format='%.6g')
print('done.')
