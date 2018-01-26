#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import gzip
import os
import itertools
from collections import defaultdict
import tempfile


def load_isoform_results(rsem_output, cols=None):
    dtype = {'transcript_id':str, 'gene_id':str,
        'length':np.int32, 'effective_length':np.float32,
        'expected_count':np.float32, 'TPM':np.float32,
        'FPKM':np.float32, 'IsoPct':np.float32}
    if cols is None:
        return pd.read_csv(rsem_output, sep='\t',index_col=0, dtype=dtype)
    else:
        return pd.read_csv(rsem_output, sep='\t',index_col=0, dtype=dtype, usecols=['transcript_id']+cols)


def load_gene_results(rsem_output, cols=None):
    dtype = {'gene_id':str, 'transcript_id(s)':str,
        'length':np.float32, 'effective_length':np.float32,
        'expected_count':np.float32, 'TPM':np.float32, 'FPKM':np.float32}
    if cols is None:
        return pd.read_csv(rsem_output, sep='\t',index_col=0, dtype=dtype)
    else:
        return pd.read_csv(rsem_output, sep='\t',index_col=0, dtype=dtype, usecols=['gene_id']+cols)


def aggregate_rsem_results(file_list, col_ids, rsem_loader):
    """
    Concatenate columns 'col_ids' from multiple RSEM output files;
    return as dict of DataFrames
    """
    df_dict = defaultdict(list)
    for k,f in enumerate(file_list):
        filename = os.path.split(f)[1]
        sample_id = filename.split('.')[0]
        print('\rProcessing RSEM output {0:d}/{1:d}'.format(k+1, len(file_list)), end='')
        rsem_df = rsem_loader(f, cols=col_ids)
        for c in col_ids:
            s = rsem_df[c]
            s.name = sample_id
            df_dict[c].append(s)

    for c in col_ids:
        df_dict[c] = pd.concat(df_dict[c], axis=1)
    return df_dict


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Aggregate RSEM expression from multiple samples.')
    parser.add_argument('rsem_output_list', help='File listing RSEM output files, with format $sample_id.rsem.{genes|isoforms}.results')
    parser.add_argument('col_ids', choices=['expected_count', 'TPM', 'FPKM', 'IsoPct'], nargs='+', help='Column header')
    parser.add_argument('output_prefix', help='Prefix for output file: ${prefix}.txt.gz')
    parser.add_argument('--chunk_size', default=500, type=int, help='Files to process simultaneously')
    parser.add_argument('--parquet', action='store_true', help='Write to parquet format instead of txt.gz')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    args = parser.parse_args()


    with open(args.rsem_output_list) as f:
        file_list = f.read().strip().split('\n')
    sample_ids = [os.path.split(i)[1].split('.')[0] for i in file_list]

    if np.all(['isoform' in os.path.split(i)[1] for i in file_list]):
        prefix = args.output_prefix+'.rsem_transcripts_'
        rsem_loader = load_isoform_results
        index_df = load_isoform_results(file_list[0], cols=['gene_id'])
    elif np.all(['gene' in os.path.split(i)[1] for i in file_list]):
        prefix = args.output_prefix+'.rsem_genes_'
        rsem_loader = load_gene_results
        index_df = load_gene_results(file_list[0], cols=['transcript_id(s)'])
    else:
        raise ValueError('Unrecognized input format.')

    # merge outputs in chunks, store to hdf
    tmp_store = tempfile.NamedTemporaryFile(dir=args.output_dir)
    nchunks = int(np.ceil(len(file_list)/args.chunk_size))
    iargs = [iter(file_list)] * args.chunk_size
    for k,sub_list in enumerate(itertools.zip_longest(*iargs)):
        sub_list = [j for j in sub_list if j is not None]  # last chunk
        print('Processing chunk {0:d}/{1:d}'.format(k+1, nchunks))
        df_dict = aggregate_rsem_results(sub_list, args.col_ids, rsem_loader)
        for c in args.col_ids:
            df_dict[c].to_hdf(tmp_store.name, '{0:s}{1:d}'.format(c,k))
        print()

    # aggregate chunks for each output type
    for c in args.col_ids:
        dfs = [index_df]
        for k in range(nchunks):
            print('\rLoading chunk {0:d}/{1:d}'.format(k+1, nchunks), end='')
            dfs.append(pd.read_hdf(tmp_store.name, '{0:s}{1:d}'.format(c,k)))
        print()
        dfs = pd.concat(dfs, axis=1)

        if args.parquet:
            fname = prefix+c.lower()+'.parquet'
            print('Writing {}'.format(fname))
            dfs.to_parquet(os.path.join(args.output_dir, fname))

        else:  # txt.gz
            fname = prefix+c.lower()+'.txt.gz'
            print('Writing {}'.format(fname))
            with gzip.open(os.path.join(args.output_dir, fname), 'wt', compresslevel=6) as f:
                dfs.to_csv(f, sep='\t', float_format='%.5g')

    tmp_store.close()
