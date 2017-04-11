#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import gzip
import os
import itertools
from collections import defaultdict
import tempfile


def load_isoform_results(rsem_output):
    return pd.read_csv(rsem_output, sep='\t',index_col=0,
        dtype={'transcript_id':str, 'gene_id':str, 'length':np.int32, 'effective_length':np.float32, 'expected_count':np.float32,
            'TPM':np.float32, 'FPKM':np.float32, 'IsoPct':np.float32})


def load_gene_results(rsem_output):
    return pd.read_csv(rsem_output, sep='\t',index_col=0,
        dtype={'gene_id':str, 'transcript_id(s)':str, 'length':np.float32, 'effective_length':np.float32, 'expected_count':np.float32,
            'TPM':np.float32, 'FPKM':np.float32})


def aggregate_rsem_results(file_list, col_ids):
    """
    Concatenate columns 'col_ids' from multiple RSEM output files;
    return as dict of DataFrames
    """
    df_dict = defaultdict(list)
    for k,f in enumerate(file_list):
        filename = os.path.split(f)[1]
        sample_id = filename.split('.')[0]
        print('\rProcessing RSEM output {0:d}/{1:d}'.format(k+1, len(file_list)), end='')
        if 'isoforms' in filename:
            rsem_df = load_isoform_results(f)
        elif 'genes' in filename:
            rsem_df = load_gene_results(f)
        else:
            raise ValueError('Unrecognized input format.')
        for i in col_ids:
            s = rsem_df[i]
            s.name = sample_id
            df_dict[i].append(s)

    for i in col_ids:
        df_dict[i] = pd.concat(df_dict[i], axis=1)
    return df_dict, rsem_df[[rsem_df.columns[0]]]


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Aggregate RSEM expression from multiple samples.')
    parser.add_argument('rsem_output_list', help='File listing RSEM output files, with format $sample_id.rsem.{genes|isoforms}.results')
    parser.add_argument('col_id', choices=['expected_count', 'TPM', 'FPKM', 'IsoPct'], nargs='+', help='Column header')
    parser.add_argument('output_prefix', help='Prefix for output file: ${prefix}.gct.gz')
    parser.add_argument('--chunk_size', default=250, type=int, help='Files to process simultaneously')
    parser.add_argument('--write_hdf', action='store_true', help='Also produce output in HDF5 format')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    args = parser.parse_args()


    with open(args.rsem_output_list) as f:
        file_list = f.read().strip().split('\n')
    sample_ids = [os.path.split(i)[1].split('.')[0] for i in file_list]

    if np.all(['isoform' in os.path.split(i)[1] for i in file_list]):
        prefix = args.output_prefix+'_transcripts_'
    elif np.all(['gene' in os.path.split(i)[1] for i in file_list]):
        prefix = args.output_prefix+'_genes_'
    else:
        raise ValueError('Unrecognized input format.')


    if len(sample_ids)<=args.chunk_size:
        df_dict, index_df = aggregate_rsem_results(file_list, args.col_id)
        print()
        for i in args.col_id:
            dfs = pd.concat([index_df, df_dict[i]], axis=1)
        
            if args.write_hdf:
                fname = prefix+i.lower()+'.hdf'
                print('Writing {}'.format(fname))
                dfs.to_hdf(os.path.join(args.output_dir, fname), i.lower())
            
            fname = prefix+i.lower()+'.txt.gz'
            print('Writing {}'.format(fname))
            with gzip.open(os.path.join(args.output_dir, fname), 'wt', compresslevel=6) as f:
                dfs.to_csv(f, sep='\t', float_format='%.4g')
    else:
        # process chunks
        tmp_store = tempfile.NamedTemporaryFile(dir=args.output_dir)
        nchunks = int(np.ceil(len(file_list)/args.chunk_size))
        iargs = [iter(file_list)] * args.chunk_size
        for k,sub_list in enumerate(itertools.zip_longest(*iargs)):
            sub_list = [j for j in sub_list if j is not None]
            print('Processing chunk {0:d}/{1:d}'.format(k+1, nchunks))
            df_dict, index_df = aggregate_rsem_results(sub_list, args.col_id)
            for i,j in df_dict.items():
                j.to_hdf(tmp_store.name, '{0:s}{1:d}'.format(i,k))
            print()

        # aggregate chunks
        for i in args.col_id:
            dfs = [index_df]
            for k in range(nchunks):
                print('\rProcessing chunk {0:d}/{1:d}'.format(k+1, nchunks), end='')
                dfs.append(pd.read_hdf(tmp_store.name, '{0:s}{1:d}'.format(i,k)))
            print()
            dfs = pd.concat(dfs, axis=1)
            
            if args.write_hdf:
                fname = prefix+i.lower()+'.hdf'
                print('Writing {}'.format(fname))
                dfs.to_hdf(os.path.join(args.output_dir, fname), i.lower())
            
            fname = prefix+i.lower()+'.txt.gz'
            print('Writing {}'.format(fname))
            with gzip.open(os.path.join(args.output_dir, fname), 'wt', compresslevel=6) as f:
                dfs.to_csv(f, sep='\t', float_format='%.4g')

        tmp_store.close()
