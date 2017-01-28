#!/usr/bin/env python3
# Author: Francois Aguet
import pandas as pd
import numpy as np
import argparse
import os
from collections import defaultdict
import gzip

def write_gct(df, gct_path):
    with gzip.open(gct_path, 'wt') as f:
        f.write('#1.2\n')
        f.write('{0:d}\t{1:d}\n'.format(df.shape[0], 1))
        df.to_csv(f, sep='\t')

parser = argparse.ArgumentParser(description='Convert read counts from RNA-SeQC to GCT')
parser.add_argument('rpkm_gct', help='RPKM GCT from RNA-SeQC')
parser.add_argument('exon_intron_report', help='exon_intron_report from RNA-SeQC')
parser.add_argument('--exon_report', nargs=2, default='', metavar=('exon_report', 'genes_gtf'), help='Generate GCT with exon-level counts')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

sample_id = os.path.split(args.exon_intron_report)[1].split('.')[0]
print('Writing read counts to GCT format ... ', end='', flush=True)
rpkm_df = pd.read_csv(args.rpkm_gct, sep='\t', skiprows=3, header=None, usecols=[0,1], index_col=0, names=['Name','Description'])
reads_df = pd.read_csv(args.exon_intron_report, sep='\t', skiprows=1, header=None, usecols=[0,2], index_col=0, names=['Name', sample_id])
rpkm_df[sample_id] = 0
rpkm_df.loc[reads_df.index, sample_id] = reads_df[sample_id]
write_gct(rpkm_df, os.path.join(args.output_dir, sample_id+'.gene_reads.gct.gz'))
print('done.', flush=True)

if args.exon_report:
    exon_df = pd.read_csv(args.exon_report[0], sep='\t', usecols=['Exon', 'Transcript', 'Gene_Name', 'Exon_Reads'], index_col=0)
    # the exon report only contains counts for genes with at least one read -> get list of exons from GTF
    exon_count = defaultdict(int)
    gene_name = {}
    print('Parsing GTF ... ', end='', flush=True)
    with open(args.exon_report[1]) as gtf:
        for line in gtf:
            if line[0]=='#': continue
            line = line.strip().split('\t')
            attr = line[8]
            if attr[-1]==';':
                attr = attr[:-1]
            attr = dict([i.split(' ') for i in attr.replace('"','').split('; ')])
            if line[2]=='exon':
                exon_count[attr['gene_id']] += 1
            elif line[2]=='gene':
                gene_name[attr['gene_id']] = attr['gene_name']
    print('done.', flush=True)
    
    print('Writing exon read counts ... ', end='', flush=True)
    output_df = pd.DataFrame(index=[i+'_'+str(j) for i in rpkm_df.index for j in np.arange(exon_count[i])], columns=['Description', sample_id])
    output_df.index.name = 'Name'
    output_df['Description'] = rpkm_df.loc[[i.split('_')[0] for i in output_df.index], 'Description'].values
    output_df[sample_id] = 0
    output_df.loc[exon_df.index, sample_id] = exon_df['Exon_Reads']
    write_gct(output_df, os.path.join(args.output_dir, sample_id+'.exon_reads.gct.gz'))
    print('done.', flush=True)
