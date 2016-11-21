#!/usr/bin/env python3
# Author: Francois Aguet
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='Convert read counts from RNA-SeQC to GCT')
parser.add_argument('rpkm_gct', help='RPKM GCT from RNA-SeQC')
parser.add_argument('exon_intron_report', help='exon_intron_report from RNA-SeQC')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

sample_id = os.path.split(args.exon_intron_report)[1].split('.')[0]
rpkm_df = pd.read_csv(args.rpkm_gct, sep='\t', skiprows=3, header=None, usecols=[0,1], index_col=0, names=['Name','Description'])
reads_df = pd.read_csv(args.exon_intron_report, sep='\t', skiprows=1, header=None, usecols=[0,2], index_col=0, names=['Name', sample_id])
rpkm_df[sample_id] = 0
rpkm_df.loc[reads_df.index, sample_id] = reads_df[sample_id]

with open(os.path.join(args.output_dir, sample_id+'.gene_reads.gct'), 'w') as f:
    f.write('#1.2\n')
    f.write('{0:d}\t{1:d}\n'.format(rpkm_df.shape[0], 1))
    rpkm_df.to_csv(f, sep='\t')
