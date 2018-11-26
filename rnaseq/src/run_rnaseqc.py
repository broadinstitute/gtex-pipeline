#!/usr/bin/env python3
# Author: Francois Aguet
import argparse
import subprocess
from datetime import datetime
import os


parser = argparse.ArgumentParser(description='Wrapper for RNA-SeQC 2')
parser.add_argument('genes_gtf', type=str, help='Gene annotation GTF')
parser.add_argument('bam_file', type=str, help='BAM file')
parser.add_argument('prefix', type=str, default='Reads', help='Prefix for output files; usually sample_id')
parser.add_argument('-o', '--output_dir', default=os.getcwd(), help='Output directory')
parser.add_argument('--stranded', default=None, choices=['rf', 'fr'], help='Strandedness for stranded libraries')
parser.add_argument('--bed', default=None, help='BED file with intervals for estimating insert size distribution')
args = parser.parse_args()

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running RNA-SeQC', flush=True)

cmd = 'rnaseqc {} {} {}'.format(args.genes_gtf, args.bam_file, args.output_dir) \
    + ' -s '+args.prefix \
    + ' -vv'
if args.stranded is not None:
    cmd += ' --stranded '+args.stranded
if args.bed is not None:
    cmd += ' --bed '+args.bed
print('  * command: "{}"'.format(cmd), flush=True)
subprocess.check_call(cmd, shell=True)

# gzip GCTs
subprocess.check_call('gzip {0}.exon_reads.gct {0}.gene_tpm.gct {0}.gene_reads.gct'.format(args.prefix), shell=True)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finished RNA-SeQC', flush=True)
