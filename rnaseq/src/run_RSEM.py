#!/usr/bin/env python3
# Author: Francois Aguet

import argparse
import os.path
import subprocess
from datetime import datetime
import contextlib

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

parser = argparse.ArgumentParser(description='Run RSEM')
parser.add_argument('rsem_ref_dir', type=str, help='Path to RSEM reference files generated with rsem-prepare-reference. The reference file prefix (for files within rsem_ref_dir) must be ''rsem_reference''')
parser.add_argument('input_file', type=str, help='BAM file or .gz.list file with paths to fastq.gz files')
parser.add_argument('prefix', help='Prefix for output file names')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
parser.add_argument('--max_frag_len', default='1000', help='Maximum fragment length')
parser.add_argument('--estimate_rspd', type=str.lower, choices=['true', 'false'], default='true', help='Set to estimate the read start position distribution from data (recommended)')
parser.add_argument('--calc_ci', type=str.lower, choices=['true', 'false'], default='false', help='Calculate 95% credibility intervals and posterior mean estimates')
parser.add_argument('--is_stranded', type=str.lower, choices=['true', 'false'], default='false', help='Stranded protocol')
parser.add_argument('--paired_end', type=str.lower, choices=['true', 'false'], default='true', help='Paired-end protocol')
parser.add_argument('-t', '--threads', default='4', help='Number of threads')
parser.add_argument('--bowtie_version', choices=['1', '2'], default='2', help='Select Bowtie version')
args = parser.parse_args()

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running RSEM', flush=True)
with cd(args.output_dir):
    cmd = 'rsem-calculate-expression --num-threads '+args.threads+' --fragment-length-max '+args.max_frag_len+' --no-bam-output'

    if args.paired_end=='true':
        cmd += ' --paired-end'

    if args.estimate_rspd=='true':
        cmd += ' --estimate-rspd'

    if args.calc_ci=='true':
        cmd += ' --calc-ci'

    if args.is_stranded=='true':
        cmd += ' --forward-prob 0.0'

    if os.path.splitext(args.input_file)[1]=='.bam':
        cmd += ' --bam '+args.input_file+' '+os.path.join(args.rsem_ref_dir,'rsem_reference')+' '+args.prefix+'.rsem'
    else:
        with open(args.input_file) as fqlist:
            fastq1 = fqlist.readline().strip()
            fastq2 = fqlist.readline().strip()
        if args.bowtie_version=='2':
            cmd += ' --bowtie2'
        cmd += ' --bowtie-chunkmbs 128 <(gunzip -c '+fastq1+') <(gunzip -c '+fastq2+') '+os.path.join(args.rsem_ref_dir,'rsem_reference')+' '+args.prefix+'.rsem'

    # run RSEM
    print('  * command: '+cmd, flush=True)
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finished RSEM', flush=True)
