#!/usr/bin/env python3
# Author: Francois Aguet

import argparse
import os
import struct
import subprocess
from datetime import datetime
import contextlib

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)


parser = argparse.ArgumentParser(description='Convert BAM/CRAM to FASTQ using Picard SamToFastq.')
parser.add_argument('bam_file', type=str, help='BAM or CRAM file')
parser.add_argument('-p', '--prefix', type=str, default='Reads', help='Prefix for output files; usually <sample_id>')
parser.add_argument('-o', '--output_dir', default=os.getcwd(), help='Directory to which FASTQs will be written')
parser.add_argument('-m', '--memory', default='8', type=str, help='Memory, in GB')
parser.add_argument('--reference_fasta', default=None, help='Path to reference sequence FASTA (required if input is CRAM)')
parser.add_argument('--jar', default='/opt/picard-tools/picard.jar', help='Path to Picard jar')
parser.add_argument('--gzip', type=str.lower, default='1', help='gzip compression level for FASTQs; see "man gzip"')
parser.add_argument('--include_non_pf_reads', type=str.lower, choices=['true', 'false'], default='true', help='Sets INCLUDE_NON_PF_READS option (PF: passed filtering). SamToFastq default: false')
parser.add_argument('--include_non_primary_alignments', type=str.lower, choices=['true', 'false'], default='false', help='Sets INCLUDE_NON_PRIMARY_ALIGNMENTS option. SamToFastq default: false')
args = parser.parse_args()

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Starting SamToFastq', flush=True)

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# Make named pipes for gzip
with cd(args.output_dir):
    fastq1 = args.prefix+'_1.fastq.gz'
    fastq2 = args.prefix+'_2.fastq.gz'
    fastq0 = args.prefix+'_unpaired.fastq.gz'

    subprocess.check_call('mkfifo read1_pipe read2_pipe read0_pipe', shell=True)

    # Set gzip streams
    subprocess.check_call('gzip -'+args.gzip+' -c < read1_pipe > '+fastq1+' &', shell=True)
    subprocess.check_call('gzip -'+args.gzip+' -c < read2_pipe > '+fastq2+' &', shell=True)
    subprocess.check_call('gzip -'+args.gzip+' -c < read0_pipe > '+fastq0+' &', shell=True)

    # SamToFastq (write to pipes)
    cmd = 'java -jar -Xmx'+str(int(args.memory))+'g '+args.jar+' SamToFastq INPUT='+args.bam_file\
        +' INCLUDE_NON_PF_READS='+args.include_non_pf_reads\
        +' INCLUDE_NON_PRIMARY_ALIGNMENTS='+args.include_non_primary_alignments\
        +' VALIDATION_STRINGENCY=SILENT FASTQ=read1_pipe SECOND_END_FASTQ=read2_pipe UNPAIRED_FASTQ=read0_pipe'
    if args.reference_fasta is not None:
        cmd += ' REFERENCE_SEQUENCE={}'.format(args.reference_fasta)
    subprocess.check_call(cmd, shell=True)

    # Delete named pipes
    subprocess.check_call('rm read1_pipe read2_pipe read0_pipe', shell=True)

    # Delete unpaired reads FASTQ if empty
    with open(fastq0, 'rb') as f0:
        f0.seek(-4,2)
        if struct.unpack('<I', f0.read(4))[0]==0:  # empty file
            os.remove(fastq0)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finished SamToFastq', flush=True)
