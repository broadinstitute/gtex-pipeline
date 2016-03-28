#!/usr/bin/env python3
# Author: Francois Aguet

import argparse
import os
import struct
import subprocess
from datetime import datetime

parser = argparse.ArgumentParser(description='Convert BAM to FASTQ using SamToFastq from Picard.')
parser.add_argument('bam_file', type=str, help='BAM file')
parser.add_argument('-p', '--prefix', type=str, default='Reads', help='Prefix for output files; usually <sample_id>')
parser.add_argument('-o', '--outputdir', default=os.getcwd(), help='Directory to which FASTQs will be written')
parser.add_argument('--gzip', type=str.lower, default='1', help='gzip compression level for FASTQs; see "man gzip"')
parser.add_argument('--include_non_pf_reads', type=str.lower, choices=['true', 'false'], default='true', help='Sets INCLUDE_NON_PF_READS option (PF: passed filtering). SamToFastq default: false')
parser.add_argument('--include_non_primary_alignments', type=str.lower, choices=['true', 'false'], default='false', help='Sets INCLUDE_NON_PRIMARY_ALIGNMENTS option. SamToFastq default: false')
args = parser.parse_args()

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Starting SamToFastq', flush=True)

fastq1 = os.path.join(args.outputdir, args.prefix+'_1.fastq.gz')
fastq2 = os.path.join(args.outputdir, args.prefix+'_2.fastq.gz')
fastq0 = os.path.join(args.outputdir, args.prefix+'_unpaired.fastq.gz')

# Make named pipes for gzip
subprocess.check_call('mkfifo /tmp/read1_pipe /tmp/read2_pipe /tmp/read0_pipe', shell=True)

# Set gzip streams
subprocess.check_call('gzip -'+args.gzip+' -c < /tmp/read1_pipe > '+fastq1+' &', shell=True)
subprocess.check_call('gzip -'+args.gzip+' -c < /tmp/read2_pipe > '+fastq2+' &', shell=True)
subprocess.check_call('gzip -'+args.gzip+' -c < /tmp/read0_pipe > '+fastq0+' &', shell=True)

# SamToFastq (write to pipes)
subprocess.check_call('java -jar -Xmx8g /picard-tools/picard.jar SamToFastq INPUT='+args.bam_file\
    +' INCLUDE_NON_PF_READS='+args.include_non_pf_reads\
    +' INCLUDE_NON_PRIMARY_ALIGNMENTS='+args.include_non_primary_alignments\
    +' VALIDATION_STRINGENCY=SILENT FASTQ=/tmp/read1_pipe SECOND_END_FASTQ=/tmp/read2_pipe UNPAIRED_FASTQ=/tmp/read0_pipe', shell=True)

# Delete named pipes
subprocess.check_call('rm /tmp/read1_pipe /tmp/read2_pipe /tmp/read0_pipe', shell=True)

# Delete unpaired reads FASTQ if empty
with open(fastq0, 'rb') as f0:
    f0.seek(-4,2)
    if struct.unpack('<I', f0.read(4))[0]==0:  # empty file
        os.remove(fastq0)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finished SamToFastq', flush=True)
