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


parser = argparse.ArgumentParser(description='Convert BAM to FASTQ using SamToFastq from Picard.')
parser.add_argument('input_bam', type=str, help='BAM file')
parser.add_argument('prefix', type=str, help='Prefix for output files; usually <sample_id>')
parser.add_argument('-o', '--output_dir', default=os.getcwd(), help='Output directory')
parser.add_argument('-m', '--memory', default='3', type=str, help='Memory, in GB')
parser.add_argument('--optical_duplicate_pixel_distance', default=100, help='Maximum offset between two duplicate clusters. 100 (default) is appropriate for unpatterned, 2500 recommended for patterned flowcells.')
parser.add_argument('--jar', default='/opt/picard-tools/picard.jar', help='Path to Picard jar')
args = parser.parse_args()

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Starting MarkDuplicates', flush=True)

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

with cd(args.output_dir):
    subprocess.check_call('java -jar -Xmx'+args.memory+'g '+args.jar\
        +' MarkDuplicates I='+args.input_bam\
        +' O='+os.path.split(args.input_bam)[1].replace('.bam', '.md.bam')\
        +' PROGRAM_RECORD_ID=null'\
        +' M='+args.prefix+'.marked_dup_metrics.txt'+' ASSUME_SORT_ORDER=coordinate OPTICAL_DUPLICATE_PIXEL_DISTANCE='+str(args.optical_duplicate_pixel_distance), shell=True)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finished MarkDuplicates', flush=True)
