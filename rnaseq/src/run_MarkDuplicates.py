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
parser.add_argument('-m', '--memory', default=3, type=int, help='Memory, in GB')
parser.add_argument('--max_records_in_ram', default=500000, type=int,
                    help='Number of records stored in RAM before spilling to disk')
parser.add_argument('--sorting_collection_size_ratio', default=0.25, type=float)
parser.add_argument('--tagging_policy', default='DontTag', choices=['All', 'OpticalOnly', 'DontTag'])
parser.add_argument('--optical_duplicate_pixel_distance', default=100,
                    help='Maximum offset between two duplicate clusters. 100 (default) is appropriate for unpatterned, 2500 recommended for patterned flowcells.')
parser.add_argument('--jar', default='/opt/picard-tools/picard.jar', help='Path to Picard jar')
args = parser.parse_args()

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Starting MarkDuplicates', flush=True)

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

with cd(args.output_dir):
    subprocess.check_call('java -jar -Xmx{}g {}'.format(args.memory, args.jar)\
        +' MarkDuplicates I={}'.format(args.input_bam)\
        +' O={}'.format(os.path.basename(args.input_bam).replace('.bam', '.md.bam'))\
        +' PROGRAM_RECORD_ID=null'\
        +' MAX_RECORDS_IN_RAM={}'.format(args.max_records_in_ram)\
        +' SORTING_COLLECTION_SIZE_RATIO={}'.format(args.sorting_collection_size_ratio)\
        +' TMP_DIR={}'.format(args.output_dir)\
        +' M={}.marked_dup_metrics.txt'.format(args.prefix)\
        +' ASSUME_SORT_ORDER=coordinate'\
        +' TAGGING_POLICY={}'.format(args.tagging_policy)\
        +' OPTICAL_DUPLICATE_PIXEL_DISTANCE={}'.format(args.optical_duplicate_pixel_distance),
    shell=True)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finished MarkDuplicates', flush=True)
