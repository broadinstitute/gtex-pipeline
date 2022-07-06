#!/usr/bin/env python3
# Author: Francois Aguet

import argparse
import os
import subprocess
from datetime import datetime

parser = argparse.ArgumentParser(description='Run GATK ASEReadCounter')
parser.add_argument('gatk_jar', help='GATK4 .jar file')
parser.add_argument('genome_fasta', help='FASTA reference genome')
parser.add_argument('het_vcf', help='VCF with het sites (biallelic only)')
parser.add_argument('bam_file', help='RNA-seq BAM file filtered for biased reads, e.g., using STAR/WASP (grep -v "vW:i:[2-7]")')
parser.add_argument('prefix', help='Prefix for output file names')
parser.add_argument('--min_depth', default=1, type=int, help='Minimum depth after filters')
parser.add_argument('--min_mapping_quality', default=255, type=int, help='Minimum read mapping quality (255 for unique mapping reads from STAR).')
parser.add_argument('--min_base_quality', default=10, type=int, help='Minimum base quality')
parser.add_argument('--disable_drf', action='store_true', help='Disable DuplicateRead filter')
args = parser.parse_args()

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running GATK ASEReadCounter', flush=True)

cmd = f"java -jar {args.gatk_jar} \
    ASEReadCounter \
    -R {args.genome_fasta} \
    -I {args.bam_file} \
    -V {args.het_vcf} \
    -O {args.prefix}.readcounts.txt \
    -min-depth {args.min_depth} \
    --min-mapping-quality {args.min_mapping_quality} \
    --min-base-quality {args.min_base_quality}"

if args.disable_drf:
    cmd += ' -DF NotDuplicateReadFilter'

subprocess.check_call(cmd, shell=True, executable='/bin/bash')
subprocess.check_call(f'gzip {args.prefix}.readcounts.txt', shell=True)
print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Done', flush=True)
