#!/usr/bin/env python3
# Author: Francois Aguet

import argparse
import os
import subprocess
from datetime import datetime

parser = argparse.ArgumentParser(description='Run GATK ASEReadCounter')
parser.add_argument('gatk_jar', help='GATK .jar file')
parser.add_argument('genome_fasta', help='FASTA reference genome')
parser.add_argument('het_vcf', help='VCF with het sites')
parser.add_argument('bam_file', help='RNA-seq BAM file')
parser.add_argument('prefix', help='Prefix for output file names')
parser.add_argument('--minDepth', default='1', type=str, help='')
parser.add_argument('--minMappingQuality', default='255', type=str, help='')
parser.add_argument('--minBaseQuality', default='10', type=str, help='')
parser.add_argument('--disable_drf', action='store_true', help='Disables DuplicateRead filter')
args = parser.parse_args()

cmd = 'java -jar '+args.gatk_jar+' \
    -R '+args.genome_fasta+' \
    -T ASEReadCounter \
    -o '+args.prefix+'.readcounts.txt \
    -I '+args.bam_file+' \
    -sites '+args.het_vcf+' \
    -U ALLOW_N_CIGAR_READS \
    -minDepth '+args.minDepth+' \
    --minMappingQuality '+args.minMappingQuality+' \
    --minBaseQuality '+args.minBaseQuality+' \
    --allow_potentially_misencoded_quality_scores \
    -L '+args.het_vcf  # unclear why this has an effect, but results differ if omitted (bug?)

if args.disable_drf:
    cmd += ' -drf DuplicateRead'

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running GATK ASEReadCounter', flush=True)
subprocess.check_call(cmd, shell=True, executable='/bin/bash')
subprocess.check_call('gzip {}.readcounts.txt'.format(args.prefix), shell=True)
print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Done', flush=True)
