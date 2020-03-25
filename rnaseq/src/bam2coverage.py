#!/usr/bin/env python3
# Author: Francois Aguet

import argparse
import subprocess
import os
from datetime import datetime
import gzip
import pyBigWig

parser = argparse.ArgumentParser(description='Computes coverage in bigWig and/or bedGraph format from input BAM.')
parser.add_argument('bam_file', help='Input BAM file')
parser.add_argument('chr_sizes', help='Chromosome sizes for the reference genome used')
parser.add_argument('prefix', help='Prefix for output file names')
parser.add_argument('--intersect', default=None, type=str, help='BED file containing intervals to calculate coverage on')
parser.add_argument('-f', '--format', default=['bigwig'], type=str.lower, nargs='+', choices=['bigwig', 'bedgraph'])
parser.add_argument('--sam_flags', default='-F 768', help='Flags for samtools. Default: filter out secondary and QC-failed reads')
parser.add_argument('-o', '--output', default='.', help='Output directory')
args = parser.parse_args()


print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Starting coverage computation', flush=True)

print('  * generating bedGraph from BAM', flush=True)
# bedGraph file must be sorted for compatibility with bedGraphToBigWig
bgpath = os.path.join(args.output, args.prefix + '.bedGraph')
# exclude reads that are: not primary alignment; failed platform/vendor quality checks (flags 8 & 9)
cmd = 'samtools view '+args.sam_flags+' -b '+args.bam_file+' | bedtools genomecov -bga -split -ibam - | sort -k1,1 -k2,2n > '+bgpath
subprocess.call(cmd, shell=True, executable='/bin/bash')

if 'bigwig' in args.format:
    print('  * generating bigWig from bedGraph', flush=True)
    bwpath = os.path.join(args.output, args.prefix + '.bigWig')
    cmd = 'bedGraphToBigWig '+bgpath+' '+args.chr_sizes+' '+bwpath
    subprocess.call(cmd, shell=True, executable='/bin/bash')

    if args.intersect is not None:
        print('  * calculating coverage on BED intervals', flush=True)
        bw = pyBigWig.open(bwpath)
        with gzip.open(bwpath.replace('.bigWig', '.coverage.gz'), 'wt') as f, gzip.open(args.intersect) as bed:
            f.write('\t'.join(['gene_id', 'chr', 'start', 'end', 'coverage'])+'\n')
            for line in bed:
                line = line.decode().strip().split()
                c = bw.values(line[0], int(line[1]), int(line[2]), numpy=True)
                f.write('\t'.join([line[3], line[0], line[1], line[2], ','.join(c.astype(int).astype(str))])+'\n')
        bw.close()

if 'bedgraph' not in args.format:
    os.remove(bgpath)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Done', flush=True)
