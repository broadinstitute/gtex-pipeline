#!/usr/bin/env python3
import numpy as np
import argparse
import subprocess
from datetime import datetime
import os

parser = argparse.ArgumentParser(description="Extract Phase Informative Reads (PIR) from a BAM or CRAM file.")
parser.add_argument('--vcf', required=True, help='VCF file')
parser.add_argument('--bam', required=True, help='BAM or CRAM file')
parser.add_argument('--participant_id', required=True, help='Participant/sample ID in the VCF corresponding to BAM')
parser.add_argument('--fasta', default=None, help='Reference genome, required for CRAM input')
parser.add_argument('--chr', default=None, help='Chromosome to process')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running bam_to_pir', flush=True)

# convert CRAM to BAM
if args.bam.endswith('.cram'):
    if args.fasta is None:
        raise ValueError('A reference FASTA must be provided for CRAM files.')
    bam_file = os.path.join(args.output_dir, os.path.split(args.bam)[1].replace('.cram', '.bam'))
    cmd = 'samtools view -bh -T '+args.fasta+' -o '+bam_file+' '+args.bam
    if args.chr is not None:
        cmd += ' '+args.chr
        print('Converting CRAM to BAM ({})'.format(args.chr), flush=True)
    else:
        print('Converting CRAM to BAM', flush=True)
    subprocess.check_call(cmd, shell=True)
    subprocess.check_call('samtools index '+bam_file, shell=True)
else:
    assert args.bam.endswith('.bam')
    bam_file = args.bam

if args.chr is None:
    # get chromosomes, check against reference sequence IDs in BAM header
    chrs = subprocess.check_output('tabix --list-chroms {}'.format(args.vcf), shell=True).decode().strip().split('\n')
    sq = set(subprocess.check_output('samtools view -H '+bam_file+' | grep "@SQ" | cut -f2 | awk -F":" \'{print $2}\'', shell=True).decode().strip().split('\n'))
    if not np.all([i in sq for i in chrs]):
        raise ValueError('Reference sequence IDs in BAM do not match VCF contig names.')
else:
    chrs = [args.chr]

# split each chromosome
chr_vcf = os.path.join(args.output_dir, 'tmp.vcf.gz')
for c in chrs:
    print('Processing chromosome {}'.format(c), flush=True)
    bam_list = os.path.join(args.output_dir, 'bam_list.txt')
    with open(bam_list, 'w') as f:
        f.write(args.participant_id+' '+bam_file+' '+c+'\n')

    # generate VCF for current chr
    print('  * subsetting VCF', flush=True)
    subprocess.check_call('tabix -h '+args.vcf+' '+c+' | bgzip > '+chr_vcf, shell=True)
    print('  * indexing VCF', flush=True)
    subprocess.check_call('tabix '+chr_vcf, shell=True)
    print('  * extracting PIRs', flush=True)
    subprocess.check_call('extractPIRs --vcf '+chr_vcf+' --bam '+bam_list+' --out '+os.path.join(args.output_dir, args.participant_id+'.'+c+'.pir'), shell=True)

os.remove(bam_list)
os.remove(chr_vcf)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] done.', flush=True)
