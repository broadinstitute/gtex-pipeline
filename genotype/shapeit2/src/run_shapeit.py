#!/usr/bin/env python3
import argparse
import sys
import subprocess
import os
import pandas as pd
from datetime import datetime

def shapeit(vcf, pir, prefix, seed=1, num_threads=1, start=None, end=None, sex=None, force=True):
    """
    Wrapper for SHAPEIT

    start, end: limit to this interval if provided
    sex: required for NONPAR region

    https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html
    """
    # shapeit assemble
    cmd = f'shapeit -assemble --input-vcf {vcf} --input-pir {pir} -O {prefix} -T {num_threads} --seed {seed}'
    if start is not None and end is not None:
        cmd += f' --input-from {start} --input-to {end}'
    if sex is not None:
        cmd += f' --chrX --input-sex {sex}'
    if force:
        cmd += ' --force'
    subprocess.check_call(cmd, shell=True)

    # shapeit convert: convert to VCF
    phased_vcf = prefix+'.phased.vcf.gz'
    subprocess.check_call(f'shapeit -convert --input-haps {prefix} --output-vcf {phased_vcf} --seed {seed}', shell=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Run SHAPEIT")
    parser.add_argument('vcf', help='VCF file, must contain one chromosome only.')
    parser.add_argument('pir', help='PIR file (output from extractPIRs)')
    parser.add_argument('prefix', help='Prefix for output files, e.g., <prefix>.phased.vcf.gz')
    parser.add_argument('--sex', default=None, help='Sex annotation for the samples. TSV with sample ID: 1 (male) or 2 (female); no header')
    parser.add_argument('--par_bed', default=None, help='BED file with PAR1, PAR2 and NONPAR regions for chrX')
    parser.add_argument('--output_dir', default='.', help='Output directory')
    parser.add_argument('--num_threads', default=1, help='Number of threads')
    parser.add_argument('--seed', default=1, help='Seed for MCMC')
    args = parser.parse_args()

    # get chromosome name from VCF
    chrom = subprocess.check_output(f'tabix --list-chroms {args.vcf}', shell=True).decode().strip().split('\n')
    if len(chrom)!=1:
        raise ValueError('Input VCF must contain a single chromosome.')
    else:
        chrom = chrom[0]

    print(f"[{datetime.now().strftime('%b %d %H:%M:%S')}] Running SHAPEIT on chromosome {chrom}", flush=True)
    prefix = os.path.join(args.output_dir, args.prefix+'.'+chrom)

    if chrom.endswith('X'):
        if args.sex is None:
            raise ValueError('Sex annotation must be provided for chrX.')
        if args.par_bed is None:
            raise ValueError('BED file with PAR1, NONPAR, and PAR2 intervals must be provided for chrX.')

        par_df = pd.read_csv(args.par_bed, sep='\t', header=None, names=['chr', 'start', 'end', 'name'], index_col='name')
        par_df['start'] += 1 # change from 0-based [) to 1-based [)
        par_df['end'] += 1

        print('  * Processing PAR1', flush=True)
        shapeit(args.vcf, args.pir, prefix+'.PAR1', start=par_df.loc['PAR1', 'start'],
                end=par_df.loc['PAR1', 'end'], num_threads=args.num_threads, seed=args.seed, force=True)
        print('  * Processing NONPAR', flush=True)
        shapeit(args.vcf, args.pir, prefix+'.NONPAR', start=par_df.loc['NONPAR', 'start'],
                end=par_df.loc['NONPAR', 'end'], sex=args.sex, num_threads=args.num_threads, seed=args.seed, force=True)
        print('  * Processing PAR2', flush=True)
        shapeit(args.vcf, args.pir, prefix+'.PAR2', start=par_df.loc['PAR2', 'start'],
                end=par_df.loc['PAR2', 'end'], num_threads=args.num_threads, seed=args.seed, force=True)

        # concatenate VCFs
        print('  * Concatenating PAR1, NONPAR, PAR2', flush=True)
        # SHAPEIT outputs gzipped VCFs --> convert to bgzip, index
        subprocess.check_call('zcat {0}.PAR1.phased.vcf.gz | bgzip -c > {0}.PAR1.vcf.gz && tabix {0}.PAR1.vcf.gz'.format(prefix), shell=True)
        subprocess.check_call('zcat {0}.NONPAR.phased.vcf.gz | bgzip -c > {0}.NONPAR.vcf.gz && tabix {0}.NONPAR.vcf.gz'.format(prefix), shell=True)
        subprocess.check_call('zcat {0}.PAR2.phased.vcf.gz | bgzip -c > {0}.PAR2.vcf.gz && tabix {0}.PAR2.vcf.gz'.format(prefix), shell=True)
        # concatenate and index
        subprocess.check_call('bcftools concat {0}.PAR1.vcf.gz {0}.NONPAR.vcf.gz {0}.PAR2.vcf.gz --output-type z --output {0}.phased.vcf.gz'.format(prefix), shell=True)
        subprocess.check_call('tabix {0}.phased.vcf.gz'.format(prefix), shell=True)
        subprocess.check_call('rm {0}.PAR1.* {0}.NONPAR.* {0}.PAR2.*'.format(prefix), shell=True)
    else:
        shapeit(args.vcf, args.pir, prefix, num_threads=args.num_threads, seed=args.seed, force=True)

    print('[' + datetime.now().strftime("%b %d %H:%M:%S") + '] done.', flush=True)
