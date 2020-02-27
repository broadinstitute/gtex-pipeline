#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import os
import subprocess
import gzip


def has_dependency(name):
    return subprocess.call('which '+name, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0


def patch_bim_ids(bim_in, bim_out):
    """Modifies long indel variant IDs in BIM file for compatibility with eigensoft"""
    with open(bim_in) as f_in, open(bim_out, 'w') as f_out:
        for line in f_in:
            line = line.strip().split('\t')
            line[1] = line[0]+':'+line[1].split('_')[1]
            if len(line[4])>20:
                line[4] = 'INDEL'
            if len(line[5])>20:
                line[5] = 'INDEL'
            f_out.write('\t'.join(line)+'\n')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculates genotype PCs')
    parser.add_argument('vcf', help="VCF file. To reproduce GTEx results, this VCF should only include biallelic sites (i.e., if the input is a GTEx analysis freeze VCF, sites marked as 'wasSplit' should be filtered out first. This can be done with the command: bcftools annotate -x INFO,QUAL,^FORMAT/GT -e 'INFO/wasSplit=1' -Oz -o $out_vcf $in_vcf)")
    parser.add_argument('--keep', action='store_true', help='Keep intermediary files')
    # parser.add_argument('--prefix', default=None, help='Prefix for output file names. Default: VCF file name.')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory. Default: current directory.')
    args = parser.parse_args()

    has_plink2 = has_dependency('plink2')
    has_plink  = has_dependency('plink')
    if not (has_plink2 or has_plink):
        raise ValueError('PLINK v1.9 or v2 must be installed and in $PATH.')
    if not has_dependency('smartpca.perl'):
        raise ValueError('EIGENSOFT must be installed and in $PATH.')

    # 1) convert to PLINK format if input is VCF
    if args.vcf.endswith('.vcf.gz') or args.vcf.endswith('.vcf'):
        plink_prefix_path = os.path.join(args.output_dir, os.path.split(args.vcf)[1]).replace('.vcf.gz','')
        if has_plink2:
            cmd = 'plink2 --vcf {} --out {} --make-bed --max-alleles 2 --output-chr chrM'.format(args.vcf, plink_prefix_path)
        else:
            cmd = 'plink --vcf {} --out {} --biallelic-only strict --output-chr chrM --keep-allele-order'.format(args.vcf, plink_prefix_path)
        subprocess.check_call(cmd, shell=True)
    else:
        if not np.all([os.path.exists(args.vcf+ext) for ext in ['.bed', '.bim', '.fam']]):
            raise ValueError('Unsupported input format. Must be VCF or PLINK BED prefix (with .bed, .bim, and .fam files).')
        plink_prefix_path = args.vcf
    plink_filtered_path = os.path.join(args.output_dir, os.path.split(plink_prefix_path)[1])+'.maf05_geno01'

    # 2) filter by minor allele frequency and missingness
    # 3) prune with --indep-pairwise 200 100 0.1
    if has_dependency('plink2'):  # use PLINK 2 if available, since much faster
        cmds = [
            'plink2 --make-bed --output-chr chrM --bfile {} --maf 0.05 --geno 0.01 --out {}'.format(plink_prefix_path, plink_filtered_path),
            'plink2 --bfile {0} --indep-pairwise 200 100 0.1 --out {0}'.format(plink_filtered_path),
            # output to pgen first since --sort-vars is not yet supported for bed in PLINK 2:
            'plink2 --bfile {0} --output-chr chrM --extract {0}.prune.in --out {0}.pruned --sort-vars --make-pgen'.format(plink_filtered_path),
            'plink2 --pfile {0} --output-chr chrM --make-bed --out {0}'.format(plink_filtered_path+'.pruned'),
        ]
    else:  # use PLINK 1.9
        cmds = [
            'plink --bfile {} --maf 0.05 --geno 0.01 --make-bed --output-chr chrM --keep-allele-order --out {}'.format(plink_prefix_path, plink_filtered_path),
            'plink --bfile {0} --indep-pairwise 200 100 0.1 --out {0}'.format(plink_filtered_path),
            'plink --bfile {0} --extract {0}.prune.in --out {0}.pruned --make-bed'.format(plink_filtered_path),
        ]
    for cmd in cmds:
        subprocess.check_call(cmd, shell=True)

    # 4) patch BIM
    bim_file = plink_filtered_path+'.pruned.bim'
    os.rename(bim_file, bim_file+'.orig')
    patch_bim_ids(bim_file+'.orig', bim_file)

    # 5) run smartpca (EIGENSOFT)
    #  -i: genotypes (plink bed)
    #  -a: snps (plink bim)
    #  -b: individual (plink fam)
    #  -k: num PCs
    #  -m: num outliers
    #  -e: eigenvalues
    #  -p: plot
    #  -l: log
    subprocess.check_call('smartpca.perl -i {0}.bed -a {0}.bim -b {0}.fam -k 20 -m 0 -o {0}.pca -e {0}.eval -p {0}.plot -l {0}.log'.format(plink_filtered_path+'.pruned'), shell=True)

    # 6) delete intermediate files
    if not args.keep:
        for i in ['prune.in', 'prune.out', 'log', 'bed', 'bim', 'fam']:
            os.remove(plink_filtered_path+'.'+i)
        os.remove(plink_filtered_path+'.pruned.bim.orig')
        if os.path.exists(plink_filtered_path+'.pruned.pgen'):
            os.remove(plink_filtered_path+'.pruned.pgen')
            os.remove(plink_filtered_path+'.pruned.psam')
            os.remove(plink_filtered_path+'.pruned.pvar')
        for i in ['bed', 'bim', 'fam']:
            os.remove(plink_filtered_path+'.pruned.'+i)
        os.remove(plink_filtered_path+'.pruned.pca.par')
        os.remove(plink_filtered_path+'.pruned.plot.ps')
        os.remove(plink_filtered_path+'.pruned.plot.xtxt')
