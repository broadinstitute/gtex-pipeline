#!/usr/bin/env python3
# Author: Francois Aguet

import numpy as np
import pandas as pd
import gzip
import subprocess
import scipy.stats as stats
import argparse
import os


def gtf2bed(annotation_gtf, feature='gene', exclude_chrs=[]):
    """
    Parse genes from GTF, create placeholder DataFrame for BED output    
    """
    chrom = []
    start = []
    end = []
    gene_id = []
    with open(annotation_gtf, 'r') as gtf:
        for row in gtf:
            row = row.strip().split('\t')
            if row[0][0]=='#' or row[2]!=feature: continue # skip header
            chrom.append(row[0])

            # TSS: gene start (0-based coordinates for BED)
            if row[6]=='+':
                start.append(np.int64(row[3])-1)
                end.append(np.int64(row[3]))
            elif row[6]=='-':
                start.append(np.int64(row[4])-1)  # last base of gene
                end.append(np.int64(row[4]))
            else:
                raise ValueError('Strand not specified.')

            gene_id.append(row[8].split(';',1)[0].split(' ')[1].replace('"',''))

    bed_df = pd.DataFrame(data={'chr':chrom, 'start':start, 'end':end, 'gene_id':gene_id}, columns=['chr', 'start', 'end', 'gene_id'], index=gene_id)
    # drop rows corresponding to excluded chromosomes
    mask = np.ones(len(chrom), dtype=bool)
    for k in exclude_chrs:
        mask = mask & (bed_df['chr']!=k)
    bed_df = bed_df[mask]
    return bed_df
    
    
def normalize_quantiles(M, inplace=False):
    """
    Note: replicates behavior of R function normalize.quantiles from library("preprocessCore")  

    Reference:
     [1] Bolstad et al., Bioinformatics 19(2), pp. 185-193, 2003
    
    Adapted from https://github.com/andrewdyates/quantile_normalize
    """
    if not inplace:
        M = M.copy()
    
    Q = M.argsort(axis=0)
    m,n = M.shape

    # compute quantile vector
    quantiles = np.zeros(m)
    for i in range(n):
        quantiles += M[Q[:,i],i]
    quantiles = quantiles / n
    
    for i in range(n):
        # Get equivalence classes; unique values == 0
        dupes = np.zeros(m, dtype=np.int)
        for j in range(m-1):
            if M[Q[j,i],i]==M[Q[j+1,i],i]:
                dupes[j+1] = dupes[j]+1
                
        # Replace column with quantile ranks
        M[Q[:,i],i] = quantiles

        # Average together equivalence classes
        j = m-1
        while j >= 0:
            if dupes[j] == 0:
                j -= 1
            else:
                idxs = Q[j-dupes[j]:j+1,i]
                M[idxs,i] = np.median(M[idxs,i])
                j -= 1 + dupes[j]
        assert j == -1
    
    if not inplace:
        return M
        

def inverse_quantile_normalization(M):
    """
    After quantile normalization of samples, standardize expression of each gene
    """
    R = stats.mstats.rankdata(M,axis=1)  # ties are averaged
    Q = stats.norm.ppf(R/(M.shape[1]+1))
    return Q
        
        
def get_donors_from_vcf(vcfpath):
    """
    Extract donor IDs from VCF
    """
    with gzip.open(vcfpath) as vcf:
        for line in vcf:
            if line.decode()[:2]=='##': continue
            break
    return line.decode().strip().split('\t')[9:]


def normalize_expression(expression_df, counts_df, expression_threshold=0.1, count_threshold=5, min_samples=10):
    """
    Genes are thresholded based on the following expression rules:
      >=min_samples with >expression_threshold expression values
      >=min_samples with >count_threshold read counts
    """
    donor_ids = ['-'.join(i.split('-')[:2]) for i in expression_df.columns]
    
    # expression thresholds
    mask = ((np.sum(expression_df>expression_threshold,axis=1)>=min_samples) & (np.sum(counts_df>count_threshold,axis=1)>=min_samples)).values
    
    # apply normalization
    M = normalize_quantiles(expression_df.loc[mask].values, inplace=False)
    R = inverse_quantile_normalization(M)

    quant_std_df = pd.DataFrame(data=R, columns=donor_ids, index=expression_df.loc[mask].index)    
    quant_df = pd.DataFrame(data=M, columns=donor_ids, index=expression_df.loc[mask].index)
    return quant_std_df, quant_df

    
def read_gct(gct_file, donor_ids):
    """
    Load GCT as DataFrame
    """    
    df = pd.read_csv(gct_file, sep='\t', skiprows=2, index_col=0)
    df.drop('Description', axis=1, inplace=True)
    df.index.name = 'gene_id'
    return df[[i for i in df.columns if '-'.join(i.split('-')[:2]) in donor_ids]]


def write_bed(bed_df, header, output_name):
    """
    Write DataFrame to BED format
    """
    bed_df.to_csv(output_name, sep='\t', index=False, header=header)
    subprocess.check_call('bgzip -f '+output_name, shell=True, executable='/bin/bash')
    subprocess.check_call('tabix -f '+output_name+'.gz', shell=True, executable='/bin/bash')


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Generate normalized expression BED files for eQTL analyses')
    parser.add_argument('expression_gct', help='GCT file with expression in normalized units, e.g., TPM or FPKM')
    parser.add_argument('counts_gct', help='GCT file with read counts')
    parser.add_argument('annotation_gtf', help='GTF annotation')
    parser.add_argument('vcf', help='VCF file with donor IDs')
    parser.add_argument('prefix', help='Prefix for output file names')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    parser.add_argument('--expression_threshold', type=np.double, default=0.1, help='Selects genes with > expression_threshold expression in at least min_samples')
    parser.add_argument('--count_threshold', type=np.int32, default=5, help='Selects genes with > count_threshold reads in at least min_samples')
    parser.add_argument('--min_samples', type=np.int32, default=10, help='Minimum number of samples that must satisfy thresholds')
    args = parser.parse_args()
    
    print('Generating normalized expression files ... ', end='', flush=True)
    donor_ids = get_donors_from_vcf(args.vcf)
    expression_df = read_gct(args.expression_gct, donor_ids)
    counts_df = read_gct(args.counts_gct, donor_ids)

    quant_std_df, quant_df = normalize_expression(expression_df, counts_df,
        expression_threshold=args.expression_threshold, count_threshold=args.count_threshold, min_samples=args.min_samples)

    # for consistency with v6/v6p pipeline results, write unsorted expression file for PEER factor calculation
    quant_std_df.to_csv(os.path.join(args.output_dir, args.prefix+'.expression.txt'), sep='\t')

    bed_df = gtf2bed(args.annotation_gtf, feature='transcript')
    quant_std_df = pd.merge(bed_df, quant_std_df, left_index=True, right_index=True)
    quant_df = pd.merge(bed_df, quant_df, left_index=True, right_index=True)

    # sort by start position
    chr_groups = quant_std_df.groupby('chr', sort=False, group_keys=False)
    quant_std_df = chr_groups.apply(lambda x: x.sort_values('start'))
    quant_df = quant_df.loc[quant_std_df.index]

    # exclude chromosomes
    chrs = subprocess.check_output('tabix --list-chroms '+args.vcf, shell=True, executable='/bin/bash')
    chrs = chrs.decode().strip().split()
    quant_std_df = quant_std_df[quant_std_df.chr.isin(chrs)]
    quant_df = quant_df[quant_df.chr.isin(chrs)]

    # header must be commented in BED format
    header_str = quant_std_df.columns.values.copy()
    header_str[0] = '#chr'
    write_bed(quant_std_df, header_str, os.path.join(args.output_dir, args.prefix+'.expression.bed'))
    write_bed(quant_df, header_str, os.path.join(args.output_dir, args.prefix+'.expression.fpkm.bed'))
    print('done.')
