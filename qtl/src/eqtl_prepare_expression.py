#!/usr/bin/env python3
# Author: Francois Aguet

import numpy as np
import pandas as pd
import gzip
import subprocess
import scipy.stats as stats
import argparse
import os
import feather

import rnaseqnorm


def gtf_to_bed(annotation_gtf, feature='gene', exclude_chrs=[]):
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


def prepare_bed(df, bed_template_df, chr_subset=None):
    bed_df = pd.merge(bed_template_df, df, left_index=True, right_index=True)
    # sort by start position
    bed_df = bed_df.groupby('chr', sort=False, group_keys=False).apply(lambda x: x.sort_values('start'))
    if chr_subset is not None:
        # subset chrs from VCF
        bed_df = bed_df[bed_df.chr.isin(chr_subset)]
    return bed_df


def write_bed(bed_df, output_name):
    """
    Write DataFrame to BED format
    """
    assert bed_df.columns[0]=='chr' and bed_df.columns[1]=='start' and bed_df.columns[2]=='end'
    # header must be commented in BED format
    header = bed_df.columns.values.copy()
    header[0] = '#chr'
    bed_df.to_csv(output_name, sep='\t', index=False, header=header)
    subprocess.check_call('bgzip -f '+output_name, shell=True, executable='/bin/bash')
    subprocess.check_call('tabix -f '+output_name+'.gz', shell=True, executable='/bin/bash')


def read_gct(gct_file, sample_ids=None, dtype=None):
    """
    Load GCT as DataFrame. The first two columns must be 'Name' and 'Description'.
    """
    if sample_ids is not None:
        sample_ids = ['Name']+list(sample_ids)

    if gct_file.endswith('.gct.gz') or gct_file.endswith('.gct'):
        if dtype is not None:
            with gzip.open(gct_file, 'rt') as gct:
                gct.readline()
                gct.readline()
                sample_ids = gct.readline().strip().split()
            dtypes = {i:dtype for i in sample_ids[2:]}
            dtypes['Name'] = str
            dtypes['Description'] = str
            df = pd.read_csv(gct_file, sep='\t', skiprows=2, usecols=sample_ids, index_col=0, dtype=dtypes)
        else:
            df = pd.read_csv(gct_file, sep='\t', skiprows=2, usecols=sample_ids, index_col=0)
    elif gct_file.endswith('.parquet'):
        df = pd.read_parquet(gct_file, columns=sample_ids)
    elif gct_file.endswith('.ft'):  # feather format
        df = feather.read_dataframe(gct_file, columns=sample_ids)
        df = df.set_index('Name')
    else:
        raise ValueError('Unsupported input format.')
    df.index.name = 'gene_id'
    if 'Description' in df.columns:
        df = df.drop('Description', axis=1)
    return df


def prepare_expression(counts_df, tpm_df, vcf_lookup_s, sample_frac_threshold=0.2, count_threshold=6, tpm_threshold=0.1, mode='tmm'):
    """
    Genes are thresholded based on the following expression rules:
      TPM >= tpm_threshold in >= sample_frac_threshold*samples
      read counts >= count_threshold in sample_frac_threshold*samples
    
    vcf_lookup: lookup table mapping sample IDs to VCF IDs
    
    Between-sample normalization modes:
      tmm: TMM from edgeR
      qn:  quantile normalization
    """

    ix = np.intersect1d(counts_df.columns, vcf_lookup_s.index)
    tpm_df = tpm_df[ix]
    counts_df = counts_df[ix]
    ns = tpm_df.shape[1]

    # expression thresholds
    mask = (
        (np.sum(tpm_df>=tpm_threshold,axis=1)>=sample_frac_threshold*ns) &
        (np.sum(counts_df>=count_threshold,axis=1)>=sample_frac_threshold*ns)
    ).values

    # apply normalization
    if mode.lower()=='tmm':
        tmm_counts_df = rnaseqnorm.edgeR_cpm(counts_df, normalized_lib_sizes=True)
        norm_df = rnaseqnorm.inverse_normal_transform(tmm_counts_df[mask])
    elif mode.lower()=='qn':
        qn_df = rnaseqnorm.normalize_quantiles(tpm_df.loc[mask])
        norm_df = rnaseqnorm.inverse_normal_transform(qn_df)
    else:
        raise ValueError('Unsupported mode {}'.format(mode))

    return norm_df



if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Generate normalized expression BED files for eQTL analyses')
    parser.add_argument('tpm_gct', help='GCT file with expression in normalized units, e.g., TPM or FPKM')
    parser.add_argument('counts_gct', help='GCT file with read counts')
    parser.add_argument('annotation_gtf', help='GTF annotation')
    parser.add_argument('sample_participant_lookup', help='Lookup table linking samples to participants')
    parser.add_argument('vcf_chr_list', help='List of chromosomes in VCF')
    parser.add_argument('prefix', help='Prefix for output file names')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    parser.add_argument('--sample_id_list', default=None, help='File listing sample IDs to include')
    parser.add_argument('--convert_tpm', action='store_true', help='Convert to TPM (in case input is in RPKM/FPKM)')
    parser.add_argument('--legacy_mode', action='store_true', help='Run in legacy mode (generates separate output for PEER factor calculation)')
    parser.add_argument('--tpm_threshold', type=np.double, default=0.1, help='Selects genes with > expression_threshold expression in at least sample_frac_threshold')
    parser.add_argument('--count_threshold', type=np.int32, default=6, help='Selects genes with >= count_threshold reads in at least sample_frac_threshold samples')
    parser.add_argument('--sample_frac_threshold', type=np.double, default=0.2, help='Minimum fraction of samples that must satisfy thresholds')
    parser.add_argument('--normalization_method', default='tmm', help='Normalization method: TMM or quantile normalization (qn)')
    args = parser.parse_args()

    print('Loading expression data', flush=True)
    sample_ids = None
    if args.sample_id_list is not None:
        with open(args.sample_id_list) as f:
            sample_ids = f.read().strip().split('\n')
            print('  * Loading {} samples'.format(len(sample_ids)), flush=True)

    counts_df = read_gct(args.counts_gct, sample_ids)
    tpm_df = read_gct(args.tpm_gct, sample_ids)
    sample_participant_lookup_s = pd.read_csv(args.sample_participant_lookup, sep='\t', index_col=0, dtype=str, squeeze=True)

    # check inputs
    if not np.all(counts_df.columns == tpm_df.columns):
        raise ValueError('Sample IDs in the TPM and read counts files must match.')
    if not np.all(counts_df.columns.isin(sample_participant_lookup_s.index)):
        raise ValueError('Sample IDs in expression files and participant lookup table must match.')

    if args.convert_tpm:
        print('  * Converting to TPM', flush=True)
        tpm_df = tpm_df/tpm_df.sum(0)*1e6

    print('Normalizing data ({})'.format(args.normalization_method), flush=True)
    norm_df = prepare_expression(counts_df, tpm_df, sample_participant_lookup_s, sample_frac_threshold=args.sample_frac_threshold,
        count_threshold=args.count_threshold, tpm_threshold=args.tpm_threshold, mode=args.normalization_method)
    print('  * {} genes in input tables.'.format(counts_df.shape[0]), flush=True)
    print('  * {} genes remain after thresholding.'.format(norm_df.shape[0]), flush=True)

    # change sample IDs to participant IDs
    norm_df.rename(columns=sample_participant_lookup_s.to_dict(), inplace=True)

    bed_template_df = gtf_to_bed(args.annotation_gtf, feature='transcript')
    with open(args.vcf_chr_list) as f:
        chr_list = f.read().strip().split('\n')
    norm_bed_df = prepare_bed(norm_df, bed_template_df, chr_subset=chr_list)
    print('  * {} genes remain after removing contigs absent from VCF.'.format(norm_bed_df.shape[0]), flush=True)
    print('Writing BED file', flush=True)
    write_bed(norm_bed_df, os.path.join(args.output_dir, args.prefix+'.expression.bed'))

    if args.legacy_mode:
        # for consistency with v6/v6p pipeline results, write unsorted expression file for PEER factor calculation
        norm_df.to_csv(os.path.join(args.output_dir, args.prefix+'.expression.txt'), sep='\t')
