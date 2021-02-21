#!/usr/bin/env python3
# Author: Francois Aguet

import numpy as np
import pandas as pd
import gzip
import subprocess
import scipy.stats as stats
import argparse
import os
import qtl.io
import qtl.norm


def prepare_bed(df, bed_template_df, chr_subset=None):
    bed_df = pd.merge(bed_template_df, df, left_index=True, right_index=True)
    # sort by start position
    bed_df = bed_df.groupby('chr', sort=False, group_keys=False).apply(lambda x: x.sort_values('start'))
    if chr_subset is not None:
        # subset chrs from VCF
        bed_df = bed_df[bed_df.chr.isin(chr_subset)]
    return bed_df


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
        tmm_counts_df = qtl.norm.edger_cpm(counts_df, normalized_lib_sizes=True)
        norm_df = qtl.norm.inverse_normal_transform(tmm_counts_df[mask])
    elif mode.lower()=='qn':
        qn_df = qtl.norm.normalize_quantiles(tpm_df.loc[mask])
        norm_df = qtl.norm.inverse_normal_transform(qn_df)
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

    counts_df = qtl.io.read_gct(args.counts_gct, sample_ids=sample_ids, load_description=False)
    tpm_df = qtl.io.read_gct(args.tpm_gct, sample_ids=sample_ids, load_description=False)

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

    bed_template_df = qtl.io.gtf_to_tss_bed(args.annotation_gtf, feature='transcript')
    with open(args.vcf_chr_list) as f:
        chr_list = f.read().strip().split('\n')
    norm_bed_df = prepare_bed(norm_df, bed_template_df, chr_subset=chr_list)
    print('  * {} genes remain after removing contigs absent from VCF.'.format(norm_bed_df.shape[0]), flush=True)
    print('Writing BED file', flush=True)
    qtl.io.write_bed(norm_bed_df, os.path.join(args.output_dir, args.prefix+'.expression.bed.gz'))

    if args.legacy_mode:
        # for consistency with v6/v6p pipeline results, write unsorted expression file for PEER factor calculation
        norm_df.to_csv(os.path.join(args.output_dir, args.prefix+'.expression.txt'), sep='\t')
