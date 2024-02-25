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


def prepare_expression(counts_df, tpm_df, sample_frac_threshold=0.2,
                       count_threshold=6, tpm_threshold=0.1, mode='tmm'):
    """
    Genes are filtered using the following expression thresholds:
      TPM >= tpm_threshold in >= sample_frac_threshold * samples
      read counts >= count_threshold in sample_frac_threshold * samples

    The filtered counts matrix is then normalized using:
      TMM (mode='tmm'; default) or
      quantile normalization (mode='qn')
    """

    # expression thresholds
    ns = tpm_df.shape[1]
    mask = (
        (np.sum(tpm_df >= tpm_threshold, axis=1) >= sample_frac_threshold * ns) &
        (np.sum(counts_df >= count_threshold, axis=1) >= sample_frac_threshold * ns)
    ).values

    # apply normalization
    if mode.lower() == 'tmm':
        tmm_counts_df = qtl.norm.edger_cpm(counts_df, normalized_lib_sizes=True)
        norm_df = qtl.norm.inverse_normal_transform(tmm_counts_df[mask])
    elif mode.lower() == 'qn':
        qn_df = qtl.norm.normalize_quantiles(tpm_df.loc[mask])
        norm_df = qtl.norm.inverse_normal_transform(qn_df)
    else:
        raise ValueError(f'Unsupported mode {mode}')

    return norm_df



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate normalized expression BED files for eQTL analyses')
    parser.add_argument('tpm_gct', help='GCT or Parquet file with expression in normalized units, e.g., TPM or FPKM')
    parser.add_argument('counts_gct', help='GCT or Parquet file with read counts')
    parser.add_argument('annotation_gtf', help='GTF annotation used for generating expression matrices')
    parser.add_argument('sample_to_participant', help='TSV linking sample IDs (columns in expression matrices) to participant IDs (VCF IDs)')
    parser.add_argument('prefix', help='Prefix for output file names')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    parser.add_argument('--sample_ids', default=None, help='File listing sample IDs to include')
    parser.add_argument('--chrs', help='File listing chromosomes to include (default: chr1-22 + chrX)')
    parser.add_argument('--convert_tpm', action='store_true', help='Convert to TPM (in case input is in RPKM/FPKM)')
    parser.add_argument('--tpm_threshold', type=np.double, default=0.1, help='Selects genes with > expression_threshold expression in at least sample_frac_threshold')
    parser.add_argument('--count_threshold', type=np.int32, default=6, help='Selects genes with >= count_threshold reads in at least sample_frac_threshold samples')
    parser.add_argument('--sample_frac_threshold', type=np.double, default=0.2, help='Minimum fraction of samples that must satisfy thresholds')
    parser.add_argument('--normalization_method', default='tmm', help='Normalization method: TMM or quantile normalization (qn)')
    parser.add_argument('--parquet', action='store_true', help='Write output in Parquet format')
    args = parser.parse_args()

    print('Loading expression data', flush=True)
    sample_ids = None
    if args.sample_ids is not None:
        with open(args.sample_ids) as f:
            sample_ids = f.read().strip().split('\n')
            print(f'  * Loading {len(sample_ids)} samples', flush=True)

    counts_df = qtl.io.read_gct(args.counts_gct, sample_ids=sample_ids, load_description=False)
    tpm_df = qtl.io.read_gct(args.tpm_gct, sample_ids=sample_ids, load_description=False)

    sample_to_participant_s = pd.read_csv(args.sample_to_participant, sep='\t', index_col=0,
                                          header=None, dtype=str).squeeze('columns')

    # check inputs
    if not counts_df.columns.equals(tpm_df.columns):
        raise ValueError('Sample IDs in the TPM and read counts files must match.')
    missing_ids = ~counts_df.columns.isin(sample_to_participant_s.index)
    if missing_ids.any():
        raise ValueError(f"Sample IDs in expression files and participant lookup table must match ({missing_ids.sum()} sample IDs missing from {os.path.basename(args.sample_to_participant)}).")

    if args.convert_tpm:
        print('  * Converting to TPM', flush=True)
        tpm_df = tpm_df / tpm_df.sum(0) * 1e6

    print(f'Normalizing data ({args.normalization_method})', flush=True)
    norm_df = prepare_expression(counts_df, tpm_df,
                                 sample_frac_threshold=args.sample_frac_threshold,
                                 count_threshold=args.count_threshold,
                                 tpm_threshold=args.tpm_threshold,
                                 mode=args.normalization_method)
    print(f'  * {counts_df.shape[0]} genes in input tables', flush=True)
    print(f'  * {norm_df.shape[0]} genes remain after thresholding', flush=True)

    # change sample IDs to participant IDs
    norm_df.rename(columns=sample_to_participant_s, inplace=True)

    if args.chrs is not None:
        with open(args.chrs) as f:
            chrs = f.read().strip().split('\n')
    else:
        chrs = [f'chr{i}' for i in range(1,23)] + ['chrX']

    # prepare BED
    bed_template_df = qtl.io.gtf_to_tss_bed(args.annotation_gtf, feature='transcript')
    bed_template_df = bed_template_df[bed_template_df['chr'].isin(chrs)]
    bed_df = pd.merge(bed_template_df, norm_df, left_index=True, right_index=True)
    qtl.io.sort_bed(bed_df, inplace=True)
    print(f'  * {bed_df.shape[0]} genes remain after selecting chromosomes', flush=True)

    print('Writing BED file', flush=True)
    if not args.parquet:
        qtl.io.write_bed(bed_df, os.path.join(args.output_dir, f'{args.prefix}.expression.bed.gz'))
    else:
        bed_df.to_parquet(os.path.join(args.output_dir, f'{args.prefix}.expression.bed.parquet'))
