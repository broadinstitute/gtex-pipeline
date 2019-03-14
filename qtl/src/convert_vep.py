#!/usr/bin/env python3
import pandas as pd
import numpy as np
import gzip
import pickle
import sys
import argparse
import os


def choose_variant_annotation(csq_string, variant_annotation_rank_dict, gene_ix, conseq_ix):
    """Parse lowest ranking consequence and gene ID from CSQ string"""
    minrank = len(variant_annotation_rank_dict)
    gene_id = 'NA'
    conseq = 'NA'

    for tr in csq_string.split(','):
        annots = tr.split('|')
        for v in annots[conseq_ix].split('&'):
            if v in variant_annotation_rank_dict:
                r = variant_annotation_rank_dict[v]
                if r<minrank:
                    minrank = r
                    gene_id = annots[gene_ix]
                    conseq = v
    return gene_id, conseq


def get_vep_format(vep_vcf):
    """Format of the CSQ string"""
    fmt = None
    with gzip.open(vep_vcf) as f:
        for line in f:
            line = line.decode().strip()
            if line[0]!='#':
                break
            elif line[:6]=='##INFO' and 'CSQ' in line:
                fmt = line.split('Format: ')[1].replace('">','').split('|')
                break
    if fmt is None:
        raise ValueError('CSQ format not found in VCF header')
    return fmt


# from http://useast.ensembl.org/info/genome/variation/predicted_data.html
variant_annotation_rank_list = [
    'transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost',
    'transcript_amplification',
    'inframe_insertion',
    'inframe_deletion',
    'missense_variant',
    'protein_altering_variant',
    'splice_region_variant',
    'incomplete_terminal_codon_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant',
    'mature_miRNA_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'non_coding_transcript_exon_variant',
    'intron_variant']

variant_annotation_rank_dict = {j:i for i,j in enumerate(variant_annotation_rank_list)}


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Extract single effect from VEP VCF')
    parser.add_argument('vcf', help='VCF')
    parser.add_argument('vep_vcf', help='Variant Effect Predictor VCF')
    parser.add_argument('-o', '--output_dir', default='.')
    args = parser.parse_args()

    #------------------------------------------------
    # Pre-processing: extract variant set from VCF
    #------------------------------------------------
    prefix = os.path.basename(args.vcf).split('.vcf')[0]
    variant_id_file = os.path.join(args.output_dir, prefix+'.variant_ids.txt.gz')
    if not os.path.exists(variant_id_file):
        print('Extracting variant IDs to {}'.format(variant_id_file))
        subprocess.check_call('zcat '+args.vcf+' | grep -v "#" | cut -f3 | gzip -c -1 > '+variant_id_file)

    print('Loading variant IDs')
    with gzip.open(variant_id_file) as f:
        v = f.read()
    variant_set = set(v.decode().strip().split('\n'))

    # get position of Gene and Consequence fields in CSQ string (changes with VEP versions)
    fmt = get_vep_format(args.vep_vcf)
    gene_ix = np.where(np.array(fmt)=='Gene')[0][0]
    conseq_ix = np.where(np.array(fmt)=='Consequence')[0][0]

    # parse VEP VCF and write variant_id, gene, consequence to file
    # use pandas parser for speed/process by chunks to reduce memory usage
    print('Parsing VEP VCF')
    vep_file = os.path.join(args.output_dir, prefix+'.vep.txt.gz')
    with gzip.open(vep_file, 'wt') as f:
        f.write('variant_id\tgene_id\tvep\n')
        for k,chunk in enumerate(pd.read_csv(args.vep_vcf, sep='\t', comment='#', header=None,
            names=['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info'], dtype=bytes, usecols=['id', 'info'],
            iterator=True, chunksize=100000)):

            print('\rProcessing chunk {0:d}'.format(k+1), end='')
            for i,c in zip(chunk['id'], chunk['info']):
                if len(c.split('|')) > 10:
                    cs = [i for i in c.split(';') if i.startswith('CSQ')][0]
                    gene_id, conseq = choose_variant_annotation(cs, variant_annotation_rank_dict, gene_ix, conseq_ix)
                    if gene_id != 'NA':
                        for j in i.split(';'):
                            if j in variant_set:
                                f.write(j+'\t'+gene_id+'\t'+conseq+'\n')
            print()

    print('Converting to dictionary')
    vep_df = pd.read_csv(vep_file, sep='\t', index_col=0)
    vep_df['combined'] = [(i,j) for i,j in zip(vep_df['gene_id'], vep_df['vep'])]
    vep_dict = vep_df['combined'].to_dict()

    print('Saving as pickle')
    with open(os.path.join(args.output_dir, prefix+'.vep_dict.pickle'), 'wb') as f:
        pickle.dump(vep_dict, f, pickle.HIGHEST_PROTOCOL)
