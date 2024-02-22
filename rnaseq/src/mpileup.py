#!/usr/bin/env python3
import numpy as np
import pandas as pd
import re
import gzip
import argparse
from collections import defaultdict
import subprocess
import os
import tempfile


def samtools_mpileup(bam_file, fasta, region=None, regions_bed=None, variant_ids=None,
                     output_extra=None, prefix=None, output_dir='.', out_file=None,
                     primary_only=True, paired_end=True, min_mq=0, min_bq=0, max_depth=100000):
    """
    Wrapper for samtools mpileup

    Inputs:
      bam_file: input BAM
      fasta: reference sequence FASTA
      region: chr:start-end string
      regions_bed: BED file defining positions to query
      variant_ids: list of variant positions to query (e.g., ['chr1_123', 'chr2_456'])
      output_extra: list tags to include in output (e.g., output_extra='CB,UB')
      prefix: prefix for output file: ${prefix}.mpileup.gz
    """
    if out_file is None:
        if prefix is None:
            out_file = bam_file.replace('.bam', '.mpileup.gz')
        else:
            out_file = os.path.join(output_dir, f"{prefix}.mpileup.gz")

    cmd = f"samtools mpileup -d {max_depth} --no-BAQ -q {min_mq} -Q {min_bq} --output-MQ"
    if primary_only:
        cmd += ' --excl-flags 3328'  # not primary, PCR or optical duplicate, supplementary
    if paired_end:
        cmd += ' --incl-flags 2'
    if output_extra is not None:
        cmd += f" --output-extra {output_extra}"

    if region is not None:
        cmd += f" --region {region}"
    elif regions_bed is not None or variant_ids is not None:
        if variant_ids is not None:
            if isinstance(variant_ids, str):
                variant_ids = [variant_ids]
            f = tempfile.NamedTemporaryFile('wt', dir=os.path.dirname(bam_file))
            regions_bed = f.name
            for variant_id in variant_ids:
                chrom, pos = variant_id.split('_')[:2]
                pos = int(pos)
                f.write(f"{chrom}\t{pos-1}\t{pos}\n")
            f.flush()
        cmd += f" -l {regions_bed}"

    cmd += f" -f {fasta} {bam_file} | gzip -1 > {out_file}"
    subprocess.check_call(cmd, shell=True, stderr=subprocess.DEVNULL)
    if variant_ids is not None:
        f.close()


def remove_indels(read_bases):
    """Remove indels from mpileup output"""

    # read_bases = ',,,-19tcttagtaagattacacat,,,+2at...'
    indel = re.compile('[-+]\d+')

    b = bytearray(read_bases.encode())

    ix = []
    # ins = []
    # dels = []
    for m in indel.finditer(read_bases):
        s = m.group()
        n = int(s[1:])
        # if s[0] == '+':
        #     ins.append([m.end(), m.end()+n])
        # else:
        #     dels.append([m.end(), m.end()+n])
        ix.append([m.start(), m.end()+n])

    # ins = [read_bases[i[0]:i[1]] for i in ins]
    # dels = [read_bases[i[0]:i[1]] for i in dels]
    # ins = '|'.join(['{}:{}'.format(j,i) for i,j in pd.Series(ins).value_counts().items()])
    # dels = '|'.join(['{}:{}'.format(j,i) for i,j in pd.Series(dels).value_counts().items()])
    # if ins == '':
    #     ins = 'NA'
    # if dels == '':
    #     dels = 'NA'

    for i in ix[::-1]:
        b[i[0]:i[1]] = b''

    return b.decode()#, ins, dels


def process_mpileup(mpileup_file, min_baseq=10, min_mapq=93, verbose=False):
    """
    Convert samtools mpileup output into base counts at each position

    https://en.wikipedia.org/wiki/Pileup_format
    http://samtools.sourceforge.net/pileup.shtml
    """
    if min_mapq > 93:
        min_mapq = 93  # maximum Phred score

    # regex for start of a read segment followed by mapping quality
    nbchar = re.compile('\\^.|\\$')

    # check if there are barcodes
    with gzip.open(mpileup_file, 'rt') as f:
        line = f.readline().strip().split('\t')
    if len(line) > 7:
        has_barcodes = True
        out = defaultdict(list)
    else:
        has_barcodes = False
        out = []

    with gzip.open(mpileup_file, 'rt') as f:
        for n,line in enumerate(f, 1):
            if verbose and n % 100 == 0:
                print(f"\r  * parsing mpileup line {n}", end='')
            line = line.strip().split('\t')
            chrom = line[0]
            pos = int(line[1])
            ref_base = line[2]
            depth = int(line[3])

            read_bases = line[4]
            base_quals = line[5]
            map_quals = line[6]
            if has_barcodes:
                barcodes = np.array(line[7].split(','))  # cell barcodes (CB tag)

            # remove indels, read starts(^) and ends($)
            read_bases = nbchar.sub('', read_bases)
            # read_bases, ins, dels = remove_indels(read_bases)
            read_bases = remove_indels(read_bases)

            assert len(read_bases) == len(base_quals) == len(map_quals) # == len(barcodes)

            # convert to Phred quality scores
            base_quals = np.array([ord(i) for i in base_quals]) - 33
            map_quals = np.array([ord(i) for i in map_quals]) - 33
            read_bases = np.array(list(read_bases))

            # filter out reference skips
            skip_set = set(['<', '>'])
            s = [i not in skip_set for i in read_bases]
            if not np.all(s):
                read_bases = read_bases[s]
                base_quals = base_quals[s]
                map_quals = map_quals[s]
                if has_barcodes:
                    barcodes = barcodes[s]

            # apply quality filters
            ix = (base_quals >= min_baseq) & (map_quals >= min_mapq)

            if not has_barcodes:
                # count
                raw_depth = len(read_bases)
                low_mapq = np.sum(map_quals < min_mapq)
                low_baseq = np.sum(base_quals < min_baseq)

                ncount = {}
                for b in ['A','T','C','G','a','t','c','g']:
                    if b == ref_base:
                        ncount[b] = np.sum(read_bases[ix] == '.')
                    elif b == ref_base.lower():
                        ncount[b] = np.sum(read_bases[ix] == ',')
                    else:
                        ncount[b] = np.sum(read_bases[ix] == b)
                total_count = np.sum(list(ncount.values()))
                other_count = np.sum(ix) - total_count

                res_s = {
                    'contig': chrom,
                    'pos': pos,
                    'ref_base': ref_base,
                    'A': ncount['A'],
                    'T': ncount['T'],
                    'C': ncount['C'],
                    'G': ncount['G'],
                    'a': ncount['a'],
                    't': ncount['t'],
                    'c': ncount['c'],
                    'g': ncount['g'],
                    # 'ins',
                    # 'del',
                    'total_count': total_count,
                    'low_mapq_depth': low_mapq,
                    'low_baseq_depth': low_baseq,
                    'raw_depth': raw_depth,
                    'other_count': other_count,
                }
                out.append(res_s)

            else:  # iterate over barcodes
                res = defaultdict(list)
                # return barcodes
                for barcode in np.unique(barcodes):
                    m = barcodes == barcode
                    ix2 = ix & m
                    # count
                    raw_depth = len(read_bases[m])
                    low_mapq = np.sum(map_quals[m] < min_mapq)
                    low_baseq = np.sum(base_quals[m] < min_baseq)

                    ncount = {}
                    for b in ['A','T','C','G','a','t','c','g']:
                        if b == ref_base:
                            ncount[b] = np.sum(read_bases[ix2] == '.')
                        elif b == ref_base.lower():
                            ncount[b] = np.sum(read_bases[ix2] == ',')
                        else:
                            ncount[b] = np.sum(read_bases[ix2] == b)
                    total_count = np.sum(list(ncount.values()))
                    other_count = np.sum(ix2) - total_count

                    res_s = {
                        'contig': chrom,
                        'pos': pos,
                        'ref_base': ref_base,
                        'A': ncount['A'],
                        'T': ncount['T'],
                        'C': ncount['C'],
                        'G': ncount['G'],
                        'a': ncount['a'],
                        't': ncount['t'],
                        'c': ncount['c'],
                        'g': ncount['g'],
                        # 'ins',
                        # 'del',
                        'total_count': total_count,
                        'low_mapq_depth': low_mapq,
                        'low_baseq_depth': low_baseq,
                        'raw_depth': raw_depth,
                        'other_count': other_count,
                    }
                    out[barcode].append(res_s)

    if verbose:
        print(f"\r  * parsing mpileup line {n}")
    if not has_barcodes:
        out = pd.DataFrame(out)
        out.index = out['contig'] + '_' + out['pos'].astype(str) + '_' + out['ref_base']
    else:
        for k in out:
            out[k] = pd.DataFrame(out[k])
            out[k].index = out[k]['contig'] + '_' + out[k]['pos'].astype(str) + '_' + out[k]['ref_base']
    return out


def collapse_strands(df):
    """Merge counts from both strands for each base"""
    if isinstance(df, pd.DataFrame):
        for n in ['a', 't', 'c', 'g']:
            df[n.upper()] += df[n]
        df.drop(['a', 't', 'c', 'g'], axis=1, inplace=True)
    elif isinstance(df, dict):
        for k in df:
            collapse_strands(df[k])
    else:
        raise ValueError(f'Input type {df.__class__} not supported')


def get_ase_counts(mpileup_df, sites_vcf, inplace=True):
    # parse variants from sites VCF
    variant_dict = {}
    with gzip.open(sites_vcf, 'rt') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                # store ID, ref, alt
                variant_dict[f"{line[0]}_{line[1]}_{line[3]}"] = line[2:5]

    # count ref/alt bases only; recompute total count
    count_dicts = {}
    for b in 'ACGT':
        count_dicts[b] = mpileup_df[b].to_dict()

    ref_count = {}
    alt_count = {}
    alt_base = {}
    for i in mpileup_df.index:
        d = variant_dict[i]
        ref_count[i] = count_dicts[d[1]][i]
        alt_count[i] = count_dicts[d[2]][i]
        alt_base[i] = d[2]
    ref_count = pd.Series(ref_count, name='ref_count')
    alt_count = pd.Series(alt_count, name='alt_count')
    alt_base = pd.Series(alt_base, name='alt_base')

    if not inplace:
        mpileup_df = mpileup_df.copy()
    mpileup_df.insert(3, 'alt_base', alt_base)
    mpileup_df.insert(4, 'ref_count', ref_count)
    mpileup_df.insert(5, 'alt_count', alt_count)
    mpileup_df['total_count_multi'] = mpileup_df['total_count']
    mpileup_df['total_count'] = mpileup_df['ref_count'] + mpileup_df['alt_count']
    mpileup_df.drop(['A', 'T', 'C', 'G'], axis=1, inplace=True)
    if not inplace:
        return mpileup_df


if __name__=='__main__':

    parser = argparse.ArgumentParser(description="Wrapper for samtools mpileup")
    parser.add_argument("bam_file")
    parser.add_argument("reference_fasta")
    parser.add_argument("prefix", default=None)
    parser.add_argument("--min_baseq", type=int, default=10)
    parser.add_argument("--min_mapq", type=int, default=93)
    parser.add_argument("--region", default=None)
    parser.add_argument("--regions_bed", default=None)
    parser.add_argument("--variant_ids", default=None)
    parser.add_argument("--output_extra", default=None)
    parser.add_argument("--max_depth", type=int, default=100000)
    parser.add_argument("-o, ", "--output_dir", type=str, default='.')
    parser.add_argument('--collapse_strands', action='store_true', help='Collapse counts from both strands')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    args = parser.parse_args()

    samtools_mpileup(args.bam_file, args.reference_fasta, region=args.region, regions_bed=args.regions_bed,
                     variant_ids=args.variant_ids, output_extra=args.output_extra,
                     prefix=args.prefix, output_dir=args.output_dir,
                     primary_only=True, paired_end=True, min_mq=0, min_bq=0, max_depth=args.max_depth)

    mpileup_file = os.path.join(args.output_dir, f"{args.prefix}.mpileup.gz")
    df = process_mpileup(mpileup_file, min_baseq=args.min_baseq, min_mapq=args.min_mapq, verbose=args.verbose)
    if args.collapse_strands:
        collapse_strands(df)
    df.to_parquet(os.path.join(args.output_dir, f"{args.prefix}.mpileup.parquet"))
