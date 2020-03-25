import numpy as np
# import pandas as pd
import re
import gzip
import argparse
from collections import defaultdict


def remove_indels(read_bases):
    """Remove indels from pileup output"""

    # read_bases = ',,,-19tcttagtaagattacacat,,,+2at...'
    indel = re.compile('[-+]\d+')

    b = bytearray(read_bases.encode())

    ix = []
    # ins = []
    # dels = []
    for m in indel.finditer(read_bases):
        s = m.group()
        n = int(s[1:])
        # if s[0]=='+':
        #     ins.append([m.end(), m.end()+n])
        # else:
        #     dels.append([m.end(), m.end()+n])
        ix.append([m.start(), m.end()+n])

    # ins = [read_bases[i[0]:i[1]] for i in ins]
    # dels = [read_bases[i[0]:i[1]] for i in dels]
    # ins = '|'.join(['{}:{}'.format(j,i) for i,j in pd.Series(ins).value_counts().items()])
    # dels = '|'.join(['{}:{}'.format(j,i) for i,j in pd.Series(dels).value_counts().items()])
    # if ins=='':
    #     ins = 'NA'
    # if dels=='':
    #     dels = 'NA'

    for i in ix[::-1]:
        b[i[0]:i[1]] = b''

    return b.decode()#, ins, dels


def process_pileup(mpileup_file, output_file, min_baseq=10, min_mapq=93):
    """
    Convert samtools mpileup output into base counts at each position
    
    https://en.wikipedia.org/wiki/Pileup_format
    http://samtools.sourceforge.net/pileup.shtml
    """
    if min_mapq>93:
        min_mapq = 93  # maximum Phred score

    # regex for start of a read segment followed by mapping quality
    nbchar = re.compile('\\^.|\\$')

    header = [
        'contig',
        'pos',
        'ref_base',
        'A',
        'T',
        'C',
        'G',
        'a',
        't',
        'c',
        'g',
        # 'ins',
        # 'del',
        'total_count',
        'low_mapq_depth',
        'low_baseq_depth',
        'raw_depth',
        'other_count',
    ]

    with gzip.open(mpileup_file, 'rt') as f, \
         gzip.open(output_file, 'wt') as out:

        out.write('\t'.join(header)+'\n')

        for n,line in enumerate(f,1):
            line = line.strip().split('\t')
            chrom = line[0]
            pos = int(line[1])
            ref_base = line[2]
            depth = int(line[3])

            read_bases = line[4]
            base_quals = line[5]
            map_quals = line[6]

            # remove indels, read starts(^) and ends($)
            read_bases = nbchar.sub('', read_bases)
            # read_bases, ins, dels = remove_indels(read_bases)
            read_bases = remove_indels(read_bases)

            assert len(read_bases)==len(base_quals)==len(map_quals)

            # convert to Phred quality scores
            base_quals = np.array([ord(i) for i in base_quals]) - 33
            map_quals = np.array([ord(i) for i in map_quals]) - 33
            read_bases = np.array(list(read_bases))

            # remove reference skips
            skip_set = set(['<', '>'])
            s = [i not in skip_set for i in read_bases]
            if not np.all(s):
                read_bases = read_bases[s]
                base_quals = base_quals[s]
                map_quals = map_quals[s]

            # count
            raw_depth = len(read_bases)
            low_mapq = np.sum(map_quals<min_mapq)
            low_baseq = np.sum(base_quals<min_baseq)
            # is_ref = (read_bases=='.') | (read_bases==',')

            # apply quality filters
            ix = (base_quals>=min_baseq) & (map_quals>=min_mapq)

            ncount = {}
            for b in ['A','T','C','G','a','t','c','g']:
                if b==ref_base:
                    ncount[b] = np.sum(read_bases[ix]=='.')
                elif b==ref_base.lower():
                    ncount[b] = np.sum(read_bases[ix]==',')
                else:
                    ncount[b] = np.sum(read_bases[ix]==b)
            total_count = np.sum(list(ncount.values()))
            other_count = np.sum(ix) - total_count

            v = [
                chrom,
                str(pos),
                ref_base,
                str(ncount['A']),
                str(ncount['T']),
                str(ncount['C']),
                str(ncount['G']),
                str(ncount['a']),
                str(ncount['t']),
                str(ncount['c']),
                str(ncount['g']),
                str(total_count),
                str(low_mapq),
                str(low_baseq),
                str(raw_depth),
                str(other_count),
            ]

            out.write('\t'.join(v)+'\n')
            if np.mod(n, 100)==0:
                print('\rProcessed {} sites'.format(n), end='', flush=True)
        print('\rProcessed {} sites'.format(n))


def process_pileup_10x(mpileup_file, variant_dict, output_file, min_baseq=10, min_mapq=93, min_depth=1):
    """
    Convert samtools mpileup output into base counts at each position
    
    https://en.wikipedia.org/wiki/Pileup_format
    http://samtools.sourceforge.net/pileup.shtml
    """
    if min_mapq>93:
        min_mapq = 93  # maximum Phred score

    # regex for start of a read segment followed by mapping quality
    nbchar = re.compile('\\^.|\\$')

    header = [
        'contig',
        'pos',
        'variant_id',
        'ref',
        'alt',
        'ref_count',
        'alt_count',
        'total_count',
        'low_mapq_depth',
        'low_baseq_depth',
        'raw_depth',
        'other_count',
        'umi_counts',
    ]

    with gzip.open(mpileup_file, 'rt') as f, \
         gzip.open(output_file, 'wt') as out:

        out.write('\t'.join(header)+'\n')

        for n,line in enumerate(f,1):
            line = line.strip().split('\t')

            chrom = line[0]
            pos = int(line[1])
            ref = variant_dict[line[0]][pos]['ref']
            alt = variant_dict[line[0]][pos]['alt']

            # ref_base = line[2]
            depth = int(line[3])

            read_bases = line[4]
            base_quals = line[5]
            map_quals = line[6]
            barcodes = np.array(line[7].split(','))

            # remove indels, read starts(^) and ends($)
            read_bases = nbchar.sub('', read_bases)
            read_bases = remove_indels(read_bases)

            assert len(read_bases)==len(base_quals)==len(map_quals)==len(barcodes)

            # convert to Phred quality scores
            base_quals = np.array([ord(i) for i in base_quals]) - 33
            map_quals = np.array([ord(i) for i in map_quals]) - 33
            read_bases = np.array(list(read_bases))

            # remove reference skips
            skip_set = set(['<', '>'])
            s = [i not in skip_set for i in read_bases]
            if not np.all(s):
                read_bases = read_bases[s]
                base_quals = base_quals[s]
                map_quals = map_quals[s]
                barcodes = barcodes[s]

            # count
            raw_depth = len(read_bases)
            low_mapq = np.sum(map_quals<min_mapq)
            low_baseq = np.sum(base_quals<min_baseq)
            is_ref = (read_bases=='.') | (read_bases==',')
            is_alt = (read_bases==alt.upper()) | (read_bases==alt.lower())
            is_other = ~is_ref & ~is_alt

            # apply quality filters
            ix = (base_quals>=min_baseq) & (map_quals>=min_mapq)
            is_ref = is_ref[ix]
            is_alt = is_alt[ix]
            is_other = is_other[ix]
            ref_count = np.sum(is_ref)
            alt_count = np.sum(is_alt)
            other_count = np.sum(is_other)




            # cell-specific counts
            read_bases = read_bases[ix]
            base_quals = base_quals[ix]
            # map_quals = map_quals[ix]
            barcodes = barcodes[ix]

            # drop duplicates
            barcode_dict = {}
            for i,j,k in zip(barcodes, base_quals, range(len(barcodes))):
                if i not in barcode_dict:
                    barcode_dict[i] = (j,k)
                elif j>barcode_dict[i][0]:
                    barcode_dict[i] = (j,k)
            ix2 = np.sort([barcode_dict[i][1] for i in barcode_dict])

            if len(ix2)>0:
                # count by barcode (ref/alt only)
                barcode_dict = defaultdict(lambda: defaultdict(int))
                for b,i,j in zip(barcodes[ix2],is_ref[ix2],is_alt[ix2]):
                    if i:
                        barcode_dict[b.split(':')[0].split('-')[0]]['ref'] += 1
                    elif j:
                        barcode_dict[b.split(':')[0].split('-')[0]]['alt'] += 1

                # cell_barcode:ref:alt
                cell_str = ','.join(['{}:{}:{}'.format(k, barcode_dict[k]['ref'], barcode_dict[k]['alt']) for k in barcode_dict])
            else:
                cell_str = ''

            out.write('\t'.join([
                line[0],
                line[1],
                variant_dict[line[0]][pos]['id'],
                ref,
                alt,
                str(ref_count),
                str(alt_count),
                str(ref_count + alt_count),
                str(low_mapq),
                str(low_baseq),
                str(raw_depth),
                str(other_count),
                cell_str,
            ])+'\n')
            if np.mod(n,1000)==0:
                print('\rProcessed {} sites'.format(n), end='')
        print()



if __name__=='__main__':

    parser = argparse.ArgumentParser(description="Parser for samtools mpileup output")
    parser.add_argument("mpileup_file")
    parser.add_argument("--min_baseq", type=int, default=10)
    parser.add_argument("--min_mapq", type=int, default=93)
    parser.add_argument("-o, ", "--out", type=str, default=None)
    args = parser.parse_args()

    if args.out is None:
        output_file = args.mpileup_file.replace('.mpileup.gz', '.mpileup.counts.txt.gz')
    else:
        output_file = args.out

    process_pileup(args.mpileup_file, output_file, min_baseq=args.min_baseq, min_mapq=args.min_mapq)
