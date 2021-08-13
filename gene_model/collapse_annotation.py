#!/usr/bin/env python3
# Author: Francois Aguet
import numpy as np
import pandas as pd
from collections import defaultdict
from bx.intervals.intersection import IntervalTree
import argparse
import os
import gzip


class Exon:
    def __init__(self, exon_id, number, transcript, start_pos, end_pos):
        self.id = exon_id
        self.number = int(number)
        self.transcript = transcript
        self.start_pos = start_pos
        self.end_pos = end_pos

class Transcript:
    def __init__(self, transcript_id, transcript_name, transcript_type, gene, start_pos, end_pos):
        self.id = transcript_id
        self.name = transcript_name
        self.type = transcript_type
        self.gene = gene
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.exons = []

class Gene:
    def __init__(self, gene_id, gene_name, gene_type, chrom, strand, start_pos, end_pos):
        self.id = gene_id
        self.name = gene_name
        self.biotype = gene_type
        self.chr = chrom
        self.strand = strand
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.transcripts = []

class Annotation:
    def __init__(self, gtfpath):
        """Parse GTF and construct gene/transcript/exon hierarchy"""

        if gtfpath.endswith('.gtf.gz'):
            opener = gzip.open(gtfpath, 'rt')
        else:
            opener = open(gtfpath, 'r')

        self.genes = []
        with opener as gtf:
            for row in gtf:
                row = row.strip().split('\t')

                if row[0][0] == '#': continue # skip header

                chrom = row[0]
                annot_type = row[2]
                start_pos = int(row[3])
                end_pos  = int(row[4])
                strand = row[6]

                attributes = defaultdict(list)
                for a in row[8].replace('"', '').replace('_biotype', '_type').split(';')[:-1]:
                    kv = a.strip().split(' ')
                    if kv[0]!='tag':
                        attributes[kv[0]] = kv[1]
                    else:
                        attributes['tags'].append(kv[1])

                if annot_type == 'gene':
                    assert 'gene_id' in attributes
                    if 'gene_name' not in attributes:
                        attributes['gene_name'] = attributes['gene_id']
                    gene_id = attributes['gene_id']
                    g = Gene(gene_id, attributes['gene_name'], attributes['gene_type'],
                             chrom, strand, start_pos, end_pos)
                    g.source = row[1]
                    g.phase = row[7]
                    g.attributes_string = row[8].replace('_biotype', '_type')
                    self.genes.append(g)

                elif annot_type == 'transcript':
                    assert 'transcript_id' in attributes
                    if 'transcript_name' not in attributes:
                        attributes['transcript_name'] = attributes['transcript_id']
                    transcript_id = attributes['transcript_id']
                    t = Transcript(attributes.pop('transcript_id'), attributes.pop('transcript_name'),
                                   attributes.pop('transcript_type'), g, start_pos, end_pos)
                    t.attributes = attributes
                    g.transcripts.append(t)

                elif annot_type == 'exon':
                    if 'exon_id' in attributes:
                        e = Exon(attributes['exon_id'], attributes['exon_number'], t, start_pos, end_pos)
                    else:
                        e = Exon(str(len(t.exons)+1), len(t.exons)+1, t, start_pos, end_pos)
                    t.exons.append(e)

                if len(self.genes) % 1000 == 0:
                    print(f'\rParsing GTF: {len(self.genes)} genes processed', end='')
            print(f'\rParsing GTF: {len(self.genes)} genes processed')

        self.genes = np.array(self.genes)


def interval_union(intervals):
    """
    Returns the union of all intervals in the input list
      intervals: list of tuples or 2-element lists
    """
    intervals.sort(key=lambda x: x[0])
    union = [intervals[0]]
    for i in intervals[1:]:
        if i[0] <= union[-1][1]:  # overlap w/ previous
            if i[1] > union[-1][1]:  # only extend if larger
                union[-1][1] = i[1]
        else:
            union.append(i)
    return union


def subtract_segment(a, b):
    """
    Subtract segment a from segment b,
      return 'a' if no overlap
    """
    if a[0]>=b[0] and a[0]<=b[1] and a[1]>b[1]:
        return (b[1]+1,a[1])
    elif a[0]<b[0] and a[1]>=b[0] and a[1]<=b[1]:
        return (a[0], b[0]-1)
    elif a[0]<b[0] and a[1]>b[1]:
        return [(a[0],b[0]-1), (b[1]+1,a[1])]
    elif a[0]>=b[0] and a[1]<=b[1]:
        return []
    else:
        return a


def add_transcript_attributes(attributes_string):
    """
    Adds transcript attributes if they were missing
    (see https://www.gencodegenes.org/pages/data_format.html)

    'status' fields were dropped in Gencode 26 and later
    """
    # GTF specification
    if 'gene_status' in attributes_string:
        attribute_order = ['gene_id', 'transcript_id', 'gene_type', 'gene_status', 'gene_name',
                           'transcript_type', 'transcript_status', 'transcript_name']
        add_list = ['transcript_id', 'transcript_type', 'transcript_status', 'transcript_name']
    else:
        attribute_order = ['gene_id', 'transcript_id', 'gene_type',
                           'gene_name', 'transcript_type', 'transcript_name']
        add_list = ['transcript_id', 'transcript_type', 'transcript_name']
    if 'level' in attributes_string:
        attribute_order += ['level']

    attr = attributes_string.strip(';').split('; ')
    req = []
    opt = []
    for k in attr:
        if k.split()[0] in attribute_order:
            req.append(k)
        else:
            opt.append(k)
    attr_dict = {i.split()[0]:i.split()[1].replace(';','') for i in req}
    if 'gene_name' not in attr_dict:
        attr_dict['gene_name'] = attr_dict['gene_id']
    if 'transcript_id' not in attr_dict:
        attr_dict['transcript_id'] = attr_dict['gene_id']
    for k in add_list:
        if k not in attr_dict:
            attr_dict[k] = attr_dict[k.replace('transcript', 'gene')]

    return '; '.join([k+' '+attr_dict[k] for k in attribute_order] + opt)+';'


def collapse_annotation(annot, transcript_gtf, collapsed_gtf, blacklist=set(),
                        collapse_only=False, stranded=False):
    """
    Collapse transcripts into a single gene model; remove overlapping intervals

    Options:
      collapse_only: only collapses transcripts of each gene, does not remove overlaps
      stranded: only considers genes on the same strand when removing overlaps
    """

    exclude = set(['retained_intron', 'readthrough_transcript'])

    # 1) collapse each gene, excluding blacklisted transcript types
    merged_coord_dict = {}
    for g in annot.genes:
        exon_coords = []
        for t in g.transcripts:
            if ((t.id not in blacklist) and
                (t.type != 'retained_intron') and
                (('tags' not in t.attributes) or len(set(t.attributes['tags']).intersection(exclude)) == 0)):
                for e in t.exons:
                    exon_coords.append([e.start_pos, e.end_pos])
        if exon_coords:
            merged_coord_dict[g.id] = interval_union(exon_coords)

    if not collapse_only:
        # 2) build interval tree with merged domains
        interval_trees = defaultdict(IntervalTree)
        for g in annot.genes:
            if g.id in merged_coord_dict:
                for i in merged_coord_dict[g.id]:
                    # half-open intervals [a,b)
                    if stranded:
                        interval_trees[g.chr, g.strand].add(i[0], i[1]+1, [i, g.id])
                    else:
                        interval_trees[g.chr].add(i[0], i[1]+1, [i, g.id])

        # 3) query intervals of each gene, remove overlaps
        new_coord_dict = {}
        for g in annot.genes:
            if g.id in merged_coord_dict:
                new_intervals = []
                for i in merged_coord_dict[g.id]:  # loop merged exons
                    if stranded:
                        ints = interval_trees[g.chr, g.strand].find(i[0], i[1]+1)
                    else:
                        ints = interval_trees[g.chr].find(i[0], i[1]+1)
                    # remove self
                    ints = [r[0] for r in ints if r[1] != g.id]
                    m = set([tuple(i)])
                    for v in ints:
                        m = [subtract_segment(mx, v) for mx in m]
                        # flatten
                        m0 = []
                        for k in m:
                            if isinstance(k, tuple):
                                m0.append(k)
                            else:
                                m0.extend(k)
                        m = m0
                    new_intervals.extend(m)
                if new_intervals:
                    new_coord_dict[g.id] = new_intervals

        # 4) remove genes containing single-base exons only
        for g in annot.genes:
            if g.id in new_coord_dict:
                exon_lengths = np.array([i[1]-i[0]+1 for i in new_coord_dict[g.id]])
                if np.all(exon_lengths == 1):
                    new_coord_dict.pop(g.id)
    else:
        new_coord_dict = merged_coord_dict

    # 5) write to GTF
    if transcript_gtf.endswith('.gtf.gz'):
        opener = gzip.open(transcript_gtf, 'rt')
    else:
        opener = open(transcript_gtf, 'r')

    with open(collapsed_gtf, 'w') as output_gtf, opener as input_gtf:
        # copy header
        for line in input_gtf:
            if line[0] == '#':
                output_gtf.write(line)
            else:
                break
        output_gtf.write('##collapsed version generated by GTEx pipeline\n')
        for g in annot.genes:
            if g.id in new_coord_dict:
                start_pos = str(np.min([i[0] for i in new_coord_dict[g.id]]))
                end_pos = str(np.max([i[1] for i in new_coord_dict[g.id]]))
                if 'transcript_id' in g.attributes_string:
                    attr = g.attributes_string
                else:
                    attr = add_transcript_attributes(g.attributes_string)
                output_gtf.write('\t'.join([g.chr, g.source, 'gene', start_pos, end_pos, '.', g.strand, g.phase, attr])+'\n')
                output_gtf.write('\t'.join([g.chr, g.source, 'transcript', start_pos, end_pos, '.', g.strand, g.phase, attr])+'\n')
                if g.strand == '-':
                    new_coord_dict[g.id] = new_coord_dict[g.id][::-1]
                for k,i in enumerate(new_coord_dict[g.id], 1):
                    output_gtf.write('\t'.join([
                        g.chr, g.source, 'exon', str(i[0]), str(i[1]), '.', g.strand, g.phase,
                        attr+f' exon_id "{g.id}_{k}; exon_number {k}";'])+'\n')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Collapse isoforms into single transcript per gene and remove overlapping intervals between genes')
    parser.add_argument('transcript_gtf', help='Transcript annotation in GTF format')
    parser.add_argument('output_gtf', help='Name of the output file')
    parser.add_argument('--transcript_blacklist', help='List of transcripts to exclude (e.g., unannotated readthroughs)')
    parser.add_argument('--collapse_only', action='store_true', help='Only collapse transcripts of each gene, do not remove overlaps.')
    parser.add_argument('--stranded', action='store_true', help='Only consider genes on the same strand when removing overlaps.')
    args = parser.parse_args()

    annotation = Annotation(args.transcript_gtf)

    if args.transcript_blacklist:
        blacklist_df = pd.read_csv(args.transcript_blacklist, sep='\t')
        blacklist = set(blacklist_df[blacklist_df.columns[0]].values)
    else:
        blacklist = set()

    print('Collapsing transcripts')
    collapse_annotation(annotation, args.transcript_gtf, args.output_gtf,
                        blacklist=blacklist, collapse_only=args.collapse_only, stranded=args.stranded)
