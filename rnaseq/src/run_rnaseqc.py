#!/usr/bin/env python3
# Author: Francois Aguet
import argparse
import pandas as pd
import subprocess
from datetime import datetime
import os, shutil
import gzip
from collections import defaultdict
import contextlib

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)


def write_gct(df, gct_path):
    """
    Write DataFrame to GCT format
    """
    assert df.index.name=='Name' and df.columns[0]=='Description'
    with gzip.open(gct_path, 'wt') as f:
        f.write('#1.2\n')
        f.write('{0:d}\t{1:d}\n'.format(df.shape[0], df.shape[1]-1))
        df.to_csv(f, sep='\t')


def convert_counts(rpkm_gct, exon_intron_report, exon_report, genes_gtf, sample_id, output_dir='.'):
    """
    """
    print('Writing gene read counts to GCT', flush=True)
    rpkm_df = pd.read_csv(rpkm_gct, sep='\t', skiprows=3, header=None, usecols=[0,1], index_col=0, names=['Name','Description'])
    reads_df = pd.read_csv(exon_intron_report, sep='\t', skiprows=1, header=None, usecols=[0,2], index_col=0, names=['Name', sample_id])
    rpkm_df[sample_id] = 0
    rpkm_df.loc[reads_df.index, sample_id] = reads_df[sample_id]
    write_gct(rpkm_df, os.path.join(output_dir, sample_id+'.gene_reads.gct.gz'))

    exon_df = pd.read_csv(exon_report, sep='\t', usecols=['Exon', 'Transcript', 'Gene_Name', 'Exon_Reads'], index_col=0)

    # the exon report only contains counts for genes with at least one read -> get list of exons from GTF
    exon_count = defaultdict(int)
    gene_name = {}
    print('Parsing GTF', flush=True)
    with open(genes_gtf) as gtf:
        for line in gtf:
            if line[0]=='#': continue
            line = line.strip().split('\t')
            attr = line[8]
            if attr[-1]==';':
                attr = attr[:-1]
            attr = dict([i.split(' ') for i in attr.replace('"','').split('; ')])
            if line[2]=='exon':
                exon_count[attr['gene_id']] += 1
            elif line[2]=='gene':
                gene_name[attr['gene_id']] = attr['gene_name']

    print('Writing exon read counts to GCT', flush=True)
    output_df = pd.DataFrame(index=[i+'_'+str(j) for i in rpkm_df.index for j in range(exon_count[i])], columns=['Description', sample_id])
    output_df.index.name = 'Name'
    output_df['Description'] = rpkm_df.loc[[i.split('_')[0] for i in output_df.index], 'Description'].values
    output_df[sample_id] = 0
    output_df.loc[exon_df.index, sample_id] = exon_df['Exon_Reads']
    write_gct(output_df, os.path.join(output_dir, sample_id+'.exon_reads.gct.gz'))


parser = argparse.ArgumentParser(description='Convert BAM to FASTQ using SamToFastq from Picard.')
parser.add_argument('bam_file', type=str, help='BAM file')
parser.add_argument('genes_gtf', type=str, help='Gene annotation GTF')
parser.add_argument('genome_fasta', type=str, help='Reference genome FASTA')
parser.add_argument('prefix', type=str, default='Reads', help='Prefix for output files; usually <sample_id>')
parser.add_argument('-o', '--output_dir', default=os.getcwd(), help='Directory to which FASTQs will be written')
parser.add_argument('-m', '--memory', default='4', type=str, help='Memory, in GB')
parser.add_argument('--jar', default='/opt/RNA-SeQC_1.1.9/RNA-SeQC.jar', help='Path to RNA-SeQC.jar')
parser.add_argument('--java_path', default=None, type=str, help='Path to Java 1.7')
parser.add_argument('--rnaseqc_flags', type=str, default=['noDoC', 'strictMode'], nargs='+', help='Optional flags for RNA-SeQC')
parser.add_argument('--gatk_flags', type=str, default=['allow_potentially_misencoded_quality_scores'], nargs='+', help='Optional flags for GATK')
args = parser.parse_args()

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running RNA-SeQC', flush=True)

if args.java_path is None:
    cmd = 'java'
else:
    cmd = args.java_path

cmd += ' -Xmx{}g'.format(args.memory)\
    + ' -jar '+args.jar\
    + ' -n 1000'\
    + ' -s {},{},{}'.format(args.prefix, args.bam_file, args.prefix)\
    + ' -t '+args.genes_gtf\
    + ' -r '+args.genome_fasta\
    + ' -o '+args.output_dir\
    + ' '+' '.join(['-'+i for i in args.rnaseqc_flags])
if args.gatk_flags:
    cmd += ' -gatkFlags ' +' '.join(['--'+i for i in args.gatk_flags])
print('  * command: "{}"'.format(cmd), flush=True)

with cd(args.output_dir):
    subprocess.check_call(cmd, shell=True)

    # remove tmp files
    os.remove('genes.rpkm.gct')
    os.remove('{0}/{0}.metrics.tmp.txt'.format(args.prefix))
    os.remove('{0}/{0}.metrics.txt'.format(args.prefix))
    shutil.move('{0}/{0}.transcripts.rpkm.gct'.format(args.prefix), '{}.gene_rpkm.gct'.format(args.prefix))
    shutil.move('metrics.tsv', '{}.metrics.tsv'.format(args.prefix))
    subprocess.check_call('gzip -f {}.gene_rpkm.gct'.format(args.prefix), shell=True)

    convert_counts('{}.gene_rpkm.gct.gz'.format(args.prefix),
        '{0}/{0}.exon_intron_report.txt'.format(args.prefix),
        '{0}/{0}.exon_report.txt'.format(args.prefix),
        args.genes_gtf, args.prefix, output_dir=args.output_dir)

    subprocess.check_call('tar -czf {0}.tar.gz {0}/*'.format(args.prefix), shell=True)
    shutil.rmtree(os.path.join(args.output_dir, args.prefix))
    for i in ['refGene.txt', 'refGene.txt.idx', 'countMetrics.html']:
        os.remove(i)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finished RNA-SeQC', flush=True)
