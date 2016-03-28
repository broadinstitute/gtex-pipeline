#!/usr/bin/env python3
# Author: Francois Aguet

import argparse
import os
import subprocess
import gzip
import shutil
from datetime import datetime
import contextlib

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)


parser = argparse.ArgumentParser(description='Run STAR')
parser.add_argument('index', help='Path to STAR index')
parser.add_argument('fastq1', help='FASTQ file containing the first read mates')
parser.add_argument('fastq2', help='FASTQ file containing the second read mates')
parser.add_argument('prefix', help='Prefix for output file names')
parser.add_argument('-o', '--output_dir', default='./', help='Output directory')
parser.add_argument('--annotation_gtf', default=None, help='Annotation in GTF format')
parser.add_argument('--outFilterMultimapNmax', default='20')
parser.add_argument('--alignSJoverhangMin', default='8')
parser.add_argument('--alignSJDBoverhangMin', default='1')
parser.add_argument('--outFilterMismatchNmax', default='999')
parser.add_argument('--outFilterMismatchNoverLmax', default='0.1')
parser.add_argument('--alignIntronMin', default='20')
parser.add_argument('--alignIntronMax', default='1000000')
parser.add_argument('--alignMatesGapMax', default='1000000')
parser.add_argument('--outFilterType', default='BySJout')
parser.add_argument('--outFilterScoreMinOverLread', default='0.33')
parser.add_argument('--outFilterMatchNminOverLread', default='0.33')
parser.add_argument('--limitSjdbInsertNsj', default='1200000')
parser.add_argument('--outSAMstrandField', default='intronMotif')
parser.add_argument('--outFilterIntronMotifs', default='None', help="Use 'RemoveNoncanonical' for Cufflinks compatibility")
parser.add_argument('--alignSoftClipAtReferenceEnds', default='Yes')
parser.add_argument('--quantMode', default='TranscriptomeSAM GeneCounts', help='Outputs read counts, and a BAM with reads in transcriptome coordinates')
parser.add_argument('--outSAMtype', default='BAM Unsorted')
parser.add_argument('--outSAMunmapped', default='Within', help='Keep unmapped reads in output BAM')
parser.add_argument('--outSAMattrRGline', default='ID:rg1 SM:sm1', help='Adds read group line to BAM header; required by GATK')
parser.add_argument('--outSAMattributes', default='NH HI AS nM NM')
parser.add_argument('--chimSegmentMin', default='15', help='Minimum chimeric segment length; switches on detection of chimeric (fusion) alignments')
parser.add_argument('--chimJunctionOverhangMin', default='15', help='Minimum overhang for a chimeric junction')
parser.add_argument('--genomeLoad', default='NoSharedMemory')
parser.add_argument('-t', '--threads', default='4', help='Number of threads')
args = parser.parse_args()

if args.fastq1.endswith('.gz'):
    with gzip.open(args.fastq1) as f:
        f.readline()
        seq = f.readline().strip()
        read_length = len(seq)
else:
    with open(args.fastq1) as f:
        f.readline()
        seq = f.readline().strip()
        read_length = len(seq)
overhang = read_length-1

if read_length<=101:
    starcmd='STAR'
else:
    starcmd='STARlong'

# set up command
cmd = starcmd+' --runMode alignReads --runThreadN '+args.threads+' --genomeDir '+args.index
if args.annotation_gtf is not None:  # only needed if genome index was built w/o annotation
    cmd += ' --sjdbGTFfile '+args.annotation_gtf
cmd += ' --twopassMode Basic --sjdbOverhang '+str(overhang)\
    +' --outFilterMultimapNmax '+args.outFilterMultimapNmax\
    +' --alignSJoverhangMin '+args.alignSJoverhangMin+' --alignSJDBoverhangMin '+args.alignSJDBoverhangMin\
    +' --outFilterMismatchNmax '+args.outFilterMismatchNmax+' --outFilterMismatchNoverLmax '+args.outFilterMismatchNoverLmax\
    +' --alignIntronMin '+args.alignIntronMin+' --alignIntronMax '+args.alignIntronMax+' --alignMatesGapMax '+args.alignMatesGapMax\
    +' --outFilterType '+args.outFilterType\
    +' --outFilterScoreMinOverLread '+args.outFilterScoreMinOverLread+' --outFilterMatchNminOverLread '+args.outFilterMatchNminOverLread\
    +' --limitSjdbInsertNsj '+args.limitSjdbInsertNsj\
    +' --readFilesIn '+args.fastq1+' '+args.fastq2
if args.fastq1.endswith('.gz'):
    cmd += ' --readFilesCommand zcat'
cmd += ' --outFileNamePrefix '+os.path.join(args.output_dir, args.prefix)+'.'\
    +' --outSAMstrandField '+args.outSAMstrandField+' --outFilterIntronMotifs '+args.outFilterIntronMotifs\
    +' --alignSoftClipAtReferenceEnds '+args.alignSoftClipAtReferenceEnds+' --quantMode '+args.quantMode\
    +' --outSAMtype '+args.outSAMtype+' --outSAMunmapped '+args.outSAMunmapped+' --genomeLoad '+args.genomeLoad
if int(args.chimSegmentMin)>0:
    cmd += ' --chimSegmentMin '+args.chimSegmentMin+' --chimJunctionOverhangMin '+args.chimJunctionOverhangMin
cmd += ' --outSAMattributes '+args.outSAMattributes+' --outSAMattrRGline '+args.outSAMattrRGline

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir )

# run STAR
subprocess.check_call(cmd, shell=True, executable='/bin/bash')

# postprocessing
with cd(args.output_dir):
    # set permissions
    for r,d,f in os.walk(args.prefix+'._STARpass1'):
        os.chmod(r, 0o755)

    # delete unneeded files
    shutil.rmtree(args.prefix+'._STARgenome')
    if os.path.exists(args.prefix+'._STARtmp'):
        shutil.rmtree(args.prefix+'._STARtmp')

    # sort BAM (use samtools to get around the memory gluttony of STAR)
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Sorting BAM', flush=True)
    cmd = 'samtools sort --threads '+args.threads+' -o '+args.prefix+'.Aligned.sortedByCoord.out.bam '+args.prefix+'.Aligned.out.bam'
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    os.remove(args.prefix+'.Aligned.out.bam')
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finished sorting BAM', flush=True)

    # index BAM
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Indexing BAM', flush=True)
    cmd = 'samtools index '+args.prefix+'.Aligned.sortedByCoord.out.bam'
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finished indexing BAM', flush=True)

    # sort and index chimeric BAM
    if int(args.chimSegmentMin)>0:
        cmd = 'samtools sort --threads '+args.threads+' -o '+args.prefix+'.Chimeric.out.sorted.bam '+args.prefix+'.Chimeric.out.sam'
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        cmd = 'samtools index '+args.prefix+'.Chimeric.out.sorted.bam'
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        os.remove(args.prefix+'.Chimeric.out.sam')
