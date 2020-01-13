from __future__ import print_function
import argparse
import os
import sys
import subprocess
import pandas
import gzip
import glob
import re

PIPELINE_PATH = os.path.join(
    os.path.abspath(os.path.dirname(__file__)),
)

chr_pattern = re.compile(r'chr(\w+)\..+')


# def filetype(arg):
#     if os.path.isfile(arg):
#         return os.path.abspath(arg)
#     raise argparse.ArgumentTypeError("No such file: "+arg)

def filetype(*exts):
    def checker(arg):
        if os.path.isfile(arg):
            extensions = os.path.basename(arg).lower().split('.')
            while len(extensions) and extensions[-1] in {'tar', 'gz', 'bz2', 'zip'}:
                extensions.pop()
            if extensions[-1] in exts:
                return os.path.abspath(arg)
            raise argparse.ArgumentTypeError("Invalid file type: %s (%s)" % (
                arg,
                extensions[-1]
            ))
        raise argparse.ArgumentTypeError("No such file: "+arg)
    return checker

def tsv_or_file(ext):
    def checker(arg):
        bamtype = filetype(ext)
        if not os.path.isfile(arg):
            raise argparse.ArgumentTypeError("No such file: "+arg)
        try:
            return [bamtype(arg)]
        except argparse.ArgumentTypeError as e:
            with open(arg) as reader:
                return [bamtype(line.strip()) for line in reader]
    return checker



def dirtype(arg):
    if os.path.isfile(arg):
        raise argparse.ArgumentTypeError(arg+" is a file")
    try:
        os.makedirs(arg)
    finally:
        if not os.path.isdir(arg):
            raise argparse.ArgumentTypeError("Unable to create directory: "+arg)
        return os.path.abspath(arg)

def run(args):
    print("Preparing...")

    args.input = [
        filepath
        for entry in args.input
        for filepath in entry
    ]
    vcf_ids = {}
    with gzip.open(args.vcf, 'r') as reader:
        for line in reader:
            if '#CHR' in line:
                cols = line.strip().split('\t')
                for col in cols[9:]:
                    sample_id = '-'.join(col.split('-')[:2])
                    vcf_ids[sample_id] = ''+col
                break

    if args.individual not in vcf_ids:
        sys.exit("ERROR: Individual %s not found in VCF" % args.individual)

    bams = []
    mapqs = []
    isizes = []
    haplo_count_excluded_bams = []

    if args.wes_bam:
        bams.append(args.wes_bam)
        mapqs.append('30')
        isizes.append('500')
        haplo_count_excluded_bams.append(len(bams))

    if args.wgs_bam:
        bams.append(args.wgs_bam)
        mapqs.append('30')
        isizes.append('1000')
        haplo_count_excluded_bams.append(len(bams))

    for bam in args.input:
        bams.append(bam)
        mapqs.append('255')
        isizes.append('1e6')

    print("Running pHASER...")
    for chromosome in (range(1,23) if args.chr is None else [args.chr]):
        command = [
            sys.executable,
            os.path.join(
                PIPELINE_PATH,
                'phaser',
                'phaser.py'
            ),
            '--chr',
            str(chromosome),
            '--temp_dir',
            '.',
            '--bam',
            ','.join(bams),
            '--vcf',
            args.vcf,
            '--sample',
            args.individual,
            '--baseq',
            '10',
            '--mapq',
            ','.join(mapqs),
            '--isize',
            ','.join(isizes),
            '--paired_end',
            '1',
            '--gw_af_field',
            args.af_field,
            '--gw_phase_method',
            '1',
            '--o',
            os.path.join(
                args.output,
                '%s.%s' % (
                    args.individual,
                    str(chromosome)
                )
            ),
            '--include_indels',
            '0',
            '--gw_phase_vcf',
            '1',
        ]
        if args.blacklist:
            command += [
                '--blacklist',
                args.blacklist
            ]
        if args.haplo_count_blacklist:
            command += [
                '--haplo_count_blacklist',
                args.haplo_count_blacklist
            ]
        if len(haplo_count_excluded_bams):
            command += [
                '--haplo_count_bam_exclude',
                ','.join(haplo_count_excluded_bams)
            ]
        command = ' '.join(command)
        print(command)
        subprocess.check_call(
            'set -euo pipefail && ' + command,
            shell=True,
        executable='/bin/bash'
        )

def post(args):
    print("Merging chromosome VCFs...")
    command = [
        'bcftools',
        'concat',
    ]
    command += sorted(
        args.vcfs, #>args.vcfs : output/*chr*.vcf.gz
        key=lambda path:int(chr_pattern.search(os.path.basename(path)).group(1))
    )
    command += [
        '|',
        'bgzip',
        '-c',
        '>',
        os.path.join(
            args.output,
            '%s.vcf.gz' % args.individual #<output/individual.vcf.gz
        )
    ]
    command = ' '.join(command)
    print(command)
    subprocess.check_call(
        'set -euo pipefail && ' + command,
        shell=True,
        executable='/bin/bash',
    )


    command = [
        'tabix',
        '-fp',
        'vcf',
        os.path.join(
            args.output,
            '%s.vcf.gz' % args.individual #<output/individual.vcf.gz
        )
    ]
    command = ' '.join(command)
    print(command)
    subprocess.check_call(
        'set -euo pipefail && ' + command,
        shell=True,
        executable='/bin/bash'
    )


    print("Generating gene AE...")
    command = [
        'head',
        '-n',
        '1',
        #we only need the header row
        args.haplotypic_counts[0], #>args.haplotypic_counts : output/*chr*.haplotypic_counts.txt
        '>',
        os.path.join(
            args.output,
            '%s.haplotypic_counts.txt' % args.individual #<output/individual.haplotypic_counts.txt
        )
    ]
    command = ' '.join(command)
    print(command)
    subprocess.check_call(
        'set -euo pipefail && ' + command,
        shell=True,
        executable='/bin/bash'
    )

    command = [
        'cat'
    ]
    command += sorted(
        args.haplotypic_counts,
        key=lambda path:int(chr_pattern.search(os.path.basename(path)).group(1))
    )
    command += [
        '|',
        'grep',
        '-v',
        'contig', #needs quotes?
        '|',
        'sed',
        's/.Aligned.sortedByCoord.out.patched//g',
        '>>',
        os.path.join(
            args.output,
            '%s.haplotypic_counts.txt' % args.individual #<output/individual.haplotypic_counts.txt
        )
    ]
    command = ' '.join(command)
    print(command)
    subprocess.check_call(
        'set -euo pipefail && ' + command,
        shell=True,
        executable='/bin/bash'
    )


    command = [
        sys.executable,
        os.path.join(
            PIPELINE_PATH,
            'phaser_gene_ae',
            'phaser_gene_ae.py'
        ),
        '--haplotypic_counts',
        os.path.join(
            args.output,
            '%s.haplotypic_counts.txt' % args.individual #<output/individual.haplotypic_counts.txt
        ),
        '--features',
        args.gene_models,
        '--o',
        os.path.join(
            args.output,
            '%s.gene_ae.txt' % args.individual #<output/individual.gene_ae.txt
        ),
        '--min_haplo_maf',
        '0.05'
    ]
    command = ' '.join(command)
    print(command)
    subprocess.check_call(
        'set -euo pipefail && ' + command,
        shell=True,
        executable='/bin/bash'
    )


    command = [
        'gzip',
        '-f',
        os.path.join(
            args.output,
            '%s.gene_ae.txt' % args.individual #<output/individual.gene_ae.txt
        )
    ]
    command = ' '.join(command)
    print(command)
    subprocess.check_call(
        'set -euo pipefail && ' + command,
        shell=True,
        executable='/bin/bash'
    )


    command = [
        'gzip',
        '-f',
        os.path.join(
            args.output,
            '%s.haplotypic_counts.txt' % args.individual #<output/individual.haplotypic_counts.txt
        )
    ]
    command = ' '.join(command)
    print(command)
    subprocess.check_call(
        'set -euo pipefail && ' + command,
        shell=True,
        executable='/bin/bash'
    )


    print("Cleaning up temporary files")
    #No need to do this.  It'll be passed in as cromwell args, so don't extract as output
    for path in glob.iglob(os.path.join(args.output, '*chr*')):
        os.remove(path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        "phASER"
    )
    subparsers = parser.add_subparsers()
    phaser_parser = subparsers.add_parser('phase')
    phaser_parser.set_defaults(func=run)
    phaser_parser.add_argument(
        'individual',
        help="Unique ID for the individual"
    )
    phaser_parser.add_argument(
        'input',
        type=tsv_or_file('bam'),
        nargs='+',
        help="Input RNA-Seq Bams.  May be supplied multiple times."
        " You may also provide a txt or single-column tsv with a filepath "
        "to a Bam on each line"
    )
    phaser_parser.add_argument(
        'vcf',
        type=filetype('vcf'),
        help="Genotype VCF. Sample IDs listed in the input VCF "
        "must appear in the VCF's columns"
    )
    phaser_parser.add_argument(
        'gene_models',
        type=filetype('bed'),
        help="BED file containing gene models for phaser_gene_ae",
    )
    phaser_parser.add_argument(
        'output',
        type=dirtype,
        help="Output directory.  Will be created if it doesn't exist"
    )
    phaser_parser.add_argument(
        '-a', '--af-field',
        help="Field in VCF to use for allele frequency",
        default='AF'
    )
    phaser_parser.add_argument(
        '-e', '--wes-bam',
        type=filetype('bam'),
        help='Path to file containing RNA-Seq WES Bam locations',
        default=None
    )
    phaser_parser.add_argument(
        '-g', '--wgs-bam',
        type=filetype('bam'),
        help='Path to file containing RNA-Seq WGS Bam locations',
        default=None
    )
    phaser_parser.add_argument(
        '-c', '--haplo-count-blacklist',
        type=filetype('bed'),
        help='BED file containing genomic intervals to be excluded from haplotypic counts.'
    )
    phaser_parser.add_argument(
        '-b', '--blacklist',
        type=filetype('bed'),
        help="BED file containing genomic intervals to be excluded from phasing"
    )
    phaser_parser.add_argument(
        '--chr',
        # type=int,
        help="Chromosome to run.  All chromosomes will be run in series if this is not specified",
        default=None
    )

    postprocess_parser = subparsers.add_parser('postprocess')
    postprocess_parser.set_defaults(func=post)
    postprocess_parser.add_argument(
        'individual',
        help="Unique ID for the individual"
    )
    postprocess_parser.add_argument(
        'vcfs',
        type=tsv_or_file('vcf'),
        help="Input gziped vcfs, or a file with a vcf filepath on each line"
    )
    postprocess_parser.add_argument(
        'haplotypic_counts',
        type=tsv_or_file('txt'),
        help="Input haplotypic counts txts, or a file with a txt filepath on each line"
    )
    postprocess_parser.add_argument(
        'gene_models',
        type=filetype('bed'),
        help="BED file containing gene models for phaser_gene_ae",
    )
    postprocess_parser.add_argument(
        'output',
        type=dirtype,
        help="Output directory.  Will be created if it doesn't exist"
    )

    args = parser.parse_args()

    args.func(args)
