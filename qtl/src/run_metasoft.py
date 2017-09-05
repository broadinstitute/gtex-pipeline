#!/usr/bin/env python3
import argparse
import subprocess
import os
import tempfile

parser = argparse.ArgumentParser(description='Run METASOFT.')
parser.add_argument('metasoft_jar', help='metasoft.jar')
parser.add_argument('metasoft_input', help='')
parser.add_argument('prefix', help='Prefix for output file: <prefix>.metasoft.txt.gz')
parser.add_argument('--pvalue_table', default=None, help='HanEskinPvalueTable.txt if not in same directory as metasoft.jar')
parser.add_argument('--seed', default=100, type=int, help='')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

if args.pvalue_table is not None:
    pvalue_table = args.pvalue_table
else:
    pvalue_table = os.path.join(os.path.split(args.metasoft_jar)[0], 'HanEskinPvalueTable.txt')
    assert os.path.exists(pvalue_table)

if args.metasoft_input.endswith('.gz'):
    input_file = os.path.splitext(os.path.split(args.metasoft_input)[1])[0]
    subprocess.check_call('gunzip -c {} > {}'.format(args.metasoft_input, input_file), shell=True, executable='/bin/bash')
else:
    input_file = args.metasoft_input

# strip header if present
with open(input_file) as f:
    header = f.readline().strip().split('\t')
    if header[0]=='pair_id':
        subprocess.check_call('tail -n+2 {0} > tmp.txt && mv tmp.txt {0}'.format(input_file), shell=True, executable='/bin/bash')

cmd = 'java -jar '+args.metasoft_jar\
    +' -input '+input_file\
    +' -pvalue_table '+pvalue_table\
    +' -output '+os.path.join(args.output_dir, args.prefix+'.metasoft.txt')\
    +' -mvalue_p_thres 1.0'\
    +' -mvalue_method mcmc'\
    +' -seed '+str(args.seed)\
    +' -log '+os.path.join(args.output_dir, args.prefix+'.metasoft.log')\
    +' -mvalue'
subprocess.check_call(cmd, shell=True, executable='/bin/bash')
subprocess.check_call('gzip -1 '+args.prefix+'.metasoft.txt', shell=True, executable='/bin/bash')
