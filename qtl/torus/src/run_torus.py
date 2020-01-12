import argparse
import os
import subprocess
from datetime import datetime

parser = argparse.ArgumentParser(description='Run TORUS')
parser.add_argument('qtl_file', help='')
parser.add_argument('annotation_file', help='')
parser.add_argument('prefix', help='Prefix for output file names')
parser.add_argument('-o', '--output_dir', default='./', help='Output directory')
args = parser.parse_args()


def run_torus(qtl_file, annotation_file, prefix, output_dir='.'):
    out_file = os.path.join(output_dir, prefix+'.torus_enrichment.txt')
    log_file = os.path.join(output_dir, prefix+'.torus_enrichment.log')

    s = '['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running Torus'
    print(s, flush=True)
    with open(log_file, 'w') as f:
        f.write(s+'\n')
        f.flush()

    cmd = 'torus \
    -d {} \
    --fastqtl \
    -est \
    -annot {} \
    1> {} \
    2> >(tee -a {} >&2)'.format(qtl_file, annotation_file, out_file, log_file)
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')

    s = '['+datetime.now().strftime("%b %d %H:%M:%S")+'] Done'
    print(s, flush=True)
    with open(log_file, 'a') as f:
        f.write(s+'\n')
        f.flush()


if __name__=='__main__':
    run_torus(args.qtl_file, args.annotation_file, args.prefix, output_dir=args.output_dir)
