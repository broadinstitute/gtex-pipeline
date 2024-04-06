import pandas as pd
import numpy as np
import argparse
import os
import gzip
import qtl.io

parser = argparse.ArgumentParser(description='Convert STAR junction output (unique mapping reads) to GCT.')
parser.add_argument('star_junction_output', help='SJ.out.tab from STAR')
parser.add_argument('reference_junctions', help='File (tsv) containing columns: chr, intron_start, intron_end, gene_id')
parser.add_argument('prefix', help='Prefix for output file: <prefix>.SJ.gct.gz')
parser.add_argument('--parquet', action='store_true', help='Write to parquet format instead of GCT')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

# STAR junction format:
#  0) chromosome
#  1) 1st base of intron (1-based)
#  2) last base of intron (1-based)
#  3) strand (0: undefined, 1: +, 2: -)
#  4) intron motif: 0: non-canonical
#                   1: GT/AG
#                   2: CT/AC
#                   3: GC/AG
#                   4: CT/GC
#                   5: AT/AC
#                   6: GT/AT
#  5) annotation status: 0: unannotated
#                        1: annotated
#  6) number of uniquely mapping reads crossing the junction
#  7) number of multi-mapping reads crossing the junction
#  8) maximum spliced alignment overhang


# read STAR output
columns = ['chr', 'intron_start', 'intron_end', 'strand', 'motif', 'status', 'n_unique', 'n_multi', 'max_overhang']
dtype = {'chr':str, 'intron_start':np.int32, 'intron_end':np.int32, 'strand':np.int32, 'motif':np.int32,
         'status':np.int32, 'n_unique':np.int32, 'n_multi':np.int32, 'max_overhang':np.int32}
junctions_df = pd.read_csv(args.star_junction_output, sep='\t', header=None, names=columns, dtype=dtype)
junctions_df['strand'] = junctions_df['strand'].map({0: '?', 1:'+', 2:'-'})
junctions_df.index = (junctions_df['chr'] + ':' + junctions_df['intron_start'].astype(str)
                      + '-' + junctions_df['intron_end'].astype(str) + ':' + junctions_df['strand'])

# read reference
reference_df = pd.read_csv(args.reference_junctions, sep='\t',
                           dtype={'chr':str, 'intron_start':np.int32, 'intron_end':np.int32, 'gene_id':str})
reference_df.index = (reference_df['chr'] + ':' + reference_df['intron_start'].astype(str)
                      + '-' + reference_df['intron_end'].astype(str) + ':' + reference_df['strand'])
assert not reference_df.index.duplicated().any(), "Junction annotation must not contain any duplicated entries"

# use unique-mapping reads
gct_df = reference_df[['gene_id']].join(junctions_df['n_unique'])
gct_df['n_unique'] = gct_df['n_unique'].fillna(0).astype(np.int32)
# write as GCT
gct_df.rename(columns={'gene_id':'Description', 'n_unique':args.prefix}, inplace=True)
gct_df.index.name = 'Name'
if args.parquet:
    gct_df.to_parquet(os.path.join(args.output_dir, f"{args.prefix}.SJ.parquet"))
else:
    qtl.io.write_gct(gct_df, os.path.join(args.output_dir, f"{args.prefix}.SJ.gct.gz"))
