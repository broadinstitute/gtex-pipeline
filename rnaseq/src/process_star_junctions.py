import pandas as pd
import numpy as np
import argparse
import os
import gzip

parser = argparse.ArgumentParser(description='Convert STAR junction output (unique mapping reads) to GCT.')
parser.add_argument('star_junction_output', help='SJ.out.tab from STAR')
parser.add_argument('reference_junctions', help='File (tsv) containing columns: chr, intron_start, intron_end, gene_id')
parser.add_argument('prefix', help='Prefix for output file: <prefix>.SJ.gct.gz')
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
dtype = {'chr':str, 'intron_start':np.int32, 'intron_end':np.int32, 'strand':np.int32, 'motif':np.int32, 'status':np.int32, 'n_unique':np.int32, 'n_multi':np.int32, 'max_overhang':np.int32}
junctions_df = pd.read_csv(args.star_junction_output, sep='\t', header=None, names=columns, dtype=dtype)
junctions_df.index = junctions_df['chr']+'_'+junctions_df['intron_start'].astype(str)+'_'+junctions_df['intron_end'].astype(str)

# read reference
reference_df = pd.read_csv(args.reference_junctions, sep='\t', dtype={'chr':str, 'intron_start':np.int32, 'intron_end':np.int32, 'gene_id':str})
reference_df.index = reference_df['chr']+'_'+reference_df['intron_start'].astype(str)+'_'+reference_df['intron_end'].astype(str)

# only consider unique-mapping reads
reference_df['n_unique'] = 0
idx = np.intersect1d(reference_df.index, junctions_df.index)
reference_df.loc[idx, 'n_unique'] = junctions_df.loc[idx, 'n_unique']

# save as GCT
gct_df = reference_df[['gene_id', 'n_unique']].rename(columns={'gene_id':'Description', 'n_unique':args.prefix})
gct_df.index.name = 'Name'
with gzip.open(os.path.join(args.output_dir, args.prefix+'.SJ.gct.gz'), 'wt') as gct:
    gct.write('#1.2\n')
    gct.write('{0:d}\t{1:d}\n'.format(gct_df.shape[0], 1))
    gct_df.to_csv(gct, sep='\t')
