from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
import subprocess
import os
import gzip
import contextlib
from datetime import datetime
import tempfile
import shutil
import glob
from sklearn.decomposition import PCA
import qtl.io


@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Run leafcutter clustering, prepare for FastQTL')
    parser.add_argument('junc_files_list', help='File with paths to ${sample_id}.regtools_junc.txt.gz files')
    parser.add_argument('exons', help='Exon definitions file, with columns: chr, start, end, strand, gene_id, gene_name')
    parser.add_argument('genes_gtf', help='Collapsed gene annotation in GTF format')
    parser.add_argument('prefix', help='Prefix for output files (sample set ID)')
    parser.add_argument('sample_participant_lookup', help='Lookup table linking samples to participants')
    parser.add_argument('--min_clu_reads', default=30, type=int, help='Minimum number of reads supporting each cluster')
    parser.add_argument('--min_clu_ratio', default='0.001', type=str, help='Minimum fraction of reads in a cluster that support a junction')
    parser.add_argument('--max_intron_len', default=500000, type=int, help='Maximum intron length')
    parser.add_argument('--num_pcs', default=10, type=int, help='Number of principal components to calculate')
    parser.add_argument('--leafcutter_dir', default='/opt/leafcutter',
                        help="leafcutter directory, containing 'clustering' directory")
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    args = parser.parse_args()

    print(f'[{datetime.now().strftime("%b %d %H:%M:%S")}] leafcutter clustering')
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    print('  * running leafcutter clustering')
    # generates ${prefix}_perind_numers.counts.gz and ${prefix}_perind.counts.gz
    cmd = f"python3 {os.path.join(args.leafcutter_dir, 'clustering', 'leafcutter_cluster_regtools.py')}" \
        + f' --juncfiles {args.junc_files_list}' \
        + f' --rundir {args.output_dir}' \
        + f' --outprefix {args.prefix}' \
        + f' --minclureads {args.min_clu_reads}' \
        + f' --mincluratio {args.min_clu_ratio}' \
        + f' --maxintronlen {args.max_intron_len}'
    subprocess.check_call(cmd, shell=True)

    # delete sorted junc files
    with open(args.junc_files_list) as f:
        junc_files = f.read().strip().split('\n')
    sorted_files = [os.path.join(args.output_dir, f'{os.path.basename(i)}.{args.prefix}.sorted.gz') for i in junc_files]
    for f in sorted_files:
        os.remove(f)
    os.remove(os.path.join(args.output_dir, f'{args.prefix}_sortedlibs'))

    with cd(args.output_dir):
        print('  * compressing outputs')
        subprocess.check_call(f'gzip -f {args.prefix}_pooled', shell=True)
        subprocess.check_call(f'gzip -f {args.prefix}_refined', shell=True)

        print('  * mapping clusters to genes')
        cmd = 'Rscript' \
            + ' '+os.path.abspath(os.path.join(os.path.dirname(__file__), 'map_clusters_to_genes.R')) \
            + ' '+os.path.join(args.output_dir, f'{args.prefix}_perind.counts.gz') \
            + f' {args.exons}' \
            + f' {args.prefix}.leafcutter.clusters_to_genes.txt'
        subprocess.check_call(cmd, shell=True)

    print('  * filtering counts')
    counts_df = pd.read_csv(os.path.join(args.output_dir, f'{args.prefix}_perind.counts.gz'), sep='\s+').set_index('chrom')
    calculate_frac = lambda x: float(x[0])/float(x[1]) if x[1] > 0 else 0
    frac_df = counts_df.applymap(lambda x: calculate_frac([int(i) for i in x.split('/')]))
    pct_zero = (frac_df == 0).sum(1) / frac_df.shape[1]  # for zero counts, frac is zero
    n_unique = frac_df.apply(lambda x: len(x.unique()), axis=1)
    zscore_df = ((frac_df.T-frac_df.mean(1)) / frac_df.std(1)).T

    # filter out introns with low counts or low complexity
    n = np.floor(frac_df.shape[1]*0.1)
    if n < 10:
        n = 10
    mask = (pct_zero <= 0.5) & (n_unique >= n)
    # additional filter for low complexity
    ns = zscore_df.shape[1]
    mask2 = ((zscore_df.abs()<0.25).sum(1) >= ns-3) & ((zscore_df.abs() > 6).sum(1) <= 3)
    if np.any(mask & mask2):
        print(f'    ** dropping {np.sum(mask & mask2)} introns with low variation')
    mask = mask & ~mask2

    filtered_counts_df = counts_df.loc[mask].copy()
    cluster_ids = np.unique(counts_df.index.map(lambda x: x.split(':')[-1]))
    filtered_cluster_ids = np.unique(filtered_counts_df.index.map(lambda x: x.split(':')[-1]))
    print('    ** dropping {} introns with counts in fewer than 50% of samples\n'
          '       {}/{} introns remain ({}/{} clusters)'.format(
               counts_df.shape[0]-filtered_counts_df.shape[0], filtered_counts_df.shape[0],
               counts_df.shape[0], len(filtered_cluster_ids), len(cluster_ids))
        )
    col_dict = {i:i.split('.')[0] for i in filtered_counts_df.columns}
    filtered_counts_df.rename(columns=col_dict, inplace=True)
    sample_participant_lookup_s = pd.read_csv(args.sample_participant_lookup,
                                              sep='\t', index_col=0, dtype=str, squeeze=True)
    assert filtered_counts_df.columns.isin(sample_participant_lookup_s.index).all()

    filtered_counts_file = os.path.join(args.output_dir, f'{args.prefix}_perind.counts.filtered.gz')
    filtered_counts_df.to_csv(filtered_counts_file, sep=' ')

    print('  * preparing phenotype table')
    subprocess.check_call(
        'python3 '+os.path.join(args.leafcutter_dir, 'scripts', 'prepare_phenotype_table.py') \
        +f' {filtered_counts_file}' \
        +f' -p {args.num_pcs}', shell=True)

    print('  * concatenating chromosome-level BED files')
    bed_files = sorted(glob.glob(os.path.join(args.output_dir, '*_perind.counts.filtered.gz.qqnorm_*')))
    bed_df = []
    for f in bed_files:
        bed_df.append(pd.read_csv(f, sep='\t', dtype=str))
    bed_df = pd.concat(bed_df, axis=0)
    bed_df['chr_ix'] = bed_df['#Chr'].str.replace('chr','').str.replace('X','23').str.replace('Y','24').astype(np.int32)
    for c in ['start', 'end']:
        bed_df[c] = bed_df[c].astype(np.int32)
    bed_df.sort_values(['chr_ix', 'start', 'end'], inplace=True)
    bed_df.drop('chr_ix', axis=1, inplace=True)
    bed_df.rename(columns={'#Chr':'#chr'}, inplace=True)
    bed_df.rename(columns=col_dict, inplace=True)
    print('    ** writing merged BED')
    bed_file = os.path.join(args.output_dir, f'{args.prefix}.perind.counts.filtered.qqnorm.bed.gz')
    qtl.io.write_bed(bed_df, bed_file)

    print('  * converting cluster coordinates to gene coordinates')
    tss_df = qtl.io.gtf_to_tss_bed(args.genes_gtf)
    cluster2gene_dict = pd.read_csv(os.path.join(args.output_dir, f'{args.prefix}.leafcutter.clusters_to_genes.txt'),
        sep='\t', index_col=0, squeeze=True).to_dict()

    print('    ** assigning introns to gene mapping(s)')
    n = 0
    gene_bed_df = []
    group_s = {}
    for _,r in bed_df.iterrows():
        s = r['ID'].split(':')
        cluster_id = s[0]+':'+s[-1]
        if cluster_id in cluster2gene_dict:
            gene_ids = cluster2gene_dict[cluster_id].split(',')
            for g in gene_ids:
                gi = r['ID']+':'+g
                gene_bed_df.append(tss_df.loc[g, ['chr', 'start', 'end']].tolist() + [gi] + r.iloc[4:].tolist())
                group_s[gi] = g
        else:
            n += 1
    if n > 0:
        print(f'    ** discarded {n} introns without a gene mapping')

    print('  * writing BED files for QTL mapping')
    gene_bed_df = pd.DataFrame(gene_bed_df, columns=bed_df.columns)
    # sort by TSS
    gene_bed_df = gene_bed_df.groupby('#chr', sort=False, group_keys=False).apply(lambda x: x.sort_values('start'))
    # change sample IDs to participant IDs
    gene_bed_df.rename(columns=sample_participant_lookup_s, inplace=True)
    qtl.io.write_bed(gene_bed_df, os.path.join(args.output_dir, f'{args.prefix}.leafcutter.bed.gz'))
    gene_bed_df[['start', 'end']] = gene_bed_df[['start', 'end']].astype(np.int32)
    gene_bed_df[gene_bed_df.columns[4:]] = gene_bed_df[gene_bed_df.columns[4:]].astype(np.float32)
    gene_bed_df.to_parquet(os.path.join(args.output_dir, f'{args.prefix}.leafcutter.bed.parquet'))
    pd.Series(group_s).sort_values().to_csv(
        os.path.join(args.output_dir, f'{args.prefix}.leafcutter.phenotype_groups.txt'),
        sep='\t', header=False)

    print('  * calculating PCs')
    pca = PCA(n_components=args.num_pcs)
    pca.fit(bed_df[bed_df.columns[4:]])
    pc_df = pd.DataFrame(pca.components_, index=[f'PC{i}' for i in range(1,11)],
        columns=['-'.join(i.split('-')[:2]) for i in bed_df.columns[4:]])
    pc_df.index.name = 'ID'
    pc_df.to_csv(os.path.join(args.output_dir, f'{args.prefix}.leafcutter.PCs.txt'), sep='\t')

    # delete intermediary files
    files = glob.glob(os.path.join(args.output_dir, f'{args.prefix}_perind.counts.filtered.gz.qqnorm_chr*')) \
          + glob.glob(os.path.join(args.output_dir, f'{args.prefix}_perind.counts.filtered.gz.phen_chr*'))
    for f in files:
        os.remove(f)
    os.remove(os.path.join(args.output_dir, f'{args.prefix}_perind.counts.filtered.gz_prepare.sh'))
    os.remove(os.path.join(args.output_dir, f'{args.prefix}.leafcutter.clusters_to_genes.txt'))
    os.remove(os.path.join(args.output_dir, f'{args.prefix}.perind.counts.filtered.qqnorm.bed.gz'))
    os.remove(os.path.join(args.output_dir, f'{args.prefix}.perind.counts.filtered.qqnorm.bed.gz.tbi'))
    os.remove(os.path.join(args.output_dir, f'{args.prefix}_perind.counts.filtered.gz.ave'))
    os.remove(os.path.join(args.output_dir, f'{args.prefix}_perind.counts.filtered.gz.PCs'))

    print(f'[{datetime.now().strftime("%b %d %H:%M:%S")}] done')
