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


def get_introns(intron_counts_df):
    intron_df = pd.DataFrame([i.split(':') for i in intron_counts_df.index],
                             columns=['chr', 'start', 'end', 'clu'])
    intron_df['start'] = intron_df['start'].astype(np.int32)
    intron_df['end'] = intron_df['end'].astype(np.int32)
    return intron_df


def map_cluster_to_genes(intron_df, exon_df):
    matches_df = []
    for c in sorted(intron_df['chr'].unique()):
        introns_chr = intron_df[intron_df['chr'] == c]
        exons_chr = exon_df[exon_df['chr'] == c]
        three_prime_matches = introns_chr.merge(exons_chr, left_on='end', right_on='start', how='inner')
        five_prime_matches = introns_chr.merge(exons_chr, left_on='start', right_on='end', how='inner')

        all_matches = pd.concat([three_prime_matches, five_prime_matches])[['clu', 'gene_id']].drop_duplicates()
        all_matches['clu'] = c + ':' + all_matches['clu']
        matches_df.append(all_matches)
    matches_df = pd.concat(matches_df).reset_index(drop=True)
    clu_s = matches_df.groupby('clu', sort=True).apply(
        lambda x: ','.join(x['gene_id']), include_groups=False).rename('genes')
    return clu_s


def run_leafcutter_clustering(junc_files_list, leafcutter_dir, prefix, output_dir='.',
                              min_clu_reads=30, min_clu_ratio=0.001, max_intron_len=500000):
    """
    Generates ${prefix}_perind_numers.counts.gz, ${prefix}_perind.counts.gz, ${prefix}.leafcutter.clusters_to_genes.txt
    """
    print('  * running leafcutter clustering')
    cmd = f"python3 {os.path.join(leafcutter_dir, 'clustering', 'leafcutter_cluster_regtools.py')} \
        --juncfiles {junc_files_list} \
        --outprefix {prefix} \
        --rundir {output_dir} \
        --minclureads {min_clu_reads} \
        --mincluratio {min_clu_ratio} \
        --maxintronlen {max_intron_len}"
    subprocess.check_call(cmd, shell=True)

    # delete sorted junc files
    with open(junc_files_list) as f:
        junc_files = f.read().strip().split('\n')
    sorted_files = [os.path.join(output_dir, f'{os.path.basename(i)}.{prefix}.sorted.gz') for i in junc_files]
    for f in sorted_files:
        os.remove(f)
    os.remove(os.path.join(output_dir, f'{prefix}_sortedlibs'))

    with cd(output_dir):
        print('  * compressing outputs')
        subprocess.check_call(f'gzip -f {prefix}_pooled', shell=True)
        subprocess.check_call(f'gzip -f {prefix}_refined', shell=True)


def filter_counts(counts_df):
    calculate_frac = lambda x: float(x[0])/float(x[1]) if x[1] > 0 else 0
    frac_df = counts_df.map(lambda x: calculate_frac([int(i) for i in x.split('/')]))
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
    # drop chrY
    filtered_counts_df = filtered_counts_df[~filtered_counts_df.index.str.startswith(('chrY:','Y:'))]
    return filtered_counts_df



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

    sample_participant_lookup_s = pd.read_csv(args.sample_participant_lookup,
                                              sep='\t', index_col=0, dtype=str).squeeze('columns')

    run_leafcutter_clustering(args.junc_files_list, args.leafcutter_dir, args.prefix, output_dir=args.output_dir,
                              min_clu_reads=args.min_clu_reads, min_clu_ratio=args.min_clu_ratio, max_intron_len=args.max_intron_len)

    counts_df = pd.read_csv(os.path.join(args.output_dir, f'{args.prefix}_perind.counts.gz'), sep='\s+', index_col=0)
    # change columns to sample IDs
    col_dict = {i:i.split('.')[0] for i in counts_df.columns}
    counts_df.rename(columns=col_dict, inplace=True)
    assert counts_df.columns.isin(sample_participant_lookup_s.index).all()

    print('  * mapping clusters to genes')
    intron_df = get_introns(counts_df)
    exon_df = pd.read_csv(args.exons, sep='\t')
    clu_s = map_cluster_to_genes(intron_df, exon_df)
    clu_s.to_csv(os.path.join(args.output_dir, f'{args.prefix}.leafcutter.clusters_to_genes.txt'), sep='\t')
    cluster2gene_dict = clu_s.to_dict()

    print('  * filtering counts')
    filtered_counts_df = filter_counts(counts_df)
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
    bed_df = pd.concat(bed_df, axis=0).rename(columns={'#Chr':'chr'})
    # sort BED
    bed_df['chr_ix'] = bed_df['chr'].str.replace('chr','').str.replace('X','23').astype(np.int32)
    for c in ['start', 'end']:
        bed_df[c] = bed_df[c].astype(np.int32)
    bed_df.sort_values(['chr_ix', 'start', 'end'], inplace=True)
    bed_df.drop('chr_ix', axis=1, inplace=True)

    print('    ** writing merged BED')
    bed_file = os.path.join(args.output_dir, f'{args.prefix}.perind.counts.filtered.qqnorm.bed.gz')
    qtl.io.write_bed(bed_df, bed_file)

    print('    ** adding TSSs')
    tss_df = qtl.io.gtf_to_tss_bed(args.genes_gtf)
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
    gene_bed_df = pd.concat([gdf.sort_values(['start', 'end']) for _,gdf in gene_bed_df.groupby('chr', sort=False, as_index=False)]).reset_index(drop=True)
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
    pc_df = pd.DataFrame(pca.components_, index=[f'PC{i}' for i in range(1, args.num_pcs+1)],
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
