<!-- Author: Francois Aguet -->
## sQTL mapping pipeline

This document describes the sQTL mapping pipeline used by the GTEx Consortium. For additional details please see Section 4.3 of the Supplementary Materials for [GTEx Consortium, Science, 2020](https://www.science.org/doi/suppl/10.1126/science.aaz1776/suppl_file/aaz1776_aguet_sm.pdf).

### Docker image
The pipeline components described below are available in a [Docker image](https://hub.docker.com/r/francois4/leafcutter/). To download the image, run:
```bash
docker pull francois4/leafcutter:latest
```

#### Image contents and pipeline components
The following tools are included in the image:

1. [regtools](https://regtools.readthedocs.io/en/latest/) for extracting exon-exon junctions from RNAseq BAM files.
2. [LeafCutter](https://davidaknowles.github.io/leafcutter/) for generating intron excision ratios across samples.

### Running the pipeline
This pipeline consists of four steps, summarized below: variant-aware alignment to correct for allelic mapping bias, extraction of exon-exon junctions, generation of phenotype tables with intron excision ratios, and QTL mapping.

#### 1) Generating WASP-corrected alignments with STAR
Allelic mapping bias can be a significant confounder for sQTL mapping and quantifying allele-specific expression ([Castel et al., 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0762-6)). This can be mitigated using variant-aware alignment, e.g., with the (WASP correction)[https://www.nature.com/articles/nmeth.3582] implemented in [STAR](https://github.com/alexdobin/STAR). Using the [STAR wrapper script](../../../rnaseq/src/run_STAR.py) in this repository, this step is applied when supplying a VCF with the participant's SNPs via the `--varVCFfile` option. The participant VCFs can be generated using [this WDL](../../../genotype/participant_vcfs.wdl).

#### 2) Generating exon-exon junction counts with regtools
Next, exon-exon junction counts are extracted from the BAM files using
```bash
samtools view -h -q 255 $bam_file | grep -v "vW:i:[2-7]" | samtools view -b > $filtered_bam
regtools junctions extract -a 8 -m 50 -M 500000 -s 0 $filtered_bam | gzip -c > ${sample_id}.regtools_junc.txt.gz
```
The first step filters out multi-mapping reads and reads that do not pass WASP filtering. With the default options in the [STAR wrapper script](../../../rnaseq/src/run_STAR.py), spliced reads are tagged with the `XS` strand attribute, which is parsed in regtools using the `-s 0` option.
A wrapper is provided in [`leafcutter_bam_to_junc.wdl`](leafcutter_bam_to_junc.wdl).

#### 3) Generating intron excision ratios with LeafCutter
The `cluster_prepare_fastqtl.py` script wraps LeafCutter's `leafcutter_cluster_regtools.py` script, applies filters to remove introns with low counts or low complexity, and generates input files formatted for QTL mappers.
```bash
python3 cluster_prepare_fastqtl.py \
    ${junc_files_list} \
    ${exons_list} \
    ${collapsed_annotation} \
    ${prefix} \
    ${sample_participant_map}
```
where `${junc_files_list}` is a text file containing paths to the `*.regtools_junc.txt.gz` files generated in the previous step. The list of exons in ${exons_list} can be generated from the [collapsed](https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model) reference annotation ${collapsed_annotation} using, e.g.,
```python
import pandas as pd
import qtl.annotation

annot = qtl.annotation.Annotation('gencode.v26.GRCh38.genes.gtf')
exon_df = pd.DataFrame([[g.chr, e.start_pos, e.end_pos, g.strand, g.id, g.name]
                        for g in annot.genes for e in g.transcripts[0].exons],
                       columns=['chr', 'start', 'end', 'strand', 'gene_id', 'gene_name'])
exon_df.to_csv('gencode.v26.GRCh38.genes.exons.txt.gz', sep='\t', index=False)
```
`${sample_participant_map}` is a headerless tab-delimited file linking sample IDs to participant IDs (used for renaming files, such that IDs in the output files match those in the VCF).

This step generates a BED file and index (with the BED start/end coordinates corresponding to the TSS), as well as a file mapping phenotype IDs (in the form `${chr}:${intron_start}:${intron_end}:${cluster_id}_${strand}:${gene_id}`) to gene IDs:
```bash
${prefix}.leafcutter.bed.gz
${prefix}.leafcutter.bed.gz.tbi
${prefix}.leafcutter.phenotype_groups.txt
```

#### 4) Mapping sQTLs with tensorQTL
The files generated above can be provided to QTL mappers such as [tensorQTL](https://github.com/broadinstitute/tensorqtl) along with covariates including PEER factors and genotype PCs (see [eQTL pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl)).

The command for running tensorQTL with the phenotype groups defined in `${prefix}.leafcutter.phenotype_groups.txt` is given below:
```python
import pandas as pd
from tensorqtl import cis, post

phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(f'{prefix}.leafcutter.bed.gz')
covariates_df = pd.read_csv(f'{prefix}.combined_covariates.txt', sep='\t', index_col=0).T
group_s = pd.read_csv(f'{prefix}.leafcutter.phenotype_groups.txt',
                      sep='\t', header=None, index_col=0, squeeze=True)
cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                     covariates_df, group_s=group_s)
post.calculate_qvalues(cis_df, fdr=0.05, qvalue_lambda=0.85)
```
