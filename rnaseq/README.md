# RNA-seq pipeline for the [GTEx Consortium](www.gtexportal.org)

This repository contains all components of the RNA-seq pipeline used by the GTEx Consortium, including alignment, expression quantification, and quality control.

## Docker image
The GTEx RNA-seq pipeline is provided as a Docker image, available at https://hub.docker.com/r/broadinstitute/gtex_rnaseq/

To download the image, run:
```bash
docker pull broadinstitute/gtex_rnaseq:V10
```

#### Image contents and pipeline components
The following tools are included in the Docker image:

* [SamToFastq](http://broadinstitute.github.io/picard/command-line-overview.html#SamToFastq): BAM to FASTQ conversion
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): sequencing quality control
* [STAR](https://github.com/alexdobin/STAR): spliced alignment of RNA-seq reads
* [Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates): mark duplicate reads
* [RSEM](http://deweylab.github.io/RSEM): transcript expression quantification
* [bamsync](bamsync): utility for transferring QC flags and re-generating read group IDs when realigning BAMs
* [RNA-SeQC](https://github.com/getzlab/rnaseqc): RNA-seq quality control (metrics and gene-level expression quantification)

Versions used across GTEx releases*:
|         | V7      | V8      | V10      | V11      |
| ------- | ------- | ------- | -------- | -------- |
| STAR    | v2.4.2a | v2.5.3a | v2.7.10a | v2.7.11b |
| RSEM    | v.1.2.22| v1.3.0  | v1.3.3   | v1.3.3   |
| RNA-SeQC| v1.1.8  | v1.1.9  | v2.4.2   | v2.4.3   |
| Genome  | GRCh37  | GRCh38  | GRCh38   | GRCh38   |
| GENCODE | [v19](https://www.gencodegenes.org/human/release_19.html) | [v26](https://www.gencodegenes.org/human/release_26.html) | [v39](https://www.gencodegenes.org/human/release_39.html) | [v47](https://www.gencodegenes.org/human/release_47.html) |

*V9 did not include any RNA-seq updates

##  Setup steps
#### Reference genome and annotation
Reference indexes for STAR and RSEM are needed to run the pipeline. All reference files are available at [gs://gtex-resources](https://console.cloud.google.com/storage/browser/gtex-resources).

GTEx releases from V8 onward are based on the GRCh38/hg38 reference genome. Please see [TOPMed_RNAseq_pipeline.md](https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md) for details and links for this reference. Releases up to V7 were based on the GRCh37/hg19 reference genome ([download](http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta)). 

For hg19-based analyses, the GENCODE annotation should be patched to use Ensembl chromosome names:
```
zcat gencode.v19.annotation.gtf.gz | \
    sed 's/chrM/chrMT/;s/chr//' > gencode.v19.annotation.patched_contigs.gtf
```
The collapsed version for RNA-SeQC was generated with:
```
python collapse_annotation.py --transcript_blacklist gencode19_unannotated_readthrough_blacklist.txt \
    gencode.v19.annotation.patched_contigs.gtf gencode.v19.annotation.patched_contigs.collapsed.gtf
```

#### Building the indexes
The STAR index should be built to match the sequencing read length, specified by the `sjdbOverhang` parameter. GTEx samples were sequenced using a 2x76 bp paired-end sequencing protocol, and the matching `sjdbOverhang` is 75.

```bash
# build the STAR index:
mkdir $path_to_references/star_index_oh75
docker run --rm -v $path_to_references:/data -t broadinstitute/gtex_rnaseq:V10 \
    /bin/bash -c "STAR \
        --runMode genomeGenerate \
        --genomeDir /data/star_index_oh75 \
        --genomeFastaFiles /data/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta \
        --sjdbGTFfile /data/gencode.v39.GRCh38.annotation.gtf \
        --sjdbOverhang 75 \
        --runThreadN 4"

# build the RSEM index:
docker run --rm -v $path_to_references:/data -t broadinstitute/gtex_rnaseq:V10 \
    /bin/bash -c "rsem-prepare-reference \
        /data/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta \
        /data/rsem_reference/rsem_reference \
        --gtf /data/gencode.v39.GRCh38.annotation.gtf \
        --num-threads 4"
```

## Running the pipeline
Individual components of the pipeline can be run using the commands below. It is assumed that the `$path_to_data` directory  contains the input data and reference indexes.

```bash
# BAM to FASTQ conversion
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq \
    /bin/bash -c "/src/run_SamToFastq.py /data/$input_bam -p ${sample_id} -o /data"

# STAR alignment
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq:V10 \
    /bin/bash -c "/src/run_STAR.py \
        /data/star_index_oh75 \
        /data/${sample_id}_1.fastq.gz \
        /data/${sample_id}_2.fastq.gz \
        ${sample_id} \
        --threads 4 \
        --output_dir /tmp/star_out && mv /tmp/star_out /data/star_out"

# sync BAMs (optional; copy QC flags and read group IDs)
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq:V10 \
    /bin/bash -c "/src/run_bamsync.sh \
        /data/$input_bam \
        /data/star_out/${sample_id}.Aligned.sortedByCoord.out.bam \
        /data/star_out/${sample_id}"

# mark duplicates (Picard)
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq:V10 \
    /bin/bash -c "/src/run_MarkDuplicates.py \
        /data/star_out/${sample_id}.Aligned.sortedByCoord.out.patched.bam \
        ${sample_id}.Aligned.sortedByCoord.out.patched.md \
        --output_dir /data"

# RNA-SeQC
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq:V10 \
    /bin/bash -c "/src/run_rnaseqc.py \
    ${sample_id}.Aligned.sortedByCoord.out.patched.md.bam \
    ${genes_gtf} \
    ${genome_fasta} \
    ${sample_id} \
    --output_dir /data"

# RSEM transcript quantification
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq:V10 \
    /bin/bash -c "/src/run_RSEM.py \
        /data/rsem_reference \
        /data/star_out/${sample_id}.Aligned.toTranscriptome.out.bam \
        /data/${sample_id} \
        --threads 4"
```

#### Aggregating outputs
Sample-level outputs in GCT format can be concatenated using `combine_GCTs.py`:
```
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq:V10 \
    /bin/bash -c "python3 /src/combine_GCTs.py \
        ${rnaseqc_tpm_gcts} ${sample_set_id}.rnaseqc_tpm"
```
