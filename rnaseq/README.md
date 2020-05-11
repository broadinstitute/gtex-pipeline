# RNA-seq pipeline for the [GTEx Consortium](www.gtexportal.org)

This repository contains all components of the RNA-seq pipeline used by the GTEx Consortium, including alignment, expression quantification, and quality control.

## Docker image
The GTEx RNA-seq pipeline is provided as a Docker image, available at https://hub.docker.com/r/broadinstitute/gtex_rnaseq/

To download the image, run:
```bash
docker pull broadinstitute/gtex_rnaseq:V8
```

#### Image contents and pipeline components
The following tools are included in the Docker image (V8 pipeline):

* [SamToFastq](http://broadinstitute.github.io/picard/command-line-overview.html#SamToFastq): BAM to FASTQ conversion
* [STAR](https://github.com/alexdobin/STAR): spliced alignment of RNA sequence reads (v2.5.3a)
* [Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates): mark duplicate reads
* [RSEM](http://deweylab.github.io/RSEM) transcript expression quantification (v1.3.0)
* bamsync: utility for transferring QC flags from the input BAM and for re-generating read group IDs
* [RNA-SeQC](https://github.com/francois-a/rnaseqc): QC metrics and gene-level expression quantification (v1.1.9)

Version in V7 pipeline:
* STAR v2.4.2a
* RSEM v1.2.22
* RNA-SeQC v1.1.8

##  Setup steps
#### Reference genome and annotation
Reference indexes for STAR and RSEM are needed to run the pipeline. All reference files are available at [gs://gtex-resources](https://console.cloud.google.com/storage/browser/gtex-resources).

GTEx releases from V8 onward are based on the GRCh38/hg38 reference genome. Please see [TOPMed_RNAseq_pipeline.md](https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md) for details and links for this reference. Releases up to V7 were based on the GRCh37/hg19 reference genome ([download](http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta)). 

Release V8 uses the [GENCODE v26](https://www.gencodegenes.org/human/release_26.html) annotation. Releases V6/V6p and V7 used [GENCODE v19](https://www.gencodegenes.org/human/release_19.html).

For hg19-based analyses, the GENCODE annotation should be patched to use Ensembl chromosome names:
```
zcat gencode.v19.annotation.gtf.gz | \
    sed 's/chrM/chrMT/;s/chr//' > gencode.v19.annotation.patched_contigs.gtf
```

#### Building the indexes
The STAR index should be built to match the sequencing read length, specified by the `sjdbOverhang` parameter. GTEx samples were sequenced using a 2x76 bp paired-end sequencing protocol, and the matching `sjdbOverhang` is 75.

```bash
# build the STAR index:
mkdir $path_to_references/star_index_oh75
docker run --rm -v $path_to_references:/data -t broadinstitute/gtex_rnaseq:V8 \
    /bin/bash -c "STAR \
        --runMode genomeGenerate \
        --genomeDir /data/star_index_oh75 \
        --genomeFastaFiles /data/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta \
        --sjdbGTFfile /data/gencode.v26.GRCh38.annotation.gtf \
        --sjdbOverhang 75 \
        --runThreadN 4"

# build the RSEM index:
docker run --rm -v $path_to_references:/data -t broadinstitute/gtex_rnaseq:V8 \
    /bin/bash -c "rsem-prepare-reference \
        /data/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta \
        /data/rsem_reference/rsem_reference \
        --gtf /data/gencode.v26.GRCh38.annotation.gtf \
        --num-threads 4"
```

## Running the pipeline
Individual components of the pipeline can be run using the commands below. It is assumed that the `$path_to_data` directory  contains the input data and reference indexes.

```bash
# BAM to FASTQ conversion
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq \
    /bin/bash -c "/src/run_SamToFastq.py /data/$input_bam -p ${sample_id} -o /data"

# STAR alignment
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq:V8 \
    /bin/bash -c "/src/run_STAR.py \
        /data/star_index_oh75 \
        /data/${sample_id}_1.fastq.gz \
        /data/${sample_id}_2.fastq.gz \
        ${sample_id} \
        --threads 4 \
        --output_dir /tmp/star_out && mv /tmp/star_out /data/star_out"

# sync BAMs (optional; copy QC flags and read group IDs)
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq:V8 \
    /bin/bash -c "/src/run_bamsync.sh \
        /data/$input_bam \
        /data/star_out/${sample_id}.Aligned.sortedByCoord.out.bam \
        /data/star_out/${sample_id}"

# mark duplicates (Picard)
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq:V8 \
    /bin/bash -c "/src/run_MarkDuplicates.py \
        /data/star_out/${sample_id}.Aligned.sortedByCoord.out.patched.bam \
        ${sample_id}.Aligned.sortedByCoord.out.patched.md \
        --output_dir /data"

# RNA-SeQC
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq:V8 \
    /bin/bash -c "/src/run_rnaseqc.py \
    ${sample_id}.Aligned.sortedByCoord.out.patched.md.bam \
    ${genes_gtf} \
    ${genome_fasta} \
    ${sample_id} \
    --output_dir /data"

# RSEM transcript quantification
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq:V8 \
    /bin/bash -c "/src/run_RSEM.py \
        /data/rsem_reference \
        /data/star_out/${sample_id}.Aligned.toTranscriptome.out.bam \
        /data/${sample_id} \
        --threads 4"
```

#### Aggregating outputs
Sample-level outputs in GCT format can be concatenated using `combine_GCTs.py`:
```
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq:V8 \
    /bin/bash -c "python3 /src/combine_GCTs.py \
        ${rnaseqc_rpkm_gcts} ${sample_set_id}.rnaseqc_rpkm"
```
