# RNA-seq pipeline for the [GTEx Consortium](www.gtexportal.org)

This repository contains all components of the RNA-seq pipeline used by the GTEx Consortium, including alignment, expression quantification, and quality control.

## Docker image
The GTEx RNA-seq pipeline is provided as a Docker image, available at https://hub.docker.com/r/broadinstitute/gtex_rnaseq/

To download the image, run:
```bash
docker pull broadinstitute/gtex_rnaseq
```

#### Image contents and pipeline components
The following tools are included in the Docker image:

1. [SamToFastq](http://broadinstitute.github.io/picard/command-line-overview.html#SamToFastq): BAM to FASTQ conversion
2. [STAR](https://github.com/alexdobin/STAR): spliced alignment of RNA sequence reads (v2.4.2a)
3. [RSEM](http://deweylab.github.io/RSEM) transcript expression quantification (v1.2.22)
4. bamsync: utility for transferring QC flags from the input BAM and for re-generating read group IDs

##  Setup steps
#### Reference genome and annotation
To run the pipeline, reference indexes must first be built for STAR and RSEM.

GTEx currently uses the GRCh37/Hg19 genome reference ([download](http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta)) and the GENCODE v19 annotation ([download](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz)).

The GENCODE annotation should be patched to use Ensembl chromosome names:
```
zcat gencode.v19.annotation.gtf.gz | sed 's/chrM/chrMT/;s/chr//' > gencode.v19.annotation.patched_contigs.gtf
```

#### Building the indexes
The STAR index should be built to match the sequencing read length, specified by the _sjdbOverhang_ parameter. GTEx samples were sequenced using a 2x76 bp paired-end sequencing protocol, and _sjdbOverhang_ should therefore be set to 75.

```bash
# To build the STAR index, run:
mkdir $path_to_references/star_index_oh75
docker run --rm -v $path_to_references:/data -t broadinstitute/gtex_rnaseq /bin/bash -c "STAR --runMode genomeGenerate --genomeDir /data/star_index_oh75 --genomeFastaFiles /data/Homo_sapiens_assembly19.fasta --sjdbGTFfile /data/gencode.v19.annotation.patched_contigs.gtf --sjdbOverhang 75 --runThreadN 4"

# To build the RSEM index, run:
docker run --rm -v $path_to_references:/data -t broadinstitute/gtex_rnaseq /bin/bash -c "rsem-prepare-reference --num-threads 4 --gtf /data/gencode.v19.annotation.patched_contigs.gtf /data/Homo_sapiens_assembly19.fasta /data/rsem_reference/rsem_reference"
```

## Running the pipeline
Individual components of the pipeline can be run using the commands below. Note that the `$path_to_data` directory should contain both the input data and the reference indexes.

```bash
# BAM to FASTQ conversion
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq /bin/bash -c "/src/run_SamToFastq.py /data/$input_bam -p $prefix -o /data"

# STAR alignment
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq /bin/bash -c "/src/run_STAR.py /data/star_index_oh75 /data/$fastq1 /data/$fastq2 $prefix --threads 4 --output_dir /tmp/star_out && mv /tmp/star_out /data/star_out"

# sync BAMs (QC flags and read group IDs)
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq /bin/bash -c "/src/run_bamsync.sh /data/$input_bam /data/star_out/$prefix.Aligned.sortedByCoord.out.bam /data/star_out/$prefix"

# RSEM transcript quantification
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_rnaseq /bin/bash -c "/src/run_RSEM.py /data/rsem_reference /data/star_out/$prefix.Aligned.toTranscriptome.out.bam /data/$prefix --threads 1"
```
