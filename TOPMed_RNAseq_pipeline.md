# TOPMed RNA-seq pipeline harmonization summary

### *Note: this is a draft of the next pipeline iteration. Software and reference versions may change.*

This document is intended to ensure reproducibility of analyses and facilitate harmonization of pipelines, defines primary file locations and accessible download locations, and formalizes the TOPMed RNA-seq pipeline components and parameters. All relevant compute parameters are assumed to be software defaults for the versions listed unless specified otherwise (with the exception of directory/file paths, number of threads, etc.).

All wrapper scripts are available from the GTEx pipeline repository: [https://github.com/broadinstitute/gtex-pipeline](https://github.com/broadinstitute/gtex-pipeline)

*Changes to this specification may be proposed based on rigorous benchmarking, pending approval of NHLBI.*

#### Previous versions of the pipeline
* TOPMed [Phase 4 studies](https://www.nhlbiwgs.org/group/project-studies): available [here](https://github.com/broadinstitute/gtex-pipeline/releases/tag/TOPMed_RNAseq_v2).
* TOPMed MESA RNA-seq pilot: available [here](https://github.com/broadinstitute/gtex-pipeline/releases/tag/TOPMed_MESA_RNAseq_pilot).

### Pipeline summary
The GTEx/TOPMed RNA-Seq pipeline generates, for each sample:

1. Aligned RNA-seq reads in BAM format.
2. Standard quality control metrics derived from the FASTQs and aligned reads.
3. Gene-level expression quantifications based on a collapsed version of a reference transcript annotation, provided as read counts and TPM.
4. Transcript-level expression quantifications, provided as TPM, expected read counts, and isoform percentages.

This document also describes the generation of the reference files required for each pipeline component.

### Pipeline components
* Sequencing quality control: [FastQC 0.11.9](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* Alignment: [STAR 2.7.10a](https://github.com/alexdobin/STAR)
  * Post-processing: [Picard 2.27.1](https://github.com/broadinstitute/picard) [MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)
* Gene quantification and quality control: [RNA-SeQC 2.4.2](https://github.com/broadinstitute/rnaseqc)
* Transcript quantification: [RSEM 1.3.3](https://deweylab.github.io/RSEM/)
* Utilities: [SAMtools 1.19.2](https://github.com/samtools/samtools/releases) and [HTSlib 1.19.1](https://github.com/samtools/htslib/releases)

### Reference files
This section describes the GRCh38 reference genome and GENCODE 39 annotation used, including the addition of ERCC spike-in annotations.

The reference files described in this section can be obtained through the following links:
* Reference genome for RNA-seq alignment:
  * [Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.dict](https://console.cloud.google.com/storage/browser/_details/gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.dict)
  * [Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta](https://console.cloud.google.com/storage/browser/_details/gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta)
  * [Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta.fai](https://console.cloud.google.com/storage/browser/_details/gtex-resources/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta.fai)
* Collapsed gene model: [gencode.v39.GRCh38.ERCC.genes.stranded.gtf](https://console.cloud.google.com/storage/browser/_details/gtex-resources/GENCODE/gencode.v39.GRCh38.ERCC.genes.stranded.gtf)
* STAR index: [STARv275a_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v39_oh100.tar.gz](https://console.cloud.google.com/storage/browser/_details/gtex-resources/STAR_genomes/STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v39_oh100.tar.gz)
* RSEM reference: [rsem_reference_GRCh38_gencode39_ercc.tar.gz](https://console.cloud.google.com/storage/browser/_details/gtex-resources/RSEM_references/rsem_reference_GRCh38_gencode39_ercc.tar.gz)

*Note: the reference genome is based on the Broad Institute's GRCh38 reference, which is also used for aligning TOPMed whole genome sequence data.*

#### Genome reference

For RNA-seq analyses, a reference FASTA excluding ALT, HLA, and Decoy contigs was generated (as of the writing of this document, none of the RNA-seq pipeline components were designed to properly handle these regions).

*Note: The reference produced after filtering out ALT, HLA, and Decoy contigs is identical to the one used by ENCODE ([FASTA file](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)). However, the ENCODE reference does not  include ERCC spike-in sequences.*

1. The Broad Institute's [GRCh38 reference](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) can be obtained from the Broad Institute's FTP server (ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/) or from [Google Cloud](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0).
2. ERCC spike-in reference annotations were downloaded from [ThermoFisher](https://tools.thermofisher.com/content/sfs/manuals/ERCC92.zip) (the archive contains two files, ERCC92.fa and ERCC92.gtf) and processed as detailed below. The patched ERCC references are also available [here](https://github.com/broadinstitute/gtex-pipeline/tree/master/rnaseq/references).<br>The ERCC FASTA was patched for compatibility with RNA-SeQC/GATK using
    ```bash
    sed 's/ERCC-/ERCC_/g' ERCC92.fa > ERCC92.patched.fa
    ```
3. ALT, HLA, and Decoy contigs were excluded from the reference genome FASTA using the following Python code:
    ```python
    with open('Homo_sapiens_assembly38.fasta', 'r') as fasta:
        contigs = fasta.read()
    contigs = contigs.split('>')
    contig_ids = [i.split(' ', 1)[0] for i in contigs]

    # exclude ALT, HLA and decoy contigs
    filtered_fasta = '>'.join([c for i,c in zip(contig_ids, contigs)
        if not (i[-4:]=='_alt' or i[:3]=='HLA' or i[-6:]=='_decoy')])

    with open('Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta', 'w') as fasta:
        fasta.write(filtered_fasta)
    ```

4. ERCC spike-in sequences were appended:
    ```bash
    cat Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta ERCC92.patched.fa \
        > Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta
    ```
5. The FASTA index and dictionary (required for Picard/GATK) were generated:
    ```bash
    samtools faidx Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta
    java -jar picard.jar \
        CreateSequenceDictionary \
        R=Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta \
        O=Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.dict
    ```

#### Reference annotation
The reference annotations were prepared as follows:
1. The comprehensive gene annotation was downloaded from [GENCODE](https://www.gencodegenes.org/human/release_39.html): [gencode.v39.annotation.gtf.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz).

2. For gene-level quantifications, the annotation was collapsed with the [script](https://github.com/broadinstitute/gtex-pipeline/blob/master/gene_model/collapse_annotation.py) used in the [GTEx pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model):
    ```bash
    python3 collapse_annotation.py \
        --stranded \
        gencode.v39.GRCh38.annotation.gtf \
        gencode.v39.GRCh38.genes.stranded.gtf
    ```
3. Gene- and transcript-level attributes were added to the ERCC GTF with the following Python code:
    ```python
    with open('ERCC92.gtf') as exon_gtf, open('ERCC92.patched.gtf', 'w') as gene_gtf:
        for line in exon_gtf:
            f = line.strip().split('\t')
            f[0] = f[0].replace('-', '_')  # required for RNA-SeQC/GATK (no '-' in contig name)
            f[5] = '.'

            attr = f[8]
            if attr[-1] == ';':
                attr = attr[:-1]
            attr = dict([i.split(' ') for i in attr.replace('"', '').split('; ')])
            # add gene_name, gene_type
            attr['gene_name'] = attr['gene_id']
            attr['gene_type'] = 'ercc_control'
            attr['gene_status'] = 'KNOWN'
            attr['level'] = 2
            for k in ['id', 'type', 'name', 'status']:
                attr[f'transcript_{k}'] = attr[f'gene_{k}']

            attr_str = []
            for k in ['gene_id', 'transcript_id', 'gene_type', 'gene_status', 'gene_name',
                'transcript_type', 'transcript_status', 'transcript_name']:
                attr_str.append(f'{k} "{attr[k]}";')
            attr_str.append(f"level {attr['level']};")
            f[8] = ' '.join(attr_str)

            # write gene, transcript, exon
            gene_gtf.write('\t'.join(f[:2] + ['gene'] + f[3:]) + '\n')
            gene_gtf.write('\t'.join(f[:2] + ['transcript'] + f[3:]) + '\n')
            f[8] = ' '.join(attr_str[:2])
            gene_gtf.write('\t'.join(f[:2] + ['exon'] + f[3:]) + '\n')
    ```

4. The ERCC annotation was appended to the reference GTFs:
    ```bash
    cat gencode.v39.GRCh38.annotation.gtf ERCC92.patched.gtf \
        > gencode.v39.GRCh38.annotation.ERCC.gtf
    cat gencode.v39.GRCh38.genes.gtf ERCC92.patched.gtf \
        > gencode.v39.GRCh38.ERCC.genes.gtf
    ```

#### STAR index
All RNA-seq samples were sequenced as 2x101bp paired-end, and the STAR index was generated accordingly:
```bash
STAR --runMode genomeGenerate \
    --genomeDir STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v39_oh100 \
    --genomeFastaFiles Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta ERCC92.patched.fa \
    --sjdbGTFfile gencode.v39.GRCh38.annotation.ERCC92.gtf \
    --sjdbOverhang 100 --runThreadN 10
```
#### RSEM reference
The RSEM references were generated using:
```bash
rsem-prepare-reference --num-threads 10 \
    --gtf gencode.v39.GRCh38.annotation.ERCC92.gtf \
    Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta,ERCC92.patched.fa \
    rsem_reference
```

### Installation of pipeline components
This section lists the source repositories and installation instructions for the pipeline components. The instruction replicate those in the pipeline [Dockerfile](https://github.com/broadinstitute/gtex-pipeline/blob/master/rnaseq/Dockerfile).

1. FastQC v0.11.9:
    ```
    cd /opt && \
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && mv FastQC FastQC-0.11.9 && cd FastQC-0.11.9 && chmod 775 fastqc
    PATH=/opt/FastQC-0.11.9:$PATH
    ```

2. STAR 2.7.10a:
    ```
    cd /opt && \
    wget --no-check-certificate https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz && \
    tar -xf 2.7.10a.tar.gz && rm 2.7.10a.tar.gz
    PATH=/opt/STAR-2.7.10a/bin/Linux_x86_64_static:$PATH
    ```

3. Picard v2.27.1 or later (for MarkDuplicates):
    ```
    mkdir /opt/picard-tools && \
    wget --no-check-certificate \
        -P /opt/picard-tools/ \
        https://github.com/broadinstitute/picard/releases/download/2.27.1/picard.jar
    ```

4. RSEM v1.3.3:
    ```
    cd /opt && \
    wget --no-check-certificate \
        https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz && \
    tar -xvf v1.3.3.tar.gz && rm v1.3.3.tar.gz && cd RSEM-1.3.3 && make
    PATH=/opt/RSEM-1.3.1:$PATH
    ```

5. RNA-SeQC v2.4.2:
    ```
    mkdir /opt/rnaseqc && cd /opt/rnaseqc && \
    wget https://github.com/getzlab/rnaseqc/releases/download/v2.4.2/rnaseqc.v2.4.2.linux.gz && \
    gunzip rnaseqc.v2.4.2.linux.gz && mv rnaseqc.v2.4.2.linux rnaseqc && chmod 775 rnaseqc
    PATH=/opt/rnaseqc:$PATH
    pip3 install rnaseqc
    ```

### Pipeline parameters
This section contains the full list of inputs and parameters provided to each method.

The following variables must be defined:
* `star_index`: path to the directory containing the STAR index
* `fastq1` and `fastq2`: paths to the two FASTQ files
* `sample_id`: sample identifier; this will be prepended to output files
* `rsem_reference`: path to the directory containing the RSEM reference
* `genome_fasta`: path to the reference genome (`Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta` as described above)
* `genes_gtf`: path to the collapsed, gene-level GTF (`gencode.v30.GRCh38.ERCC.genes.gtf` as described above)
* `star_bam_file`: name/file path of the BAM generated by the STAR aligner, by default `${sample_id}.Aligned.sortedByCoord.out.bam`.
* `md_bam_file`: name of the BAM generated by Picard MarkDuplicates.

1. STAR
    ```bash
    STAR --runMode alignReads \
        --runThreadN 8\
        --genomeDir ${star_index} \
        --twopassMode Basic \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFilterType BySJout \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --limitSjdbInsertNsj 1200000 \
        --readFilesIn ${fastq1} ${fastq2} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${sample_id} \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs None \
        --alignSoftClipAtReferenceEnds Yes \
        --quantMode TranscriptomeSAM GeneCounts \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --genomeLoad NoSharedMemory \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --chimOutType Junctions WithinBAM SoftClip \
        --chimMainSegmentMultNmax 1 \
        --outSAMattributes NH HI AS nM NM ch \
        --outSAMattrRGline ID:rg1 SM:sm1
    ```
2. MarkDuplicates
    ```bash
    java -jar picard.jar \
        MarkDuplicates I=${star_bam_file} \
        O=${md_bam_file} \
        M=${sample_id}.marked_dup_metrics.txt \
        ASSUME_SORT_ORDER=coordinate
    ```
3. RSEM
    ```bash
    rsem-calculate-expression \
        --num-threads 2 \
        --fragment-length-max 1000 \
        --no-bam-output \
        --paired-end \
        --estimate-rspd \
        --forward-prob 0.0 \
        --bam ${sample_id}.Aligned.toTranscriptome.out.bam \
        ${rsem_reference} ${sample_id}
    ```
4. RNA-SeQC
    ```bash
    rnaseqc ${genes_gtf} ${bam_file} . \
        -s ${sample_id} --stranded rf -vv
    ```

### Appendix: wrapper scripts from the GTEx pipeline
This section provides the commands used to run each step of the pipeline based on the wrapper scripts from the [GTEx pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/rnaseq).

The following variables must be defined:
* `star_index`: path to the directory containing the STAR index
* `fastq1` and `fastq2`: paths to the two FASTQ files
* `sample_id`: sample identifier; this will be prepended to output files
* `rsem_reference`: path to the directory containing the RSEM reference
* `genome_fasta`: path to the reference genome (`Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta` as described above)
* `genes_gtf`: path to the collapsed, gene-level GTF (`gencode.v39.GRCh38.ERCC.genes.gtf` as described above)

1. STAR ([run_STAR.py](rnaseq/src/run_STAR.py))
    ```bash
    python3 run_STAR.py \
        ${star_index} ${fastq1} ${fastq2} ${sample_id} \
        --output_dir star_out \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFilterType BySJout \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --limitSjdbInsertNsj 1200000 \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs None \
        --alignSoftClipAtReferenceEnds Yes \
        --quantMode TranscriptomeSAM GeneCounts \
        --outSAMattrRGline ID:rg1 SM:sm1 \
        --outSAMattributes NH HI AS nM NM ch \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --chimOutType Junctions WithinBAM SoftClip \
        --chimMainSegmentMultNmax 1 \
        --threads 8
    ```
2. MarkDuplicates ([run_MarkDuplicates.py](rnaseq/src/run_MarkDuplicates.py))
    ```bash
    python3 -u run_MarkDuplicates.py ${sample_id}.Aligned.sortedByCoord.out.bam ${sample_id}
    ```
3. RSEM ([run_RSEM.py](rnaseq/src/run_RSEM.py))
    ```bash
    python3 run_RSEM.py \
        --max_frag_len 1000 \
        --estimate_rspd true \
        --is_stranded true \
        --threads 2 \
        ${rsem_reference} ${sample_id}.Aligned.toTranscriptome.out.bam ${sample_id}
    ```
4. RNA-SeQC ([run_rnaseqc.py](rnaseq/src/run_rnaseqc.py))
    ```bash
    python3 run_rnaseqc.py \
        ${genes_gtf}
        ${sample_id}.Aligned.sortedByCoord.out.md.bam \
        ${sample_id} \
        --stranded rf
    ```
