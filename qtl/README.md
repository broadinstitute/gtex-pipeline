<!-- Author: Francois Aguet -->
# eQTL discovery pipeline for the [GTEx Consortium](www.gtexportal.org)

This repository contains all components of the eQTL discovery pipeline used by the GTEx Consortium, including data normalization, QTL mapping, and annotation steps. This document describes the pipeline used for the V7 and V8 data releases; for settings specific to the V6p analyses presented in [[GTEx Consortium, 2017](https://www.nature.com/articles/nature24277)], please see the last section.

## Docker image
The GTEx eQTL pipeline components are provided in a Docker image, available at https://hub.docker.com/r/broadinstitute/gtex_eqtl/

To download the image, run:
```bash
docker pull broadinstitute/gtex_eqtl:V8
```

#### Image contents and pipeline components
The following tools are included in the Docker image:

1. [FastQTL](https://github.com/francois-a/fastqtl): QTL mapping software ([Ongen et al., Bioinformatics, 2016](http://bioinformatics.oxfordjournals.org/content/32/10/1479.abstract))
2. R 3.2
3. Python 3.5

## Prerequisites
The following input files are needed:

* VCF file with genotype information. Must be bgzip compressed and indexed with tabix.
* Expression tables in GCT format. Two tables are needed: read counts and normalized (FPKM or TPM).
* Gene annotation in GTF format.


## Running the pipeline
Additional [documentation](http://gtexportal.org/home/documentationPage#staticTextAnalysisMethods) and details about parameter choices are provided on the [GTEx Portal](gtexportal.org).

This pipeline requires gene-level expression data. A collapsed reference GTF can be generated for this purpose using the [`collapse_annotation.py`](https://github.com/broadinstitute/gtex-pipeline/blob/master/gene_model/collapse_annotation.py) script available in the [gene model](https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model) directory. In the code below, it is assumed that `${annotation_gtf}` was generated using this script.

#### 1) Generate normalized expression in BED format
The expression data are normalized as follows: 
1. Read counts are normalized between samples using TMM ([Robinson & Oshlack, Genome Biology, 2010](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25))
2. Genes are selected based on the following expression thresholds: 
   - ≥0.1 TPM in ≥20% samples AND
   - ≥6 reads (unnormalized) in ≥20% samples
3. Each gene is inverse normal transformed across samples.
```bash
eqtl_prepare_expression.py ${tpm_gct} ${counts_gct} ${annotation_gtf} \
    ${sample_participant_lookup} ${vcf_chr_list} ${prefix} \
    --tpm_threshold 0.1 \
    --count_threshold 6 \
    --sample_frac_threshold 0.2 \
    --normalization_method tmm
```
The file `${vcf_chr_list}` lists the chromosomes in the VCF, and can be generated using
```
tabix --list-chroms ${vcf} > ${vcf_chr_list}
```
The file `${sample_participant_lookup}` must contain two columns, `sample_id` and `participant_id`, mapping IDs in the expression files to IDs in the VCF (these can be the same).

This step generates the following BED file and index:
```bash
${prefix}.expression.bed.gz
${prefix}.expression.bed.gz.tbi
```

#### 2) Calculate PEER factors
```bash
Rscript run_PEER.R ${prefix}.expression.bed.gz ${prefix} ${num_peer}
```
The number of PEER factors was selected as function of sample size (N):
- 15 factors for N < 150
- 30 factors for 150 ≤ N < 250
- 45 factors for 250 ≤ N < 350
- 60 factors for N ≥ 350

For information on how these thresholds were determined, please see the [Supplementary Information](https://media.nature.com/original/nature-assets/nature/journal/v550/n7675/extref/nature24277-s1.pdf) of [[GTEx Consortium, 2017](https://www.nature.com/articles/nature24277)].

This step will generate 3 files:
```bash
${prefix}.PEER_residuals.txt
${prefix}.PEER_alpha.txt
${prefix}.PEER_covariates.txt
```

#### 3) Combine covariates
This step generates a combined covariates file, containing genotype PCs, PEER factors, and additional explicit covariates (e.g., genotyping platform).
```bash
combine_covariates.py ${prefix}.PEER_covariates.txt ${prefix} \
    --genotype_pcs ${genotype_pcs} \
    --add_covariates ${add_covariates}
```
The covariate files should have one covariate per row, with an identifier in the first column, and a header line with sample identifiers. This step will generate the file `${prefix}.combined_covariates.txt`

#### 4) Run FastQTL
A wrapper script for multithreaded execution is provided in the docker image (`/opt/fastqtl/python/run_FastQTL_threaded.py`) and at https://github.com/francois-a/fastqtl
```bash
# nominal pass
run_FastQTL_threaded.py ${vcf} ${prefix}.expression.bed.gz ${prefix} \
    --covariates ${prefix}.combined_covariates.txt \
    --window 1e6 --chunks 100 --threads 16

# permutation pass
run_FastQTL_threaded.py ${vcf} ${prefix}.expression.bed.gz ${prefix} \
    --covariates ${prefix}.combined_covariates.txt \
    --window 1e6 --chunks 100 --threads 16 \
    --permute 1000 10000 
```
The following files will be generated:
```bash
${prefix}.allpairs.txt.gz
${prefix}.egenes.txt.gz
```

### Using docker
The steps described above can be run using docker. This assumes that the `$path_to_data` directory contains all required input files.
```bash
# Docker command for step 1:
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_eqtl:V8 /bin/bash \
    -c "/src/eqtl_prepare_expression.py /data/${tpm_gct} /data/${counts_gct} \
        /data/${annotation_gtf} /data/${sample_participant_lookup} /data/${vcf_chr_list} ${prefix} \
        --tpm_threshold 0.1 --count_threshold 6 --sample_frac_threshold 0.2 --normalization_method tmm"        
```

### V6p pipeline settings

#### Expression normalization
The expression data were normalized as follows: 
1. Genes were selected based on the following exression thresholds: 
   * &gt;0.1 RPKM in ≥10 samples AND
   * ≥6 reads (unnormalized) in ≥10 samples
2. RPKMs were normalized between samples using quantile normalization
3. Each gene was inverse normal transformed across samples.
