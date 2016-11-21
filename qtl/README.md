<!-- Author: Francois Aguet -->
# eQTL discovery pipeline for the [GTEx Consortium](www.gtexportal.org)

This repository contains all components of the eQTL discovery pipeline used by the GTEx Consortium, including data normalization, QTL mapping, and annotation steps.

## Docker image
The GTEx eQTL pipeline components are provided in a Docker image, available at https://hub.docker.com/r/broadinstitute/gtex_eqtl/

To download the image, run:
```bash
docker pull broadinstitute/gtex_eqtl
```

#### Image contents and pipeline components
The following tools are included in the Docker image:

1. [FastQTL](https://github.com/francois-a/fastqtl): QTL mapping software ([Ongen et al., Bioinformatics, 2016](http://bioinformatics.oxfordjournals.org/content/32/10/1479.abstract))
2. R 3.2
3. Python 3.5

## Prerequisites
The following input files are needed:

* VCF file with genotype information. Must be bgzip compressed and tabix indexed.
* Expression tables in GCT format. Two tables are needed: read counts and normalized (e.g., FPKM or TPM).
* Gene annotation in GTF format.


## Running the pipeline
Additional [documentation](http://gtexportal.org/home/documentationPage#staticTextAnalysisMethods) and details about parameter choices are provided on the [GTEx Portal](gtexportal.org).

#### 1) Generate normalized expression in BED format
The expression data are normalized as follows: (i) expression values are quantile normalized to the average empirical distribution observed across samples; (ii) for each gene, expression values are inverse quantile normalized to a standard normal distribution across samples.
```bash
normalize_expression.py ${normalized_gct} ${counts_gct} ${annotation_gtf} ${vcf} ${prefix} --expression_threshold 0.1 --count_threshold 5 --min_samples 10
```
Using these settings, genes are selected based on expression thresholds of >0.1 RPKM in ≥10 samples and >5 reads in ≥10 samples. This step will generate 5 files:

```bash
${prefix}.expression.bed.gz
${prefix}.expression.bed.gz.tbi
${prefix}.expression.fpkm.bed.gz
${prefix}.expression.fpkm.bed.gz.tbi
${prefix}.expression.txt
```

#### 2) Calculate PEER factors
```bash
Rscript run_PEER.R ${prefix}.expression.txt ${prefix} ${num_peer}
```
This will generate 3 files:
```bash
${prefix}_PEER_residuals.txt
${prefix}_PEER_alpha.txt
${prefix}_PEER_covariates.txt
```

#### 3) Combine covariates
This step generates a combined covariates file, containing genotype PCs, PEER factors, and additional explicit covariates (e.g., genotyping platform).
```bash
combine_covariates.py ${genotype_pcs} ${prefix}_PEER_covariates.txt ${prefix} --add_covariates ${explicit_cov}
```
The covariate files should have one covariate per row, with an identifier in the first column, and a header line with sample identifiers. This step will generate the file `${prefix}.combined_covariates.txt`

#### 4) Run FastQTL
A wrapper script for multithreaded execution is provided in the docker image (`/opt/fastqtl/python/run_FastQTL_threaded.py`) and at https://github.com/francois-a/fastqtl
```bash
# nominal pass
run_FastQTL_threaded.py ${vcf} ${prefix}.expression.bed.gz ${prefix} --covariates ${prefix}.combined_covariates.txt --window 1e6 --ma_sample_threshold 10 --maf_threshold 0.01 --chunks 100 --threads 12

# permutation pass
run_FastQTL_threaded.py ${vcf} ${prefix}.expression.bed.gz ${prefix} --covariates ${prefix}.combined_covariates.txt --permute 1000 10000 --window 1e6 --ma_sample_threshold 10 --maf_threshold 0.01 --chunks 100 --threads 16
```
The following files will be generated:
```bash
${prefix}.allpairs.txt.gz
${prefix}.egenes.txt.gz
```

### Using docker
The steps described above can be run using docker. Note that the `$path_to_data` directory should the required input files.

```bash
# Docker command for step 1:
docker run --rm -v $path_to_data:/data -t broadinstitute/gtex_eqtl /bin/bash -c "/src/normalize_expression.py /data/${normalized_gct} /data/${counts_gct} /data/${annotation_gtf} /data/${vcf} ${prefix} --expression_threshold 0.1 --count_threshold 5 --min_samples 10"
```
