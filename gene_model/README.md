<!-- Author: Francois Aguet -->

This repository contains utilities for the generation of gene models and annotations used in the RNA-seq and eQTL pipelines.

## Collapsed gene model

Gene-level expression and eQTLs from the GTEx project are calculated based on a collapsed gene model (i.e., combining all isoforms of a gene into a single transcript), according to the following rules:

1. Transcripts annotated as “retained_intron” or “read_through” are excluded. Additionally, transcripts that overlap with annotated read-through transcripts may be blacklisted (blacklists for GENCODE v19, 24 & 25 are provided in this repository; no transcripts were blacklisted for v26).
2. The union of all exon intervals of each gene is calculated.
3. Overlapping intervals between genes are excluded from all genes.

Command:
```bash
python3 collapse_annotation.py gencode.v26.GRCh38.annotation.gtf gencode.v26.GRCh38.genes.gtf
```
where `gencode.v26.GRCh38.annotation.gtf` is the GTF from [GENCODE](https://www.gencodegenes.org/releases/current.html).

Further documentation is available on the [GTEx Portal](https://gtexportal.org/home/documentationPage#staticTextAnalysisMethods).
