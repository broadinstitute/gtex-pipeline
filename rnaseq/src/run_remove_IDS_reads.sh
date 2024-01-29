#!/bin/bash
# Author: Francois Aguet

set -euo pipefail
IFS=$'\n\t'

input_bam=$1
prefix=$2

output_bam=${prefix}.Aligned.toTranscriptome_noIDS.out.bam

samtools view ${input_bam} | awk '{print $1}' | sort | uniq > readids_all
samtools view ${input_bam} | awk '$6 ~ "I"' | awk '{print $1}' | sort | uniq > readids_flagI
samtools view ${input_bam} | awk '$6 ~ "D"' | awk '{print $1}' | sort | uniq > readids_flagD
samtools view ${input_bam} | awk '$6 ~ "S"' | awk '{print $1}' | sort | uniq > readids_flagS
cat readids_flagI readids_flagD readids_flagS | sort | uniq > readids_flagIDS
comm -23 readids_all readids_flagIDS > readids_flagnonIDS

echo readids all = $(cat readids_all | wc -l)
echo readids for I flag = $(cat readids_flagI | wc -l)
echo readids for D flag = $(cat readids_flagD | wc -l)
echo readids for S flag = $(cat readids_flagS | wc -l)
echo readids for I/D/S flag = $(cat readids_flagIDS | wc -l)
echo readids for non I/D/S flag = $(cat readids_flagnonIDS | wc -l)

samtools view -N readids_flagnonIDS -hb -o ${output_bam} ${input_bam}
samtools index ${output_bam}

