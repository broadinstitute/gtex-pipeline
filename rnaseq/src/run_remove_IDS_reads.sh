#!/bin/bash
# Author: Abhishek Choudhary

set -euo pipefail
IFS=$'\n\t'

input_bam=$1
prefix=$2

output_bam=${prefix}.Aligned.toTranscriptome_noIDS.out.bam

samtools view ${input_bam} | awk '{print $1}' | sort | uniq > readids_all
samtools view ${input_bam} | awk '$6 ~ "I|D|S"' | awk '{print $1}' | sort | uniq > readids_flagIDS
comm -23 readids_all readids_flagIDS > readids_flagnonIDS

echo readids all = $(cat readids_all | wc -l)
echo readids for I/D/S flag = $(cat readids_flagIDS | wc -l)
echo readids for non I/D/S flag = $(cat readids_flagnonIDS | wc -l)

samtools view -N readids_flagnonIDS -hb -o ${output_bam} ${input_bam}
