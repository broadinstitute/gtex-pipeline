#!/bin/bash
# Author: Francois Aguet

set -euo pipefail
IFS=$'\n\t'

source_bam=$1
target_bam=$2
prefix=$3

output_bam=${prefix}.Aligned.sortedByCoord.out.patched.bam
bamsync ${source_bam} ${target_bam} -o ${output_bam}
samtools index ${output_bam}

# verify that both BAMs have same number of reads
nreads_in=$(samtools idxstats ${target_bam} | cut -f3-4 | awk '{s+=$1+$2}END{print s}')
nreads_out=$(samtools idxstats ${output_bam} | cut -f3-4 | awk '{s+=$1+$2}END{print s}')
if (( ${nreads_in}!=${nreads_out} )); then
    echo "Read numbers do not match: ${nreads_in} vs ${nreads_out}"
	exit 1
fi
