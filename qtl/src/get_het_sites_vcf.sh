#!/bin/bash
# Author: Francois Aguet

# inputs
vcf=$1
sample_id=$2

line_num=$(zcat $vcf | grep -n -m1 CHROM | cut -d':' -f1)

# get column headers (last line of VCF header)
col_header=$(zcat $vcf | head -$line_num | tail -1)

col_num=$(echo $col_header | awk -v sample_id=$sample_id '{for (i=1;i<=NF;++i) {if ($i==sample_id) print i}; exit}')

date +"[%b %d %H:%M:%S] Extracting het sites"
echo "VCF column for $sample_id: $col_num"

# extract het sites for this sample; filter for SNPs and exclude chr. X
cat <( zcat $vcf | head -n $((line_num-2)) )\
    <( echo "${col_header}" | cut -f1-9,$col_num )\
    <( zcat $vcf | awk -F"\t" -v start=${line_num} -v col=${col_num} 'NR>start && ($col=="0|1"||$col=="1|0") {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$col} {OFS="\t"}' ) |\
    bcftools view -v snps -t ^X -t ^chrX | bgzip -c > $sample_id.hets.vcf.gz

tabix $sample_id.hets.vcf.gz
date +"[%b %d %H:%M:%S] Done"
